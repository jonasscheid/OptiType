"""
Core HLA typing logic for OptiType.

This module handles coverage matrix construction, statistical calculations,
allele pruning, and visualization.
"""

import logging
import warnings
from collections.abc import Callable

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from optitype._logging import elapsed

logger = logging.getLogger(__name__)


def feature_order(feature: tuple[str, int]) -> int:
    """
    Return sort order for feature tuples.

    Returns:
        0 for ('UTR', 1)
        2 for exon, 1
        3 for intron, 1
        4 for exon, 2
        ...
        999 for UTR, 2
    """
    feat_type, feat_number = feature
    if feat_type not in ("intron", "exon", "UTR"):
        raise ValueError(f"Unknown feature type: {feat_type!r}")
    if not isinstance(feat_number, int):
        raise TypeError(f"Feature number must be int, got {type(feat_number).__name__}")
    if feature == ("UTR", 1):
        return 0
    elif feature == ("UTR", 2):
        return 999
    else:
        return feat_number * 2 + (feat_type == "intron")


def get_compact_model(
    hit_df: pd.DataFrame,
    weak_hit_df: pd.DataFrame | None = None,
    weight: float | None = None,
) -> tuple[pd.DataFrame, dict]:
    """
    Turn a binary hit matrix into a compact model.

    Removes duplicate rows and creates an occurrence vector with the number
    of rows each representative read represents.

    Args:
        hit_df: Binary hit matrix DataFrame.
        weak_hit_df: Optional DataFrame of "weak" hits (e.g., unpaired reads).
        weight: Weight for weak hits (0 < weight <= 1).

    Returns:
        Tuple of (unique_matrix, occurrence_dict).
    """
    hit_df = hit_df.loc[hit_df.any(axis=1)]  # Remove all-zero rows
    occurrence = {
        r[0]: len(r) for r in hit_df.groupby(hit_df.columns.tolist()).groups.values()
    }

    if weak_hit_df is not None:
        weak_hit_df = weak_hit_df.loc[weak_hit_df.any(axis=1)]
        if not (0 < weight <= 1):
            raise ValueError(f"weak hit weight must be in (0, 1], got {weight}")
        weak_occ = {
            r[0]: len(r) * weight
            for r in weak_hit_df.groupby(weak_hit_df.columns.tolist()).groups.values()
        }
        occurrence.update(weak_occ)
        unique_mtx = pd.concat([hit_df.drop_duplicates(), weak_hit_df.drop_duplicates()])
    else:
        unique_mtx = hit_df.drop_duplicates()

    return unique_mtx, occurrence


def mtx_to_sparse_dict(hit_df: pd.DataFrame) -> dict:
    """
    Convert a hit matrix to a sparse dictionary.

    Creates dictionary of (read, allele): 1 tuples for hits.

    Args:
        hit_df: Hit matrix DataFrame.

    Returns:
        Dictionary mapping (read_id, allele_id) to 1 for each hit.
    """
    all_hits = {}
    for read_id, alleles in hit_df.iterrows():
        hit_alleles = alleles[alleles != 0].index
        for hit_allele in hit_alleles:
            all_hits[(read_id, hit_allele)] = 1
    return all_hits


def prune_identical_alleles(
    binary_mtx: pd.DataFrame, report_groups: bool = False
) -> pd.DataFrame | tuple[pd.DataFrame, dict]:
    """
    Remove identical allele columns from the binary matrix.

    Args:
        binary_mtx: Binary hit matrix.
        report_groups: If True, also return a dict mapping representatives to their groups.

    Returns:
        Pruned matrix, or tuple of (pruned_matrix, groups_dict) if report_groups=True.
    """
    hash_columns = binary_mtx.transpose().dot(np.random.rand(binary_mtx.shape[0]))

    if report_groups:
        grouper = hash_columns.groupby(hash_columns)
        groups = {g[1].index[0]: g[1].index.tolist() for g in grouper}

    alleles_to_keep = hash_columns.drop_duplicates().index

    if report_groups:
        return binary_mtx[alleles_to_keep], groups
    return binary_mtx[alleles_to_keep]


def prune_identical_reads(binary_mtx: pd.DataFrame) -> pd.DataFrame:
    """
    Remove identical read rows from the binary matrix.

    This should only be used to compactify a matrix before passing it to
    the function finding overshadowed alleles. Final compactifying should
    be done on the original binary matrix.

    Args:
        binary_mtx: Binary hit matrix.

    Returns:
        Pruned matrix with unique rows.
    """
    reads_to_keep = (
        binary_mtx.dot(np.random.rand(binary_mtx.shape[1])).drop_duplicates().index
    )
    return binary_mtx.loc[reads_to_keep]


def prune_overshadowed_alleles(binary_mtx: pd.DataFrame) -> pd.Index:
    """
    Find and remove alleles that are overshadowed by others.

    An allele is overshadowed if another allele has all of its hits plus more.

    For optimal performance, first prune identical columns and rows:
    prune_overshadowed_alleles(prune_identical_alleles(prune_identical_reads(binary_mtx)))

    Args:
        binary_mtx: Binary hit matrix.

    Returns:
        Index of non-overshadowed alleles.
    """
    # Handle potential overflow in dot product
    if (binary_mtx.shape[0] < np.iinfo(np.uint16).max) or (
        binary_mtx.sum(axis=0).max() < np.iinfo(np.uint16).max
    ):
        bb = (
            binary_mtx
            if all(binary_mtx.dtypes == np.uint16)
            else binary_mtx.astype(np.uint16)
        )
    else:
        bb = binary_mtx.astype(np.uint32)

    covariance = bb.transpose().dot(bb)
    diagonal = pd.Series(
        [covariance[ii][ii] for ii in covariance.columns], index=covariance.columns
    )

    new_covariance = covariance[covariance.columns].copy()
    for ii in new_covariance.columns:
        new_covariance.loc[ii, ii] = 0

    overshadowed = []
    for ii in new_covariance.columns:
        potential_superiors = new_covariance[ii][new_covariance[ii] == diagonal[ii]].index
        if any(diagonal[potential_superiors] > diagonal[ii]):
            overshadowed.append(ii)

    non_overshadowed = covariance.columns.difference(overshadowed)
    return non_overshadowed


def create_paired_matrix(
    binary_1: pd.DataFrame,
    binary_2: pd.DataFrame,
    id_cleaning: Callable | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Create paired, mispaired, and unpaired matrices from two binary hit matrices.

    Args:
        binary_1: Binary hit matrix for first read end.
        binary_2: Binary hit matrix for second read end.
        id_cleaning: Optional function to clean read IDs for matching.

    Returns:
        Tuple of (paired, mispaired, unpaired) DataFrames.
    """
    if id_cleaning is not None:
        binary_1.index = list(map(id_cleaning, binary_1.index))
        binary_2.index = list(map(id_cleaning, binary_2.index))

    common_read_ids = binary_1.index.intersection(binary_2.index)
    only_1 = binary_1.index.difference(binary_2.index)
    only_2 = binary_2.index.difference(binary_1.index)

    b_1 = binary_1.loc[common_read_ids]
    b_2 = binary_2.loc[common_read_ids]
    b_12 = b_1 * b_2  # Elementwise AND

    b_ispaired = b_12.any(axis=1)
    b_paired = b_12.loc[b_ispaired]
    b_mispaired = b_1.loc[~b_ispaired] + b_2.loc[~b_ispaired]
    b_unpaired = pd.concat([binary_1.loc[only_1], binary_2.loc[only_2]])

    logger.debug(
        "%s Alignment pairing completed. %d paired, %d unpaired, %d discordant",
        elapsed(), b_paired.shape[0], b_unpaired.shape[0], b_mispaired.shape[0],
    )

    return b_paired, b_mispaired, b_unpaired


def get_features(
    allele_id: str, features: pd.DataFrame, feature_list: list[tuple[str, int]]
) -> pd.DataFrame:
    """
    Get features for an allele from the features DataFrame.

    Args:
        allele_id: Allele ID (can be HLA12345 or HLA12345_HLA67890 for reconstructed).
        features: Features DataFrame.
        feature_list: List of (feature_type, feature_number) tuples to include.

    Returns:
        DataFrame of features to include.
    """
    if "_" in allele_id:
        partial_allele, complete_allele = allele_id.split("_")
    else:
        complete_allele = allele_id
        partial_allele = None

    feats_complete = {
        (of["feature"], of["number"]): of
        for _, of in features.loc[features["id"] == complete_allele].iterrows()
    }

    if partial_allele:
        feats_partial = {
            (of["feature"], of["number"]): of
            for _, of in features.loc[features["id"] == partial_allele].iterrows()
        }
    else:
        feats_partial = feats_complete

    feats_to_include = []
    for feat in sorted(feature_list, key=feature_order):
        if feat in feats_partial:
            feats_to_include.append(feats_partial[feat])
        elif feat in feats_complete:
            feats_to_include.append(feats_complete[feat])
        else:
            warnings.warn(f"Feature {feat} not found for allele {allele_id}")

    return pd.DataFrame(feats_to_include)


def calculate_coverage(
    alignment: tuple,
    features: pd.DataFrame,
    alleles_to_plot: pd.Series,
    features_used: list[tuple[str, int]],
) -> list[tuple[str, np.ndarray]]:
    """
    Calculate coverage matrices for alleles.

    Args:
        alignment: Tuple of alignment data (pos, read_details, ...).
        features: Features DataFrame.
        alleles_to_plot: Series of allele IDs to plot.
        features_used: List of features to use.

    Returns:
        List of (allele_id, coverage_matrix) tuples.
    """
    if len(alignment) not in (2, 4, 5):
        raise ValueError(f"Alignment tuple must have 2, 4 or 5 elements, got {len(alignment)}")
    has_pairing_info = len(alignment) == 5

    if len(alignment) == 2:
        matrix_pos, read_details = alignment[:2]
        pos = matrix_pos[alleles_to_plot]
        pairing = np.sign(pos) * 2
        hit_counts = np.sign(pos).sum(axis=1)
        max_ambiguity = max(1, hit_counts.max())
        to_process = [(pos, read_details, hit_counts, pairing)]
    elif len(alignment) == 4:
        matrix_pos1, read_details1, matrix_pos2, read_details2 = alignment
        pos1 = matrix_pos1[alleles_to_plot]
        pos2 = matrix_pos2[alleles_to_plot]
        pairing1 = np.sign(pos1) * 2
        pairing2 = np.sign(pos2) * 2
        hit_counts1 = np.sign(pos1).sum(axis=1)
        hit_counts2 = np.sign(pos2).sum(axis=1)
        max_ambiguity = max(1, hit_counts1.max(), hit_counts2.max())
        to_process = [
            (pos1, read_details1, hit_counts1, pairing1),
            (pos2, read_details2, hit_counts2, pairing2),
        ]
    else:
        matrix_pos1, read_details1, matrix_pos2, read_details2, pairing_binaries = alignment
        bin_p, bin_u, bin_m = pairing_binaries
        pos1 = matrix_pos1[alleles_to_plot]
        pos2 = matrix_pos2[alleles_to_plot]
        pairing = pd.concat([
            bin_p[alleles_to_plot],
            bin_u[alleles_to_plot] * 2,
            bin_m[alleles_to_plot] * 3,
        ])
        pairing1 = pairing.loc[pos1.index]
        pairing2 = pairing.loc[pos2.index]
        hit_counts1 = np.sign(pairing1).sum(axis=1)
        hit_counts2 = np.sign(pairing2).sum(axis=1)
        max_ambiguity = max(1, hit_counts1.max(), hit_counts2.max())
        to_process = [
            (pos1, read_details1, hit_counts1, pairing1),
            (pos2, read_details2, hit_counts2, pairing2),
        ]

    coverage_matrices = []

    for allele in alleles_to_plot:
        allele_features = get_features(allele, features, features_used)
        allele_length = allele_features["length"].sum()

        coverage = np.zeros((2, 3, max_ambiguity, allele_length), dtype=int)

        for pos, read_details, hit_counts, pairing_info in to_process:
            reads = pos[pos[allele] != 0].index
            for i_pos, i_read_length, i_mismatches, i_hitcount, i_pairing in zip(
                pos.loc[reads][allele],
                read_details.loc[reads]["read_length"],
                read_details.loc[reads]["mismatches"],
                hit_counts[reads],
                pairing_info.loc[reads][allele],
            ):
                if not i_pairing:
                    continue
                coverage[int(bool(i_mismatches))][i_pairing - 1][i_hitcount - 1][
                    i_pos - 1 : i_pos - 1 + i_read_length
                ] += 1

        coverage_matrices.append((allele, coverage))

    return coverage_matrices


def plot_coverage(
    outfile: str,
    coverage_matrices: list[tuple[str, np.ndarray]],
    allele_data: pd.DataFrame,
    features: pd.DataFrame,
    features_used: list[tuple[str, int]],
    columns: int = 2,
) -> None:
    """
    Generate coverage plot PDF.

    Args:
        outfile: Output PDF path.
        coverage_matrices: List of (allele_id, coverage_matrix) tuples.
        allele_data: Allele table DataFrame.
        features: Features DataFrame.
        features_used: List of features used.
        columns: Number of columns in the plot.
    """
    def start_end_zeros(cov_array):
        return np.append(np.append([0], cov_array), [0])

    def allele_sorter(allele_cov_mtx_tuple):
        allele, _ = allele_cov_mtx_tuple
        return allele_data.loc[allele.split("_")[0]]["type"]

    def get_allele_locus(allele):
        return allele_data.loc[allele.split("_")[0]]["locus"]

    number_of_loci = len(set(get_allele_locus(allele) for allele, _ in coverage_matrices))

    dpi = 50
    box_size = (7, 1)
    subplot_rows = 3 * number_of_loci + 1

    area_colors = [
        (0.26, 0.76, 0.26),  # perfect, paired, unique
        (0.40, 0.84, 0.40),  # perfect, paired, shared
        (0.99, 0.75, 0.20),  # perfect, unpaired, unique
        (0.99, 0.75, 0.20),  # perfect, mispaired, unique
        (0.99, 0.85, 0.35),  # perfect, unpaired, shared
        (0.99, 0.85, 0.35),  # perfect, mispaired, shared
        (0.99, 0.23, 0.23),  # mismatch, paired, unique
        (0.99, 0.49, 0.49),  # mismatch, paired, shared
        (0.14, 0.55, 0.72),  # mismatch, unpaired, unique
        (0.14, 0.55, 0.72),  # mismatch, mispaired, unique
        (0.33, 0.70, 0.88),  # mismatch, unpaired, shared
        (0.33, 0.70, 0.88),  # mismatch, mispaired, shared
    ]

    figure = plt.figure(figsize=(box_size[0] * columns, box_size[1] * subplot_rows), dpi=dpi)

    coverage_matrices = sorted(coverage_matrices, key=allele_sorter)
    prev_locus = ""
    i_locus = -1

    for allele, coverage in coverage_matrices:
        if "_" in allele:
            partial, complete = allele.split("_")
            plot_title = f"{allele_data.loc[partial]['type']} (introns from {allele_data.loc[complete]['type']})"
        else:
            plot_title = allele_data.loc[allele]["type"]

        if prev_locus != get_allele_locus(allele):
            i_locus += 1
            i_allele_in_locus = 0
        else:
            i_allele_in_locus = 1

        prev_locus = get_allele_locus(allele)

        plot = plt.subplot2grid(
            (subplot_rows, columns),
            (3 * i_locus, i_allele_in_locus),
            rowspan=3,
            adjustable="box",
        )

        _, _, max_ambig, seq_length = coverage.shape
        shared_weighting = np.reciprocal(np.arange(max_ambig) + 1.0)
        shared_weighting[0] = 0

        perfect_paired_unique = start_end_zeros(coverage[0][0][0])
        mismatch_paired_unique = start_end_zeros(coverage[1][0][0])
        perfect_unpaired_unique = start_end_zeros(coverage[0][1][0])
        mismatch_unpaired_unique = start_end_zeros(coverage[1][1][0])
        perfect_mispaired_unique = start_end_zeros(coverage[0][2][0])
        mismatch_mispaired_unique = start_end_zeros(coverage[1][2][0])

        perfect_paired_shared = start_end_zeros(shared_weighting.dot(coverage[0][0]))
        mismatch_paired_shared = start_end_zeros(shared_weighting.dot(coverage[1][0]))
        perfect_unpaired_shared = start_end_zeros(shared_weighting.dot(coverage[0][1]))
        mismatch_unpaired_shared = start_end_zeros(shared_weighting.dot(coverage[1][1]))
        perfect_mispaired_shared = start_end_zeros(shared_weighting.dot(coverage[0][2]))
        mismatch_mispaired_shared = start_end_zeros(shared_weighting.dot(coverage[1][2]))

        # Exon annotation
        i_start = 1
        for _, ft in get_features(allele, features, features_used).iterrows():
            if ft["feature"] == "exon":
                plot.axvspan(i_start, i_start + ft["length"], facecolor="black", alpha=0.1, linewidth=0, zorder=1)
            i_start += ft["length"]

        areas = plot.stackplot(
            np.arange(seq_length + 2),
            perfect_paired_unique + 0.001,
            perfect_paired_shared,
            perfect_unpaired_unique,
            perfect_mispaired_unique,
            perfect_unpaired_shared,
            perfect_mispaired_shared,
            mismatch_paired_unique,
            mismatch_paired_shared,
            mismatch_unpaired_unique,
            mismatch_mispaired_unique,
            mismatch_unpaired_shared,
            mismatch_mispaired_shared,
            linewidth=0,
            colors=area_colors,
            zorder=5,
        )

        for aa in areas:
            aa.set_edgecolor(aa.get_facecolor())

        plot.tick_params(axis="both", labelsize=10, direction="out", which="both", top=False)
        plot.text(
            0.015, 0.97, plot_title,
            horizontalalignment="left",
            verticalalignment="top",
            transform=plot.transAxes,
            fontsize=10,
            zorder=6,
        )
        _, _, _, y2 = plot.axis()
        plot.axis((0, seq_length, 1, y2))
        plot.set_yscale("log")
        plot.set_ylim(bottom=0.5)

    # Legend
    legend = plt.subplot2grid((subplot_rows, columns), (subplot_rows - 1, 0), colspan=2, adjustable="box")
    legend.add_patch(mpatches.Rectangle((0, 2), 2, 2, color=area_colors[0]))
    legend.add_patch(mpatches.Rectangle((0, 0), 2, 2, color=area_colors[1]))
    legend.add_patch(mpatches.Rectangle((25, 2), 2, 2, color=area_colors[2]))
    legend.add_patch(mpatches.Rectangle((25, 0), 2, 2, color=area_colors[4]))
    legend.add_patch(mpatches.Rectangle((50, 2), 2, 2, color=area_colors[6]))
    legend.add_patch(mpatches.Rectangle((50, 0), 2, 2, color=area_colors[7]))
    legend.add_patch(mpatches.Rectangle((75, 2), 2, 2, color=area_colors[8]))
    legend.add_patch(mpatches.Rectangle((75, 0), 2, 2, color=area_colors[10]))
    legend.text(2.5, 3, "paired, no mismatches, unique", va="center", size="smaller")
    legend.text(2.5, 1, "paired, no mismatches, ambiguous", va="center", size="smaller")
    legend.text(27.5, 3, "unpaired, no mismatches, unique", va="center", size="smaller")
    legend.text(27.5, 1, "unpaired, no mismatches, ambiguous", va="center", size="smaller")
    legend.text(52.5, 3, "paired, mismatched, unique", va="center", size="smaller")
    legend.text(52.5, 1, "paired, mismatched, ambiguous", va="center", size="smaller")
    legend.text(77.5, 3, "unpaired, mismatched, unique", va="center", size="smaller")
    legend.text(77.5, 1, "unpaired, mismatched, ambiguous", va="center", size="smaller")
    legend.set_xlim(0, 100)
    legend.set_ylim(0, 4)
    legend.axison = False

    figure.tight_layout()
    figure.savefig(outfile)
    plt.close(figure)
