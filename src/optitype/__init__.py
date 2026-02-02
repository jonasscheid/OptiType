"""
OptiType: Precision HLA typing from next-generation sequencing data.

OptiType is a novel HLA genotyping algorithm based on integer linear
programming, capable of producing accurate 4-digit HLA genotyping predictions
from NGS data by simultaneously selecting all minor and major HLA-I alleles.
"""

__version__ = "2.0.0"
__author__ = "Andras Szolek, Benjamin Schubert, Christopher Mohr, Jonas Scheid"

from optitype.api import HLATypingConfig, HLATypingResult, run_hla_typing

__all__ = [
    "__version__",
    "run_hla_typing",
    "HLATypingResult",
    "HLATypingConfig",
]
