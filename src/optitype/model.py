"""
OptiType ILP model for HLA typing.

This module implements the integer linear programming model for optimal
HLA allele selection using Pyomo.
"""

import itertools
from collections import defaultdict

import pandas as pd
from pyomo.environ import (
    Binary,
    ConcreteModel,
    Constraint,
    ConstraintList,
    Objective,
    Param,
    Set,
    Var,
    maximize,
)
from pyomo.opt import SolverFactory, TerminationCondition


class OptiType:
    """
    OptiType ILP model for HLA genotyping.

    Uses integer linear programming to select the optimal combination
    of HLA alleles that best explain the observed read mappings.
    """

    def __init__(
        self,
        cov: dict,
        occ: dict,
        groups_4digit: dict,
        allele_table: pd.DataFrame,
        beta: float,
        t_max_allele: int = 2,
        solver: str = "glpk",
        threads: int = 1,
        verbosity: bool = False,
    ):
        """
        Initialize the OptiType model.

        Args:
            cov: Coverage dictionary mapping (read, allele) to 1 for hits.
            occ: Occurrence dictionary mapping read patterns to counts.
            groups_4digit: Dictionary mapping 4-digit types to allele lists.
            allele_table: DataFrame with allele information.
            beta: Homozygosity detection parameter (0.0 to 0.1).
            t_max_allele: Maximum alleles per locus (default 2).
            solver: ILP solver to use ('glpk', 'cplex', 'cbc').
            threads: Number of solver threads.
            verbosity: Whether to print verbose output.
        """
        self._allele_table = allele_table
        self._beta = float(beta)
        self._t_max_allele = t_max_allele
        self._solver = SolverFactory(solver)
        self._threads = threads
        self._opts = {"threads": threads} if threads > 1 else {}
        self._verbosity = verbosity
        self._changed = True
        self._ks = 1
        self._groups_4digit = groups_4digit

        # Group alleles by locus
        loci_alleles = defaultdict(list)
        for type_4digit, group_alleles in groups_4digit.items():
            loci_alleles[type_4digit.split("*")[0]].extend(group_alleles)
        loci = loci_alleles

        self._allele_to_4digit = {
            allele: type_4digit
            for type_4digit, group in groups_4digit.items()
            for allele in group
        }

        # Build the ILP model
        model = ConcreteModel()

        # Initialize Sets
        model.LociNames = Set(initialize=list(loci.keys()))
        model.Loci = Set(model.LociNames, initialize=lambda m, l: loci[l])

        L = list(itertools.chain(*list(loci.values())))
        reconst = {allele_id: 0.01 for allele_id in L if "_" in allele_id}
        R = set(r for (r, _) in cov.keys())
        model.L = Set(initialize=L)
        model.R = Set(initialize=R)

        # Initialize Parameters
        model.cov = Param(model.R, model.L, initialize=lambda model, r, a: cov.get((r, a), 0))
        model.reconst = Param(model.L, initialize=lambda model, a: reconst.get(a, 0))
        model.occ = Param(model.R, initialize=occ)
        model.t_allele = Param(initialize=self._t_max_allele, mutable=True)
        model.beta = Param(
            initialize=self._beta,
            validate=lambda val, model: 0.0 <= float(self._beta) <= 0.999,
            mutable=True,
        )
        model.nof_loci = Param(initialize=len(loci))

        # Initialize Variables
        model.x = Var(model.L, domain=Binary)
        model.y = Var(model.R, domain=Binary)
        model.re = Var(model.R, bounds=(0.0, None))
        model.hetero = Var(bounds=(0.0, model.nof_loci))

        # Objective function
        model.read_cov = Objective(
            rule=lambda model: sum(
                model.occ[r] * (model.y[r] - model.beta * model.re[r]) for r in model.R
            )
            - sum(model.reconst[a] * model.x[a] for a in model.L),
            sense=maximize,
        )

        # Constraints
        model.max_allel_selection = Constraint(
            model.LociNames,
            rule=lambda model, l: sum(model.x[a] for a in model.Loci[l]) <= model.t_allele,
        )
        model.min_allel_selection = Constraint(
            model.LociNames,
            rule=lambda model, l: sum(model.x[a] for a in model.Loci[l]) >= 1,
        )
        model.is_read_cov = Constraint(
            model.R,
            rule=lambda model, r: sum(model.cov[r, a] * model.x[a] for a in model.L)
            >= model.y[r],
        )
        model.heterozygot_count = Constraint(
            rule=lambda model: model.hetero >= sum(model.x[a] for a in model.L) - model.nof_loci
        )

        # Regularization constraints
        model.reg1 = Constraint(
            model.R, rule=lambda model, r: model.re[r] <= model.nof_loci * model.y[r]
        )
        model.reg2 = Constraint(model.R, rule=lambda model, r: model.re[r] <= model.hetero)
        model.reg3 = Constraint(
            model.R,
            rule=lambda model, r: model.re[r]
            >= model.hetero - model.nof_loci * (1 - model.y[r]),
        )

        # Constraint list for solution enumeration
        model.c = ConstraintList()

        self._instance = model

    def set_beta(self, beta: float) -> None:
        """Set the beta parameter for homozygosity detection."""
        self._changed = True
        getattr(self._instance, str(self._instance.beta)).set_value(float(beta))

    def set_t_max_allele(self, t_max_allele: int) -> None:
        """Set the upper bound of alleles selected per locus."""
        self._changed = True
        getattr(self._instance, str(self._instance.t_allele)).set_value(t_max_allele)

    def solve(self, ks: int) -> pd.DataFrame:
        """
        Solve the ILP model and enumerate top k solutions.

        Args:
            ks: Number of solutions to enumerate.

        Returns:
            DataFrame with typing results.
        """
        d = defaultdict(list)

        if self._changed or self._ks != ks:
            self._ks = ks

            for k in range(ks):
                expr = 0

                try:
                    res = self._solver.solve(
                        self._instance, options=self._opts, tee=self._verbosity
                    )
                except Exception:
                    print(
                        "WARNING: Solver does not support multi-threading. "
                        "Please change the config file accordingly. "
                        "Falling back to single-threading."
                    )
                    res = self._solver.solve(self._instance, options={}, tee=self._verbosity)

                self._instance.solutions.load_from(res)

                if res.solver.termination_condition != TerminationCondition.optimal:
                    print("Optimal solution hasn't been obtained. This is a terminal problem.")
                    break

                selected = []
                indices = []
                encountered_4digit = []

                for j in self._instance.x:
                    if self._allele_to_4digit[j][0] in "HJG":
                        if 0.99 <= self._instance.x[j].value <= 1.01:
                            selected.append(j)
                        indices.append(j)
                        continue

                    if 0.99 <= self._instance.x[j].value <= 1.01:
                        selected.append(j)
                        exp_i = 0
                        exp_i += self._instance.x[j]

                        if self._allele_to_4digit[j] in encountered_4digit:
                            continue

                        encountered_4digit.append(self._allele_to_4digit[j])
                        for i_allele in self._groups_4digit[self._allele_to_4digit[j]]:
                            if self._instance.x[i_allele].value <= 0:
                                exp_i += self._instance.x[i_allele]
                            indices.append(i_allele)
                        expr += 1.0 - exp_i

                zero_indices = set(j for j in self._instance.x).difference(set(indices))
                for j in zero_indices:
                    expr += self._instance.x[j]

                self._instance.c.add(expr >= 1)

                # Build result dictionary
                aas = [self._allele_to_4digit[x].split("*")[0] for x in selected]
                c = dict.fromkeys(aas, 1)

                for i in range(len(aas)):
                    if aas.count(aas[i]) < 2:
                        d[aas[i] + "1"].append(selected[i])
                        d[aas[i] + "2"].append(selected[i])
                    else:
                        d[aas[i] + str(c[aas[i]])].append(selected[i])
                        c[aas[i]] += 1

                nof_reads = sum(
                    self._instance.occ[j] * self._instance.y[j].value for j in self._instance.y
                )
                d["obj"].append(self._instance.read_cov())
                d["nof_reads"].append(nof_reads)

            self._instance.c.clear()
            self._changed = False
            self._enumeration = pd.DataFrame(d)

            return self._enumeration
        else:
            return self._enumeration

    def solve_for_k_alleles(self, k: int, ks: int = 1) -> pd.DataFrame:
        """
        Solve with exactly k alleles (experimental).

        Args:
            k: Exact number of alleles to select.
            ks: Number of solutions to enumerate.

        Returns:
            DataFrame with typing results.
        """
        if k < int(self._instance.nof_loci.value) or k > int(
            self._instance.nof_loci.value * self._t_max_allele
        ):
            raise ValueError(
                f"k {k} is out of range [{self._instance.nof_loci.value}, "
                f"{self._instance.nof_loci.value * self._t_max_allele}]"
            )

        inst = self._instance.clone()
        getattr(inst, str(inst.beta)).set_value(float(0.0))

        inst.del_component("heterozygot_count")
        inst.del_component("reg1")
        inst.del_component("reg2")
        inst.del_component("reg3")

        expr1 = sum(inst.x[j] for j in inst.x)
        inst.c.add(expr1 == k)

        d = defaultdict(list)

        for _ in range(ks):
            try:
                res = self._solver.solve(inst, options=self._opts, tee=self._verbosity)
            except Exception:
                print(
                    "WARNING: Solver does not support multi-threading. "
                    "Falling back to single-threading."
                )
                res = self._solver.solve(inst, options={}, tee=self._verbosity)

            inst.solutions.load_from(res)

            if self._verbosity:
                res.write(num=1)

            if res.solver.termination_condition != TerminationCondition.optimal:
                print("Optimal solution hasn't been obtained.")
                break

            selected = []
            expr = 0
            indices = []
            encountered_4digit = []

            for j in inst.x:
                if 0.99 <= inst.x[j].value <= 1.01:
                    exp_i = 0
                    selected.append(j)
                    exp_i += inst.x[j]

                    if self._allele_to_4digit[j] in encountered_4digit:
                        continue

                    encountered_4digit.append(self._allele_to_4digit[j])
                    for i_allele in self._groups_4digit[self._allele_to_4digit[j]]:
                        if inst.x[i_allele].value <= 0:
                            exp_i += inst.x[i_allele]
                        indices.append(i_allele)
                    expr += 1 - exp_i

            zero_indices = set(j for j in inst.x).difference(set(indices))
            for j in zero_indices:
                expr += inst.x[j]

            inst.c.add(expr >= 1)

            if self._verbosity:
                print(selected)

            aas = [self._allele_to_4digit[x].split("*")[0] for x in selected]
            c = dict.fromkeys(aas, 1)

            for i in range(len(aas)):
                if aas.count(aas[i]) < 2:
                    d[aas[i] + "1"].append(selected[i])
                    d[aas[i] + "2"].append(selected[i])
                else:
                    d[aas[i] + str(c[aas[i]])].append(selected[i])
                    c[aas[i]] += 1

            nof_reads = sum(inst.occ[j] * inst.y[j].value for j in inst.y)
            d["obj"].append(inst.read_cov())
            d["nof_reads"].append(nof_reads)

        return pd.DataFrame(d)

    def solve_enforced_zygosity(self, gosity_dict: dict) -> pd.DataFrame:
        """
        Solve with enforced homo/heterozygosity per locus (experimental).

        Args:
            gosity_dict: Dictionary mapping locus to number of alleles.

        Returns:
            DataFrame with typing results.
        """
        inst = self._instance.clone()
        getattr(inst, str(inst.beta)).set_value(float(0.0))

        inst.del_component("heterozygot_count")
        inst.del_component("reg1")
        inst.del_component("reg2")
        inst.del_component("reg3")
        inst.del_component("max_allel_selection")

        for locus in inst.LociNames:
            cons = sum(inst.x[a] for a in inst.Loci[locus])
            inst.c.add(cons <= gosity_dict.get(locus, 2))

        d = defaultdict(list)

        try:
            res = self._solver.solve(inst, options=self._opts, tee=self._verbosity)
        except Exception:
            print(
                "WARNING: Solver does not support multi-threading. "
                "Falling back to single-threading."
            )
            res = self._solver.solve(inst, options={}, tee=self._verbosity)

        inst.solutions.load_from(res)

        selected = [al for al in inst.x if 0.99 <= inst.x[al].value <= 1.01]
        aas = [self._allele_to_4digit[x].split("*")[0] for x in selected]
        c = dict.fromkeys(aas, 1)

        for q in range(len(aas)):
            if aas.count(aas[q]) < 2:
                d[aas[q] + "1"].append(selected[q])
                d[aas[q] + "2"].append(selected[q])
            else:
                d[aas[q] + str(c[aas[q]])].append(selected[q])
                c[aas[q]] += 1

        nof_reads = sum(inst.occ[h] * inst.y[h].value for h in inst.y)
        d["obj"].append(inst.read_cov())
        d["nof_reads"].append(nof_reads)

        return pd.DataFrame(d)
