"""
Group/class FDR and local FDR (lFDR) for reference and non-canonical peptides.

Dr James Wright 2026 - The Institute of Cancer Research -

Decoy assignment
------------------------

Three strategies for decoy_assignment:

  "proportional"  (default)
      Randomly split decoys between groups in proportion to the target counts. 

  "score_kde"
      Assign each decoy to group based on PEP distribution using a KDE likelihood ratio. 
      Best for the lFDR calculation since the same KDE is used.

  "all_to_each"
      Give every group the all decoys. Over-estimates FDR making it the most conservative method.
      Requires zero assumptions about decoy composition.

Results written to PeptideInfo:

  group_fdr  (Dict[str, float])
    "qvalue"           — per-group target-decoy q-value  (conservative)
    "running_fdr"      — instantaneous running FDR at this peptide's rank
    "global_qvalue"    — pooled q-value, diagnostic only

  group_lfdr  (Dict[str, float])
    "lfdr"             — local FDR from per-group KDE density ratio
    "qvalue"           — q-value derived from per-group lFDR
    "global_lfdr"      — pooled KDE lFDR, diagnostic only
    "global_qvalue"    — q-value from pooled lFDR, diagnostic only

References
----------
Elias & Gygi (2007) Nat Methods 4:207-214.
Kall  et al. (2008) J Proteome Res 7:40-44.
Chong et al. (2020) Nat Commun 11:1293.       [NewAnce group FDR]
"""

from __future__ import annotations

import logging, sys
import warnings
import numpy as np
import argparse
from dataclasses import dataclass, field
from typing import Dict, Literal, Tuple

from Results_Annotator import PeptideInfo

logger = logging.getLogger(__name__)

try:
    from scipy.stats import gaussian_kde
    _SCIPY_AVAILABLE = True
except ImportError:
    _SCIPY_AVAILABLE = False

DecoyAssignment = Literal["proportional", "score_kde", "all_to_each"]

# ---------------------------------------------------------------------------
# Parse User Arguments
# ---------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description='''
        Compute group FDR and lFDR for reference and non-canonical peptides.
        - Requires numpy and optionally scipy for KDE-based methods.
        - Imports peptide dataclasses from Results_Annotator.py.
        - Input is a TSV peptide table with columns: peptide, type, and a score column (default=bestPSMScore).
        - Output is TSV table with global and group qvalues for each peptide.
        ---- James Wright (2026) - The Institute of Cancer Research, London ----'''
    )

    parser.add_argument(
        "--peptidetable", "-p", dest='input_file', type=str, required=True, 
        help="Input peptide table (TSV) with columns: peptide, type, psm_pep."
    )

    parser.add_argument(
        "--fdr_threshold", "-ft", type=float, dest='fdr_threshold', default=0.01,
        help="Q-value threshold for summary logging (default: 0.01)"
    )
    
    parser.add_argument(
        "--score_higher_better", '-sh', dest='score_higher_better', action="store_true",
        help="Whether higher scores are better (default: False, e.g. for PEP)"
    )

    parser.add_argument(
        "--decoy_ratio", '-dr', type=float, dest='decoy_ratio', default=1.0,
        help="Decoy-to-target ratio in the database (default: 1.0)"
    )

    parser.add_argument(
        "--score_field", "-f", dest='score_field', type=str, default="bestPSMScore",
        help="PeptideInfo attribute to use as score (default: 'bestPSMScore')"
    )

    parser.add_argument(
        "--decoy_assignment", '-da', dest='decoy_assignment', type=str, default="proportional",
        choices=["proportional", "score_kde", "all_to_each"],
        help=(
            "Strategy to assign untyped decoys to groups "
            "(default: 'proportional', choices: 'proportional', 'score_kde', 'all_to_each')"
        )
    )
    parser.add_argument(
        "--lfdr_bandwidth", "-lfb", dest='lfdr_band', type=str, default="scott",
        help=(
            "KDE bandwidth method for lFDR ('scott', 'silverman', or a float; "
            "default: 'scott')"
        )
    )
    parser.add_argument(
        "--min_decoys_for_lfdr", '-md', dest='min_decoys', type=int, default=50,
        help=(
            "Minimum number of decoys in a group to attempt KDE lFDR calculation "
            "(default: 50). Groups with fewer decoys will have NaN lFDR values."
        )
    )
    parser.add_argument(
        "--random_seed", '-rs', dest='seed', type=int, default=42,
        help="Random seed for reproducibility of proportional split (default: 42)"
    )
    parser.add_argument(
        "--output_file", '-o', dest='output_file', type=str, default="group_fdr_results.tsv",
        help="Output TSV file to write peptides with group FDR and lFDR results (default: 'group_fdr_results.tsv')"
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Function to compute FDR and qvalues from target-decoy scores
# ---------------------------------------------------------------------------

def td_fdr_and_qvalues(
    scores: np.ndarray,
    is_decoy: np.ndarray,
    *,
    decoy_ratio: float = 1.0,
    higher_better: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Target-decoy running FDR and monotone-minimum q-values (Käll et al. 2008).

    FDR(i) = (cum_decoys(i) + 1) / max(cum_targets(i), 1) / decoy_ratio
    q(i)   = min(FDR(j) for j >= i)
    """
    n = len(scores)
    if n == 0:
        return np.empty(0), np.empty(0)

    idx  = np.argsort(scores)[::-1] if higher_better else np.argsort(scores)
    sd   = is_decoy[idx].astype(bool)
    st   = ~sd
    cd   = np.cumsum(sd).astype(float)
    ct   = np.cumsum(st).astype(float)

    fdr_s  = np.clip((cd + 1.0) / (np.maximum(ct, 1.0) * decoy_ratio), 0.0, 1.0)
    qval_s = np.minimum.accumulate(fdr_s[::-1])[::-1]

    fdr_out  = np.empty(n)
    qval_out = np.empty(n)
    fdr_out[idx]  = fdr_s
    qval_out[idx] = qval_s
    return fdr_out, qval_out

# ---------------------------------------------------------------------------
# Function to compute KDE lFDR from target-decoy scores
# ---------------------------------------------------------------------------

def kde_lfdr(
    scores: np.ndarray,
    is_decoy: np.ndarray,
    *,
    decoy_ratio: float = 1.0,
    bandwidth = "scott",
    higher_better: bool = False,
) -> np.ndarray:
    """
    KDE density-ratio local FDR.

    lFDR(s) = [n_d * f_decoy(s)]  /
              [n_d * f_decoy(s) + (n_t * decoy_ratio) * f_target(s)]
    """
    t_scores = scores[~is_decoy]
    d_scores = scores[is_decoy]
    n_t, n_d = len(t_scores), len(d_scores)

    if n_t < 2 or n_d < 2:
        return np.full(len(scores), np.nan)
    try:
        f_t = gaussian_kde(t_scores, bw_method=bandwidth).evaluate(scores)
        f_d = gaussian_kde(d_scores, bw_method=bandwidth).evaluate(scores)
        w_t = float(n_t) * decoy_ratio
        w_d = float(n_d)
        denom = np.maximum(w_d * f_d + w_t * f_t, 1e-300)
        return np.clip(w_d * f_d / denom, 0.0, 1.0)
    except Exception as exc:
        logger.warning("KDE lFDR failed: %s", exc)
        return np.full(len(scores), np.nan)

# ---------------------------------------------------------------------------
# Function to compute q-values from lFDR values
# ---------------------------------------------------------------------------

def qvalues_from_lfdr(
    scores: np.ndarray,
    lfdr: np.ndarray,
    is_decoy: np.ndarray,
    *,
    higher_better: bool = False,
) -> np.ndarray:
    """
    Q-values from lFDR: running mean of lFDR over targets sorted best-first,
    with monotone-minimum correction.

    q(i) = min over j>=i of [mean(lFDR_target(rank 0 .. j))]
    """
    n    = len(scores)
    idx  = np.argsort(scores)[::-1] if higher_better else np.argsort(scores)
    sl   = lfdr[idx]
    st   = ~is_decoy[idx]
    cum  = np.cumsum(sl * st)
    cnt  = np.cumsum(st).astype(float)

    with np.errstate(invalid="ignore", divide="ignore"):
        qv_s = np.where(cnt > 0, cum / cnt, 1.0)

    qv_s = np.minimum.accumulate(qv_s[::-1])[::-1]
    out  = np.empty(n)
    out[idx] = qv_s
    return out

# ---------------------------------------------------------------------------
# Function to extract score and decoy arrays for a given group of sequences
# ---------------------------------------------------------------------------

def extract_arrays(
    peptides: Dict[str, PeptideInfo],
    seqs: list[str],
    score_field: str,
) -> Tuple[np.ndarray, np.ndarray]:
    scores   = np.fromiter((getattr(peptides[s], score_field) for s in seqs),
                           dtype=float, count=len(seqs))
    is_decoy = np.fromiter((peptides[s].decoy for s in seqs),
                           dtype=bool,  count=len(seqs))
    return scores, is_decoy

# ---------------------------------------------------------------------------
# Function to log group statistics at a given FDR threshold
# ---------------------------------------------------------------------------

def log_group_stats(
    peptides, groups, threshold, fdr_key, qval_key, label,
) -> None:
    for g, seqs in groups.items():
        if not seqs:
            continue
        n_t    = sum(1 for s in seqs if not peptides[s].decoy)
        n_d    = sum(1 for s in seqs if peptides[s].decoy)
        n_pass = sum(
            1 for s in seqs
            if not peptides[s].decoy
            and getattr(peptides[s], fdr_key).get(qval_key, 99.0) <= threshold
        )
        logger.info(
            "  [%-16s] %4d targets  %4d decoys  → %4d pass at %.1f%%  (%s)",
            g, n_t, n_d, n_pass, threshold * 100, label,
        )

# ---------------------------------------------------------------------------
# Function to log combined summary comparing group FDR and lFDR results at a given threshold
# ---------------------------------------------------------------------------

def log_combined_summary(peptides, groups, threshold) -> None:
    lines = [
        "",
        "=" * 72,
        f"  Group FDR Summary  |  threshold = {threshold * 100:.1f}%",
        f"  {'Group':<18}  {'Targets':>7}  {'Decoys':>7}  "
        f"{'q-val pass':>10}  {'lFDR pass':>9}  {'vs global':>9}",
        "=" * 72,
    ]
    for g, seqs in groups.items():
        targets  = [s for s in seqs if not peptides[s].decoy]
        decoys   = [s for s in seqs if peptides[s].decoy]
        pass_qv  = sum(1 for s in targets
                       if peptides[s].group_fdr.get("qvalue",        99) <= threshold)
        pass_lf  = sum(1 for s in targets
                       if peptides[s].group_lfdr.get("qvalue",       99) <= threshold)
        pass_glb = sum(1 for s in targets
                       if peptides[s].group_fdr.get("global_qvalue", 99) <= threshold)
        vs       = pass_qv - pass_glb
        lines.append(
            f"  {g:<18}  {len(targets):>7}  {len(decoys):>7}  "
            f"{pass_qv:>10}  {pass_lf:>9}  {vs:>+9}"
        )
    lines.append("=" * 72)
    logger.info("\n".join(lines))

# ---------------------------------------------------------------------------
# Function to assign untyped decoys to groups based on the specified strategy
# ---------------------------------------------------------------------------

def assign_untyped_decoys(
    peptides: Dict[str, PeptideInfo],
    typed_targets: Dict[str, list[str]],
    typed_decoys: Dict[str, list[str]],
    untyped_decoys: list[str],
    strategy: DecoyAssignment,
    score_field: str,
    decoy_ratio: float,
    higher_better: bool,
    bandwidth,
    rng: np.random.Generator,
) -> Tuple[Dict[str, list[str]], Dict[str, float]]:
    """
    Assign untyped decoys to groups and compute effective decoy ratios.

    Returns
    -------
    groups : Dict[str, list[str]]
        Each group's full sequence list (targets + all decoys for that group).
    effective_ratios : Dict[str, float]
        Effective decoy_ratio per group, adjusted for the assignment strategy.
    """
    group_names = list(typed_targets.keys())

    # Count typed targets per group (determines proportions)
    n_targets = {g: len(typed_targets[g]) for g in group_names}
    n_total_targets = max(sum(n_targets.values()), 1)
    fractions = {g: n_targets[g] / n_total_targets for g in group_names}

    # Default: each group starts with its typed targets + typed decoys
    groups: Dict[str, list[str]] = {
        g: typed_targets[g] + typed_decoys[g] for g in group_names
    }
    # Default effective ratio = provided decoy_ratio
    effective_ratios: Dict[str, float] = {g: decoy_ratio for g in group_names}

    # If there are no untyped decoys, we're done.
    if not untyped_decoys:
        return groups, effective_ratios

    if strategy == "proportional":
        assigned = _split_proportional(untyped_decoys, fractions, rng)
        for g in group_names:
            groups[g].extend(assigned[g])

    elif strategy == "score_kde":
        assigned = _split_by_score_kde(
            peptides, untyped_decoys, typed_targets,
            score_field, higher_better, bandwidth,
        )
        for g in group_names:
            groups[g].extend(assigned[g])

    elif strategy == "all_to_each":
        # Every group gets all untyped decoys.
        # Scale effective_ratio so FDR formula remains calibrated:
        #   We observe N_d decoys but expect f_g * N_d for group g.
        #   Divide by f_g to correct: effective_ratio(g) = decoy_ratio / f_g
        for g in group_names:
            groups[g].extend(untyped_decoys)
            f = fractions[g] if fractions[g] > 0 else 1.0
            effective_ratios[g] = decoy_ratio / f

    else:
        raise ValueError(
            f"Unknown decoy_assignment '{strategy}'. "
            "Choose: 'proportional', 'score_kde', 'all_to_each'."
        )

    return groups, effective_ratios

# ---------------------------------------------------------------------------
# Function to randomly split untyped decoys proportional to group target fractions
# ---------------------------------------------------------------------------

def _split_proportional(
    untyped_decoys: list[str],
    fractions: Dict[str, float],
    rng: np.random.Generator,
) -> Dict[str, list[str]]:
    """
    Randomly split untyped decoys proportional to group target fractions.

    The split is deterministic given random_seed and preserves the exact
    counts as closely as possible (integer rounding on all but the last group).
    """
    group_names = list(fractions.keys())
    shuffled    = list(untyped_decoys)
    rng.shuffle(shuffled)
    n = len(shuffled)

    assigned: Dict[str, list[str]] = {g: [] for g in group_names}
    cursor = 0
    for i, g in enumerate(group_names[:-1]):
        k = round(fractions[g] * n)
        assigned[g] = shuffled[cursor: cursor + k]
        cursor += k
    # Last group gets the remainder (avoids off-by-one from rounding)
    assigned[group_names[-1]] = shuffled[cursor:]

    for g, seqs in assigned.items():
        logger.info(
            "  Proportional split → '%s': %d untyped decoys assigned",
            g, len(seqs),
        )
    return assigned

# ---------------------------------------------------------------------------
# Function to assign untyped decoys to groups based on KDE likelihood ratio of their scores
# ---------------------------------------------------------------------------

def _split_by_score_kde(
    peptides: Dict[str, PeptideInfo],
    untyped_decoys: list[str],
    typed_targets: Dict[str, list[str]],
    score_field: str,
    higher_better: bool,
    bandwidth,
) -> Dict[str, list[str]]:
    """
    Assign each untyped decoy to the group whose target PEP distribution
    it is most consistent with (KDE log-likelihood ratio).

    For each decoy score s, computes:
        log p(s | group_g targets)   for each group g
    and assigns the decoy to the group with the highest log-likelihood.

    If a group has fewer than 2 targets (KDE not possible), falls back
    to proportional split based on target counts.
    """
    group_names  = list(typed_targets.keys())
    assigned: Dict[str, list[str]] = {g: [] for g in group_names}

    # Build KDE per group from target scores
    kdes: Dict[str, gaussian_kde | None] = {}
    for g in group_names:
        t_scores = np.array([
            getattr(peptides[s], score_field) for s in typed_targets[g]
        ])
        if len(t_scores) >= 2:
            try:
                kdes[g] = gaussian_kde(t_scores, bw_method=bandwidth)
            except Exception:
                kdes[g] = None
        else:
            kdes[g] = None

    if all(k is None for k in kdes.values()):
        logger.warning(
            "score_kde: all group KDEs failed — falling back to proportional."
        )
        fractions = {
            g: len(typed_targets[g]) / max(sum(len(v) for v in typed_targets.values()), 1)
            for g in group_names
        }
        rng = np.random.default_rng(0)
        return _split_proportional(untyped_decoys, fractions, rng)

    for seq in untyped_decoys:
        s = getattr(peptides[seq], score_field)

        log_likelihoods: Dict[str, float] = {}
        for g in group_names:
            if kdes[g] is not None:
                ll = float(np.log(np.maximum(kdes[g].evaluate([s])[0], 1e-300)))
            else:
                # No KDE: assign neutral log-likelihood so this group is
                # only chosen if the others are also neutral or worse
                ll = float("-inf")
            log_likelihoods[g] = ll

        best_group = max(log_likelihoods, key=log_likelihoods.__getitem__)
        assigned[best_group].append(seq)

    for g, seqs in assigned.items():
        logger.info(
            "  score_kde split → '%s': %d untyped decoys assigned",
            g, len(seqs),
        )
    return assigned

# ---------------------------------------------------------------------------
# Function to filter peptides based on a specified FDR method and threshold
# ---------------------------------------------------------------------------

def filter_peptides(
    peptides: Dict[str, PeptideInfo],
    fdr_threshold: float = 0.01,
    method: str = "group_qvalue",
    targets_only: bool = True,
) -> Dict[str, PeptideInfo]:
    """
    Return peptides passing the FDR threshold.

    Parameters
    ----------
    method : str
        "group_qvalue"        — per-group target-decoy q-value (conservative)
        "lfdr_qvalue"         — per-group lFDR-derived q-value
        "global_qvalue"       — pooled target-decoy q-value (diagnostic only)
        "global_lfdr_qvalue"  — pooled lFDR q-value (diagnostic only)
    """
    _KEY_MAP = {
        "group_qvalue":       ("group_fdr",  "qvalue"),
        "lfdr_qvalue":        ("group_lfdr", "qvalue"),
        "global_qvalue":      ("group_fdr",  "global_qvalue"),
        "global_lfdr_qvalue": ("group_lfdr", "global_qvalue"),
    }
    if method not in _KEY_MAP:
        raise ValueError(f"method must be one of {list(_KEY_MAP)}")
    attr, key = _KEY_MAP[method]

    return {
        seq: info
        for seq, info in peptides.items()
        if getattr(info, attr, {}).get(key, float("inf")) <= fdr_threshold
        and (not targets_only or not info.decoy)
    }

# ---------------------------------------------------------------------------
# Function to summarize group FDR results by type, counting targets passing each FDR method at the threshold
# ---------------------------------------------------------------------------

def group_fdr_summary(
    peptides: Dict[str, PeptideInfo],
    fdr_threshold: float = 0.01,
) -> Dict[str, Dict]:
    """Structured per-type summary with counts for both FDR approaches."""
    summary: Dict[str, Dict] = {}
    for info in peptides.values():
        if not info.group_fdr:
            continue
        t = info.type
        if t not in summary:
            summary[t] = {
                "n_targets": 0, "n_decoys": 0,
                "n_pass_group_qvalue": 0,
                "n_pass_lfdr_qvalue": 0,
                "n_pass_global_qvalue": 0,
                "_sc_qv": [], "_sc_lf": [],
            }
        b = summary[t]
        if info.decoy:
            b["n_decoys"] += 1
            continue
        b["n_targets"] += 1
        if info.group_fdr.get("qvalue", 99.0)  <= fdr_threshold:
            b["n_pass_group_qvalue"] += 1
            b["_sc_qv"].append(info.psm_score)
        if info.group_lfdr.get("qvalue", 99.0) <= fdr_threshold:
            b["n_pass_lfdr_qvalue"] += 1
            b["_sc_lf"].append(info.psm_score)
        if info.group_fdr.get("global_qvalue", 99.0) <= fdr_threshold:
            b["n_pass_global_qvalue"] += 1

    for stats in summary.values():
        # PEP: higher cutoff = more peptides pass (score_higher_better=False)
        stats["score_cutoff_group_qvalue"] = (
            float(max(stats["_sc_qv"])) if stats["_sc_qv"] else float("nan")
        )
        stats["score_cutoff_lfdr_qvalue"] = (
            float(max(stats["_sc_lf"])) if stats["_sc_lf"] else float("nan")
        )
        del stats["_sc_qv"], stats["_sc_lf"]

    return summary



# ---------------------------------------------------------------------------
# Function to compute group FDR and lFDR, main entry point
# ---------------------------------------------------------------------------

def compute_group_fdr(
    peptides: Dict[str, PeptideInfo],
    *,
    fdr_threshold: float = 0.01,
    use_top_rank: bool = True,
    score_higher_better: bool = False,      # PEP: lower = better
    decoy_ratio: float = 1.0,
    score_field: str = "psm_score",
    decoy_assignment: DecoyAssignment = "proportional",
    lfdr_bandwidth: str | float = "scott",
    min_decoys_for_lfdr: int = 50,
    random_seed: int = 42,
) -> Dict[str, PeptideInfo]:
    
    """
    Compute per-group q-values and per-group lFDR.
    Untyped decoys are distributed between groups according to decoy_assignment.
    Updates group_fdr and group_lfdr in PeptideInfo Dataclass. 

    Parameters
    ----------
    peptides : Dict[str, PeptideInfo]
        Peptide dictionary keyed by sequence. 

    fdr_threshold : float
        Q-value threshold used for summary logging only. Default 1%.

    use_top_rank : bool
        If True (default), only peptides with PeptideInfo.top_rank=True are
        used for FDR calculations. Avoids issue where conflicting peptides or lower ranked peptides are included in the input. 

    score_higher_better : bool
        False (default) — lower score is better (PEP, p-value, E-value).
        True            — higher score is better.

    decoy_ratio : float
        Decoy-to-target DB size ratio (1.0 for a 1:1 combined FASTA).
        Only relevant when decoy_assignment is "all_to_each", where it is
        additionally scaled by the group target fraction.

    score_field : str
        PeptideInfo attribute to use as ranking score. Default: "psm_score".

    decoy_assignment : {"proportional", "score_kde", "all_to_each"}
        Strategy for distributing untyped decoys into groups:

        "proportional"  — Random split proportional to group target sizes. (Default)

        "score_kde"     — KDE likelihood ratio assigns each decoy to the
                          group whose target score distribution it best fits.

        "all_to_each"   — All untyped decoys go to every group. Conservative
                          (over-estimates FDR). 

    lfdr_bandwidth : str or float
        KDE bandwidth method. "scott" (default), "silverman", or a float.

    min_decoys_for_lfdr : int
        Minimum decoys per group to attempt KDE lFDR. Default 50.

    random_seed : int
        Seed for the proportional random split. Default 42.
    """

    if decoy_assignment == "score_kde" and not _SCIPY_AVAILABLE:
        warnings.warn(
            "decoy_assignment='score_kde' requires scipy. "
            "Falling back to 'proportional'.",
            ImportWarning, stacklevel=2,
        )
        decoy_assignment = "proportional"

    if not _SCIPY_AVAILABLE:
        warnings.warn(
            "scipy not found — lFDR (KDE) will be skipped. "
            "Install with: pip install scipy",
            ImportWarning, stacklevel=2,
        )

    rng = np.random.default_rng(random_seed)

    # ------------------------------------------------------------------
    # 1. Partition into typed targets, typed decoys, and untyped decoys
    # ------------------------------------------------------------------
    typed_targets:  Dict[str, list[str]] = {}
    typed_decoys:   Dict[str, list[str]] = {}
    untyped_decoys: list[str] = []

    for seq, info in peptides.items():

        # Skip non-top-rank peptides if requested. 
        # This is relevant for multi-database/engine searches 
        if use_top_rank and not info.top_rank:
            continue

        if info.decoy:
            
            if info.type == "unknown":
                untyped_decoys.append(seq)   # ← the common case
            else:
                if info.type not in typed_decoys:
                    typed_decoys[info.type] = []

                typed_decoys[info.type].append(seq)
               
        else:
            if info.type not in typed_targets:
                typed_targets[info.type] = []
            
            if info.type not in typed_decoys:
                typed_decoys[info.type] = []

            typed_targets[info.type].append(seq)    

    n_typed_d   = sum(len(v) for v in typed_decoys.values())
    n_untyped_d = len(untyped_decoys)

    #Report on inventory and assignment strategy
    logger.info(
        "Decoy inventory: %d typed, %d untyped → assignment='%s'",
        n_typed_d,
        n_untyped_d,
        decoy_assignment,
    )

    for g, seqs in typed_decoys.items():
        logger.info("  Typed decoys for group '%s': %d", g, len(seqs))

    logger.info(
        "Target inventory: %d typed → groups: %s",
        sum(len(v) for v in typed_targets.values()),
        ", ".join(f"'{g}': {len(s)}" for g, s in typed_targets.items()),
    )

    for g in typed_targets:
        logger.info("  Typed targets for group '%s': %d", g, len(typed_targets[g]))

    # ------------------------------------------------------------------
    # 2. Assign untyped decoys to groups
    # ------------------------------------------------------------------
    # groups[g] = targets + typed decoys + assigned untyped decoys for group g
    # effective_ratio[g] = adjusted decoy_ratio for the "all_to_each" case
    groups, effective_ratios = assign_untyped_decoys(
        peptides=peptides,
        typed_targets=typed_targets,
        typed_decoys=typed_decoys,
        untyped_decoys=untyped_decoys,
        strategy=decoy_assignment,
        score_field=score_field,
        decoy_ratio=decoy_ratio,
        higher_better=score_higher_better,
        bandwidth=lfdr_bandwidth,
        rng=rng,
    )

    # Log group composition
    for g, seqs in groups.items():
        n_t = sum(1 for s in seqs if not peptides[s].decoy)
        n_d = sum(1 for s in seqs if peptides[s].decoy)
        logger.info("  Group '%-16s': %d targets, %d decoys (ratio=%.2f)",
                    g, n_t, n_d, effective_ratios[g])

    # ------------------------------------------------------------------
    # 3. Per-group target-decoy q-values - group_fdr!
    # ------------------------------------------------------------------
    all_seqs = []
    for group_name, seqs in groups.items():
        if not seqs:
            logger.warning("Group '%s' is empty — skipping.", group_name)
            continue
        all_seqs.extend(seqs)

        scores, is_decoy = extract_arrays(peptides, seqs, score_field)
        run_fdr, qvals   = td_fdr_and_qvalues(
            scores, is_decoy,
            decoy_ratio=effective_ratios[group_name],
            higher_better=score_higher_better,
        )
        for seq, rf, qv in zip(seqs, run_fdr, qvals):
            peptides[seq].group_fdr["running_fdr"] = float(rf)
            peptides[seq].group_fdr["qvalue"]      = float(qv)

    # Deduplicate all_seqs (relevant for "all_to_each" where decoys appear twice)
    seen: set[str] = set()
    unique_all = [s for s in all_seqs if not (s in seen or seen.add(s))]  # type: ignore[func-returns-value]

    # Global pooled target-decoy q-value — diagnostic
    if unique_all:
        scores, is_decoy = extract_arrays(peptides, unique_all, score_field)
        _, gqvals = td_fdr_and_qvalues(
            scores, is_decoy,
            decoy_ratio=decoy_ratio,
            higher_better=score_higher_better,
        )
        for seq, qv in zip(unique_all, gqvals):
            peptides[seq].group_fdr["global_qvalue"] = float(qv)

    log_group_stats(peptides, groups, fdr_threshold,
                     fdr_key="group_fdr", qval_key="qvalue",
                     label="target-decoy q-value")

    # ------------------------------------------------------------------
    # 4. Per-group KDE lFDR  →  group_lfdr
    # ------------------------------------------------------------------
    if _SCIPY_AVAILABLE:
        for group_name, seqs in groups.items():
            if not seqs:
                continue
            scores, is_decoy = extract_arrays(peptides, seqs, score_field)
            n_d = int(is_decoy.sum())

            if n_d < min_decoys_for_lfdr:
                logger.warning(
                    "Group '%s': %d decoys < min_decoys_for_lfdr=%d — "
                    "lFDR set to NaN.", group_name, n_d, min_decoys_for_lfdr,
                )
                for seq in seqs:
                    peptides[seq].group_lfdr.setdefault("lfdr",   float("nan"))
                    peptides[seq].group_lfdr.setdefault("qvalue", float("nan"))
                continue

            lfdr_vals  = kde_lfdr(scores, is_decoy,
                                   decoy_ratio=effective_ratios[group_name],
                                   bandwidth=lfdr_bandwidth,
                                   higher_better=score_higher_better)
            lfdr_qvals = qvalues_from_lfdr(scores, lfdr_vals, is_decoy,
                                            higher_better=score_higher_better)
            for seq, lf, lq in zip(seqs, lfdr_vals, lfdr_qvals):
                peptides[seq].group_lfdr["lfdr"]   = float(lf)
                peptides[seq].group_lfdr["qvalue"] = float(lq)

        # Global pooled lFDR — diagnostic
        if unique_all:
            scores, is_decoy = extract_arrays(peptides, unique_all, score_field)
            g_lfdr  = kde_lfdr(scores, is_decoy, decoy_ratio=decoy_ratio,
                                bandwidth=lfdr_bandwidth,
                                higher_better=score_higher_better)
            g_lqval = qvalues_from_lfdr(scores, g_lfdr, is_decoy,
                                         higher_better=score_higher_better)
            for seq, lf, lq in zip(unique_all, g_lfdr, g_lqval):
                peptides[seq].group_lfdr["global_lfdr"]   = float(lf)
                peptides[seq].group_lfdr["global_qvalue"] = float(lq)

        log_group_stats(peptides, groups, fdr_threshold,
                         fdr_key="group_lfdr", qval_key="qvalue",
                         label="lFDR q-value")
        
    else:
        logger.warning("scipy not available — skipping lFDR computation.")

    log_combined_summary(peptides, groups, fdr_threshold)
    return peptides


# ---------------------------------------------------------------------------
# Peptide Table Parsing Function 
# ---------------------------------------------------------------------------
def parse_peptide_table(file_path: str, score_field: str) -> Dict[str, PeptideInfo]:

    logger.info("Parsing peptide table from '%s' with score field '%s'", file_path, score_field)

    peptides: Dict[str, PeptideInfo] = {}
    with open(file_path, 'r') as f:
        header = f.readline().strip().split('\t')
        col_idx = {col: idx for idx, col in enumerate(header)}

        for line in f:
            fields = line.strip().split('\t')
            seq = None
            if 'PeptideSequence' in col_idx:
                seq = fields[col_idx['PeptideSequence']]
            elif 'Sequence' in col_idx:
                seq = fields[col_idx['Sequence']]
            elif 'Peptide' in col_idx:
                seq = fields[col_idx['Peptide']]
            else:
                raise ValueError(
                    "Input file must contain a 'PeptideSequence', 'Sequence', or 'Peptide' column."
                )

            decoy = False
            type = None
            if 'Type' in col_idx:
                type = fields[col_idx['Type']] 
            elif 'Peptide_Type' in col_idx:
                type = fields[col_idx['Peptide_Type']]
            else:
                type = "unknown"
                logger.warning(
                    "No 'Type' or 'Peptide_Type' column found. All peptides will be assigned type='unknown'."
                )

            if type == "DECOY":
                    decoy=True
                    type = "unknown"

            if 'Decoy' in col_idx:
                if fields[col_idx['Decoy']].lower() == 'true':
                    decoy = True
                    type = "unknown"

            if score_field not in col_idx:
                logger.warning(
                    f"Score field '{score_field}' not found in input columns. "
                )
                raise ValueError(f"Score field '{score_field}' not found in input columns.")
            
            top_rank = True
            if 'top_rank' in col_idx:
                top_rank = fields[col_idx['top_rank']].lower() == 'true'

            peptides[seq] = PeptideInfo(
                type=type,
                decoy=decoy,
                psm_score=float(fields[col_idx[score_field]]),
                top_rank=top_rank,
            )

    logger.info("Parsed %d peptides from input file.", len(peptides))

    return peptides

# ---------------------------------------------------------------------------
# Peptide Table Writing Function
# ---------------------------------------------------------------------------
def write_peptide_table(peptides: Dict[str, PeptideInfo], output_file: str) -> None:

    logger.info("Writing peptide table with group FDR results to '%s'", output_file)

    with open(output_file, 'w') as f:
        header = [
            "PeptideSequence", "Type", "Decoy", "psm_score", "top_rank",
            "group_fdr_qvalue", "group_lfdr_qvalue", "global_qvalue", "global_lfdr_qvalue"
        ]
        f.write('\t'.join(header) + '\n')

        for seq, info in peptides.items():
            row = [
                seq,
                info.type,
                str(info.decoy),
                str(info.psm_score),
                str(info.top_rank),
                str(info.group_fdr.get("qvalue", "na")),
                str(info.group_lfdr.get("qvalue", "na")),
                str(info.group_fdr.get("global_qvalue", "na")),
                str(info.group_lfdr.get("global_lfdr", "na")),
            ]
            f.write('\t'.join(row) + '\n')


# ---------------------------------------------------------------------------
# Main block to apply the group FDR computation to an existing peptide table
# ---------------------------------------------------------------------------
def main() -> int:

    logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                        format="[%(levelname)s] %(message)s")

    args = parse_args()
    peptides = parse_peptide_table(args.input_file, args.score_field)  

    compute_group_fdr(
        peptides = peptides,
        fdr_threshold = args.fdr_threshold,
        score_higher_better = args.score_higher_better,
        use_top_rank=True,
        decoy_ratio = args.decoy_ratio,
        decoy_assignment = args.decoy_assignment,
        lfdr_bandwidth = args.lfdr_band,
        min_decoys_for_lfdr = args.min_decoys,
        random_seed = args.seed,
    ) 

    write_peptide_table(peptides, args.output_file)

    logger.info("Group FDR computation completed successfully.")


if __name__ == "__main__":
    raise SystemExit(main())