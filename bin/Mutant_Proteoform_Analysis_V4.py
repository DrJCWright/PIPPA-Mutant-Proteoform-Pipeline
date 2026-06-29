#!/usr/bin/env python3
"""
Mutant Proteoform Analysis
==========================

Estimate mutant-proteoform abundance from TMT proteomics data using
"wild-type peptide dropout" at the sites of mapped somatic mutations.

For every channel (cell line) and every mutation that the cell line carries,
the script locates the *wild-type* tryptic peptide(s) that report on the
mutated residue (either covering the residue, or immediately flanking a
cleavage-site mutation that the variant destroys/creates). It then measures how
far the TMT signal of that WT peptide has dropped relative to the rest of the
gene's peptides. A drop indicates that part of the protein pool in that cell
line is the mutant proteoform (which no longer yields the WT peptide), so the
dropout is a proxy for mutant-proteoform fraction / bias.

Two complementary references are computed for each WT-at-site peptide:

  1. Bridge/control reference  -- ratio of the WT peptide to the same peptide in
     the designated control/bridge channel, compared against the gene-level
     control ratio. This is the classic single-reference dropout.

  2. Cross-sample reference     -- ratio of the WT peptide in the mutant channel
     to a robust summary (trimmed mean / median) of the *same peptide across all
     the other channels in the plex*. This needs no dedicated bridge channel and
     is more robust when the control itself is atypical.

Inputs (one experiment / TMT plex per run):
  --samplefile           Sample table (ExperimentID, TMTLabel, SampleID)
  --pepfile              Peptide quant TSV (one experiment)
  --mutfile              Mapped mutation TSV
  --fasta                GENCODE protein-translation FASTA
  --spectrafile          PSM spectral TSV (one experiment)
  --control-channel-name SampleID of the bridge/control cell line (e.g. SW48)

Dr James C Wright (2025), The Institute of Cancer Research, London.
Revised 2026.
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import random
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from statistics import median as _median
from statistics import mean as _mean
from typing import Dict, Iterable, List, Optional, Set, Tuple


# =============================================================================
# Constants
# =============================================================================

# Canonical TMT-11 plex order. The trailing integer on the abundance columns
# (TMT_Abundance_Raw_<N>) and on the spectral intensity columns
# (Raw_TMT_Intensity_<N>) is this 1-based plex position.
TMT_PLEX_ORDER: Dict[int, str] = {
    1: "126", 2: "127N", 3: "127C", 4: "128N", 5: "128C",
    6: "129N", 7: "129C", 8: "130N", 9: "130C", 10: "131", 11: "131C",
}

# Peptide-level quality thresholds (FDR / PEP / delta).
MAX_PEPTIDE_FDR = 0.01
MAX_BEST_PSM_PEP = 0.01
MIN_DELTA = 0.02
# PSM-level PEP threshold for spectral-table aggregation.
MAX_PSM_PEP = 0.01
# A gene needs more than this many quantified peptides in a channel for a
# meaningful gene-level reference.
MIN_GENE_PEPTIDES = 3
# Minimum sub-peptide length considered for usability / matching.
MIN_SUBPEP_LEN = 7
# Fraction trimmed from each tail when building the cross-sample reference.
CROSS_SAMPLE_TRIM = 0.1
# Minimum number of PSMs required to attempt a bootstrap CI.
MIN_BOOTSTRAP_PSMS = 3
# Cross-sample filtering: a reference channel must carry at least this fraction
# of the gene's median across-channel expression to qualify (guards against
# channels where the gene is essentially absent).
MIN_GENE_EXPR_FRACTION = 0.1
# Cross-sample filtering: the WT peptide in a reference channel must reach at
# least this fraction of the across-channel median peptide signal to qualify
# (guards against near-noise-floor peptide values).
MIN_PEP_SIGNAL_FRACTION = 0.05
# Minimum qualifying reference channels for a cross-sample estimate to be
# reported as solid (fewer than this -> value still emitted but flagged).
MIN_CROSS_SAMPLE_CHANNELS = 3


LOGGER = logging.getLogger("proteoform")


# =============================================================================
# Argument parsing & logging
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate mutant-proteoform abundance via wild-type peptide "
            "dropout at mutation sites from TMT proteomics data. "
            "Dr J C Wright - The Institute of Cancer Research."
        )
    )
    parser.add_argument("--samplefile", "-e", required=True,
                        help="Sample table TSV (ExperimentID, TMTLabel, SampleID).")
    parser.add_argument("--pepfile", "-i", required=True,
                        help="Peptide quant TSV file (single experiment).")
    parser.add_argument("--mutfile", "-m", required=True,
                        help="Mapped mutation TSV file.")
    parser.add_argument("--fasta", "-f", required=True,
                        help="GENCODE protein-translation FASTA file.")
    parser.add_argument("--spectrafile", "-s", required=True,
                        help="PSM spectral TSV table (single experiment).")
    parser.add_argument("--control-channel-name", "-c", required=True,
                        help="SampleID of the control/bridge cell line, e.g. SW48.")
    parser.add_argument("--output", "-o", default="mutant_proteoform_results.txt",
                        help="Output TSV file.")
    parser.add_argument("--bootstrap", type=int, default=1000,
                        help="Bootstrap iterations for the cross-sample CI "
                             "(0 disables; default 1000).")
    parser.add_argument("--seed", type=int, default=1234,
                        help="Random seed for the bootstrap (default 1234).")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable debug logging.")
    return parser.parse_args()


def setup_logging(verbose: bool = False) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# =============================================================================
# Small numeric / string helpers
# =============================================================================

def safe_float(value: object, default: float = 0.0) -> float:
    try:
        if value is None:
            return default
        s = str(value).strip()
        if s == "" or s == "-":
            return default
        return float(s)
    except (TypeError, ValueError):
        return default


def safe_int(value: object, default: int = 0) -> int:
    try:
        if value is None:
            return default
        s = str(value).strip()
        if s == "" or s == "-":
            return default
        return int(float(s))
    except (TypeError, ValueError):
        return default


def mean(values: Iterable[float]) -> float:
    vals = list(values)
    return _mean(vals) if vals else 0.0


def median(values: Iterable[float]) -> float:
    vals = list(values)
    return float(_median(vals)) if vals else 0.0


def variance(values: Iterable[float]) -> float:
    vals = list(values)
    if len(vals) < 2:
        return 0.0
    m = _mean(vals)
    return sum((x - m) ** 2 for x in vals) / len(vals)


def trimmed_mean(values: Iterable[float], trim: float = CROSS_SAMPLE_TRIM) -> float:
    """Symmetric trimmed mean; falls back to the plain mean for short lists."""
    vals = sorted(values)
    n = len(vals)
    if n == 0:
        return 0.0
    k = int(math.floor(n * trim))
    if 2 * k >= n:
        return float(_median(vals))
    return _mean(vals[k:n - k])


def weighted_median(pairs: List[Tuple[float, float]]) -> Optional[float]:
    """Weighted median of (value, weight) pairs. Weights <= 0 are ignored.

    Returns the value at which the cumulative weight first reaches half the
    total weight; averages the two straddling values on an exact tie.
    """
    valid = [(v, w) for v, w in pairs if w > 0]
    if not valid:
        return None
    valid.sort(key=lambda x: x[0])
    total = sum(w for _, w in valid)
    half = total / 2.0
    cum = 0.0
    for i, (v, w) in enumerate(valid):
        cum += w
        if cum > half:
            return v
        if cum == half:
            # exact tie: average with the next distinct value if present
            if i + 1 < len(valid):
                return (v + valid[i + 1][0]) / 2.0
            return v
    return valid[-1][0]


def log2(value: float) -> Optional[float]:
    try:
        return math.log2(value) if value > 0 else None
    except (TypeError, ValueError):
        return None


def fmt(value: object) -> str:
    """Format a cell for TSV output, using '-' for None / missing."""
    if value is None:
        return "-"
    if isinstance(value, float):
        if math.isnan(value):
            return "-"
        return f"{value:.6g}"
    return str(value)


def split_tryptic_peptide(seq: str) -> List[str]:
    """Split after K or R unless followed by P (fully-cleaved fragments)."""
    return re.sub(r"(K|R)([^P])", r"\1_\2", seq).split("_")


def replace_i_with_l(seq: str) -> str:
    """I/L are isobaric and indistinguishable by MS; collapse to L."""
    return seq.replace("I", "L")


def normalise_name(name: str) -> str:
    """Normalise a cell-line / sample label for robust matching.

    Removes hyphens, spaces and punctuation and lower-cases, so that the
    sample-table 'CL-40' and the mutation-table 'CL40' compare equal.
    """
    return re.sub(r"[^a-z0-9]+", "", str(name).lower())


# =============================================================================
# Data classes
# =============================================================================

@dataclass
class MutationRecord:
    """A mutation aggregated across rows (one COSMIC id within a gene)."""
    gene_id: str
    cosmic_id: str
    mutation_type: str = "-"
    transcript_positions: Dict[str, int] = field(default_factory=dict)  # ENST -> AA pos
    unique_ids: Set[str] = field(default_factory=set)
    aa_positions: Set[str] = field(default_factory=set)
    aa_changes: Set[str] = field(default_factory=set)
    na_changes: Set[str] = field(default_factory=set)
    tissues_by_cellline: Dict[str, str] = field(default_factory=dict)
    zygosity_by_cellline: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))
    reference_peptides: Set[str] = field(default_factory=set)
    mutant_peptides: Set[str] = field(default_factory=set)


@dataclass
class ChannelQuant:
    """Per-channel quant for one peptide."""
    typex: str = "reference"          # reference | nonmutant (WT-at-site) | mutant
    raw: float = 0.0
    norm: float = 0.0
    scaled: float = 0.0
    quant_count: int = 0
    noq: float = 0.0                  # depth-normalised abundance
    ratio_to_control: Optional[float] = None  # None == bridge-zero ("B0")


@dataclass
class PSMStats:
    base_peak_intensity_sum: float = 0.0
    precursor_intensity_sum: float = 0.0
    max_control_ratios: List[float] = field(default_factory=list)
    per_channel_intensities: Dict[int, List[float]] = field(default_factory=lambda: defaultdict(list))
    # Aligned per-PSM channel vectors {channel_id: intensity}, one dict per PSM.
    # Keeps the channel pairing within a PSM intact, which per_channel_intensities
    # loses, so that bootstrap resampling can resample whole PSMs.
    psm_vectors: List[Dict[int, float]] = field(default_factory=list)


@dataclass
class PeptideRecord:
    peptide: str
    peptide_type: str
    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_count: int          # number of transcript IDs in `transcript_id`
    protein_count: str             # ProteinCount column (search-level grouping count)
    peptide_fdr: float
    best_psm_fdr: float
    best_psm_pep: float
    delta: float
    all_psm_count: float
    parent_sequences: Set[str] = field(default_factory=set)
    channels: Dict[int, ChannelQuant] = field(default_factory=dict)
    psm_stats: PSMStats = field(default_factory=PSMStats)
    mean_control_noq: float = 0.0     # mean of this peptide's control NOQ across experiments (single exp here)

    def channel(self, cid: int) -> ChannelQuant:
        if cid not in self.channels:
            self.channels[cid] = ChannelQuant()
        return self.channels[cid]


@dataclass
class GeneChannelRecord:
    raw_sum: float = 0.0
    peptide_count: int = 0
    norm_sum: float = 0.0
    ratio_to_control: Optional[float] = None
    scaled: Optional[float] = None
    median_ratio: Optional[float] = None
    weighted_median_ratio: Optional[float] = None  # PSM-count-weighted median of ref-peptide ratios
    topn_norm_sum: float = 0.0
    topn_ratio: Optional[float] = None
    mean_control_norm: float = 0.0


@dataclass
class PeptideUsability:
    """Usability info for a fully-cleaved sub-peptide (first-pass)."""
    ambiguity: int
    missed_cleavages: int
    passed: int
    peptide_type: str
    parent_by_missed: Dict[int, Dict[str, str]] = field(default_factory=lambda: defaultdict(dict))


@dataclass
class ChannelHeaders:
    """Resolved peptide-quant column names for one channel."""
    channel_id: int
    raw: str
    norm: str
    scaled: str


@dataclass
class Experiment:
    experiment_id: str
    channel_label: Dict[int, str] = field(default_factory=dict)        # channel -> TMT label
    channel_sample: Dict[int, str] = field(default_factory=dict)       # channel -> SampleID (cell line)
    control_channel_id: Optional[int] = None
    total_intensity_by_channel: Dict[int, float] = field(default_factory=lambda: defaultdict(float))


# =============================================================================
# TSV reading
# =============================================================================

def read_tsv(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
    """Read a TSV (tolerating CRLF) into a header list and list of row dicts."""
    with open(path, "r", newline="") as handle:
        reader = csv.reader((line.replace("\r\n", "\n").replace("\r", "\n") for line in handle),
                            delimiter="\t")
        headers = next(reader, None)
        if headers is None:
            raise ValueError(f"No header found in TSV: {path}")
        headers = [h.strip() for h in headers]
        rows = [dict(zip(headers, [c for c in row])) for row in reader if row]
    return headers, rows


def require_columns(headers: List[str], required: Iterable[str], context: str) -> None:
    missing = [c for c in required if c not in headers]
    if missing:
        raise KeyError(f"{context}: missing required column(s): {missing}\nFound: {headers}")


# =============================================================================
# Sample table
# =============================================================================

def load_sample_table(path: str) -> Dict[str, Dict[str, str]]:
    """Return {ExperimentID: {TMTLabel: SampleID}}."""
    LOGGER.info("Reading sample table: %s", path)
    headers, rows = read_tsv(path)
    require_columns(headers, ["ExperimentID", "TMTLabel", "SampleID"], "Sample table")
    table: Dict[str, Dict[str, str]] = defaultdict(dict)
    for row in rows:
        table[row["ExperimentID"]][row["TMTLabel"].strip()] = row["SampleID"].strip()
    LOGGER.info("Loaded sample table for %d experiments", len(table))
    return table


# =============================================================================
# FASTA
# =============================================================================

def parse_fasta(path: str) -> Dict[str, Dict[str, str]]:
    """Parse GENCODE translations into {ENSG: {ENST: sequence}} (I->L applied)."""
    LOGGER.info("Reading FASTA: %s", path)
    gene_sequences: Dict[str, Dict[str, str]] = defaultdict(dict)
    header_re = re.compile(r"(ENST\d+)\.\d+.*?(ENSG\d+)\.")
    enst = ensg = None
    buf: List[str] = []

    def flush() -> None:
        if enst and ensg:
            gene_sequences[ensg][enst] = replace_i_with_l("".join(buf))

    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                flush()
                m = header_re.search(line)
                if m:
                    enst, ensg = m.group(1), m.group(2)
                    buf = []
                else:
                    enst = ensg = None
                    buf = []
            elif enst and ensg:
                buf.append(line)
        flush()

    LOGGER.info("Loaded sequences for %d genes", len(gene_sequences))
    return gene_sequences


def peptide_positions_in_gene(
    gene_id: str, peptide: str, gene_sequences: Dict[str, Dict[str, str]],
) -> List[Tuple[str, int, int]]:
    """Return (transcript, start, end) 1-based spans where peptide occurs."""
    out: List[Tuple[str, int, int]] = []
    pep = replace_i_with_l(peptide)
    for enst, seq in gene_sequences.get(gene_id, {}).items():
        idx = seq.find(pep)
        if idx >= 0:
            start = idx + 1
            out.append((enst, start, start + len(pep) - 1))
    return out


# =============================================================================
# Mutation table
# =============================================================================

MUT_COLUMNS = [
    "UNIQUE_MUTATION_ID", "GENE_ID", "TRANSCRIPT_ID", "GENE_NAME", "SAMPLE_NAME",
    "MUTATION_TYPE", "SAMPLE_TISSUE", "MUTATION_ZYGOSITY", "COSMIC_ID",
    "AA_POSITION", "AA_CHANGE", "NA_CHANGE",
    "OBSERVABLE_WT_PEPTIDES", "OBSERVABLE_MUTANT_PEPTIDES",
]


def parse_mutation_file(
    mutfile: str,
    gene_sequences: Dict[str, Dict[str, str]],
    control_name: str,
) -> Tuple[Dict[str, Dict[str, MutationRecord]], Dict[str, Dict[str, Dict[str, Dict[str, str]]]]]:
    """Parse the mapped mutation table.

    Returns
    -------
    mutations_by_gene : {ENSG: {COSMIC_ID: MutationRecord}}
    peptide_mutation_map : {peptide: {cellline_norm: {ENSG: {COSMIC_ID: 'R'|'M'}}}}
        'R' = wild-type (reference) peptide reporting on the site,
        'M' = mutant peptide. Also stores under a 'control' pseudo-cellline.
    """
    LOGGER.info("Reading mutation table: %s", mutfile)
    headers, rows = read_tsv(mutfile)
    require_columns(headers, MUT_COLUMNS, "Mutation table")

    mutations_by_gene: Dict[str, Dict[str, MutationRecord]] = defaultdict(dict)
    peptide_mutation_map: Dict[str, Dict[str, Dict[str, Dict[str, str]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(dict))
    )
    control_norm = normalise_name(control_name)

    for row in rows:
        cosmic_id = row["COSMIC_ID"].strip()
        if "COSV" not in cosmic_id:
            continue

        gene_id = row["GENE_ID"].strip()
        transcript_id = row["TRANSCRIPT_ID"].strip()
        cellline_norm = normalise_name(row["SAMPLE_NAME"])
        aa_pos = safe_int(row["AA_POSITION"], 0)

        zygosity = row["MUTATION_ZYGOSITY"].replace(";", "").strip()
        if zygosity not in {"hom", "het"}:
            zygosity = "Unknown"

        rec = mutations_by_gene[gene_id].get(cosmic_id)
        if rec is None:
            rec = MutationRecord(gene_id=gene_id, cosmic_id=cosmic_id)
            mutations_by_gene[gene_id][cosmic_id] = rec

        rec.mutation_type = row["MUTATION_TYPE"].strip() or "-"
        if aa_pos > 0:
            rec.transcript_positions[transcript_id] = aa_pos
        rec.unique_ids.add(row["UNIQUE_MUTATION_ID"].strip())
        rec.aa_positions.add(row["AA_POSITION"].strip())
        rec.aa_changes.add(row["AA_CHANGE"].strip())
        rec.na_changes.add(row["NA_CHANGE"].strip())
        rec.tissues_by_cellline[cellline_norm] = row["SAMPLE_TISSUE"].strip()
        rec.zygosity_by_cellline[cellline_norm].add(zygosity)

        is_control = cellline_norm == control_norm
        if is_control:
            rec.tissues_by_cellline["control"] = row["SAMPLE_TISSUE"].strip()
            rec.zygosity_by_cellline["control"].add(zygosity)

        # Mutant peptides: definitionally at the site.
        mut_peptides = [p for p in row["OBSERVABLE_MUTANT_PEPTIDES"].split(";") if p]
        for pep in mut_peptides:
            rec.mutant_peptides.add(pep)
            peptide_mutation_map[pep][cellline_norm][gene_id][cosmic_id] = "M"
            if is_control:
                peptide_mutation_map[pep]["control"][gene_id][cosmic_id] = "M"

        # WT peptides: keep only those that actually report on the site, i.e.
        # they either cover the residue, or immediately flank a cleavage-site
        # mutation that the variant destroys/creates (so the WT-only peptide
        # disappears in the mutant proteoform).
        wt_peptides = [p for p in row["OBSERVABLE_WT_PEPTIDES"].split(";") if p]
        seq = gene_sequences.get(gene_id, {}).get(transcript_id)
        if seq is not None and aa_pos > 0:
            seq_il = seq  # already I->L
            for pep in wt_peptides:
                pep_il = replace_i_with_l(pep)
                idx = seq_il.find(pep_il)
                if idx < 0:
                    continue
                start = idx + 1
                end = start + len(pep_il) - 1
                # cover the residue, or sit immediately N/C-terminal to it
                if start - 1 <= aa_pos <= end + 1:
                    rec.reference_peptides.add(pep)
                    peptide_mutation_map[pep][cellline_norm][gene_id][cosmic_id] = "R"
                    if is_control:
                        peptide_mutation_map[pep]["control"][gene_id][cosmic_id] = "R"

    LOGGER.info(
        "Loaded %d mutations across %d genes",
        sum(len(v) for v in mutations_by_gene.values()), len(mutations_by_gene),
    )
    return mutations_by_gene, peptide_mutation_map


# =============================================================================
# Peptide-quant column resolution
# =============================================================================

# Fixed (canonical) column names in the peptide-quant file.
PQ = {
    "peptide": "PeptideSequence",
    "type": "Type",
    "decoy": "Decoy",
    "gene_ambiguity": "Gene_Ambiguity",
    "peptide_fdr": "Peptide_FDR",
    "best_psm_pep": "bestPSMScore(PEP)",
    "best_psm_fdr": "bestPSMFDR",
    "delta": "DeltaScore",
    "gene_name": "GeneNames",
    "gene_id": "Genes",
    "transcript_id": "Transcripts",
    "transcript_count": "ProteinCount",
    "missed": None,           # derived from sequence (see count_missed_cleavages)
    "all_psm_count": "PSMCount",
    "quant_psms": "QuantPSMs",
}

RAW_RE = re.compile(r"^TMT_Abundance_Raw_(\d+)$")
NORM_RE = re.compile(r"^TMT_Abundance_Normalised_(\d+)$")
SCALED_RE = re.compile(r"^TMT_Abundance_Scaled_(\d+)$")
PSM_TMT_RE = re.compile(r"^Raw_TMT_Intensity_(\d+)$")


def resolve_channel_headers(headers: List[str]) -> List[ChannelHeaders]:
    raw: Dict[int, str] = {}
    norm: Dict[int, str] = {}
    scaled: Dict[int, str] = {}
    for h in headers:
        for rx, store in ((RAW_RE, raw), (NORM_RE, norm), (SCALED_RE, scaled)):
            m = rx.match(h)
            if m:
                store[int(m.group(1))] = h
                break
    cids = sorted(set(raw) & set(norm) & set(scaled))
    if not cids:
        raise ValueError("No TMT channels resolved from peptide-quant headers.")
    out = [ChannelHeaders(c, raw[c], norm[c], scaled[c]) for c in cids]
    LOGGER.info("Resolved %d peptide-quant channels: %s", len(out), cids)
    return out


def count_missed_cleavages(peptide: str) -> int:
    """Count internal K/R not followed by P."""
    return len(split_tryptic_peptide(peptide)) - 1


def count_transcripts(transcript_field: str) -> int:
    """Number of distinct transcript IDs in a comma/semicolon-separated field."""
    parts = re.split(r"[,;]", transcript_field)
    return len({p.strip() for p in parts if p.strip()})


def peptide_passes(row: Dict[str, str]) -> bool:
    return (
        safe_float(row.get(PQ["peptide_fdr"]), 1.0) <= MAX_PEPTIDE_FDR
        and safe_float(row.get(PQ["best_psm_pep"]), 1.0) <= MAX_BEST_PSM_PEP
        and safe_float(row.get(PQ["delta"]), 0.0) > MIN_DELTA
    )


# =============================================================================
# First pass: peptide usability map
# =============================================================================

def build_usability_map(pepfile: str) -> Dict[str, PeptideUsability]:
    """First pass: map each fully-cleaved sub-peptide to usability info."""
    LOGGER.info("Building peptide usability map from: %s", pepfile)
    headers, rows = read_tsv(pepfile)
    require_columns(
        headers,
        [PQ["peptide"], PQ["type"], PQ["gene_ambiguity"], PQ["peptide_fdr"],
         PQ["best_psm_pep"], PQ["delta"]],
        "Peptide-quant (usability pass)",
    )

    usability: Dict[str, PeptideUsability] = {}
    for row in rows:
        if row.get(PQ["decoy"], "False").strip().lower() == "true":
            continue
        peptide = row[PQ["peptide"]].strip()
        peptide_type = row[PQ["type"]].strip()
        ambiguity = 0 if safe_int(row[PQ["gene_ambiguity"]], 0) == 1 else 1
        missed = count_missed_cleavages(peptide)
        passed = 1 if peptide_passes(row) else 0

        for subpep in split_tryptic_peptide(peptide):
            if len(subpep) < 6:
                continue
            entry = usability.get(subpep)
            if entry is None:
                usability[subpep] = PeptideUsability(
                    ambiguity=ambiguity, missed_cleavages=missed,
                    passed=passed, peptide_type=peptide_type,
                )
                entry = usability[subpep]
            else:
                if entry.ambiguity < ambiguity:
                    entry.ambiguity = ambiguity
                if entry.missed_cleavages > missed:
                    entry.missed_cleavages = missed
                    entry.peptide_type = peptide_type
                elif entry.missed_cleavages == missed and entry.peptide_type != peptide_type:
                    entry.peptide_type = "Un"
                if passed == 0:
                    entry.passed = 0
            if missed > 0 and passed == 1:
                usability[subpep].parent_by_missed[missed][peptide] = peptide_type

    LOGGER.info("Usability map: %d sub-peptides", len(usability))
    return usability


# =============================================================================
# Second pass: quantitative peptide data
# =============================================================================

def determine_typex(
    peptide: str, cellline_norm: str, gene_id: str,
    peptide_mutation_map: Dict[str, Dict[str, Dict[str, Dict[str, str]]]],
) -> str:
    """reference | nonmutant (WT-at-site) | mutant for this cell line."""
    typex = "reference"
    for subpep in split_tryptic_peptide(peptide):
        labels = (
            peptide_mutation_map.get(subpep, {})
            .get(cellline_norm, {})
            .get(gene_id, {})
        )
        for label in labels.values():
            if label == "M":
                return "mutant"
            if label == "R":
                typex = "nonmutant"
    return typex


def parse_peptide_data(
    pepfile: str,
    usability: Dict[str, PeptideUsability],
    peptide_mutation_map: Dict[str, Dict[str, Dict[str, Dict[str, str]]]],
    experiment: Experiment,
    channel_headers: List[ChannelHeaders],
) -> Tuple[
    Dict[str, PeptideRecord],
    Dict[str, Dict[int, GeneChannelRecord]],
    Dict[str, Set[str]],
    Dict[str, Set[str]],
]:
    """Second pass: build peptide records, gene aggregates and parent->stored map."""
    LOGGER.info("Reading peptide quant data: %s", pepfile)
    headers, rows = read_tsv(pepfile)

    peptide_records: Dict[str, PeptideRecord] = {}
    gene_records: Dict[str, Dict[int, GeneChannelRecord]] = defaultdict(dict)
    gene_peptide_map: Dict[str, Set[str]] = defaultdict(set)
    parent_peptide_map: Dict[str, Set[str]] = defaultdict(set)

    processed = 0
    for row in rows:
        if row.get(PQ["decoy"], "False").strip().lower() == "true":
            continue
        if not peptide_passes(row):
            continue
        if safe_int(row[PQ["gene_ambiguity"]], 0) != 1:
            continue
        if safe_int(row.get(PQ["quant_psms"], "0"), 0) < 1:
            continue

        processed += 1
        parent_peptide = row[PQ["peptide"]].strip()
        gene_id = row[PQ["gene_id"]].strip()
        gene_name = row[PQ["gene_name"]].strip()
        transcript_id = row[PQ["transcript_id"]].strip()
        protein_count = row[PQ["transcript_count"]].strip()
        missed = count_missed_cleavages(parent_peptide)

        # Channel depth totals use ALL passing signal (correct for normalisation).
        # Gene raw sums are NOT accumulated here; they are rebuilt after typing
        # from mutation-clean peptides only (see rebuild_clean_gene_records).
        for ch in channel_headers:
            raw_val = safe_float(row.get(ch.raw), 0.0)
            experiment.total_intensity_by_channel[ch.channel_id] += raw_val
            if ch.channel_id not in gene_records[gene_id]:
                gene_records[gene_id][ch.channel_id] = GeneChannelRecord()

        # determine the set of usable, fully-cleaved stored peptides
        pepset: Dict[str, str] = {}
        match = 0
        subpeptides = split_tryptic_peptide(parent_peptide)
        for subpep in subpeptides:
            if len(subpep) < MIN_SUBPEP_LEN:
                continue
            entry = usability.get(subpep)
            if entry is None:
                continue
            if entry.passed == 1 and entry.ambiguity == 0:
                if entry.missed_cleavages == 0:
                    pepset[subpep] = entry.peptide_type
                    match = 1
            else:
                match = -1

        current_mc = 1
        while match == 0 and current_mc <= missed:
            for subpep in subpeptides:
                if len(subpep) < MIN_SUBPEP_LEN:
                    continue
                entry = usability.get(subpep)
                if entry is None:
                    continue
                if entry.passed == 1 and entry.ambiguity == 0:
                    if current_mc in entry.parent_by_missed:
                        for parent_seq, ptype in entry.parent_by_missed[current_mc].items():
                            if parent_seq not in pepset:
                                pepset[parent_seq] = ptype
                            elif pepset[parent_seq] != ptype:
                                pepset[parent_seq] = "Un"
                            match = 1
                else:
                    match = -1
            current_mc += 1

        if match != 1:
            continue

        for stored_peptide, peptide_type in pepset.items():
            if peptide_type == "Un":
                continue

            rec = peptide_records.get(stored_peptide)
            if rec is None:
                rec = PeptideRecord(
                    peptide=stored_peptide,
                    peptide_type=peptide_type,
                    gene_id=gene_id,
                    gene_name=gene_name,
                    transcript_id=transcript_id,
                    transcript_count=count_transcripts(transcript_id),
                    protein_count=protein_count,
                    peptide_fdr=safe_float(row[PQ["peptide_fdr"]]),
                    best_psm_fdr=safe_float(row[PQ["best_psm_fdr"]]),
                    best_psm_pep=safe_float(row[PQ["best_psm_pep"]]),
                    delta=safe_float(row[PQ["delta"]]),
                    all_psm_count=safe_float(row.get(PQ["all_psm_count"], 0.0)),
                )
                for ch in channel_headers:
                    cq = rec.channel(ch.channel_id)
                    cq.typex = determine_typex(
                        stored_peptide,
                        normalise_name(experiment.channel_sample.get(ch.channel_id, "")),
                        gene_id, peptide_mutation_map,
                    )
                peptide_records[stored_peptide] = rec
            else:
                rec.peptide_fdr = max(rec.peptide_fdr, safe_float(row[PQ["peptide_fdr"]]))
                rec.best_psm_pep = max(rec.best_psm_pep, safe_float(row[PQ["best_psm_pep"]]))
                rec.best_psm_fdr = max(rec.best_psm_fdr, safe_float(row[PQ["best_psm_fdr"]]))
                rec.delta = max(rec.delta, safe_float(row[PQ["delta"]]))
                rec.all_psm_count += safe_float(row.get(PQ["all_psm_count"], 0.0))

            rec.parent_sequences.add(parent_peptide)
            parent_peptide_map[parent_peptide].add(stored_peptide)

            for ch in channel_headers:
                cq = rec.channel(ch.channel_id)
                cq.raw += safe_float(row.get(ch.raw), 0.0)
                cq.norm += safe_float(row.get(ch.norm), 0.0)
                cq.scaled += safe_float(row.get(ch.scaled), 0.0)
                cq.quant_count += 1

            gene_peptide_map[gene_id].add(stored_peptide)

    LOGGER.info("Parsed %d quant rows -> %d stored peptides", processed, len(peptide_records))
    return peptide_records, gene_records, gene_peptide_map, parent_peptide_map


# =============================================================================
# PSM spectral table
# =============================================================================

PSM_COLUMNS = ["Base_Peak_Intensity", "Precursor_Intensity", "Peptide", "PEP",
               "Experiment_Sample_ID"]


def parse_spectral_file(
    spectrafile: str,
    peptide_records: Dict[str, PeptideRecord],
    parent_peptide_map: Dict[str, Set[str]],
    experiment: Experiment,
) -> None:
    """Aggregate confident-PSM statistics into peptide records."""
    LOGGER.info("Reading spectral PSM table: %s", spectrafile)
    headers, rows = read_tsv(spectrafile)
    require_columns(headers, PSM_COLUMNS, "Spectral table")

    psm_channel_headers: Dict[int, str] = {}
    for h in headers:
        m = PSM_TMT_RE.match(h)
        if m:
            psm_channel_headers[int(m.group(1))] = h
    if not psm_channel_headers:
        raise ValueError("No Raw_TMT_Intensity_<N> columns in spectral table.")

    control_cid = experiment.control_channel_id
    if control_cid is None or control_cid not in psm_channel_headers:
        LOGGER.warning("Control channel %s absent from spectral table; PSM control "
                       "ratios will be skipped.", control_cid)
        control_header = None
    else:
        control_header = psm_channel_headers[control_cid]

    used = 0
    for row in rows:
        peptide = row.get("Peptide", "").strip()
        if peptide in {"", "-"}:
            continue
        if safe_float(row.get("PEP"), 1.0) > MAX_PSM_PEP:
            continue
        stored_set = parent_peptide_map.get(peptide)
        if not stored_set:
            continue

        control_int = safe_float(row.get(control_header), 0.0) if control_header else 0.0
        base_peak = safe_float(row.get("Base_Peak_Intensity"), 0.0)
        precursor = safe_float(row.get("Precursor_Intensity"), 0.0)

        # per-channel intensities once per row
        channel_ints = {cid: safe_float(row.get(hdr), 0.0)
                        for cid, hdr in psm_channel_headers.items()}
        max_int = max(channel_ints.values()) if channel_ints else 0.0

        for stored_peptide in stored_set:
            rec = peptide_records.get(stored_peptide)
            if rec is None:
                continue
            rec.psm_stats.base_peak_intensity_sum += base_peak
            rec.psm_stats.precursor_intensity_sum += precursor
            for cid, val in channel_ints.items():
                rec.psm_stats.per_channel_intensities[cid].append(val)
            rec.psm_stats.psm_vectors.append(dict(channel_ints))
            if control_int > 0 and max_int > 0:
                rec.psm_stats.max_control_ratios.append(control_int / max_int)
        used += 1

    LOGGER.info("Aggregated %d confident PSMs", used)


# =============================================================================
# Normalisation & gene-level calculations
# =============================================================================

def rebuild_clean_gene_records(
    gene_records: Dict[str, Dict[int, GeneChannelRecord]],
    gene_peptide_map: Dict[str, Set[str]],
    peptide_records: Dict[str, PeptideRecord],
    channel_headers: List[ChannelHeaders],
) -> None:
    """Rebuild gene raw sums from mutation-clean peptides only.

    A peptide is excluded from the gene-level reference if it is typed
    nonmutant (WT-at-site) or mutant in ANY channel - i.e. it overlaps a
    mutation in the mutated cell line, the control, or any cross-sample
    channel. Excluding it from every channel keeps the gene baseline
    internally consistent and free of mutation-driven dropout, so the gene
    reference reflects only the bulk protein level.
    """
    LOGGER.info("Rebuilding gene reference from mutation-clean peptides")
    channel_ids = [ch.channel_id for ch in channel_headers]
    n_excluded = 0
    for gene_id, peptides in gene_peptide_map.items():
        grecs = gene_records.get(gene_id)
        if grecs is None:
            continue
        for cid in channel_ids:
            grec = grecs.get(cid)
            if grec is None:
                grec = GeneChannelRecord()
                grecs[cid] = grec
            grec.raw_sum = 0.0
            grec.peptide_count = 0
        for peptide in peptides:
            rec = peptide_records.get(peptide)
            if rec is None:
                continue
            # exclude if affected by a mutation in any channel
            if any(cq.typex in ("nonmutant", "mutant") for cq in rec.channels.values()):
                n_excluded += 1
                continue
            for cid in channel_ids:
                raw_val = rec.channel(cid).raw
                if raw_val > 0:
                    grecs[cid].raw_sum += raw_val
                    grecs[cid].peptide_count += 1
    LOGGER.info("Excluded %d mutation-affected peptide(s) from gene references", n_excluded)


def calculate_gene_normalisation(
    gene_records: Dict[str, Dict[int, GeneChannelRecord]],
    experiment: Experiment,
) -> None:
    """Depth-normalise gene raw sums and compute ratio-to-control + scaled."""
    LOGGER.info("Calculating gene-level normalisation")
    control_cid = experiment.control_channel_id
    if control_cid is None:
        return
    control_total = experiment.total_intensity_by_channel[control_cid]

    for gene_id, channels in gene_records.items():
        control_raw = channels[control_cid].raw_sum if control_cid in channels else 0.0
        for cid, grec in channels.items():
            denom = experiment.total_intensity_by_channel[cid]
            grec.norm_sum = (grec.raw_sum / denom) * control_total if denom > 0 else 0.0
            grec.ratio_to_control = (grec.norm_sum / control_raw) if control_raw > 0 else None

        mean_control = channels[control_cid].norm_sum if control_cid in channels else 0.0
        for grec in channels.values():
            grec.mean_control_norm = mean_control
            grec.scaled = (grec.ratio_to_control * mean_control
                           if grec.ratio_to_control is not None else None)


def calculate_peptide_normalisation(
    peptide_records: Dict[str, PeptideRecord],
    experiment: Experiment,
) -> None:
    """Depth-normalise peptide raw values; compute per-channel ratio-to-control."""
    LOGGER.info("Calculating peptide-level normalisation")
    control_cid = experiment.control_channel_id
    if control_cid is None:
        return
    control_total = experiment.total_intensity_by_channel[control_cid]

    for rec in peptide_records.values():
        control_raw = rec.channel(control_cid).raw
        for cid, cq in rec.channels.items():
            denom = experiment.total_intensity_by_channel[cid]
            cq.noq = (cq.raw / denom) * control_total if denom > 0 else 0.0
            cq.ratio_to_control = (cq.noq / control_raw) if control_raw > 0 else None
        rec.mean_control_noq = rec.channel(control_cid).noq


def calculate_gene_topn_and_median(
    gene_records: Dict[str, Dict[int, GeneChannelRecord]],
    gene_peptide_map: Dict[str, Set[str]],
    peptide_records: Dict[str, PeptideRecord],
    experiment: Experiment,
) -> None:
    """Per gene/channel: median peptide ratio and top-3 ratio from reference peptides."""
    LOGGER.info("Calculating gene top-3 and median peptide ratios")
    control_cid = experiment.control_channel_id
    if control_cid is None:
        return

    for gene_id, channels in gene_records.items():
        peptides = gene_peptide_map.get(gene_id, set())
        for cid, grec in channels.items():
            pep_values: Dict[str, float] = {}
            for peptide in peptides:
                rec = peptide_records.get(peptide)
                if rec is None:
                    continue
                cq = rec.channel(cid)
                cc = rec.channel(control_cid)
                if cq.noq > 0 and cq.typex == "reference" and cc.typex == "reference":
                    pep_values[peptide] = cq.noq

            ratios: List[float] = []
            weighted_pairs: List[Tuple[float, float]] = []
            for peptide, value in pep_values.items():
                control_noq = peptide_records[peptide].channel(control_cid).noq
                if control_noq > 0:
                    ratio = value / control_noq
                    ratios.append(ratio)
                    # weight by the peptide's PSM evidence so low-evidence
                    # peptides do not count equally with well-sampled ones
                    weight = peptide_records[peptide].all_psm_count
                    weighted_pairs.append((ratio, weight if weight > 0 else 1.0))
            grec.median_ratio = median(ratios) if ratios else None
            grec.weighted_median_ratio = weighted_median(weighted_pairs)

            top = sorted(pep_values, key=lambda p: pep_values[p], reverse=True)[:3]
            sum_channel = sum(pep_values[p] for p in top)
            sum_control = sum(peptide_records[p].channel(control_cid).noq for p in top)
            grec.topn_norm_sum = sum_channel
            grec.topn_ratio = (sum_channel / sum_control) if sum_control > 0 else None


# =============================================================================
# Cross-sample reference (absolute and within-gene-fraction)
# =============================================================================

@dataclass
class CrossSampleResult:
    """Result of a cross-sample reference computation for one peptide/channel."""
    test_value: Optional[float] = None      # the test channel's value (noq or fraction)
    reference: Optional[float] = None       # trimmed-mean reference from other channels
    ratio: Optional[float] = None           # test_value / reference
    n_channels: int = 0                     # qualifying reference channels
    robust_sd: Optional[float] = None       # MAD-scaled SD of reference channels
    robust_z: Optional[float] = None        # (test - reference) / robust_sd
    sufficient: bool = False                # >= MIN_CROSS_SAMPLE_CHANNELS qualifying


def _qualifying_reference_channels(
    rec: PeptideRecord,
    channel_id: int,
    control_cid: Optional[int],
    gene_records: Dict[str, Dict[int, GeneChannelRecord]],
    peptide_records: Dict[str, PeptideRecord],
) -> List[int]:
    """Reference channels passing the gene-expression and peptide-signal floors.

    A channel qualifies as a reference for this peptide only if:
      * it is not the test channel or the control,
      * this peptide is not itself mutation-affected (nonmutant/mutant) there,
      * the gene is expressed there at >= MIN_GENE_EXPR_FRACTION of the gene's
        across-channel median expression, and
      * the WT peptide signal there is >= MIN_PEP_SIGNAL_FRACTION of the
        peptide's across-channel median signal.
    """
    gene_id = rec.gene_id
    grecs = gene_records.get(gene_id, {})

    # gene across-channel median expression (norm_sum), over channels with signal
    gene_norms = [g.norm_sum for cid, g in grecs.items()
                  if cid not in (channel_id, control_cid) and g.norm_sum > 0]
    gene_med = median(gene_norms) if gene_norms else 0.0

    # peptide across-channel median signal (noq), over other channels with signal
    pep_noqs = [cq.noq for cid, cq in rec.channels.items()
                if cid not in (channel_id, control_cid) and cq.noq > 0]
    pep_med = median(pep_noqs) if pep_noqs else 0.0

    qualifying: List[int] = []
    for cid, cq in rec.channels.items():
        if cid == channel_id or cid == control_cid:
            continue
        if cq.noq <= 0:
            continue
        if cq.typex in ("nonmutant", "mutant"):
            continue
        # gene-expression floor
        g = grecs.get(cid)
        if gene_med > 0 and (g is None or g.norm_sum < MIN_GENE_EXPR_FRACTION * gene_med):
            continue
        # peptide-signal floor
        if pep_med > 0 and cq.noq < MIN_PEP_SIGNAL_FRACTION * pep_med:
            continue
        qualifying.append(cid)
    return qualifying


def cross_sample_reference_absolute(
    rec: PeptideRecord,
    channel_id: int,
    control_cid: Optional[int],
    qualifying: List[int],
) -> CrossSampleResult:
    """Absolute-NOQ cross-sample reference over qualifying channels.

    Compares the peptide's depth-normalised abundance in the test channel to
    the trimmed mean of the same peptide across qualifying reference channels.
    Sensitive to genuine between-sample differences in total gene abundance.
    """
    result = CrossSampleResult()
    test = rec.channel(channel_id).noq
    result.test_value = test if test > 0 else None
    others = [rec.channel(cid).noq for cid in qualifying if rec.channel(cid).noq > 0]
    result.n_channels = len(others)
    result.sufficient = len(others) >= MIN_CROSS_SAMPLE_CHANNELS
    if not others or test <= 0:
        return result
    ref = trimmed_mean(others)
    result.reference = ref
    if ref > 0:
        result.ratio = test / ref
    if len(others) >= 3:
        med = median(others)
        mad = median([abs(x - med) for x in others])
        if mad > 0:
            result.robust_sd = 1.4826 * mad
            result.robust_z = (test - ref) / result.robust_sd
    return result


def cross_sample_reference_genefraction(
    rec: PeptideRecord,
    channel_id: int,
    control_cid: Optional[int],
    qualifying: List[int],
    gene_records: Dict[str, Dict[int, GeneChannelRecord]],
) -> CrossSampleResult:
    """Within-gene-fraction cross-sample reference.

    Expresses the peptide as a fraction of its gene's signal within each
    channel (peptide noq / gene norm_sum), then compares the test-channel
    fraction to the trimmed mean of that fraction across qualifying reference
    channels. Removes the between-sample gene-expression confound by
    construction; the trade-off is that mutations also changing total gene
    abundance (e.g. NMD) are partly cancelled.
    """
    result = CrossSampleResult()
    grecs = gene_records.get(rec.gene_id, {})

    def fraction(cid: int) -> Optional[float]:
        g = grecs.get(cid)
        if g is None or g.norm_sum <= 0:
            return None
        noq = rec.channel(cid).noq
        return (noq / g.norm_sum) if noq > 0 else 0.0

    test_frac = fraction(channel_id)
    result.test_value = test_frac
    others = [f for f in (fraction(cid) for cid in qualifying) if f is not None and f > 0]
    result.n_channels = len(others)
    result.sufficient = len(others) >= MIN_CROSS_SAMPLE_CHANNELS
    if test_frac is None or not others:
        return result
    ref = trimmed_mean(others)
    result.reference = ref
    if ref > 0:
        result.ratio = test_frac / ref
    if len(others) >= 3:
        med = median(others)
        mad = median([abs(x - med) for x in others])
        if mad > 0:
            result.robust_sd = 1.4826 * mad
            result.robust_z = (test_frac - ref) / result.robust_sd
    return result


# =============================================================================
# Mutant-fraction estimate, bootstrap CI and orthogonal mutant/WT check (new)
# =============================================================================

def mutant_fraction_from_ratio(ratio: Optional[float]) -> Optional[float]:
    """Convert a WT-peptide observed:expected ratio into a mutant-proteoform
    fraction estimate, clamped to [0, 1].

    The WT-at-site peptide is produced only by the wild-type proteoform, so if
    the gene's total protein level is unchanged, the WT-peptide signal scales
    with the wild-type fraction: observed = expected * (1 - mutant_fraction).
    Hence mutant_fraction ~= 1 - ratio. A heterozygous, equally expressed
    mutation is expected near 0.5; a homozygous one near 1.0.
    """
    if ratio is None:
        return None
    frac = 1.0 - ratio
    return min(1.0, max(0.0, frac))


def _channel_depth_factors(experiment: Experiment) -> Dict[int, float]:
    """Per-channel multiplicative factor that maps a raw PSM intensity onto the
    depth-normalised (control-equivalent) scale, matching the peptide NOQ
    normalisation: factor = control_total / channel_total."""
    control_cid = experiment.control_channel_id
    factors: Dict[int, float] = {}
    if control_cid is None:
        return factors
    control_total = experiment.total_intensity_by_channel.get(control_cid, 0.0)
    for cid, total in experiment.total_intensity_by_channel.items():
        factors[cid] = (control_total / total) if total > 0 else 0.0
    return factors


def bootstrap_cross_sample_ci(
    rec: PeptideRecord,
    channel_id: int,
    control_cid: Optional[int],
    depth_factors: Dict[int, float],
    other_channels: List[int],
    n_boot: int = 1000,
    rng: Optional["random.Random"] = None,
) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """Bootstrap a 95% CI for the cross-sample ratio by resampling PSMs.

    Each PSM contributes an aligned channel-intensity vector. We resample whole
    PSMs (rows) with replacement; for each resample we depth-normalise, sum per
    channel, form the test-channel value and the trimmed-mean reference over the
    other channels, and take their ratio. Returns (point_estimate, lo95, hi95).
    Requires >= MIN_BOOTSTRAP_PSMS PSMs, otherwise returns (None, None, None).
    """
    vectors = rec.psm_stats.psm_vectors
    n = len(vectors)
    if n < MIN_BOOTSTRAP_PSMS or not other_channels:
        return None, None, None
    if rng is None:
        rng = random

    def ratio_for(sample_idx: List[int]) -> Optional[float]:
        test_sum = 0.0
        other_sums: Dict[int, float] = {c: 0.0 for c in other_channels}
        for i in sample_idx:
            vec = vectors[i]
            test_sum += vec.get(channel_id, 0.0) * depth_factors.get(channel_id, 0.0)
            for c in other_channels:
                other_sums[c] += vec.get(c, 0.0) * depth_factors.get(c, 0.0)
        others = [v for v in other_sums.values() if v > 0]
        if not others or test_sum <= 0:
            return None
        ref = trimmed_mean(others)
        return (test_sum / ref) if ref > 0 else None

    point = ratio_for(list(range(n)))
    boot: List[float] = []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        r = ratio_for(idx)
        if r is not None:
            boot.append(r)
    if len(boot) < n_boot * 0.5:
        return point, None, None
    boot.sort()
    lo = boot[int(0.025 * len(boot))]
    hi = boot[min(len(boot) - 1, int(0.975 * len(boot)))]
    return point, lo, hi


def orthogonal_mut_wt_ratio(
    rec_wt: PeptideRecord,
    mutant_peptides: Set[str],
    peptide_records: Dict[str, PeptideRecord],
    channel_id: int,
    gene_id: str,
    cellline_norm: str,
    peptide_mutation_map: Dict[str, Dict[str, Dict[str, Dict[str, str]]]],
) -> Tuple[Optional[float], int]:
    """Direct mutant-fraction estimate from co-quantified mutant peptides.

    Where a mutant peptide for the same site is also quantified in this channel,
    mutant_fraction ~= mut_noq / (mut_noq + wt_noq). Averaged over available
    mutant peptides. Returns (fraction, n_mutant_peptides_used).
    """
    wt_noq = rec_wt.channel(channel_id).noq
    fractions: List[float] = []
    for mpep in mutant_peptides:
        mrec = peptide_records.get(mpep)
        if mrec is None:
            continue
        # confirm this mutant peptide is flagged 'M' for this cell line/gene
        labels = (
            peptide_mutation_map.get(mpep, {})
            .get(cellline_norm, {})
            .get(gene_id, {})
        )
        if "M" not in labels.values():
            continue
        mut_noq = mrec.channel(channel_id).noq
        denom = mut_noq + wt_noq
        if denom > 0:
            fractions.append(mut_noq / denom)
    if not fractions:
        return None, 0
    return mean(fractions), len(fractions)


# =============================================================================
# Output
# =============================================================================

OUTPUT_HEADER = [
    "Experiment", "Peptide", "Type", "TypeX", "ENSG", "GENE", "Transcript_Count",
    "Protein_Count",
    "Transcripts", "Channel", "TMT_Label", "Cell_Line", "Peptide_FDR",
    "Best_PSM_FDR", "Best_PSM_PEP", "Delta", "All_PSM_Count",
    "Pep_Raw_Int", "Control_Raw_Int", "Pep_Norm_Int", "Control_Norm_Int",
    "Pep_NOQ", "Control_NOQ", "Pep_Control_Ratio",
    "Retained_PSM_Count", "PSM_Mean_Ratio_to_Control", "PSM_Ratio_Variance",
    "Sum_Precursor_Int", "Sum_Base_Peak_Int", "PSM_Control_Mean_MaxRatio",
    "Gene_Raw_Sum", "Gene_Control_Raw_Sum", "Gene_Norm_Sum", "Gene_Control_Norm_Sum",
    "Gene_Norm_Ratio", "Gene_Log2_FC", "Gene_Peptide_Count",
    "Gene_Top3_Norm_Sum", "Gene_Top3_Ratio", "Gene_Top3_Log2FC", "Gene_Median_Ratio",
    "Gene_Weighted_Median_Ratio",
    # --- bridge/control dropout ---
    "DropOut_vs_GeneRatio", "PctDropOut_vs_GeneRatio",
    "DropOut_vs_GeneMedian", "PctDropOut_vs_GeneMedian",
    # --- cross-sample dropout: absolute-NOQ approach ---
    "CrossSample_Ref_NOQ", "CrossSample_N_Channels", "CrossSample_Ratio",
    "CrossSample_Log2_FC", "DropOut_CrossSample", "PctDropOut_CrossSample",
    "CrossSample_Robust_Z", "CrossSample_Sufficient",
    # --- cross-sample dropout: within-gene-fraction approach ---
    "CrossSampleGF_N_Channels", "CrossSampleGF_Ratio", "CrossSampleGF_Log2_FC",
    "PctDropOut_CrossSampleGF", "CrossSampleGF_Robust_Z", "CrossSampleGF_Sufficient",
    # --- mutant-fraction estimates and bootstrap CI ---
    "MutantFraction_vs_GeneRatio", "MutantFraction_CrossSample",
    "MutantFraction_CrossSampleGF",
    "MutantFraction_Expected", "CrossSample_Ratio_Boot",
    "CrossSample_Ratio_CI_Lo", "CrossSample_Ratio_CI_Hi",
    "MutantFraction_CrossSample_CI_Lo", "MutantFraction_CrossSample_CI_Hi",
    "CrossSample_Significant",
    # --- orthogonal mutant/WT check ---
    "Orthogonal_MutantFraction", "Orthogonal_N_MutPeptides",
    # --- mutation annotation ---
    "CosmicID", "Zygosity", "Tissue", "MutationType", "TranscriptIDs",
    "Unique_Mutation_IDs", "AA_Pos", "AA_Change", "NA_Change",
    "RefPeptides", "MutPeptides", "Control_Mutation",
    "Peptide_Position", "Peptide_is_at_site",
]


def select_cosmic_id(
    peptide: str, cellline_norm: str, gene_id: str,
    peptide_mutation_map: Dict[str, Dict[str, Dict[str, Dict[str, str]]]],
) -> str:
    ids = list(
        peptide_mutation_map.get(peptide, {})
        .get(cellline_norm, {})
        .get(gene_id, {})
        .keys()
    )
    if len(ids) > 1:
        return "multiple_mutations"
    if len(ids) == 1:
        return ids[0]
    return "-"


def join_set(values: Set[str]) -> str:
    return ";".join(sorted(v for v in values if v)) or "-"


def write_output(
    output_file: str,
    experiment: Experiment,
    peptide_records: Dict[str, PeptideRecord],
    gene_records: Dict[str, Dict[int, GeneChannelRecord]],
    mutations_by_gene: Dict[str, Dict[str, MutationRecord]],
    peptide_mutation_map: Dict[str, Dict[str, Dict[str, Dict[str, str]]]],
    gene_sequences: Dict[str, Dict[str, str]],
    n_bootstrap: int = 1000,
    seed: int = 1234,
) -> None:
    LOGGER.info("Writing output: %s", output_file)
    control_cid = experiment.control_channel_id
    depth_factors = _channel_depth_factors(experiment)
    rng = random.Random(seed)
    n_written = n_nonmutant = 0
    genes_seen: Set[str] = set()
    celllines_seen: Set[str] = set()

    with open(output_file, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(OUTPUT_HEADER)

        for peptide, rec in peptide_records.items():
            gene_id = rec.gene_id
            psm_control_maxratio = mean(rec.psm_stats.max_control_ratios)

            for channel_id, cq in rec.channels.items():
                if channel_id == control_cid:
                    continue
                if cq.typex != "nonmutant":       # only WT-at-site peptides
                    continue

                n_nonmutant += 1
                cellline = experiment.channel_sample.get(channel_id, "-")
                cellline_norm = normalise_name(cellline)
                tmt_label = experiment.channel_label.get(channel_id, "-")
                genes_seen.add(gene_id)
                celllines_seen.add(cellline)

                grec = gene_records.get(gene_id, {}).get(channel_id)
                gctrl = gene_records.get(gene_id, {}).get(control_cid) if control_cid is not None else None
                if grec is None or grec.peptide_count <= MIN_GENE_PEPTIDES:
                    continue

                cosmic_id = select_cosmic_id(peptide, cellline_norm, gene_id, peptide_mutation_map)
                if cosmic_id == "multiple_mutations":
                    continue

                retained_psms = len(rec.psm_stats.max_control_ratios)
                psm_vector_count = len(rec.psm_stats.psm_vectors)

                # --- bridge/control dropout (only when PSM evidence exists) ---
                drop_gene = pct_gene = None
                drop_median = pct_median = None
                if psm_vector_count > 0:
                    if (cq.ratio_to_control is not None and grec.ratio_to_control is not None
                            and grec.ratio_to_control > 0):
                        drop_gene = ((100 * (cq.ratio_to_control / grec.ratio_to_control)) - 50) / 50
                        pct_gene = 50 - (50 * drop_gene)
                    if (cq.ratio_to_control is not None and grec.median_ratio is not None
                            and grec.median_ratio > 0):
                        drop_median = ((100 * (cq.ratio_to_control / grec.median_ratio)) - 50) / 50
                        pct_median = 50 - (50 * drop_median)

                # --- qualifying reference channels (shared by both approaches) ---
                qualifying = _qualifying_reference_channels(
                    rec, channel_id, control_cid, gene_records, peptide_records,
                )

                # --- cross-sample dropout: absolute-NOQ approach ---
                cs_abs = cross_sample_reference_absolute(rec, channel_id, control_cid, qualifying)
                cs_ref = cs_abs.reference
                cs_n = cs_abs.n_channels
                cs_ratio = cs_abs.ratio
                cs_z = cs_abs.robust_z
                cs_log2 = log2(cs_ratio) if cs_ratio else None
                drop_cs = pct_cs = None
                if cs_ratio is not None:
                    drop_cs = (100 * cs_ratio - 50) / 50
                    pct_cs = 50 - (50 * drop_cs)
                cs_qual = "YES" if cs_abs.sufficient else "NO"

                # --- cross-sample dropout: within-gene-fraction approach ---
                cs_gf = cross_sample_reference_genefraction(
                    rec, channel_id, control_cid, qualifying, gene_records,
                )
                gf_ratio = cs_gf.ratio
                gf_n = cs_gf.n_channels
                gf_z = cs_gf.robust_z
                gf_log2 = log2(gf_ratio) if gf_ratio else None
                pct_gf = None
                if gf_ratio is not None:
                    pct_gf = 50 - (50 * ((100 * gf_ratio - 50) / 50))
                gf_qual = "YES" if cs_gf.sufficient else "NO"

                # --- mutation annotation ---
                zygosity = tissue = mut_type = "-"
                transcript_ids = unique_ids = aa_pos = aa_change = na_change = "-"
                ref_peps = mut_peps = "-"
                mrec = mutations_by_gene.get(gene_id, {}).get(cosmic_id) if cosmic_id != "-" else None
                if mrec is not None:
                    zset = mrec.zygosity_by_cellline.get(cellline_norm, set())
                    zjoin = "".join(sorted(zset))
                    zygosity = zjoin if zjoin in {"het", "hom"} else "-"
                    tissue = mrec.tissues_by_cellline.get(cellline_norm, "-")
                    mut_type = mrec.mutation_type
                    transcript_ids = join_set(set(mrec.transcript_positions))
                    unique_ids = join_set(mrec.unique_ids)
                    aa_pos = join_set(mrec.aa_positions)
                    aa_change = join_set(mrec.aa_changes)
                    na_change = join_set(mrec.na_changes)
                    ref_peps = join_set(mrec.reference_peptides)
                    mut_peps = join_set(mrec.mutant_peptides)

                control_mutation = "YES" if (
                    "control" in peptide_mutation_map.get(peptide, {})
                    and gene_id in peptide_mutation_map[peptide]["control"]
                ) else "NO"

                positions = peptide_positions_in_gene(gene_id, peptide, gene_sequences)
                pos_str = ";".join(f"{t}[{s}-{e}]" for t, s, e in positions) if positions else "-"

                at_site = "NO"
                if mrec is not None:
                    for enst, start, end in positions:
                        if enst in mrec.transcript_positions:
                            p = mrec.transcript_positions[enst]
                            if start - 1 <= p <= end + 1:
                                at_site = "YES"
                                break

                if at_site != "YES" or control_mutation == "YES":
                    continue

                # --- mutant-fraction estimates ---
                # WT-peptide ratio -> mutant fraction (clamped to [0,1]).
                mf_gene = None
                if psm_vector_count > 0:
                    mf_gene = mutant_fraction_from_ratio(
                        (cq.ratio_to_control / grec.ratio_to_control)
                        if (cq.ratio_to_control is not None
                            and grec.ratio_to_control is not None
                            and grec.ratio_to_control > 0) else None
                    )
                mf_cross = mutant_fraction_from_ratio(cs_ratio)
                mf_genefrac = mutant_fraction_from_ratio(gf_ratio)
                mf_expected = {"het": 0.5, "hom": 1.0}.get(zygosity)

                # --- bootstrap CI on the cross-sample (absolute) ratio ---
                cs_boot = cs_lo = cs_hi = mf_lo = mf_hi = cs_sig = None
                if n_bootstrap > 0 and qualifying:
                    cs_boot, cs_lo, cs_hi = bootstrap_cross_sample_ci(
                        rec, channel_id, control_cid, depth_factors,
                        qualifying, n_boot=n_bootstrap, rng=rng,
                    )
                    if cs_lo is not None and cs_hi is not None:
                        # mutant-fraction CI is the inverted ratio CI
                        mf_lo = mutant_fraction_from_ratio(cs_hi)
                        mf_hi = mutant_fraction_from_ratio(cs_lo)
                        # significant dropout if the whole ratio CI sits below 1
                        cs_sig = "YES" if cs_hi < 1.0 else "NO"

                # --- orthogonal mutant/WT fraction (new) ---
                ortho_frac = ortho_n = None
                if mrec is not None and mrec.mutant_peptides:
                    ortho_frac, ortho_n = orthogonal_mut_wt_ratio(
                        rec, mrec.mutant_peptides, peptide_records,
                        channel_id, gene_id, cellline_norm, peptide_mutation_map,
                    )

                # PSM-level mean channel:control intensity ratio
                psm_chan = rec.psm_stats.per_channel_intensities.get(channel_id, [])
                psm_ctrl = (rec.psm_stats.per_channel_intensities.get(control_cid, [])
                            if control_cid is not None else [])
                psm_mean_ratio: Optional[float] = None
                if psm_chan and psm_ctrl:
                    ctrl_mean = mean(psm_ctrl)
                    if ctrl_mean > 0:
                        psm_mean_ratio = mean(psm_chan) / ctrl_mean

                gctrl_raw = gctrl.raw_sum if gctrl else "-"
                gctrl_norm = gctrl.norm_sum if gctrl else "-"
                control_raw = rec.channel(control_cid).raw if control_cid is not None else "-"
                control_noq = rec.channel(control_cid).noq if control_cid is not None else "-"

                row = [
                    experiment.experiment_id, peptide, rec.peptide_type, cq.typex,
                    gene_id, rec.gene_name, rec.transcript_count, rec.protein_count,
                    rec.transcript_id,
                    channel_id, tmt_label, cellline,
                    rec.peptide_fdr, rec.best_psm_fdr, rec.best_psm_pep, rec.delta,
                    rec.all_psm_count,
                    cq.raw, control_raw, cq.norm,
                    rec.channel(control_cid).norm if control_cid is not None else "-",
                    cq.noq, control_noq, cq.ratio_to_control,
                    len(rec.psm_stats.max_control_ratios), psm_mean_ratio,
                    variance(psm_chan),
                    rec.psm_stats.precursor_intensity_sum, rec.psm_stats.base_peak_intensity_sum,
                    psm_control_maxratio,
                    grec.raw_sum, gctrl_raw, grec.norm_sum, gctrl_norm,
                    grec.ratio_to_control, log2(grec.ratio_to_control) if grec.ratio_to_control else None,
                    grec.peptide_count, grec.topn_norm_sum, grec.topn_ratio,
                    log2(grec.topn_ratio) if grec.topn_ratio else None, grec.median_ratio,
                    grec.weighted_median_ratio,
                    drop_gene, pct_gene, drop_median, pct_median,
                    cs_ref, cs_n, cs_ratio, cs_log2, drop_cs, pct_cs, cs_z, cs_qual,
                    gf_n, gf_ratio, gf_log2, pct_gf, gf_z, gf_qual,
                    mf_gene, mf_cross, mf_genefrac, mf_expected, cs_boot, cs_lo, cs_hi,
                    mf_lo, mf_hi, cs_sig,
                    ortho_frac, ortho_n,
                    cosmic_id, zygosity, tissue, mut_type, transcript_ids, unique_ids,
                    aa_pos, aa_change, na_change, ref_peps, mut_peps, control_mutation,
                    pos_str, at_site,
                ]
                writer.writerow([fmt(c) for c in row])
                n_written += 1

    LOGGER.info("Output complete: %s", output_file)
    LOGGER.info("  WT-at-site peptide/channel observations: %d", n_nonmutant)
    LOGGER.info("  Rows written: %d", n_written)
    LOGGER.info("  Genes with mutant proteoforms: %d", len(genes_seen))
    LOGGER.info("  Cell lines represented: %d", len(celllines_seen))


# =============================================================================
# Experiment assembly
# =============================================================================

def build_experiment(
    spectrafile: str,
    sample_table: Dict[str, Dict[str, str]],
    control_name: str,
    channel_headers: List[ChannelHeaders],
) -> Experiment:
    """Determine the experiment id (from the spectral table) and wire up the
    channel -> TMT label -> cell line mapping using the sample table."""
    headers, rows = read_tsv(spectrafile)
    require_columns(headers, ["Experiment_Sample_ID"], "Spectral table")
    exp_ids = {r["Experiment_Sample_ID"].strip() for r in rows
               if r.get("Experiment_Sample_ID", "").strip()}
    if len(exp_ids) != 1:
        raise ValueError(f"Expected exactly one Experiment_Sample_ID in spectral "
                         f"table, found: {sorted(exp_ids)}")
    experiment_id = next(iter(exp_ids))
    LOGGER.info("Experiment: %s", experiment_id)

    if experiment_id not in sample_table:
        raise KeyError(f"Experiment '{experiment_id}' not present in sample table. "
                       f"Available: {sorted(sample_table)}")
    label_to_sample = sample_table[experiment_id]
    control_norm = normalise_name(control_name)

    experiment = Experiment(experiment_id=experiment_id)
    present_channels = {ch.channel_id for ch in channel_headers}
    for cid in sorted(present_channels):
        label = TMT_PLEX_ORDER.get(cid)
        if label is None:
            LOGGER.warning("Channel %d has no TMT-plex label; skipping.", cid)
            continue
        sample = label_to_sample.get(label)
        if sample is None:
            LOGGER.warning("TMT label %s (channel %d) not in sample table for %s.",
                           label, cid, experiment_id)
            sample = "-"
        experiment.channel_label[cid] = label
        experiment.channel_sample[cid] = sample
        if normalise_name(sample) == control_norm:
            experiment.control_channel_id = cid

    if experiment.control_channel_id is None:
        raise ValueError(f"Control sample '{control_name}' not found among channels "
                         f"for experiment {experiment_id}.")
    LOGGER.info("Control channel: %d (%s = %s)", experiment.control_channel_id,
                experiment.channel_label[experiment.control_channel_id],
                experiment.channel_sample[experiment.control_channel_id])
    return experiment


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args = parse_args()
    setup_logging(args.verbose)
    LOGGER.info("Starting mutant-proteoform analysis")

    sample_table = load_sample_table(args.samplefile)
    gene_sequences = parse_fasta(args.fasta)

    # Resolve channel headers once from the peptide-quant header.
    with open(args.pepfile) as fh:
        header_line = fh.readline().rstrip("\n").rstrip("\r").split("\t")
    channel_headers = resolve_channel_headers([h.strip() for h in header_line])

    experiment = build_experiment(args.spectrafile, sample_table,
                                  args.control_channel_name, channel_headers)

    mutations_by_gene, peptide_mutation_map = parse_mutation_file(
        args.mutfile, gene_sequences, args.control_channel_name,
    )

    usability = build_usability_map(args.pepfile)

    peptide_records, gene_records, gene_peptide_map, parent_peptide_map = parse_peptide_data(
        args.pepfile, usability, peptide_mutation_map, experiment, channel_headers,
    )

    parse_spectral_file(args.spectrafile, peptide_records, parent_peptide_map, experiment)

    rebuild_clean_gene_records(gene_records, gene_peptide_map, peptide_records, channel_headers)
    calculate_gene_normalisation(gene_records, experiment)
    calculate_peptide_normalisation(peptide_records, experiment)
    calculate_gene_topn_and_median(gene_records, gene_peptide_map, peptide_records, experiment)

    write_output(
        args.output, experiment, peptide_records, gene_records,
        mutations_by_gene, peptide_mutation_map, gene_sequences,
        n_bootstrap=args.bootstrap, seed=args.seed,
    )
    LOGGER.info("Done.")


if __name__ == "__main__":
    main()
