#!/usr/bin/env python3
"""
Parse and Annotate Peptide Results with Sample and Mutation Information

Dr James Wright (2026),
The Institute of Cancer Research, London
"""

import argparse
import csv
import logging, sys
import re
import os
import glob
import numpy as np
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple


# =============================================================================
# Define data classes
# =============================================================================


#PSM information data class
@dataclass
class PSMInfo:
    score: float = 99
    fdr: float = 99
    modifications: set[str] = field(default_factory=set)
    peptidoform: Optional[str] = None
    rt: Optional[float] = None
    mz: Optional[float] = None
    charge: Optional[int] = None
    conflicts: Dict[str, List[Tuple[str, float]]] = field(default_factory=dict)  # database -> list of (peptide, score) for conflicting assignments
    delta_score: Optional[float] = None  # difference in score between this and next best peptide assignment (maybe negative if this is not the best assignment)
    best_peptide: bool = False  # flag to indicate if this PSM has the best score for this spectrum

#Peptide information data class
@dataclass
class PeptideInfo:
    type: str = "unknown"
    decoy: bool = True
    top_rank: bool = False
    psm_score: float = 99
    psm_fdr: float = 99
    delta_score: float = -1
    best_psm: Optional[str] = None
    peptide_fdr: float = 99
    database_fdr: Dict[str, float] = field(default_factory=dict)
    database_score: Dict[str, float] = field(default_factory=dict)
    group_fdr: Dict[str, float] = field(default_factory=dict)
    group_lfdr: Dict[str, float] = field(default_factory=dict)
    psms: set[str] = field(default_factory=set)
    protein_ids: Dict[str, set[str]] = field(default_factory=lambda: defaultdict(set))
    genes: set[str] = field(default_factory=set)
    transcripts: set[str] = field(default_factory=set)
    experiments: Dict[str, float] = field(default_factory=dict)
    search_engines: Dict[str, float] = field(default_factory=dict)
    preaa: set[str] = field(default_factory=set)
    postaa: set[str] = field(default_factory=set)
    modifications: set[str] = field(default_factory=set)
    peptidoforms: set[str] = field(default_factory=set)
    minSAAVdistance: int = 0
    canonical_sequence: Optional[str] = None
    quantification: Dict[str, Dict[str, float]] = field(default_factory=dict)
    conflicts: Dict[str, float] = field(default_factory=dict)  

#Protein information data class
@dataclass
class ProteinInfo:
    type: str = "unknown"
    decoy: bool = True
    score: float = 0
    gene_ids: set[str] = field(default_factory=set)
    transcript_ids: set[str] = field(default_factory=set)
    peptides: set[str] = field(default_factory=set)
    proteotypic_peptides: set[str] = field(default_factory=set)
    experiments: set[str] = field(default_factory=set)
    search_engines: set[str] = field(default_factory=set)
    databases: set[str] = field(default_factory=set)
    length: Optional[int] = 0

@dataclass
class ProteinGroupInfo:
    master_protein: Optional[str] = None
    info : ProteinInfo = ProteinInfo()
    same_set_protein_ids: set[str] = field(default_factory=set)
    subset_protein_ids: set[str] = field(default_factory=set)

#Modification information data class
@dataclass
class ModificationInfo:
    site : str = ""
    pos: str = ""
    info: str = ""

#MZTAB Meta information data class
@dataclass
class MzTabMetaInfo:
    variable_modifications: Dict[str, ModificationInfo] = field(default_factory=dict)
    fixed_modifications: Dict[str, ModificationInfo] = field(default_factory=dict)
    quantification_method: Optional[str] = None
    mzml_file: Optional[str] = None

# =============================================================================
# Logging
# =============================================================================

#Setup Logging Handle
LOG = logging.getLogger("results_annotator")

# Function to configure logging based on the specified log level
def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
    )

# =============================================================================
# Get user parameters
# =============================================================================

# Function to load command line arguments
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='''
            Parse and Annotate Peptide and PSM Results with Sample and optionally Mutation Information
            ---- James Wright (2025) - The Institute of Cancer Research, London ----'''
    )

    parser.add_argument('--resultsfolder', '-d', dest='results_folder', help='Directory Containing Results *.mzTab and *.tab', required=True)
    parser.add_argument('--sampleinfo', '-s', dest='sample_file', default=None, help='Sample Information File if available (TSV with ExperimentID, TMTLabel, and SampleID Columns). Default=None')
    parser.add_argument('--referencedb', '-r', dest='ref_fasta', default='reference.fa', help='Reference Protein Sequences, Default=reference.fa')
    parser.add_argument('--peptidefdrthreshold', '-f', type=float, dest='pep_fdr_threshold', default=0.05, help='PEPTIDE FDR Threshold to Filter Results.  Default=0.05')
    parser.add_argument('--psmpepthreshold', '-pm', type=float, dest='psm_score_threshold', default=0.05, help='Medium PSM PEP Threshold. Default=0.05')
    parser.add_argument('--psmfdrthreshold', '-fm', type=float, dest='psm_fdr_threshold', default=0.05, help='Medium PSM FDR Threshold. Default=0.05')
    parser.add_argument('--prefix', '-o', dest='out_prefix', default='', help='Set output files prefix. Default=')
    parser.add_argument('--ms3', dest='ms3', default=False, action='store_true', help='Flag to indicate if the data is MS3 TMT data, which requires special handling for quantification. Default=false')
    parser.add_argument('--mztabout', '-mz', dest='mztab_out', default=False, action='store_true', help='Flag to indicate whether to output a new mzTab file with the updated annotations. Default=false')
    parser.add_argument('--log_level', '-l', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], help='Logging level (Default: INFO)')
    return parser.parse_args()

# =============================================================================
# Utility helpers
# =============================================================================

# Function to clean and standardize sample names by removing non-alphanumeric characters and converting to uppercase
def clean_sample_name(name: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "", name.upper())

# Function to strip version numbers from identifiers (e.g., ENSG00000354587.3 -> ENSG00000354587)
def strip_version(identifier: str) -> str:
    return identifier.split(".", 1)[0]

# =============================================================================
# Load Sample and Reference Information
# =============================================================================

#Function to read sample information from a TSV file and return a dictionary mapping TMT labels to sample IDs for each ExperimentID.
# Example:
# ExperimentID	TMTLabel	SampleID
# Exp1	126	SampleA
def load_sample_info(path: Path) -> Dict[str, Dict[str, str]]:
    sample_info: Dict[str, Dict[str, str]] = defaultdict(dict)

    #Check if file exists
    if not path.is_file():
        LOG.warning(f"Sample information file {path} not found. Sample information will not be loaded.")
        return sample_info

    with path.open("r") as sinfo:
        reader = csv.DictReader(sinfo, delimiter="\t")
        if not {"ExperimentID", "TMTLabel", "SampleID"}.issubset(reader.fieldnames):
            LOG.warning(f"Sample information file {path} is missing required columns. Expected columns: ExperimentID, TMTLabel, SampleID. Found columns: {reader.fieldnames}. Sample information will not be loaded.")
            return sample_info
        for row in reader:
            exp_id = row["ExperimentID"]
            tmt_label = row["TMTLabel"]
            sample_id = clean_sample_name(row["SampleID"])
            sample_info[exp_id][tmt_label] = sample_id
    return sample_info

# Function to read FASTA files and yield header-sequence pairs as an iterator
# Make I/L equivalent
def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header: Optional[str] = None
    sequence: List[str] = []
    with path.open("r") as FASTA:
        for line in FASTA:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line
                sequence = []
            else:
                sequence.append(line.replace('I', 'L'))
        if header is not None:
            yield header, "".join(sequence)

#Function to parse reference FASTA and store the sequences and also the gene to transcript mappings
def getReferenceProteins(ref_path: Path) -> Tuple [ Dict[str, Dict[str, str]], Dict[str, Dict[str, set[str]]], Dict[str, Dict[str, set[str]]] ]:
    reference_proteins: Dict[str, Dict[str, str]] = {}
    gene_transcript_map: Dict[str, Dict[str, set[str]]] = {}
    transcript_gene_map: Dict[str, Dict[str, set[str]]] = {}
    
    for header, sequence in read_fasta(ref_path):

        protein_id = header.split()[0][1:]  # Using the first word after '>' as the protein_id
        gene_id = ""
        gene_name = ""
        transcript_id = ""
        
        if header.startswith(">ENSP"):
            parts = header.split("|")
            transcript_id = strip_version(parts[1])
            gene_id = strip_version(parts[2])
            gene_name = parts[6]
            
        reference_proteins[protein_id] = {
            "gene_id": gene_id,
            "transcript_id": transcript_id,
            "sequence": sequence
        }

        if gene_id.startswith("ENSG"):
            # Ensure the gene_id exists in the gene_transcript_map
            if gene_id not in gene_transcript_map:
                gene_transcript_map[gene_id] = {"transcript_ids": set(), "gene_name": set()}

            # Add protein_id and gene_name to dict of sets in gene_transcript_map
            gene_transcript_map[gene_id]["transcript_ids"].add(transcript_id)
            gene_transcript_map[gene_id]["gene_name"].add(gene_name)

        if transcript_id.startswith("ENST"):
            transcript_gene_map[transcript_id] = {"gene_ids": set(), "protein_ids": set()}
            transcript_gene_map[transcript_id]["gene_ids"].add(gene_id)
            transcript_gene_map[transcript_id]["protein_ids"].add(protein_id)
    
    return reference_proteins, gene_transcript_map, transcript_gene_map 


# =============================================================================
# Results Parsing
# =============================================================================

#Function to parse mzTab files and extract PSM and Peptide information, while also updating the gene mapping and additional headers with search information.
def parse_mztab(file: str, psms: Dict[str, Dict[str, Dict[str, Dict[str, PSMInfo]]]], peptides: Dict[str, PeptideInfo], heads: Dict[str, set[str]], metadata: MzTabMetaInfo, ms3: bool, psm_score_threshold: float, psm_fdr_threshold: float) -> None:

    #Strip path from filename
    fname = os.path.basename(file)

    #Parse search information from mzTab filename
    m = re.search('(V[^_]+)_(.*)_DBx(.*)_TD[^_]*_([a-z]+).*.mzTab', fname)
    if m:

        search_id = m.group(1)
        experiment = m.group(2)
        database = m.group(3)
        engine = m.group(4)

        heads["Databases"].add(database)
        heads["SearchEngines"].add(engine)
        heads["Experiments"].add(experiment)
        heads["Searchids"].add(search_id)
        header_index: Dict[str, int] = {}

        LOG.info(f"Parsing mzTab file: {fname} with SearchID: {search_id}, Experiment: {experiment}, Database: {database}, Engine: {engine}")

        modlookup: Dict[str, str] = {}  # Dictionary to map modification ids to their names

        with open(file, "r") as f:
            for line in f:

                #Parse MTD line for recreating new mztab files in future
                if line.startswith("MTD"):
                    mtd_fields = line.strip().split("\t")
                    if len(mtd_fields) >= 3:
                        key = mtd_fields[1]
                        value = mtd_fields[2]

                        #MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
                        #MTD	fixed_mod[1]-site	C
                        #MTD	fixed_mod[1]-position	Anywhere

                        match_mod = re.match(r"(variable_mod|fixed_mod)\[(\d+)\](?:-(site|pos))?", key)
                        if match_mod:
                            
                            mod_type = match_mod.group(1)
                            mod_number = match_mod.group(2)
                            mod_attr = match_mod.group(3)

                            modkey = f"{mod_type}[{mod_number}]"

                            LOG.debug(f"Parsing modification metadata from mzTab: {modkey} with attribute '{mod_attr}' and value '{value}'")

                            if modkey not in metadata.variable_modifications and mod_type == "variable_mod":
                                metadata.variable_modifications[modkey] = ModificationInfo()

                            if modkey not in metadata.fixed_modifications and mod_type == "fixed_mod":
                                metadata.fixed_modifications[modkey] = ModificationInfo()

                            if mod_attr is None:
                                if mod_type == "variable_mod":
                                    metadata.variable_modifications[modkey].info = value
                                else:
                                    metadata.fixed_modifications[modkey].info = value
                                #Extract modsouce, modid and modname from value : example =  [UNIMOD, UNIMOD:35, Oxidation, ]
                                mod_parts = value.strip("[]").split(",")
                                if len(mod_parts) >= 3:
                                    modlookup[mod_parts[1].strip()] = mod_parts[2].strip()  # Map modid to modname for later use
                                    LOG.info(f"Parsed modification from mzTab metadata: {modkey} with source '{mod_parts[0].strip()}', id '{mod_parts[1].strip()}', name '{mod_parts[2].strip()}'")
                                else:
                                    LOG.warning(f"Unexpected modification format in mzTab metadata for {modkey}: '{value}' : '{key}'. Expected format: [source, modid, modname, ...]. This modification may not be correctly mapped in PSM parsing.")

                            elif mod_attr == "site":
                                if mod_type == "variable_mod":
                                    metadata.variable_modifications[modkey].site = value
                                else:
                                    metadata.fixed_modifications[modkey].site = value
                            
                            elif mod_attr == "pos":
                                if mod_type == "variable_mod":
                                    metadata.variable_modifications[modkey].pos = value
                                else:
                                    metadata.fixed_modifications[modkey].pos = value

                        elif "quantification_method" in key:
                            metadata.quantification_method = value

                        elif "mzml_file" in key:
                            metadata.mzml_file = value.replace("file://", "").replace('_filtered.mzML', '.mzMl')  # Remove file:// prefix if present and standardize filename if needed

                #PSH / PEH lines contain the headers for the PSM and PEP sections, we use this to create a mapping of header names to their column indices for later parsing of PSM / PEP lines
                elif line.startswith("PSH") or line.startswith("PEH"):
                    headers = line.strip().split("\t")
                    header_index = {h: i for i, h in enumerate(headers)}

                #Parse PEP lines
                #PEP lines contain the actual peptide information, we use the header_index mapping to extract relevant fields and filter based on FDR thresholds, storing the information in the peptides dictionary and updating the corresponding PSM information.
                elif line.startswith("PEP") and not ms3:  # If this is MS3 data we skip the peptide lines as the quantification is incorrect 
                    fields = line.strip().split("\t")
                    
                    sequence = fields[header_index["sequence"]]
                    best_psm_id = experiment + "_" + fields[header_index['spectra_ref']].split('=', 1)[-1]
                    
                    if sequence not in peptides:
                        peptides[sequence] = PeptideInfo()

                    #get abundance value for each TMT channel
                    #NOTE: OpenMS mzTab uses only the TMT intensity for the bestmatching PSM for each peptide.
                    #NOTE: If this is MS3 data the abundance is incorrect hence peptide line is skipped earlier.
                    for key, index in header_index.items():
                        m = re.match(r"peptide_abundance_study_variable\[(\d+)\]", key)
                        if m:
                            channel = int(m.group(1))
                            if fields[index] != "null" and fields[index] != "":
                                try:
                                    abundance = float(fields[index]) if fields[index] else 0
                                    if best_psm_id not in peptides[sequence].quantification:
                                        peptides[sequence].quantification[best_psm_id] = {}
                                    peptides[sequence].quantification[best_psm_id][channel] = abundance
                                    heads["TMTChannels"].add(channel)
                                except ValueError:
                                    LOG.warning(f"Invalid abundance value in file {fname} for peptide {sequence} in channel {channel}: '{fields[index]}', ignoring this value.")
                            
                #Parse PSM lines
                #PSM lines contain the actual PSM information, we use the header_index mapping to extract relevant fields and filter based on score and FDR thresholds, storing the information in the psms dictionary and updating the corresponding peptide information.
                elif line.startswith("PSM"):
                    fields = line.strip().split("\t")

                    sequence = fields[header_index["sequence"]]
                    spectrum_id = experiment + "_" + fields[header_index['spectra_ref']].split('=', 1)[-1]
                    
                    score_str = "99"
                    fdr_str = "99"

                    #Epifany overwrites mzTab scores incorrectly - Ignore scores for Epifany results as they are not comparable to the scores from other search engines and would cause incorrect filtering and annotation.
                    if "epifany" not in fname.lower():  

                        #check if search_engine_score[1] column is present, if not use 'search_engine_score[1]_ms_run[1]'
                        if "search_engine_score[1]" in header_index:
                            score_str = fields[header_index["search_engine_score[1]"]]
                        elif "search_engine_score[1]_ms_run[1]" in header_index:
                            score_str = fields[header_index["search_engine_score[1]_ms_run[1]"]]
                    
                        fdr_str = fields[header_index.get("PSM-level FDR", "")] if header_index.get("PSM-level FDR") is not None else "99"

                    try:
                        score = float(score_str) if score_str else 99
                        fdr = float(fdr_str) if fdr_str else 99
                    except ValueError:
                        LOG.warning(f"Invalid score or FDR value in file {fname} for PSM {spectrum_id}: score='{score_str}', fdr='{fdr_str}'")
                        score = 99
                        fdr = 99

                    mods = set(fields[header_index.get("modifications", "")].split(",")) if header_index.get("modifications") is not None else set()
                    modifications = set()
                    for mod in mods:
                        if mod == "null" or mod == "":
                            continue
                        mod_parts = mod.split('-')
                        if mod_parts[1] in modlookup:
                            mod = mod_parts[0] + "-" + modlookup[mod_parts[1]]
                            
                        modifications.add(mod)

                    peptidoform = fields[header_index.get("opt_global_cv_MS:1000889_peptidoform_sequence", "")] if header_index.get("opt_global_cv_MS:1000889_peptidoform_sequence") is not None else sequence

                    if spectrum_id not in psms:
                        psms[spectrum_id] = {}

                    if sequence not in psms[spectrum_id]:
                        psms[spectrum_id][sequence] = {}

                    if database not in psms[spectrum_id][sequence]:
                        psms[spectrum_id][sequence][database] = {}

                    if engine not in psms[spectrum_id][sequence][database]:
                        psms[spectrum_id][sequence][database][engine] = PSMInfo(
                            score=score,
                            fdr=fdr,
                            peptidoform=peptidoform,
                            modifications=modifications,
                            rt=float(fields[header_index.get("retention_time", "0")]) if header_index.get("retention_time") is not None else None,
                            mz=float(fields[header_index.get("exp_mass_to_charge", "0")]) if header_index.get("exp_mass_to_charge") is not None else None,
                            charge=int(fields[header_index.get("charge", "0")]) if header_index.get("charge") is not None else None,
                        ) 
                    else:
                        if score < psms[spectrum_id][sequence][database][engine].score:
                            psms[spectrum_id][sequence][database][engine].score = score
                        
                        if fdr < psms[spectrum_id][sequence][database][engine].fdr:
                            psms[spectrum_id][sequence][database][engine].fdr = fdr
                    

                    if (score <= psm_score_threshold and (fdr <= psm_fdr_threshold or fdr == 99)) and sequence != "":

                        #Update peptide information with this PSM
                        if sequence not in peptides:
                            peptides[sequence] = PeptideInfo()

                        peptides[sequence].psms.add(spectrum_id)
                        peptides[sequence].modifications.update(modifications)
                        peptides[sequence].preaa.update(fields[header_index.get("pre", "")].split(",") if header_index.get("pre") is not None else [])
                        peptides[sequence].postaa.update(fields[header_index.get("post", "")].split(",") if header_index.get("post") is not None else [])
                        peptides[sequence].peptidoforms.add(peptidoform)

                        if engine not in peptides[sequence].search_engines:
                            peptides[sequence].search_engines[engine] = score
                        elif score < peptides[sequence].search_engines[engine]:
                            peptides[sequence].search_engines[engine] = score

                        if experiment not in peptides[sequence].experiments:
                            peptides[sequence].experiments[experiment] = score
                        elif score < peptides[sequence].experiments[experiment]:
                            peptides[sequence].experiments[experiment] = score 
                        
                        if fdr < peptides[sequence].psm_fdr:
                            peptides[sequence].psm_fdr = fdr
                        
                        if score < peptides[sequence].psm_score:
                            peptides[sequence].psm_score = score
                            peptides[sequence].best_psm = spectrum_id

                        if database in peptides[sequence].database_score:
                            if score < peptides[sequence].database_score[database]:
                                peptides[sequence].database_score[database] = score
                        else:
                            peptides[sequence].database_score[database] = score

                        protein_list = fields[header_index.get("protein_accs", "")].split(",") if header_index.get("protein_accs") is not None else []
                        protein_list += fields[header_index.get("opt_global_cv_MS:1000889_protein_accessions", "")].split(",") if header_index.get("opt_global_cv_MS:1000889_protein_accessions") is not None else []
                        protein_list += fields[header_index.get("accession", "")].split(",") if header_index.get("accession") is not None else []
                        protein_list = [p.strip() for p in protein_list if p.strip()]  # Remove empty strings and strip whitespace

                        for protein in protein_list:
                            peptides[sequence].protein_ids[database].add(protein)

    else:
        LOG.warning(f"Filename does not match expected pattern for mzTab parsing: {file}")
        return

#Function to parse Percolator tab files and update PSM and Peptide information with scores and FDRs, filtering based on specified thresholds.
def parse_percolator_tab(file: str, psms: Dict[str, Dict[str, Dict[str, Dict[str, PSMInfo]]]], peptides: Dict[str, PeptideInfo], heads: Dict[str, set[str]], psm_score_threshold: float, psm_fdr_threshold: float) -> None:

    #Strip path from filename
    fname = os.path.basename(file)

    m = re.search('(V[^_]+)_(.*)_DBx(.*)_TD[^_]*_([a-z]+)_.*percolator.tab', fname)
    if m:

        search_id = m.group(1)
        experiment = m.group(2)
        database = m.group(3)
        engine = m.group(4)

        heads["Databases"].add(database)
        heads["SearchEngines"].add(engine)
        heads["Experiments"].add(experiment)
        heads["Searchids"].add(search_id)

        LOG.info(f"Parsing Percolator tab file: {fname} with SearchID: {search_id}, Experiment: {experiment}, Database: {database}, Engine: {engine}")

        with open(file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t", restkey="extra_proteins")
            for row in reader:
                spectrum_id = experiment + "_" + row["PSMId"].split('=', 1)[-1]
                #spectrum_id = experiment + "_" + fields[header_index['spectra_ref']].split('=', 1)[-1]
                pep_str = row.get("posterior_error_prob")
                fdr_str = row.get("q-value")
                peptide_str = row.get("peptide")

                LOG.debug(f"Processing Percolator PSM {spectrum_id} with peptide '{peptide_str}', score='{pep_str}', FDR='{fdr_str}'")

                p = re.search(r"(\D+)\.(.*)\.(\D+)", peptide_str)
                if p:
                    preaa = p.group(1)
                    pepform = p.group(2)
                    postaa = p.group(3)
                    peptide = re.sub('\[[^A-Z]+\]', '', pepform.replace('I', 'L').replace('n', ''))  # Remove modifications to get the base peptide sequence

                    try:
                        score = float(pep_str) 
                        fdr = float(fdr_str) 
                    except ValueError:
                        LOG.warning(f"Invalid score or FDR value in file {fname} for PSM {spectrum_id}: score='{pep_str}', fdr='{fdr_str}'")
                        score = 99
                        fdr = 99

                    if spectrum_id not in psms:
                        psms[spectrum_id] = {}

                    if peptide not in psms[spectrum_id]:
                        psms[spectrum_id][peptide] = {}

                    if database not in psms[spectrum_id][peptide]:
                        psms[spectrum_id][peptide][database] = {}

                    if engine not in psms[spectrum_id][peptide][database]:
                        psms[spectrum_id][peptide][database][engine] = PSMInfo(
                                score=score,
                                fdr=fdr,
                                peptidoform=pepform,
                        )
                    else:
                        if score < psms[spectrum_id][peptide][database][engine].score:
                            psms[spectrum_id][peptide][database][engine].score = score
                        
                        if fdr < psms[spectrum_id][peptide][database][engine].fdr:
                            psms[spectrum_id][peptide][database][engine].fdr = fdr

                    if (score <= psm_score_threshold or fdr <= psm_fdr_threshold) and peptide != "":

                        if peptide not in peptides:
                            LOG.debug(f"Creating new PeptideInfo for peptide '{peptide}' with initial PSM score {score} and FDR {fdr}")
                            peptides[peptide] = PeptideInfo(
                                preaa={preaa},
                                postaa={postaa},
                                peptidoforms={pepform},
                                psm_score=score,
                                psm_fdr=fdr,
                                best_psm=spectrum_id,
                                search_engines={engine: score},
                                experiments={experiment: score},
                                database_score={database: score}
                            )
                        else:
                            
                            peptides[peptide].preaa.add(preaa)
                            peptides[peptide].postaa.add(postaa)    

                            if engine not in peptides[peptide].search_engines:
                                peptides[peptide].search_engines[engine] = score
                            elif score < peptides[peptide].search_engines[engine]:
                                peptides[peptide].search_engines[engine] = score

                            if experiment not in peptides[peptide].experiments:
                                peptides[peptide].experiments[experiment] = score
                            elif score < peptides[peptide].experiments[experiment]:
                                peptides[peptide].experiments[experiment] = score 

                            peptides[peptide].peptidoforms.add(pepform)

                            if score < peptides[peptide].psm_score:
                                peptides[peptide].psm_score = score
                                peptides[peptide].best_psm = spectrum_id

                            if fdr < peptides[peptide].psm_fdr:
                                peptides[peptide].psm_fdr = fdr
                                LOG.debug(f"Updated PSM FDR for peptide '{peptide}' to {fdr} based on PSM {spectrum_id}")

                            if database not in peptides[peptide].database_score or score < peptides[peptide].database_score[database]:
                                peptides[peptide].database_score[database] = score

                        #Get protein accessions from "protein_accessions" column and also from "extra_proteins" if present, and update the peptide information with these proteins
                        protein_list = row.get("proteinIds", "").split(",") if row.get("proteinIds") is not None else []
                        extra_proteins = row.get("extra_proteins", "") if row.get("extra_proteins") is not None else []
                        for protein in protein_list + extra_proteins:
                            peptides[peptide].protein_ids[database].add(protein)

                else:
                    LOG.warning(f"Peptide string does not match expected pattern for percolator parsing in file {fname}: '{peptide_str}'. Expected format: 'preAA.peptidoform.postAA' where peptidoform may contain modifications in brackets. This PSM will be skipped.")

    else:
        LOG.warning(f"Filename does not match expected pattern for percolator parsing: {fname}")
        return


# =============================================================================
# Peptide and PSM Annotation
# =============================================================================

#Function to assess conflicting PSM assignments for each spectrum and annotate the PSM and peptide information with the best assignment, delta scores, and conflicting assignments.
# For each spectrum, we check if there are multiple peptide assignments with similar scores (e.g., within a certain threshold) and annotate the PSM with the best peptide assignment, the next best assignment, 
# the delta score between them, and a flag indicating whether there are conflicting assignments with similar scores. 
# This provides insight into the confidence of the peptide assignment for each PSM.
def psm_conflict_assessment(psms: Dict[str, Dict[str, Dict[str, Dict[str, PSMInfo]]]], peptides: Dict[str, PeptideInfo], threshold: float) -> None:

#1. Conflicting PSM assignments and delta score between them - For each PSM, we check if there are multiple peptide assignments with similar scores (e.g., within a certain threshold) 
    # and annotate the PSM with the best peptide assignment, the next best assignment, the delta score between them, 
    # and a flag indicating whether there are conflicting assignments with similar scores. 
    # This provides insight into the confidence of the peptide assignment for each PSM.
    for spectrum_id, peptide_dict in psms.items():

        best_score = float('100.00')  
        second_best_score = float('100.00') 
        best_peptide = None
        conflicting_assignments: Dict[str, float] = {}

        for peptide, database_dict in peptide_dict.items():
            for database, engine_dict in database_dict.items():
                for engine, psm_info in engine_dict.items():

                    #Record conflicting assignments and best score for each
                    if peptide not in conflicting_assignments:
                        conflicting_assignments[peptide] = psm_info.score
                    elif psm_info.score < conflicting_assignments[peptide]:
                        conflicting_assignments[peptide] = psm_info.score 

                    #Set best and second best scores and peptides for this spectrum based on the PSM scores, 
                    if psm_info.score < best_score:
                        if peptide != best_peptide and best_peptide is not None:
                            second_best_score = best_score
                        best_score = psm_info.score
                        best_peptide = peptide
                    elif psm_info.score < second_best_score and psm_info.score > best_score:
                        if peptide != best_peptide:
                            second_best_score = psm_info.score
            
        #Calculate delta score between best and second best assignment.
        #If no second best assignment is present we set delta score to the difference between the best score and a default threshold 
        # (e.g., 2 * score threshold) to provide some insight into the confidence of the assignment even when there is only one assignment.
        delta_score = second_best_score - best_score

        if best_peptide is not None and best_peptide in peptides:
            #If delta score is greater than 0, it means the best peptide is lower than the 2*score threshold and there are no other assignments with same score.
            #so we set the top_rank flag for the best peptide to indicate that it has a top ranked spectrum assignment in atleast 1 spectrum.
            if delta_score > 0:
                peptides[best_peptide].top_rank = True

            #Store the delta score for the best peptide assignment across all its PSMs.
            if delta_score > peptides[best_peptide].delta_score:
                peptides[best_peptide].delta_score = delta_score

        #Reloop through the PSMs for this spectrum and annotate the PEPTIDEInfo and PSMInfo for each peptide assignment 
        # with the delta score to the best assignment and any conflicting assignments with similar scores,
        for peptide, database_dict in peptide_dict.items():

            #Skip peptides that did not score high enough to be in peptide list
            if peptide not in peptides:
                continue

            #Annotate the peptide information with any conflicting assignments only store smallest delta.
            for pep, score in conflicting_assignments.items():
                if pep != peptide:
                    if pep not in peptides[peptide].conflicts or score < peptides[peptide].conflicts[pep]:
                        peptides[peptide].conflicts[pep] = score


            for database, engine_dict in database_dict.items():
                for engine, psm_info in engine_dict.items():

                    #Calculate Delta Score for this PSM
                    #either to the best peptide or second best if this is the same peptide as the best
                    #Can be negative if this PSm matches the same peptide as the best but has a worse score than second best peptide 
                    
                    if peptide != best_peptide:
                        psmdelta = psm_info.score - best_score
                    else:
                        psmdelta = second_best_score - psm_info.score

                    psms[spectrum_id][peptide][database][engine].delta_score = psmdelta

                    #Set best_peptide flag for this PSM if it has the best score for this spectrum
                    if psms[spectrum_id][peptide][database][engine].score == best_score:
                        psms[spectrum_id][peptide][database][engine].best_peptide = True

                    #Store conflicting peptides assignments for this PSM
                    for pep, score in conflicting_assignments.items():
                        if pep != peptide:
                            psms[spectrum_id][peptide][database][engine].conflicts[pep] = score


#Function to calculate q-values for a list of target and decoy scores using the target-decoy approach, enforcing monotonicity of q-values with respect to score thresholds.
def calculate_qvalues(target_scores: List[float], decoy_scores: List[float]) -> List[Tuple[float, float]]:
    qval_lookup: Dict[float, float] = {}
    current_qval = 1.0

    test_fdr = len(decoy_scores) / len(target_scores) if len(target_scores) > 0 else 1.0
    LOG.info(f"Calculating q-values using target-decoy approach with {len(target_scores)} target scores and {len(decoy_scores)} decoy scores. Initial FDR estimate: {test_fdr:.4f}")

    #loop concatenated and sorted target and decoy scores in descending order.
    for score in sorted(set(target_scores + decoy_scores), reverse=True):
        num_targets = sum(s <= score for s in target_scores)
        num_decoys = sum(s <= score for s in decoy_scores)

        raw_fdr = num_decoys / num_targets if num_targets > 0 else 1.0

        # Enforce monotonicity: q-value = minimum FDR seen at this score or higher
        current_qval = min(current_qval, raw_fdr)
        qval_lookup[score] = current_qval

    return qval_lookup

#Basic function to calculate the minimum SAAV distance between a given peptide sequence and all reference peptides in the reference database, returning the minimum distance and the best matching reference peptide sequence.
#Should not be applied to large peptide sets as not very efficient
#Limited to a maximum distance threshold (e.g., 1 or 2) to avoid unnecessary calculations for peptides that are very different from any reference peptide, which is often the case for non-canonical peptides.
def calculate_min_saav_distance(sequence: str, reference_proteins: Dict[str, Dict[str, str]], max_distance: int) -> Tuple[int, str]:
    min_distance = max_distance + 1
    best_ref_peptide = None

    peptide_len = len(sequence)
    peptide_array = np.array(list(sequence)) #numpy array of the peptide sequence for faster comparison

    for ref_protein in reference_proteins.values():
        
        protein_len = len(ref_protein["sequence"])
        protein_array = np.array(list(ref_protein["sequence"])) #numpy array of the reference protein sequence for faster comparison

        for i in range(protein_len - peptide_len + 1):

            window = protein_array[i:i+peptide_len]
            substitution_count = np.sum(peptide_array != window)

            if substitution_count < min_distance:
                min_distance = substitution_count
                best_ref_peptide = "".join(window)
                if min_distance == 1:
                    return min_distance, best_ref_peptide
    
    return min_distance, best_ref_peptide

#Function to compute k-mers for all reference peptides in the reference database and create an index mapping each k-mer to the set of 
#reference protein IDs that contain that k-mer.
def compute_kmers(reference_proteins: Dict[str, Dict[str, str]], k: int) -> Dict[str, Dict[str, Set[int]]]:
    kmers: Dict[str, Dict[str, Set[int]]] = defaultdict()
    for id, protein in reference_proteins.items():
        sequence = protein["sequence"]
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if kmer not in kmers:
                kmers[kmer] = defaultdict(set)

            if id not in kmers[kmer]:
                kmers[kmer][id] = set()

            #Store the kmer sequence, the matching protein and the position of the kmer
            kmers[kmer][id].add(i)

            kmers[kmer]
    return kmers

#Function to calculate the minimum SAAV distance between a given peptide sequence and reference peptides in the reference database 
#using a k-mer based approach to first identify candidate reference peptides that share k-mers with the query peptide, and then 
#calculating the SAAV distance only for those candidates. This can significantly speed up the calculation for larger peptide sets 
#and reference databases, while still providing accurate distance calculations for peptides that are similar to reference peptides.
def calculate_min_saav_distance_kmer_accelerated(sequence: str, reference_proteins: Dict[str, Dict[str, str]], k: int, kmer_index: Dict[str, Set[str]], max_distance: int) -> Tuple[int, str]:

    peptide_len = len(sequence)
    #If the peptide is shorter than k * max_distance, we cannot reliably use the k-mer approach to find candidate reference peptides, as even a single substitution could affect multiple k-mers.
    if peptide_len < k * max_distance:
        return calculate_min_saav_distance(sequence, reference_proteins, max_distance)
    
    min_distance = max_distance + 1
    best_ref_peptide = None

    candidate_proteins = set()
    for i in range(peptide_len - k + 1):
        kmer = sequence[i:i+k]
        if kmer in kmer_index:
            for id in kmer_index[kmer]:
                for pos in kmer_index[kmer][id]:
                    start = pos - i  # Calculate the potential start position of the peptide in the reference protein based on the position of the k-mer
                    if start >= 0:
                        candidate_proteins.add((id, start))

    if candidate_proteins:
        for protein_id, start_pos in candidate_proteins:
            ref_protein = reference_proteins[protein_id]
            protein_len = len(ref_protein["sequence"])

            if start_pos + peptide_len <= protein_len:

                protein_sequence = ref_protein["sequence"]
                
                window = protein_sequence[start_pos:start_pos+peptide_len]
                substitution_count = 0
                for a, b in zip(sequence, window): 
                    if a != b:
                        substitution_count += 1
                        if substitution_count > max_distance:
                            break

                if substitution_count < min_distance and substitution_count <= max_distance:
                    min_distance = substitution_count
                    best_ref_peptide = window
                    if min_distance == 1:
                        return min_distance, best_ref_peptide
    
    return min_distance, best_ref_peptide

   
#Function to annotate peptide data with:
# 1. Conflicting PSM assignments and delta score between them
# 2. Decoy status based on whether the associated proteins are decoys
# 3. Type based on whether the associated proteins are in the reference database
# 4. Gene and transcript mappings based on the associated proteins and the reference database
# 5. Peptide level FDR based on PSM scores
# 6. Database level FDR based on the best PSM score for each database
# 7. Group / Class FDR based on peptide type
# 8. Minimum SAAV distance to reference peptides for variant and non-canonical peptides

def peptide_annotation(psms: Dict[str, PSMInfo], peptides: Dict[str, PeptideInfo], proteins: Dict[str, ProteinInfo], reference_proteins: Dict[str, Dict[str, str]], heads: Dict[str, set[str]], gene_map: Dict[str, Dict[str, set[str]]], transcript_map: Dict[str, Dict[str, set[str]]], score_threshold: float, fdr_threshold: float) -> None:

    
    #1. Conflicting PSM assignments and delta score between them - For peptides that are associated with multiple PSMs, 
    # identify cases where there are conflicting assignments (e.g., one PSM indicates a canonical peptide while another 
    # indicates a non-canonical peptide) and calculate the delta score between the best PSM and the next best PSM to provide insight 
    # into the confidence of the assignment.  
    LOG.info("  ...Assessing conflicting PSM assignments for each spectrum.")              
    psm_conflict_assessment(psms, peptides, score_threshold)   
    
    target_score_list: List[float] = []
    decoy_score_list: List[float] = []
    target_database_score_list: Dict[str, List[float]] = {db: [] for db in heads["Databases"]}
    decoy_database_score_list: Dict[str, List[float]] = {db: [] for db in heads["Databases"]}

    #2. Decoy status | 3. type annotation | 4. Gene and transcript mapping
    LOG.info("  ...Annotating peptides and proteins with decoy status, type, and gene/transcript mappings.")
    for sequence, info in peptides.items():

        protein_set_no_decoys = set()
        decoy_proteins = set()

        for database, protein_ids in info.protein_ids.items():
            for protein in protein_ids:

                if protein not in proteins:
                    proteins[protein] = ProteinInfo()

                    # Determine decoy status based on protein
                    if not re.match(r"DECOY|REVERSED|RANDOM", protein, re.IGNORECASE):

                        proteins[protein].decoy = False
                        protein_set_no_decoys.add(protein)

                        # Determine protein type
                        if protein in reference_proteins:
                            ref_protein = reference_proteins[protein]
                            proteins[protein].type = "canonical"
                            proteins[protein].gene_ids.add(ref_protein["gene_id"])
                            proteins[protein].transcript_ids.add(ref_protein["transcript_id"])
                            proteins[protein].length = len(ref_protein["sequence"])
                        else:
                            gene_match = re.match(r"(ENSG\d+)", protein)
                            if gene_match:
                                gene_id = gene_match.group(1)
                                proteins[protein].gene_ids.add(gene_id)
                                if gene_id in gene_map:
                                    for transcript in gene_map[gene_id]["transcript_ids"]:
                                        for prot_id in transcript_map[transcript]["protein_ids"]:
                                            if prot_id in reference_proteins:
                                                ref_protein = reference_proteins[prot_id]
                                                if sequence in ref_protein["sequence"]:
                                                    proteins[protein].type = "canonical"
                                                    proteins[protein].transcript_ids.add(transcript)
                                                    proteins[protein].length = len(ref_protein["sequence"])

                            transcript_match = re.match(r"(ENST\d+)", protein)
                            if transcript_match:
                                transcript_id = transcript_match.group(1)
                                if transcript_id not in proteins[protein].transcript_ids:
                                    proteins[protein].transcript_ids.add(transcript_id)
                                    if transcript_id in transcript_map:
                                        proteins[protein].gene_ids.update(transcript_map[transcript_id]["gene_ids"])
                                        for prot_id in transcript_map[transcript_id]["protein_ids"]:
                                            if prot_id in reference_proteins:
                                                ref_protein = reference_proteins[prot_id]
                                                if sequence in ref_protein["sequence"]:
                                                    proteins[protein].type = "canonical"
                                                    proteins[protein].length = len(ref_protein["sequence"])

                            if proteins[protein].type != "canonical":
                                if re.match(r"^CON|sp\|cRAP\d+", protein, re.IGNORECASE) and proteins[protein].type == "unknown":
                                    proteins[protein].type = "contaminant"
                                elif re.match(r"^VAR", protein, re.IGNORECASE):
                                    proteins[protein].type = "variant"
                                elif proteins[protein].type == "unknown":
                                    proteins[protein].type = "non_canonical"
                    else:
                        decoy_proteins.add(protein)

                # Annotate proteins with peptide information
                proteins[protein].peptides.add(sequence)
                proteins[protein].databases.add(database)
                proteins[protein].search_engines.update(info.search_engines)
                proteins[protein].experiments.update(info.experiments)

                if not proteins[protein].decoy:
                    # Assign peptide type and update peptide with protein annotations if protein not a decoy
                    info.decoy = False
                    if info.type != "canonical":
                        if proteins[protein].type == "canonical":
                            info.type = "canonical"
                        elif proteins[protein].type == "variant":
                            info.type = "variant"
                        elif proteins[protein].type == "contaminant" and info.type == "unknown":
                            info.type = "contaminant"
                        elif proteins[protein].type == "non_canonical" and info.type != "variant":
                            info.type = "non_canonical"

                    info.genes.update(proteins[protein].gene_ids)
                    info.transcripts.update(proteins[protein].transcript_ids)


            if info.top_rank:
                if info.decoy:
                    decoy_database_score_list[database].append(info.database_score[database])
                else:
                    target_database_score_list[database].append(info.database_score[database])

        # If not a known canonical, conduct final check to make sure peptide not in reference proteins
        if info.type != "canonical" and not info.decoy:
            for ref_protein in reference_proteins.values():
                if sequence in ref_protein["sequence"]:
                    info.type = "canonical"
                    info.genes.add(ref_protein["gene_id"])
                    info.transcripts.add(ref_protein["transcript_id"])

        # If not a decoy peptide, remove all decoy proteins from protein list
        if not info.decoy:
            for database, protein_ids in info.protein_ids.items():
                target_proteins = set()
                for protein in protein_ids:
                    if protein not in decoy_proteins:
                        if not re.match(r"DECOY__\d+", protein, re.IGNORECASE):
                            target_proteins.add(protein)
                        else:
                            LOG.warning(f"Protein {protein} is in decoy list based on name but was not previously annotated as decoy. This protein will be treated as decoy for peptide {sequence}. Consider updating decoy annotation for this protein.")
                info.protein_ids[database] = target_proteins

                        
        # Update target and decoy score lists if this peptide is top spectrum match
        # Only include top match peptides as proteotypic for protein
        if info.top_rank:
            if info.decoy:
                decoy_score_list.append(info.psm_score)
                # Add peptides that only map to decoy proteins as proteotypic and update the protein score
                for protein in decoy_proteins:
                    proteins[protein].proteotypic_peptides.add(sequence)
                    proteins[protein].score += -10 * np.log10(info.psm_score) if info.psm_score > 0 else 200
            else:
                target_score_list.append(info.psm_score)
                # Add peptides that map to only 1 non-decoy protein as proteotypic and update the protein score
                if len(protein_set_no_decoys) == 1:
                    for protein in protein_set_no_decoys:
                        proteins[protein].proteotypic_peptides.add(sequence)
                        proteins[protein].score += -10 * np.log10(info.psm_score) if info.psm_score > 0 else 200
                    
    #5. Peptide level q-value annotations total and 6. database specific based on target and decoy score distributions using the best PSM score for each peptide
    LOG.info(f"  ...Calculating peptide level q-values based on target-decoy score distributions. Number of target peptides: {len(target_score_list)}, number of decoy peptides: {len(decoy_score_list)}.")
    qvalues = calculate_qvalues(target_score_list, decoy_score_list)
    database_qvalues: Dict[str, Dict[float, float]] = {}
    for database in heads["Databases"]:
        database_qvalues[database] = calculate_qvalues(target_database_score_list[database], decoy_database_score_list[database])

    for sequence, info in peptides.items():
        if info.psm_score in qvalues:
            peptides[sequence].peptide_fdr = qvalues[info.psm_score]
        else:
            peptides[sequence].peptide_fdr = 1.0

        for database in heads["Databases"]:
            if database in info.database_score and info.database_score[database] in database_qvalues[database]:
                peptides[sequence].database_fdr[database] = database_qvalues[database][info.database_score[database]]
            else:
                peptides[sequence].database_fdr[database] = 1.0
    
        
    #7. Group / Class FDR based of peptide type - Use groupFDR.py to calculate FDRs for each peptide type (canonical, non-canonical, decoy) based on the PSM scores, and annotate the peptides with these group FDRs.
    LOG.info("  ...Calculating group FDR based on peptide type.")
    import groupFDR
    groupFDR.compute_group_fdr(
        peptides,
        decoy_assignment="proportional",   # recommended for PEP scores
        score_higher_better=False,
        fdr_threshold=fdr_threshold,
        use_top_rank=True,
    ) 

    # 8. Minimum SAAV distance to reference peptides for variant and non-canonical peptides
    # For peptides annotated as variant or non-canonical, calculate the minimum SAAV distance to any reference peptide sequence in the reference database.
    #count non-canonical and variant peptides to get an estimate of how long this will take
    nc_count = sum(1 for info in peptides.values() if info.type in ["variant", "non_canonical"] and not info.decoy)
    if nc_count > 0 and nc_count < 1000:
        LOG.info(f"Calculating minimum SAAV distance for {nc_count} non-canonical and variant peptides. This may take some time...")
        
        LOG.info("  ...computing k-mer index for reference peptides to accelerate SAAV distance calculation.")
        kmer_length = 3
        kmer_index = compute_kmers(reference_proteins, k=kmer_length)

        LOG.info("  ...calculating minimum SAAV distance to reference peptides for non-canonical and variant peptides.")
        for sequence, info in peptides.items():
            if info.type in ["variant", "non_canonical"] and not info.decoy:
                #min_distance, best_ref_peptide = calculate_min_saav_distance(sequence, reference_proteins, max_distance=2)
                min_distance, best_ref_peptide = calculate_min_saav_distance_kmer_accelerated(sequence, reference_proteins, k=kmer_length, kmer_index=kmer_index, max_distance=2)
                info.minSAAVdistance = min_distance
                info.canonical_sequence = best_ref_peptide
    elif nc_count >= 1000:
        LOG.warning(f"Number of non-canonical and variant peptides is {nc_count}, skipping min SAAV disctance calculation to save time.")
        info.minSAAVdistance = -1

# ---------------------------------------------------------------------------
# Protein inference
# ---------------------------------------------------------------------------

def run_protein_inference(
    peptide_dict: dict[str, PeptideInfo],
    fdr_threshold: float = 0.01,
    gene_inference: bool = False,   
) -> dict[str, ProteinGroupInfo]:
    """
    Perform protein inference using the parsimony / Occam's Razor principle.

    1. Filter peptides by FDR threshold.
    2. Build protein→peptide and peptide→protein maps.
    3. Cluster all proteins by identical peptide sets (same-set clusters).
    4. Greedy set-cover on cluster representatives to get the parsimonious set.
    5. Flag non-parsimonious proteins whose peptides are a strict subset of a
        selected cluster's peptides.
    6. Compute unambiguous peptides, same-set members, and protein scores.

    Parameters
    ----------
    peptide_dict   : mapping of peptide_sequence -> PeptideInfo
    fdr_threshold  : only peptides with peptide_fdr <= this value are used

    Returns
    -------
    dict of master_protein_id -> ProteinInfo  (one entry per parsimonious group)
    """

    # ------------------------------------------------------------------
    # 1.  Filter peptides by FDR and Top Rank
    # ------------------------------------------------------------------
    passing: dict[str, PeptideInfo] = {
        seq: info
        for seq, info in peptide_dict.items()
        if info.peptide_fdr <= fdr_threshold and info.top_rank
    }

    if not passing:
        return {}

    LOG.info(f"Running protein inference on {len(passing)} peptides passing FDR threshold of {fdr_threshold}.")

    # ------------------------------------------------------------------
    # 2.  Build protein <-> peptide maps  (flatten source DBs)
    # ------------------------------------------------------------------

    LOG.info("  ...building protein-peptide maps for inference.")
    protein_to_peptides: dict[str, frozenset[str]] = {}
    peptide_to_proteins: dict[str, set[str]] = {}

    peptides_for_removal: set[str] = set()

    for seq, info in passing.items():
        proteins_for_seq: set[str] = set()

        if gene_inference:
            # For gene-level inference, we treat all proteins mapping to the same gene(s) as a single entity.
            # ignore any empy gene ids and any peptides that don't have gene annotations
            if info.genes:
                gene_match = False
                for gene in info.genes:
                    if gene != "":
                        proteins_for_seq.add(gene)
                        gene_match = True

                if gene_match:
                    peptide_to_proteins[seq] = proteins_for_seq
                    for prot in proteins_for_seq:
                        # accumulate mutable first, freeze later
                        protein_to_peptides.setdefault(prot, None)   # placeholder
                
        else:
            for db_proteins in info.protein_ids.values():
                proteins_for_seq.update(db_proteins)

            peptide_to_proteins[seq] = proteins_for_seq
            for prot in proteins_for_seq:
                # accumulate mutable first, freeze later
                protein_to_peptides.setdefault(prot, None)   # placeholder

    # build as sets first, then freeze
    prot_pep_mutable: dict[str, set[str]] = {}
    for seq, prots in peptide_to_proteins.items():
        for prot in prots:
            prot_pep_mutable.setdefault(prot, set()).add(seq)
    protein_to_peptides: dict[str, frozenset[str]] = {
        p: frozenset(s) for p, s in prot_pep_mutable.items()
    }

    LOG.info(f"  ...{len(protein_to_peptides)} proteins map to {len(peptide_to_proteins)} peptides for inference.")

    # ------------------------------------------------------------------
    # 3.  Cluster proteins by identical peptide sets  (same-set groups)
    # ------------------------------------------------------------------

    LOG.info("  ...clustering proteins by identical peptide sets.")

    pep_set_to_proteins: dict[frozenset[str], list[str]] = {}
    for prot, pep_set in protein_to_peptides.items():
        pep_set_to_proteins.setdefault(pep_set, []).append(prot)

    # Sort deterministically; first protein = cluster representative / master
    clusters: list[tuple[frozenset[str], list[str]]] = []
    for pep_set, prots in pep_set_to_proteins.items():
        clusters.append((pep_set, sorted(prots)))

    # representative -> full peptide set
    rep_to_pep_set: dict[str, frozenset[str]] = {
        prots[0]: pep_set for pep_set, prots in clusters
    }

    LOG.info(f"  ...{len(clusters)} clusters of proteins with identical peptide sets.")

    # ------------------------------------------------------------------
    # 4.  Greedy set-cover on cluster representatives  (parsimony)
    # ------------------------------------------------------------------

    LOG.info("  ...performing greedy set cover to select parsimonious representatives.")

    uncovered: set[str] = set(peptide_to_proteins.keys())
    selected_reps: set[str] = set()

    # work with representatives only
    reps = list(rep_to_pep_set.keys())

    while uncovered:
        best = max(
            reps,
            key=lambda r: (len(rep_to_pep_set[r] & uncovered), r),
        )
        selected_reps.add(best)
        uncovered -= rep_to_pep_set[best]

    LOG.info(f"  ...selected {len(selected_reps)} representative proteins covering all peptides.")

    # ------------------------------------------------------------------
    # 5.  Subset detection
    #     A non-selected representative is a subset if its peptides are a
    #     strict subset of any selected representative's peptides.
    # ------------------------------------------------------------------

    LOG.info("  ...detecting subset relationships among non-selected representatives.")

    selected_master_to_pep_set = {r: rep_to_pep_set[r] for r in selected_reps}

    # For each selected master, collect non-selected reps that are subsets,
    # plus ALL proteins in those non-selected clusters.
    master_subset_map: dict[str, list[str]] = {r: [] for r in selected_reps}

    for pep_set, prots in clusters:
        rep = prots[0]
        if rep in selected_reps:
            continue
        for master, master_peps in selected_master_to_pep_set.items():
            if pep_set < master_peps:           # strict subset
                master_subset_map[master].extend(prots)
                break                           # assign to first matching master

    LOG.info(f"  ...identified {sum(len(s) for s in master_subset_map.values())} subset proteins across {len(selected_reps)} selected representatives.")

    # ------------------------------------------------------------------
    # 6.  Unambiguous peptides
    #     A peptide is unambiguous for a master if *every* protein it maps to
    #     belongs to the same selected cluster.
    # ------------------------------------------------------------------
    # Build protein -> master lookup (selected clusters only)

    LOG.info("  ...determining unambiguous peptides for each selected representative.")

    protein_to_master: dict[str, str] = {}
    for pep_set, prots in clusters:
        master = prots[0]
        if master in selected_reps:
            for p in prots:
                protein_to_master[p] = master

    def _unambiguous(master: str, pep_set: frozenset[str]) -> list[str]:
        cluster_members = {p for p, m in protein_to_master.items() if m == master}
        return sorted(
            seq for seq in pep_set
            if peptide_to_proteins[seq].issubset(cluster_members)
        )

    LOG.info(f"  ...unambiguous peptides determined for each representative.")

    # ------------------------------------------------------------------
    # 7.  Assemble ProteinInfo objects
    # ------------------------------------------------------------------

    LOG.info("  ...assembling ProteinGroupInfo objects for each representative.")

    result: dict[str, ProteinGroupInfo] = {}

    for pep_set, prots in clusters:
        master = prots[0]
        if master not in selected_reps:
            continue

        same_set     = prots[1:]                            # rest of the cluster
        peptides     = sorted(pep_set)
        unambiguous  = _unambiguous(master, pep_set)
        subsets      = sorted(master_subset_map.get(master, []))
        # Protein score: sum of -10*log10(psm_score) for all peptides in the cluster if peptide score is 0 use 200 in sum
        score = sum(
            -10 * np.log10(passing[seq].psm_score if passing[seq].psm_score > 0 else 1e-20)
            for seq in unambiguous
        )
        
        result[master] = ProteinGroupInfo(
            master_protein      = master,
            info = ProteinInfo(),
            same_set_protein_ids= same_set,
            subset_protein_ids  = subsets,
        )

        result[master].info.peptides.update(peptides)
        result[master].info.proteotypic_peptides.update(unambiguous)
        result[master].info.score = score

    LOG.info(f"  ...assembled ProteinGroupInfo for {len(result)} representative proteins.")

    return result


def protein_annotation (peptides: Dict[str, PeptideInfo], proteins: Dict[str, ProteinInfo], reference_proteins: Dict[str, Dict[str, str]], gene_map: Dict[str, Dict[str, set[str]]], transcript_map: Dict[str, Dict[str, set[str]]], fdr_threshold: float ) -> Tuple[ Dict[str, ProteinGroupInfo], Dict[str, ProteinGroupInfo] ]:

    #Conduct protein / gene clustering and inference based on the peptide mappings
    
    LOG.info("  ...running protein inference and clustering.")
    protein_groups = run_protein_inference(peptides, fdr_threshold, False)

    for master_protein, group in protein_groups.items():

        #Add master protein info to protein cluster
        group.info.type = proteins[master_protein].type
        group.info.decoy = proteins[master_protein].decoy
        group.info.gene_ids = proteins[master_protein].gene_ids
        group.info.transcript_ids = proteins[master_protein].transcript_ids
        group.info.experiments = proteins[master_protein].experiments
        group.info.search_engines = proteins[master_protein].search_engines
        group.info.databases = proteins[master_protein].databases
        group.info.length = proteins[master_protein].length

        #loop same set proteins and add their info to the cluster
        for prot in group.same_set_protein_ids:
            if prot in proteins:
                if group.info.type != "canonical":
                    if proteins[prot].type == "canonical":
                        group.info.type = "canonical"
                    elif proteins[prot].type == "variant" and group.info.type != "canonical":
                        group.info.type = "variant"
                    elif proteins[prot].type == "contaminant" and group.info.type == "unknown":
                        group.info.type = "contaminant"
                    elif proteins[prot].type == "non_canonical" and group.info.type not in ["canonical", "variant"]:
                        group.info.type = "non_canonical"

                group.info.decoy = group.info.decoy or proteins[prot].decoy
                group.info.gene_ids.update(proteins[prot].gene_ids)
                group.info.transcript_ids.update(proteins[prot].transcript_ids)
                group.info.experiments.update(proteins[prot].experiments)
                group.info.search_engines.update(proteins[prot].search_engines)
                group.info.databases.update(proteins[prot].databases)
                group.info.length = max(group.info.length, proteins[prot].length)
            else:
                LOG.warning(f"Protein {prot} in same set cluster for master protein {master_protein} not found in proteins dictionary.")
    
        #Run a final check for reference match
        if group.info.type != "canonical":
            ptypes = set()
            for peptide in group.info.proteotypic_peptides:
                if peptide in peptides:
                    ptypes.add(peptides[peptide].type)
            
            if "canonical" in ptypes and len(ptypes) == 1:
                group.info.type = "canonical"
            elif "non_canonical" in ptypes:
                group.info.type = "non_canonical"
            elif "variant" in ptypes:
                group.info.type = "variant"
            elif "contaminant" in ptypes:
                group.info.type = "contaminant"

            if "unknown" in ptypes: 
                if len(ptypes) == 1:
                    group.info.type = "unknown"
                else:
                    LOG.warning(f"Master protein {master_protein} has proteotypic peptides with conflicting types: {ptypes}. Keeping type as {group.info.type}.")

    LOG.info(f"  ...{len(protein_groups)} protein groups identified after inference and clustering.")

    LOG.info("  ...collating gene level inference information.")

    gene_groups = run_protein_inference(peptides, fdr_threshold, True)

    #Add master gene info to gene cluster
    for master_gene, group in gene_groups.items():
        group.info.gene_ids.add(master_gene)
        if master_gene in gene_map:
            for transcript in gene_map[master_gene]["transcript_ids"]:
                for prot_id in transcript_map[transcript]["protein_ids"]:
                    if prot_id in reference_proteins:
                        group.info.type = "canonical"
            group.info.transcript_ids.update(gene_map[master_gene]["transcript_ids"])
            group.info.length = max(len(reference_proteins[prot_id]["sequence"]) for transcript in gene_map[master_gene]["transcript_ids"] for prot_id in transcript_map[transcript]["protein_ids"] if prot_id in reference_proteins)
        else:
            group.info.type = "non_canonical"

        for gene in group.same_set_protein_ids:
            group.info.gene_ids.add(gene)
            if gene in gene_map:
                for transcript in gene_map[gene]["transcript_ids"]:
                    for prot_id in transcript_map[transcript]["protein_ids"]:
                        if prot_id in reference_proteins:
                            group.info.type = "canonical"
                group.info.transcript_ids.update(gene_map[gene]["transcript_ids"])
                gene_length = max(len(reference_proteins[prot_id]["sequence"]) for transcript in gene_map[gene]["transcript_ids"] for prot_id in transcript_map[transcript]["protein_ids"] if prot_id in reference_proteins)
                group.info.length = max(group.info.length, gene_length)
        
        for peptide in group.info.proteotypic_peptides:
            if peptide in peptides:
                group.info.experiments.update(peptides[peptide].experiments)
                group.info.search_engines.update(peptides[peptide].search_engines)
            else:
                LOG.warning(f"Peptide {peptide} in proteotypic peptide set for master gene {master_gene} not found in peptides dictionary.")

    LOG.info(f"  ...{len(gene_groups)} gene groups identified after inference and clustering.")

    return protein_groups, gene_groups

# =============================================================================
# Output File Writing
# =============================================================================

def write_output_files(
    psms: Dict[str, PSMInfo], 
    peptides: Dict[str, PeptideInfo], 
    proteins: Dict[str, ProteinGroupInfo], 
    genes: Dict[str, ProteinGroupInfo], 
    heads: Dict[str, set[str]], 
    gene_map: Dict[str, Dict[str, set[str]]], 
    sample_info: Dict[str, Dict[str, str]], 
    peptide_fdr: float, 
    psm_fdr: float,
    psm_score: float,
    prefix: str ) -> None:

    #Write annotated PSMs to file
    psm_output_file = os.path.join(f"{prefix}PSMs.tsv")
    with open(psm_output_file, "w") as f:
        psm_header = [
            "Experiment",
            "SpectrumID",
            "PeptideSequence",
            "Database",
            "SearchEngine",
            "Score(PEP)",
            "FDR",
            "DeltaScoreToBest",
            "BestPeptide",
            "ConflictingPeptidesAndScores",
            "RetentionTime",
            "MZ",
            "Charge"
        ]
        f.write("\t".join(psm_header) + "\n")

        for spectrum_id, peptide_dict in psms.items():
            for peptide, database_dict in peptide_dict.items():
                for database, engine_dict in database_dict.items():
                    for engine, psm_info in engine_dict.items():

                        conflicting_str = ";".join(f"{pep}:{score:.2f}" for pep, score in psm_info.conflicts.items())

                        line = [
                            spectrum_id.split("_")[0],  # Experiment
                            spectrum_id,
                            peptide,
                            database,
                            engine,
                            f"{psm_info.score}",
                            f"{psm_info.fdr:.4f}",
                            f"{psm_info.delta_score:.4f}" if psm_info.delta_score is not None else "",
                            str(psm_info.best_peptide),
                            conflicting_str,
                            f"{psm_info.rt:.2f}" if psm_info.rt is not None else "",
                            f"{psm_info.mz:.4f}" if psm_info.mz is not None else "",
                            str(psm_info.charge) if psm_info.charge is not None else ""
                        ]
                        f.write("\t".join(line) + "\n")
        
    #Write annotated peptides to file
    peptide_output_file = os.path.join(f"{prefix}Peptides.tsv")
    with open(peptide_output_file, "w") as f:
        peptide_header = [
            "PeptideSequence",
            "Type",
            "Decoy",
            "Gene_Ambiguity",
            "Length",
            "Peptide_FDR",
            "Group_FDR",
            "Group_lFDR",
            "Global_FDR",
            "PSMCount",
            "bestPSM",
            "bestPSMScore(PEP)",
            "bestPSMFDR",
            "DeltaScore",
            "GeneNames",
            "Genes",
            "Transcripts",
            "ProteinCount",
            "Modifications",
            "PTM_Peptiforms",
            "Conflicts",
            "DistanceToReference",
            "CanonicalSequence",
            "PreAA",
            "PostAA"
        ]

        for db in sorted(heads["Databases"]):
            peptide_header.append(f"{db}_ProteinIds")
            peptide_header.append(f"{db}_BestScore(PEP)")
            peptide_header.append(f"{db}_PeptideFDR")

        for engine in sorted(heads["SearchEngines"]):
            peptide_header.append(f"{engine}_BestScore(PEP)")

        for experiment in sorted(heads["Experiments"]):
            peptide_header.append(f"{experiment}_BestScore(PEP)")

        for experiment in sorted(heads["Experiments"]):
            for tmt in sorted(heads["TMTChannels"]):
                sample = sample_info.get(experiment, {}).get(tmt, "")
                peptide_header.append(f"{experiment}_{tmt}_{sample}_RawOMSAbundance")

        f.write("\t".join(peptide_header) + "\n")

        for sequence, info in peptides.items():

            if info.peptide_fdr <= peptide_fdr and info.top_rank and (info.psm_fdr <= psm_fdr or info.psm_score <= psm_score):

                database_strs = []
                for db in sorted(heads["Databases"]):
                    if db in info.protein_ids:
                        database_strs.append(",".join(info.protein_ids[db]))
                        database_strs.append(f"{info.database_score.get(db, '99')}")
                        database_strs.append(f"{info.database_fdr.get(db, '99'):.4f}")
                    else:
                        database_strs.append("-")
                        database_strs.append("-")
                        database_strs.append("-")

                engine_score_strs = []
                for engine in sorted(heads["SearchEngines"]):
                    if engine in info.search_engines:
                        engine_score_strs.append(f"{info.search_engines.get(engine, '')}")
                    else:
                        engine_score_strs.append("")
                    
                experiment_score_strs = []
                for experiment in sorted(heads["Experiments"]):
                    if experiment in info.experiments:
                        experiment_score_strs.append(f"{info.experiments.get(experiment, '')}")
                    else:
                        experiment_score_strs.append("")

                tmt_abundance_strs = []
                for experiment in sorted(heads["Experiments"]):
                    for tmt in sorted(heads["TMTChannels"]):
                        abundance = ""
                        for psm_id in peptides[sequence].quantification:
                            if experiment in psm_id:
                                if tmt in peptides[sequence].quantification[psm_id]:
                                    abundance = f"{peptides[sequence].quantification[psm_id][tmt]:.2f}"
                                    break
                        tmt_abundance_strs.append(abundance)

                conflict_str = ";".join(f"{pep}:{score:.2f}" for pep, score in info.conflicts.items())
                group_fdr_qvalue = info.group_fdr.get("qvalue", 99.0)
                group_lfdr_qvalue = info.group_lfdr.get("qvalue", 99.0)
                global_fdr_qvalue = info.group_fdr.get("global_qvalue", 99.0)

                #remove empty gene id from info.transcripts if present
                if "" in info.genes:
                    info.genes.remove("")

                #Gather GeneNames
                genenames = set()
                for gene in info.genes:
                    if gene in gene_map:
                        genenames.update(gene_map[gene]["gene_name"])

                #remove empty transcript id from info.transcripts if present
                if "" in info.transcripts:
                    info.transcripts.remove("")

                if info.delta_score == 'inf':
                    info.delta_score = None
                
                line = [
                    sequence,
                    info.type,
                    str(info.decoy),
                    str(len(info.genes)),
                    str(len(sequence)),
                    f"{info.peptide_fdr:.4f}",
                    f"{group_fdr_qvalue:.4f}",
                    f"{group_lfdr_qvalue:.4f}",
                    f"{global_fdr_qvalue:.4f}",
                    str(len(info.psms)),
                    info.best_psm,
                    f"{info.psm_score}" if info.psm_score is not None else "",
                    f"{info.psm_fdr:.4f}" if info.psm_fdr is not None else "",
                    f"{info.delta_score:.4f}" if info.delta_score is not None else "",
                    ",".join(sorted(genenames)),
                    ",".join(sorted(info.genes)),
                    ",".join(sorted(info.transcripts)),
                    str(len(info.protein_ids)),
                    ",".join(info.modifications),
                    ",".join(info.peptidoforms),
                    conflict_str,
                    str(info.minSAAVdistance) if info.minSAAVdistance is not None else "",
                    info.canonical_sequence if info.canonical_sequence is not None else "",
                    ",".join(info.preaa),
                    ",".join(info.postaa)
                ] + database_strs + engine_score_strs + experiment_score_strs + tmt_abundance_strs
                f.write("\t".join(line) + "\n")

    #Write annotated proteins to file

    protein_output_file = os.path.join(f"{prefix}Proteins.tsv")
    with open(protein_output_file, "w") as f:
        protein_header = [
            "ProteinClusterID",
            "GroupSize",
            "MasterProteinID",
            "SameSetProteins",
            "GeneNames",
            "Type",
            "Decoy",
            "Protein_Score",
            "Gene_Ambiguity",
            "Length",
            "Genes",
            "Transcripts",
            "Peptide_Count",
            "Pepides",
            "ProteotypicPeptide_Count",
            "ProteotypicPeptides",
            "Experiments",
            "SearchEngines",
            "SubsetProteins"
        ]
        f.write("\t".join(protein_header) + "\n")

        counter = 1

        #sort proteins by score
        for protein_id, group in sorted(proteins.items(), key=lambda x: x[1].info.score, reverse=True):

            #Skip proteins with no proteotypic peptides 
            if len(group.info.proteotypic_peptides) == 0:
                continue

            #remove empty gene/transcript id from info if present
            if "" in group.info.gene_ids:
                group.info.gene_ids.remove("")
                if "" in info.genes:
                    info.genes.remove("")

            if "" in group.info.transcript_ids:
                group.info.transcript_ids.remove("")
                if "" in info.transcripts:
                    info.transcripts.remove("")

            #Gather GeneNames
            genenames = set()
            for gene in group.info.gene_ids:
                if gene in gene_map:
                    genenames.update(gene_map[gene]["gene_name"])

            line = [
                str(counter),
                str(len(group.same_set_protein_ids) + 1),
                protein_id,
                ",".join(group.same_set_protein_ids),
                ",".join(sorted(genenames)),
                group.info.type,
                str(group.info.decoy),
                f"{group.info.score:.2f}",
                str(len(group.info.gene_ids)),
                str(group.info.length) if group.info.length is not None else "",
                ",".join(group.info.gene_ids),
                ",".join(group.info.transcript_ids),
                str(len(group.info.peptides)),
                ",".join(group.info.peptides),
                str(len(group.info.proteotypic_peptides)),
                ",".join(group.info.proteotypic_peptides),
                ",".join(group.info.experiments),
                ",".join(group.info.search_engines),
                ",".join(group.subset_protein_ids)
            ]
            f.write("\t".join(line) + "\n")

            counter += 1

    #Write annotated genes to file
    gene_output_file = os.path.join(f"{prefix}Genes.tsv")
    with open(gene_output_file, "w") as f:
        gene_header = [
            "GeneClusterID",
            "GroupSize",
            "MasterGeneID",
            "SameSetGenes",
            "GeneNames",
            "Type",
            "Score",
            "Length",
            "Transcripts",
            "Peptide_Count",
            "Peptides",
            "ProteotypicPeptide_Count",
            "ProteotypicPeptides",
            "Experiments",
            "SearchEngines",
            "SubsetGenes"
        ]
        f.write("\t".join(gene_header) + "\n")

        counter = 1

        for gene_id, group in sorted(genes.items(), key=lambda x: x[1].info.score, reverse=True):

            if gene_id == "":
                continue

            #Gather GeneNames
            genenames = set()
            for gene in group.info.gene_ids:
                if gene in gene_map:
                    genenames.update(gene_map[gene]["gene_name"])

            line = [
                str(counter),
                str(len(group.same_set_protein_ids) + 1),
                gene_id,
                ",".join(group.same_set_protein_ids),
                ",".join(sorted(genenames)),
                group.info.type,
                f"{group.info.score:.2f}",
                str(group.info.length) if group.info.length is not None else "",
                ",".join(group.info.transcript_ids),
                str(len(group.info.peptides)),
                ",".join(group.info.peptides),
                str(len(group.info.proteotypic_peptides)),
                ",".join(group.info.proteotypic_peptides),
                ",".join(group.info.experiments),
                ",".join(group.info.search_engines),
                ",".join(group.subset_protein_ids)
            ]
            f.write("\t".join(line) + "\n")

            counter += 1

def write_mztab_output(
    psms: Dict[str, PSMInfo], 
    peptides: Dict[str, PeptideInfo], 
    proteins: Dict[str, ProteinInfo], 
    genes: Dict[str, ProteinInfo], 
    metadata: Dict[str, MzTabMetaInfo] , 
    heads: Dict[str, set[str]], 
    psm_fdr: float, 
    prefix: str
    ) -> None:

    #Function to write annotated results to mzTab format, including PSM, Peptide, and Protein sections with appropriate metadata and annotations based on the input data structures.

    mztab_output_path = os.path.join(f"{prefix}Annotated_Results.mzTab")
    with open(mztab_output_path, "w") as mzOUT:
        mzOUT.write("MTD\tmzTab-version\t1.0.0\n")
        mzOUT.write("MTD\tmzTab-mode\tSummary\n")

        if len(heads["TMTChannels"]) > 0:
            mzOUT.write("MTD\tmzTab-type\tQuantification\n")
        else:
            mzOUT.write("MTD\tmzTab-type\tIdentification\n")

        mzOUT.write("MTD\ttitle\t" + prefix + " PPP NF-OpenMS Final Results\n")
        mzOUT.write("MTD\tdescription\tPersonal Proteomics Pipeline Final Results Export\n")
        mzOUT.write("MTD\tprotein_search_engine_score[1]\t[, , one-peptide-rule, ]\n")
        mzOUT.write("MTD\tpeptide_search_engine_score[1]\t[MS, MS:1003114, OpenMS:Best PSM Score, ]\n")
        mzOUT.write("MTD\tpsm_search_engine_score[1]\t[, , OpenMS/ConsensusID_best Posterior Probability, ]\n")

        mzOUT.write("MTD\tsoftware[1]\t[MS, MS:1000752, TOPP software, 2.6.0-pre-exported-20201001]\n")
        mzOUT.write("MTD\tsoftware[2]\t[MS, MS:1003118, EPIFANY, 2.6.0-pre-exported-20201001]\n")
        mzOUT.write("MTD\tsoftware[3]\t[MS, MS:1002188, TOPP ConsensusID, 2.6.0-pre-exported-20201001]\n")
        mzOUT.write("MTD\tsoftware[4]\t[, , PIPPA Pipeline - NF-OpenMS (Institute of Cancer Research), ]\n")
        mzOUT.write("MTD\tsoftware[5]\t[MS, MS:1001490, Percolator, ]\n")

        so = 6
        for se in sorted(heads["SearchEngines"]):

            if se == 'mergedresults':
                soft = "[MS, MS:1002188, TOPP ConsensusID, ]"
            elif se == 'comet':
                soft = "[MS, MS:1002251, Comet, ]"
            elif se == 'msgfplus':
                soft = "[MS, MS:1002048, MS-GF+, ]"
            elif se == 'msfragger':
                soft = "[MS, MS:1003014, MSFragger, ]"
            elif se == 'sage':
                soft = "[, , SAGE, ]"
            else:
                soft = "[, , " + se + ", ]"

            mzOUT.write("MTD\tsoftware["+str(so)+"]\t"+soft+"\n")
            so += 1

        pcount = 1
        for modkey in metadata.variable_modifications:
            mzOUT.write("MTD\tvariable_mod[" + str(pcount) + "]\t" + metadata.variable_modifications[modkey].info + "\n")
            mzOUT.write("MTD\tvariable_mod[" + str(pcount) + "]-site\t" + metadata.variable_modifications[modkey].site +"\n")
            mzOUT.write("MTD\tvariable_mod[" + str(pcount) + "]-position\t" + metadata.variable_modifications[modkey].pos +"\n")
            pcount += 1
        
        pcount = 1
        for modkey in metadata.fixed_modifications:
            mzOUT.write("MTD\tfixed_mod[" + str(pcount) + "]\t" + metadata.fixed_modifications[modkey].info + "\n")
            mzOUT.write("MTD\tfixed_mod[" + str(pcount) + "]-site\t" + metadata.fixed_modifications[modkey].site +"\n")
            mzOUT.write("MTD\tfixed_mod[" + str(pcount) + "]-position\t" + metadata.fixed_modifications[modkey].pos +"\n")
            pcount += 1

        if len(heads["TMTChannels"]) > 0:
            mzOUT.write("MTD\tquantification_method\t" + metadata.quantification_method +"\n")
            sv = 1
            sa = 1
            for tag in sorted(heads["TMTChannels"]):
                mzOUT.write("MTD\tstudy_variable["+str(sv)+"]-description\tTMT "+str(tag)+" Sample "+str(sa)+"\n")
                sa += 1
                #mzOUT.write("MTD\tstudy_variable["+str(sv)+"]-description\tTMT "+str(tag)+" Channel not used!\n")
                sv += 1

        mzOUT.write("MTD\tms_run[1]-format\t[MS, MS:1000584, mzML file, ]\n")

        ## EXTRACT FROM mzTab input but need to drop _filtered for the quant as this is not used in original input mzMl file name
        mzOUT.write("MTD\tms_run[1]-location\tfile://" + metadata.mzml_file + "\n")

        mzOUT.write("MTD\tms_run[1]-id_format\t[MS, MS:1000777, spectrum identifier nativeID format, ]\n\n")

        
        ##PROTEIN RESULTS SECTION
        ##PRH - Protein header line
        mzOUT.write("PRH\taccession\tdescription\ttaxid\tspecies\tdatabase\tdatabase_version\tsearch_engine\tbest_search_engine_score[1]\tambiguity_members\tmodifications\tprotein_coverage")
        
        '''
        #Protein Quant columns - not currently implemented 
        if len(heads["TMTChannels"]) > 0:
            sv = 1
            for tag in sorted(heads["TMTChannels"]):
                mzOUT.write("\tprotein_abundance_study_variable["+str(sv)+"]\tprotein_abundance_stdev_study_variable["+str(sv)+"]\tprotein_abundance_std_error_study_variable["+str(sv)+"]")
                sv += 1
        '''

        mzOUT.write("\tnum_psms_ms_run[1]\tnum_peptides_distinct_ms_run[1]\tnum_peptides_unique_ms_run[1]")
        mzOUT.write("\topt_global_cv_PRIDE:0000303_decoy_hit\topt_global_result_type\topt_global_protein_FDR\topt_global_gene\topt_global_transcript\n")    


        #PRT WRITE PROTEINS
        for protein_id, info in proteins.items():
            
            mzOUT.write(f"PRT\t{protein_id}\t{info.type}\tnull\tnull\t{'|'.join(info.databases)}\tnull\t[MS, MS:1002188, TOPP ConsensusID, ]\t{info.score}\tnull\tnull\tnull")
            mzOUT.write(f"\tnull\t{len(info.peptides)}\t{len(info.proteotypic_peptides)}")
            mzOUT.write(f"\t{str(info.decoy)}\t{info.type}\tnull\t{'|'.join(info.gene_ids)}\t{'|'.join(info.transcript_ids)}\n")

        ##PEPTIDE RESULTS SECTION
        ##PEH - Peptide header line
        mzOUT.write("PEH\tsequence\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\tbest_search_engine_score[1]\tsearch_engine_score[1]_ms_run[1]\tmodifications\tretention_time\tretention_time_window\tcharge\tmass_to_charge\tspectra_ref")
        if len(heads["TMTChannels"]) > 0:
            sv = 1
            for tag in sorted(heads["TMTChannels"]):
                mzOUT.write("\tpeptide_abundance_study_variable["+str(sv)+"]\tpeptide_abundance_stdev_study_variable["+str(sv)+"]\tpeptide_abundance_std_error_study_variable["+str(sv)+"]")
                sv += 1
        mzOUT.write("\topt_global_cv_MS:1000889_peptidoform_sequence\topt_global_cv_MS:1002217_decoy_peptide\topt_global_PSMs\topt_global_Peptide_Type\topt_global_Peptide_FDR\topt_global_ENSG\topt_global_ENST\n")
    
        peptide_spectrum_map: Dict[str, Dict[int, Set[str]]] = defaultdict(lambda: defaultdict(set))
        for spectrum_id in psms:
            for peptide in psms[spectrum_id]:
                for database in psms[spectrum_id][peptide]:
                    for engine in psms[spectrum_id][peptide][database]:
                        charge = psms[spectrum_id][peptide][database][engine].charge
                        peptide_spectrum_map[peptide][charge].add(spectrum_id)

        #PEP WRITE PEPTIDES - Each peptide in each charge and for each protein 
        for peptide, info in peptides.items():

            score = info.psm_score if info.psm_score is not None else "null"
            rt = "null"
            mz = "null"
            modifications = ",".join(info.modifications) if info.modifications else "null"
            peptiforms = ",".join(info.peptidoforms) if info.peptidoforms else "null"
            decoy_str = str(info.decoy)

            
            type_str = info.type if info.type else "null"
            peptide_fdr_str = f"{info.peptide_fdr:.4f}" if info.peptide_fdr is not None else "null"



            for charge in peptide_spectrum_map[peptide]:

                specs = "|".join(peptide_spectrum_map[peptide][charge])
                psms = set()

                rt_values = []
                mz_values = []
                for spec in peptide_spectrum_map[peptide][charge]:
                    for database in psms[spec][peptide]:
                        for engine in psms[spec][peptide][database]:
                            if psms[spec][peptide][database][engine].charge == charge:

                                if psms[spec][peptide][database][engine].fdr < psm_fdr:
                                    psms.add(spec)

                                if psms[spec][peptide][database][engine].rt is not None:
                                    rt_values.append(psms[spec][peptide][database][engine].rt)
                                if psms[spec][peptide][database][engine].mz is not None:
                                    mz_values.append(psms[spec][peptide][database][engine].mz)

                if len(rt_values) > 0:
                    rt = sum(rt_values) / len(rt_values)

                if len(mz_values) > 0:
                    mz = sum(mz_values) / len(mz_values)

                psm_count = len(psms)

                for database in info.protein_ids:
                    for protein in info.protein_ids[database]:

                        genes = "|".join(proteins[protein].gene_ids) if proteins[protein].gene_ids else "null"
                        transcripts = "|".join(proteins[protein].transcript_ids) if proteins[protein].transcript_ids else "null"

                        unique = 0
                        if peptide in proteins[protein].proteotypic_peptides:
                            unique = 1

                        mzOUT.write(f"PEP\t{peptide}\t{protein}\t{unique}\t{database}\tnull\t[MS, MS:1002188, TOPP ConsensusID, ]\t{score}\t{score}\t{modifications}\t{rt}\tnull\t{charge}\t{mz}\t{specs}")
                        if len(heads["TMTChannels"]) > 0:
                            for tag in sorted(heads["TMTChannels"]):
                                abundance  = info.quantification[info.best_psm][tag] if info.best_psm in info.quantification and tag in info.quantification[info.best_psm] else 'null'
                                mzOUT.write(f"\t{abundance}\tnull\tnull")
                        mzOUT.write(f"\t{peptiforms}\t{decoy_str}\t{psm_count}\t{type_str}\t{peptide_fdr_str}\t{genes}\t{transcripts}\n")

        ##PSM RESULTS SECTION
        ##PSH - PSM header line
        mzOUT.write("PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\tsearch_engine_score[1]\tmodifications\tretention_time\tcharge\texp_mass_to_charge\tcalc_mass_to_charge\tspectra_ref\tpre\tpost\tstart\tend\topt_global_cv_MS:1002217_decoy_peptide\topt_global_cv_MS:1000889_peptidoform_sequence\topt_global_Peptide_Type\topt_global_PSM_FDR\topt_global_ENSG\topt_global_ENST\n")

        #PSM WRITE PSMs
        for spec in psms:
            for peptide in psms[spec]:

                decoy_str = str(peptides[peptide].decoy) if peptide in peptides else "null"
                type_str = peptides[peptide].type if peptide in peptides else "null"

                for database in psms[spec][peptide]:
                    for se, psm_info in psms[spec][peptide][database].items():

                        if psm_info.delta_score > 0:

                            if se == 'mergedresults':
                                soft = "[MS, MS:1002188, TOPP ConsensusID, ]"
                            elif se == 'comet':
                                soft = "[MS, MS:1002251, Comet, ]"
                            elif se == 'msgfplus':
                                soft = "[MS, MS:1002048, MS-GF+, ]"
                            elif se == 'msfragger':
                                soft = "[MS, MS:1003014, MSFragger, ]"
                            elif se == 'sage':
                                soft = "[, , SAGE, ]"
                            else:
                                soft = "[, , " + se + ", ]"

                            modifications = ",".join(psm_info.modifications) if psm_info.modifications else "null"
                            rt = f"{psm_info.rt:.2f}" if psm_info.rt is not None else "null"
                            charge = str(psm_info.charge) if psm_info.charge is not None else "null"
                            mz = f"{psm_info.mz:.4f}" if psm_info.mz is not None else "null"
                            psm_fdr_str = f"{psm_info.fdr:.4f}" if psm_info.fdr is not None else "null"
                            psm_score_str = f"{psm_info.score:.2f}" if psm_info.score is not None else "null"
                            peptiform = psm_info.peptidoform if psm_info.peptidoform else "null"

                            for protein in peptides[peptide].protein_ids[database]:

                                genes = "|".join(proteins[protein].gene_ids) if proteins[protein].gene_ids else "null"
                                transcripts = "|".join(proteins[protein].transcript_ids) if proteins[protein].transcript_ids else "null"

                                unique = 0
                                if peptide in proteins[protein].proteotypic_peptides:
                                    unique = 1

                                mzOUT.write(f"PSM\t{peptide}\t{spec}\t{protein}\t{unique}\t{database}\tnull\t{soft}\t{psm_score_str}\t{modifications}\t{rt}\t{charge}\t{mz}\tnull\t{spec}\tnull\tnull\tnull\tnull\t{decoy_str}\t{peptiform}\t{type_str}\t{psm_fdr_str}\t{genes}\t{transcripts}\n")


    return


# =============================================================================
# Main Function
# =============================================================================

def main() -> int:

    #Parse command-line arguments and configure logging based on the specified log level.
    args = parse_args()
    configure_logging(args.log_level)

    LOG.info("Starting Results Annotation Process.")

    sample_info = {}
    if args.sample_file:
        #Load Sample Information
        LOG.info(f"Reading Sample data from {args.sample_file}")
        sample_info = load_sample_info(Path(args.sample_file))

    #Load Reference Protein Sequences
    LOG.info(f"Reading Reference Protein Sequences from {args.ref_fasta}")
    reference_proteins, gene_transcript_map, transcript_gene_map = getReferenceProteins(Path(args.ref_fasta))

    #initialise dictionaries to hold PSM, Peptide, and Protein information
    psms: Dict[str, PSMInfo] = {}
    peptides: Dict[str, PeptideInfo] = {}
    proteins: Dict[str, ProteinInfo] = {}
    genes: Dict[str, ProteinInfo] = {}
    mztab_metadata = MzTabMetaInfo()

    additional_headers = {
        "Databases": set(),
        "SearchEngines": set(),
        "Experiments": set(),
        "Searchids": set(),
        "TMTChannels": set()
    }

    LOG.info(f"Processing Results Files in Directory: {args.results_folder}")
    #Loop through each mzTab file in results directory
    for mztab in glob.glob(args.results_folder + '/*.mzTab'):

        tabname = os.path.basename(mztab)
        LOG.info(f"Processing OpenMS mzTab file: {tabname}")
        #Parse mzTab file and extract PSM and Peptide information
        parse_mztab(
            file = mztab, 
            psms = psms,
            peptides = peptides,
            heads = additional_headers,
            metadata = mztab_metadata,
            ms3 = args.ms3,
            psm_score_threshold= args.psm_score_threshold,
            psm_fdr_threshold = args.psm_fdr_threshold
        )
    
    for percolator_tab in glob.glob(args.results_folder + '/*percolator*.tab'):
        tabname = os.path.basename(percolator_tab)
        LOG.info(f"Processing Percolator tab file: {tabname}")
        #Parse Percolator tab file and update PSM and Peptide information with scores and FDRs
        parse_percolator_tab(
            file = percolator_tab,
            psms = psms,
            peptides = peptides,
            heads = additional_headers,
            psm_score_threshold= args.psm_score_threshold,
            psm_fdr_threshold = args.psm_fdr_threshold
        )

    if len(psms) == 0:
        LOG.error("No PSMs were parsed from the results files. Please check the input files.")
        return 1
    
    psm_count = sum(len(peptide_dict) for spectrum_id, peptide_dict in psms.items())

    LOG.info(f"Parsed {psm_count} PSMs and {len(peptides)} unique peptides from the results files.")
    
    #Annotate peptides
    LOG.info("Annotating Peptides:")
    peptide_annotation(
        psms = psms,
        peptides = peptides,
        proteins = proteins,
        reference_proteins = reference_proteins,
        heads = additional_headers,
        gene_map = gene_transcript_map,
        transcript_map = transcript_gene_map,
        score_threshold = args.psm_score_threshold,
        fdr_threshold = args.pep_fdr_threshold
    )

    #Annotate proteins and genes based on peptide information and perform protein inference
    LOG.info("Annotating Proteins and Genes with Protein Inference:")
    protein_groups, gene_groups = protein_annotation(
        peptides = peptides,
        proteins = proteins,
        reference_proteins = reference_proteins,
        gene_map = gene_transcript_map,
        transcript_map = transcript_gene_map,
        fdr_threshold = args.pep_fdr_threshold
    )

    #Write output files
    LOG.info("Writing output files.")
    write_output_files(
        psms = psms,
        peptides = peptides,
        proteins = protein_groups,
        genes = gene_groups,
        heads = additional_headers,
        gene_map = gene_transcript_map,
        sample_info = sample_info,
        peptide_fdr = args.pep_fdr_threshold,
        prefix = args.out_prefix,
        psm_score = args.psm_score_threshold,
        psm_fdr = args.psm_fdr_threshold
    )

    if args.mztab_out:
        LOG.info("Writing mzTab output file.")
        write_mztab_output(
            psms= psms,
            peptides= peptides,
            proteins= proteins,
            genes= genes,
            metadata= mztab_metadata,
            heads= additional_headers,
            psm_fdr= args.psm_fdr_threshold,
            prefix= args.out_prefix
        )
    
    LOG.info("Results Annotation Process Completed Successfully.")


# Entry point of the script
if __name__ == "__main__":
    raise SystemExit(main())