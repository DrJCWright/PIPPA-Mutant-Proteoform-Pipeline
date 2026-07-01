#!/usr/bin/env python3
"""
Build cell-line-specific mutant protein FASTA databases from a COSMIC TSV and
GENCODE protein/transcript FASTA files.

James Wright (2025),
The Institute of Cancer Research, London
"""

import argparse
import csv
import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import DefaultDict, Dict, Iterator, List, Optional, Sequence, Tuple, IO


#Setup Logging Handle
LOG = logging.getLogger("cosmic_build_seq_db")

# Standard genetic code mapping codons to amino acids
GENETIC_CODE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "_", "TAG": "_",
    "TGC": "C", "TGT": "C", "TGA": "_", "TGG": "W",
}

# COSMIC variants data class to store relevant information for each mutation
@dataclass(frozen=True)
class VariantRecord:
    aa: str
    na: str
    zygosity: str
    cosmic_id: str
    mutation_type: str

# Gene information class to store ENSG and gene name for each transcript
@dataclass
class GeneInfo:
    ensg: str = ""
    gene_name: str = ""

# Sequence combination class to represent a specific combination of amino acid and nucleotide sequences, along with a label for the variant combination
@dataclass
class SequenceCombination:
    aa: List[str]
    na: List[str]
    var_label: str = ""
    vlist: List[VariantRecord] = field(default_factory=list)
    poslist: List[int] = field(default_factory=list)

    # Method to create a copy of the SequenceCombination instance, allowing for branching of variant combinations when needed
    def clone(self) -> "SequenceCombination":
        return SequenceCombination(self.aa.copy(), self.na.copy(), self.var_label, self.vlist.copy() if self.vlist else [], self.poslist.copy() if self.poslist else [])

# Nested dictionary structure to store variants by sample, protein ID, and amino acid position
VariantMap = DefaultDict[str, DefaultDict[str, DefaultDict[int, Dict[str, VariantRecord]]]]

# Error writer class to handle logging of warnings and errors during variant processing, as well as tracking counts of total and failed variants
class ErrorWriter:
    def __init__(self, path: Path):
        self.path = path
        self.handle = path.open("w", encoding="utf-8")
        self.fail_count = 0
        self.total_count = 0

    def warn(self, message: str) -> None:
        self.handle.write(message.rstrip() + "\n")

    def fail(self, message: str) -> None:
        self.fail_count += 1
        self.warn(message)

    def seen(self) -> None:
        self.total_count += 1

    def close(self) -> None:
        self.warn(f"FAILED VARS=={self.fail_count} out of {self.total_count}")
        self.handle.close()


# Function to load command line arguments
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='''
            Tool to extract mutations from a COSMIC TSV and map them into GENCODE proteins 
            to create cell-line-specific mutant FASTA databases.
            Supports missense, nonsense, in-frame indels, frameshift indels, and nonstop mutations. 
            Unrecognized or complex mutations are logged to per-sample .err files for review.
            ---- James Wright (2025) - The Institute of Cancer Research, London ----'''
    )
    parser.add_argument('--tsv', '-i', dest='cosmic_tsv', type=Path, help='COSMIC mutation export TSV, can be downloaded from cosmic website (https://cancer.sanger.ac.uk/cosmic) (i.e. CellLinesProject_GenomeScreensMutant_v103_GRCh38.tsv)', required=True)
    parser.add_argument('--aafasta', '-p', dest='protein_fasta', type=Path, help='GENCODE protein translated AA FASTA downloaded from GENCODE webpage (i.e. gencode.v36.pc_translations.fa)', required=True)
    parser.add_argument('--nafasta', '-n', dest='transcript_fasta', type=Path, help='GENCODE protein transcript NA FASTA downloaded from GENCODE webpage (i.e. gencode.v36.pc_transcripts.fa)', required=True)
    parser.add_argument('--prefix', '-o', dest='output_prefix', type=Path, help='Output file prefix')
    parser.add_argument('--sample_info', '-s', dest='sample_info', type=Path, help='Sample table TSV (ExperimentID, TMTLabel, SampleID).')
    parser.add_argument('--variant_distance', '-d', dest='variant_distance', type=int, default=30, help='AA Distance threshold for variant mutant combinations (default: 30 aa)')
    parser.add_argument('--log_level', '-l', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], help='Logging level (Default: INFO)')
    parser.add_argument('--locktranscriptversion', '-v', dest='lock_transcript_version', action='store_true', help='Lock transcript versions when matching COSMIC variants to GENCODE transcripts. By default, transcript versions are stripped to allow matching across different GENCODE releases. Enabling this option will require exact version matches between COSMIC and GENCODE transcript IDs which reduces mapping errors but some transcripts will be lost if using newer GENCODE release than COSMIC.')
    parser.add_argument('--oldcosmic', '-c', dest='old_cosmic', action='store_true', help='Use this flag if using an older COSMIC TSV export with different column headers (pre-v100). This will switch to the old column headers for parsing the TSV file.')
    return parser.parse_args()


# Function to configure logging based on the specified log level
def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="[%(levelname)s] %(message)s",
    )

# Function to clean and standardize sample names by removing non-alphanumeric characters and converting to uppercase
def clean_sample_name(name: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "", name.upper())

# Function to read sample information from a TSV file and create a mapping of experiments to sample IDs and TMT labels
def parse_sample_info(sample_info_path: Path) -> Dict[str, Dict[str, str]]:
    sample_map: Dict[str, Dict[str, str]] = defaultdict(dict)
    with sample_info_path.open("r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            experiment_id = row["ExperimentID"]
            sample_id = clean_sample_name(row["SampleID"])
            tmt_label = row["TMTLabel"]
            sample_map[experiment_id][sample_id] = tmt_label
    return sample_map

# Function to clean and standardize tissue names by replacing whitespace and slashes with underscores
def clean_tissue_name(name: str) -> str:
    return re.sub(r"\s+|/", "_", name)

# Function to strip version numbers from identifiers (e.g., ENSG00000354587.3 -> ENSG00000354587)
def strip_version(identifier: str) -> str:
    return identifier.split(".", 1)[0]

# Function to read FASTA files and yield header-sequence pairs as an iterator
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
                sequence.append(line)
        if header is not None:
            yield header, "".join(sequence)

# Function to translate a nucleotide sequence into an amino acid sequence
# Handles stop codons by terminating translation and logs warnings for unrecognized codons, using 'X' for unknown amino acids
def translate_na(naseq: str) -> str:
    aaseq: List[str] = []
    for i in range(0, len(naseq) - 2, 3):
        codon = (naseq[i:i + 3]).upper()
        aa = GENETIC_CODE.get(codon)
        if aa is None:
            LOG.warning("TRANSLATION WARNING CODON NOT FOUND::%s::%s::", codon, naseq)
            aaseq.append("X")
            continue
        if aa == "_":
            break
        aaseq.append(aa)
    return "".join(aaseq)

# Function to load protein sequences from a GENCODE PC TRANSLATION FASTA file. 
# filtering for selected transcript IDs and extracting gene information.
def load_protein_sequences( fasta: Path, keepversion: bool, selected_ids: set[str] ) -> Tuple[Dict[str, str], Dict[str, GeneInfo]]:
    proteins: Dict[str, str] = {}
    gene_info: Dict[str, GeneInfo] = {}

    for header, seq in read_fasta(fasta):

        # Only process headers that start with ">ENSP" to ensure we're looking at protein entries
        if not header.startswith(">ENSP"):
            continue

        # The GENCODE PC TRANSLATION FASTA header format is expected to have specific fields separated by "|".
        # Check if the header has enough fields to extract transcript ID, ENSG, and gene name
        parts = header[1:].split("|")
        if len(parts) < 7:
            continue
        
        transcript_id = parts[1] if keepversion else strip_version(parts[1])
       
        if transcript_id in selected_ids:

            # Clean the sequence by removing any non-alphanumeric characters and store it in the proteins dictionary
            proteins[transcript_id] = re.sub(r"\W+", "", seq)

            # Extract ENSG and gene name from the header parts
            ensg = strip_version(parts[2])
            gene_name = parts[6]
            gene_info[transcript_id] = GeneInfo(ensg=ensg, gene_name=gene_name)
        

    return proteins, gene_info

# Function to load transcript sequences from a GENCODE FASTA file.
# Filter for selected transcript IDs and extract CDS and 3'UTR positions.
def load_transcript_sequences( fasta: Path, keepversion: bool, selected_ids: set[str], require_utr: set[str] ) -> Tuple[Dict[str, str], Dict[str, str]]:
    transcripts: Dict[str, str] = {}
    utr: Dict[str, str] = {}

    #Define regular expressions to extract ENST id, CDS and 3'UTR positions from the FASTA headers
    cds_re = re.compile(r"CDS:(\d+)-(\d+)")
    utr_re = re.compile(r"UTR3:(\d+)-(\d+)")
    id_re = re.compile(r"(ENSTR?\d+.\d+)") if keepversion else re.compile(r"(ENSTR?\d+)")

    # Iterate through the FASTA records and process those that match the selected transcript IDs
    for header, seq in read_fasta(fasta):
        seq = list(re.sub(r"\W+", "", seq))
        id = id_re.search(header)
        if not id:
            continue
        transcript_id = id.group(1)
        if transcript_id not in selected_ids:
            continue

        # Extract CDS positions from the header using the defined regular expression. If no CDS is found, skip this record.
        cds = cds_re.search(header)
        if not cds:
            continue
        cds_start = int(cds.group(1)) - 1 # Convert to 0-based index
        cds_end = int(cds.group(2))
        transcripts[transcript_id] = "".join(seq[cds_start:cds_end])

        # If this transcript is in the set of transcripts that require 3'UTR sequences, extract the UTR positions and store the UTR sequence.
        if transcript_id in require_utr:
            utrcoords = utr_re.search(header)
            if utrcoords:
                utr_start = int(utrcoords.group(1)) - 1
                utr_end = int(utrcoords.group(2))
                utr[transcript_id] = "".join(seq[utr_start:utr_end])

    return transcripts, utr


# Function to extract the amino acid position from a COSMIC AA change string using regular expressions.
def extract_aa_position(aa_change: str) -> Optional[int]:
    match = re.search(r"p\.\D+(\d+)\D+", aa_change)
    return int(match.group(1)) if match else None

# Function to load variants from a COSMIC TSV file, filtering for requested samples and valid mutation types.
def load_cosmic_variants( cosmic_tsv: Path, requested_samples: set[str], keepversion: bool, oldcosmic: bool) -> Tuple[VariantMap, Dict[str, str], set[str], set[str]]:
    # Define a nested dictionary structure to store variants by sample, protein ID, and amino acid position. 
    sampledata: VariantMap = defaultdict(
        lambda: defaultdict(lambda: defaultdict(dict))
    )

    #dictionaries of old and new headers to handle changes in COSMIC TSV format over time
    old_header_map = {
        "SAMPLE_NAME": "Sample name",
        "TRANSCRIPT_ACCESSION": "Accession Number",
        "MUTATION_DESCRIPTION": "Mutation Description",
        "MUTATION_AA": "Mutation AA",
        "MUTATION_CDS": "Mutation CDS",
        "GENOMIC_MUTATION_ID": "GENOMIC_MUTATION_ID",
        "MUTATION_ZYGOSITY": "Mutation zygosity",
        "COSMIC_SAMPLE_ID": "Primary site"
    }

    header_map = {
        "SAMPLE_NAME": "SAMPLE_NAME",
        "TRANSCRIPT_ACCESSION": "TRANSCRIPT_ACCESSION",
        "MUTATION_DESCRIPTION": "MUTATION_DESCRIPTION",
        "MUTATION_AA": "MUTATION_AA",
        "MUTATION_CDS": "MUTATION_CDS",
        "GENOMIC_MUTATION_ID": "GENOMIC_MUTATION_ID",
        "MUTATION_ZYGOSITY": "MUTATION_ZYGOSITY",
        "COSMIC_SAMPLE_ID": "COSMIC_SAMPLE_ID"
    }

    if oldcosmic:
        header_map = old_header_map


    #Define dictionaries to store tissue information and sets to track needed transcripts and those requiring 3'UTR sequences.
    tissues: Dict[str, str] = {}
    selected_transcripts: set[str] = set()
    require_utr: set[str] = set()
    
    #Create dictionary to count total unique cosmic_ids (as a set) and a count of unique cosmic_ids for each mutation type for each sample to log at the end of processing
    sample_mutation_counts: Dict[str, Dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    

    exclude_count = 0
    mutant_count = 0

    # Read the COSMIC TSV file and process each row, filtering for requested samples and valid mutation types.
    #NOTE: UTF-8 encoding is not sufficient to read the COSMIC TSV file, which contains some Latin-1 encoded characters
    with cosmic_tsv.open("r", encoding="latin-1", newline="") as COSMIC:
        reader = csv.DictReader(COSMIC, delimiter="\t")

        #Check if all required headers are present in the TSV file, if not log an error and exit
        for required_header in header_map.values():
            if required_header not in reader.fieldnames:
                LOG.error(f"Required header '{required_header}' not found in COSMIC TSV file. Please check the file format and headers.")
                exit(1)

        for row in reader:
            sample = clean_sample_name(row[header_map["SAMPLE_NAME"]])
            
            # If a specific set of samples was requested, skip any samples that are not in that set.
            if requested_samples and sample not in requested_samples:
                continue

            transcript_id = strip_version(row[header_map["TRANSCRIPT_ACCESSION"]]) if not keepversion else row[header_map["TRANSCRIPT_ACCESSION"]]
            # Only process variants that are associated with valid transcript IDs
            if not transcript_id.startswith("ENST"):
                continue

            mutation_type = row[header_map["MUTATION_DESCRIPTION"]].replace(" ", "")
            # Skip variants that are classified as silent, whole gene, or unknown mutations, as these are not relevant for building mutant protein databases.
            if re.search(r"silent|Whole|Unknown|synonymous", mutation_type, flags=re.IGNORECASE):
                continue
            
            aa_change = row[header_map["MUTATION_AA"]]
            na_change = row[header_map["MUTATION_CDS"]]
            cosmic_id = row[header_map["GENOMIC_MUTATION_ID"]]
            zygosity = row[header_map["MUTATION_ZYGOSITY"]]
            if zygosity == "":
                zygosity = "Unknown"
            tissue = clean_tissue_name(row[header_map["COSMIC_SAMPLE_ID"]])
            tissues[sample] = tissue

            # Extract the amino acid position from the AA change string. 
            # If the position cannot be extracted, log a warning and skip this variant.
            aapos = extract_aa_position(aa_change)
            if aapos is None:
                LOG.debug(
                    "EXCLUDING1::%s\t%s\t%s\t%s\t%s",
                    sample, transcript_id, -1, mutation_type, aa_change,
                )
                exclude_count += 1
                continue

            # Store the variant information in the nested dictionary structure, organized by sample, transcript ID, and amino acid position.
            sampledata[sample][transcript_id][aapos][cosmic_id] = VariantRecord(
                aa=aa_change,
                na=na_change,
                zygosity=zygosity,
                cosmic_id=cosmic_id,
                mutation_type=mutation_type,
            )

            sample_mutation_counts[sample]["TOTAL"].add(cosmic_id)
            sample_mutation_counts[sample][mutation_type].add(cosmic_id)
            sample_mutation_counts[sample][zygosity].add(cosmic_id)

            mutant_count += 1

            # Add the transcript ID to the set of selected transcripts
            selected_transcripts.add(transcript_id)

            # If the mutation type indicates a nonstop mutation, add the transcript ID to the set of transcripts that require 3'UTR sequences.
            if re.search(r"Nonstop", mutation_type, flags=re.IGNORECASE) or re.search(r"Frameshift", mutation_type, flags=re.IGNORECASE):
                require_utr.add(transcript_id)

    LOG.info(f"Finished loading {mutant_count} COSMIC mutations in {len(selected_transcripts)} protein coding transcripts.")
    LOG.info(f"Excluded {exclude_count} mutations outside CDS regions.")

    LOG.info("Unique genomic mutantion (COSMIC ID) counts by sample and type:")
    for sample, sets in sample_mutation_counts.items():
        LOG.info(f"  Sample: {sample}")
        for mtype, ids in sets.items():
            LOG.info(f"    {mtype}: {len(ids)} ")

            

    return sampledata, tissues, selected_transcripts, require_utr

# Function to build a unique identifier for a mutation sample
def build_mutation_sample_id( sample: str, transcript_id: str, pos: int, variant: VariantRecord) -> str:
    return "|".join(
        [
            sample,
            transcript_id,
            variant.cosmic_id,
            variant.aa
        ]
    )

# Function to format the output FASTA header for a given variant
def format_header( transcript_id: str, sample: str, suffix: str, gene: GeneInfo, tissue: str, variant: List[VariantRecord], pos: List[int] ) -> str:
    
    cid_str = ";".join(sorted(set(v.cosmic_id for v in variant)))
    pos_str = ";".join(sorted(set(str(p) for p in pos)))
    zygosity_str = ";".join(sorted(set(v.zygosity for v in variant)))
    aa_str = ";".join(sorted(set(v.aa for v in variant)))
    na_str = ";".join(sorted(set(v.na for v in variant)))
    mtype_str = ";".join(sorted(set(v.mutation_type for v in variant)))
    
    return (
        f">COSMICMutant_{sample}_{transcript_id}{suffix}\t{gene.ensg}\tgene={gene.gene_name}\t"
        f"tissue={tissue}\tZygosity={zygosity_str}\tCOSMIC_ID={cid_str}\t"
        f"type={mtype_str}\tpos={pos_str}\tcellline={sample}\t"
        f"AAanno={aa_str}\tNAanno={na_str}"
    )

# Function to write FASTA entry.
def write_record(out_handle, header: str, sequence: str, min_length: int = 10 ) -> None:
    if len(sequence) >= min_length:
        out_handle.write(header + "\n")
        out_handle.write(sequence + "\n")

# Function to parse the amino acid change annotation and extract the reference wild-type amino acid and position.
def parse_variant_annotation(aa_change: str) -> Optional[Tuple[str, int]]:
    match = re.match(r"p\.(\D+)(\d+)", aa_change)
    if not match:
        return None
    return match.group(1), int(match.group(2)) - 1

# Function to handel missense variants 
# by parsing the AA change annotation, finding the corresponding position in the sequence combinations, and applying the amino acid change.
def handle_missense( combos: List[SequenceCombination], variant: VariantRecord, pos: int, err: ErrorWriter, sample: str, transcript_id: str) -> str:
    
    # Parse the amino acid annotation to extract the reference amino acid, position, and variant amino acid.
    match = re.match(r"p\.(\D+)(\d+)(\D+)", variant.aa)
    if not match:
        err.fail(f"[ERROR] - {variant.mutation_type} - AA format unknown:: {variant.aa} :: {variant.na} :: {variant.cosmic_id}")
        return

    ref_aa, aa_pos_s, var_aa = match.groups()
    aa_pos = int(aa_pos_s) - 1 # Convert to 0-based index

    # Iterate through the sequence combinations and apply the missense change to the appropriate position, while also handling branching for close variants.
    o_count = len(combos) 
    for i in range(o_count):
        combo = combos[i]

        # Check if the amino acid position is valid and matches the reference amino acid in the sequence combination. 
        if aa_pos >= len(combo.aa) or combo.aa[aa_pos] != ref_aa:
            continue

        # Apply the missense change by replacing the reference amino acid with the variant amino acid in the sequence combination
        # Update the variant label for the FASTA header accordingly.
        combo.aa[aa_pos] = var_aa.upper()
        combo.var_label += f"_Mis:{ref_aa}{aa_pos_s}{var_aa}"
        combo.vlist.append(variant)
        combo.poslist.append(aa_pos_s)

        # Parse the nucleotide annotation to extract the reference nucleotide, position, and variant nucleotide, and apply the change to the sequence combination.
        na_match = re.match(r"c\.(\d+)(\D+)>(\D+)", variant.na)
        if not na_match:
            err.fail(f"[ERROR] - {variant.mutation_type} - NA format unknown :: {variant.na} :: {variant.aa} :: {variant.cosmic_id}")
            continue
        na_pos = int(na_match.group(1)) - 1
        ref_na = na_match.group(2)
        var_na = na_match.group(3)
        if na_pos >= len(combo.na):
            err.fail(
                f"[WARNING] - MISSENSE NA :: {sample} :: {transcript_id} :: position {na_pos} out of range"
            )
            continue
        if combo.na[na_pos] == ref_na:
            combo.na[na_pos] = var_na.upper()
        else:
            err.fail(
                f"[WARNING] - MISSENSE NA:: {sample} :: {transcript_id} :: Positions: {pos} == {aa_pos} and {na_pos} "
                f":: does not match {ref_na} to {combo.na[na_pos]}"
            )

    return "".join(combos[0].aa)

# Function to handle in-frame deletions 
# by parsing the AA and NA change annotations, determining the affected positions, and applying the deletions to the sequence combinations.
def handle_inframe_deletion( combos: List[SequenceCombination], variant: VariantRecord, err: ErrorWriter) -> str:
    
    # Parse the amino acid annotation to determine the start and end positions of the deletion.
    aa_single = re.match(r"p\.\D+(\d+)del$", variant.aa)
    aa_range = re.match(r"p\.\D+(\d+)_\D+(\d+)del$", variant.aa)
    na_range = re.match(r"c\.(\d+)_(\d+)del$", variant.na)

    # Determine the start and end positions for the amino acid deletion based on the parsed annotation.
    if aa_single:
        aa_start = aa_end = int(aa_single.group(1)) - 1 # Convert to 0-based index
    elif aa_range:
        aa_start = int(aa_range.group(1)) - 1 
        aa_end = int(aa_range.group(2)) - 1
    else:
        err.fail(f"[ERROR] - {variant.mutation_type} - AA format unknown:: {variant.aa} :: {variant.na} :: {variant.cosmic_id}")
        return None

    # Determine the start and end positions for the nucleotide deletion based on the parsed annotation.
    if not na_range:
        err.fail(f"[ERROR] - {variant.mutation_type} - NA format unknown:: {variant.na} :: {variant.aa} :: {variant.cosmic_id}")
        return None
    
    na_start = int(na_range.group(1)) - 1
    na_end = int(na_range.group(2)) - 1

    # Iterate through the sequence combinations and apply the in-frame deletion by removing the specified amino acids and nucleotides, 
    # while also handling branching for close variants.

    LOG.debug(f"Handling in-frame deletion: {variant.aa} :: {variant.na} :: AA positions {aa_start}-{aa_end} :: NA positions {na_start}-{na_end} :: Mutation ID {variant.cosmic_id}")
    LOG.debug(f"Current sequence combination before deletion: {''.join(combos[0].aa)}")

    for combo in combos:
        for p in range(aa_start, aa_end + 1):
            if p < len(combo.aa):
                combo.aa[p] = ""
        for p in range(na_start, na_end + 1):
            if p < len(combo.na):
                combo.na[p] = ""
        combo.var_label += f"_Del:{aa_start}-{aa_end}"
        combo.vlist.append(variant)
        combo.poslist.append(aa_start)

    LOG.debug(f"Current sequence combination after deletion: {''.join(combos[0].aa)}")
    return ''.join(combos[0].aa)
            

# Function to handle in-frame insertions 
# by parsing the AA and NA change annotations, determining the affected positions, and applying the insertions to the sequence combinations.
def handle_inframe_insertion( combos: List[SequenceCombination], variant: VariantRecord, pos: int, out_handle, sample: str, transcript_id: str, gene: GeneInfo, tissue: str, err: ErrorWriter) -> Tuple[bool, str]:
    
    # First, check for a special case of a stop-gain insertion, which is indicated by an annotation like "p.X123_*124ins*". 
    # If this pattern is detected, extract the position of the stop codon and write the truncated sequence to the output FASTA file.
    stop_insert = re.match(r"p\.\D+(\d+)_\D+\d+ins\*$", variant.aa)
    if stop_insert:
        stop_pos = int(stop_insert.group(1)) - 1
        outseq = ""
        for combo in combos:
            outseq = "".join(combo.aa[:stop_pos + 1])
            write_record(
                out_handle,
                format_header(transcript_id, sample, f"{combo.var_label}_Ins:{stop_pos}*", gene, tissue, combo.vlist + [variant], combo.poslist + [stop_pos]),
                outseq,
            )

        return True, outseq

    for combo in combos:

        # Handle different patterns of in-frame insertions based on the AA change annotation.
        # The patterns include single amino acid duplications (e.g., "p.A123dup"), range duplications (e.g., "p.A123_B125dup"), and insertions of specific amino acids (e.g., "p.A123_B125insXYZ").
        if m := re.match(r"p\.(\D+)(\d+)dup$", variant.aa):
            ref_aa, aa_pos_s = m.groups()
            aa_pos = int(aa_pos_s) - 1
            if aa_pos < len(combo.aa):
                combo.aa[aa_pos] += ref_aa
                combo.var_label += f"_Ins:{aa_pos}{ref_aa}"
                combo.vlist.append(variant)
                combo.poslist.append(aa_pos)
            else:
                err.fail(
                    f"[ERROR] - AA - {variant.mutation_type} :: {sample} :: {transcript_id} :: position {aa_pos} out of range :: {variant.aa} :: {variant.na} :: {variant.cosmic_id}"
                )

        elif m := re.match(r"p\.(\D+)(\d+)_(\D+)(\d+)dup$", variant.aa):
            _, aa_start_s, _, aa_end_s = m.groups()
            aa_start = int(aa_start_s) - 1
            aa_end = int(aa_end_s) - 1
            ins = "".join(combo.aa[aa_start:aa_end + 1])
            if aa_end < len(combo.aa):
                combo.aa[aa_end] += ins
                combo.var_label += f"_Ins:{aa_end}{ins}"
                combo.vlist.append(variant)
                combo.poslist.append(aa_end)
            else:
                err.fail(
                    f"[ERROR] - AA - {variant.mutation_type} :: {sample} :: {transcript_id} :: position {aa_end} out of range :: {variant.aa} :: {variant.na} :: {variant.cosmic_id}"
                )

        elif m := re.match(r"p\.(\D+)(\d+)_(\D+)(\d+)ins(\D+)$", variant.aa):
            _, _, _, aa_end_s, ins = m.groups()
            aa_end = int(aa_end_s) - 1
            if aa_end < len(combo.aa):
                combo.aa[aa_end] += ins
                combo.var_label += f"_Ins:{aa_end}{ins}" 
                combo.vlist.append(variant)
                combo.poslist.append(aa_end)
            else:
                err.fail(
                    f"[ERROR] - AA - {variant.mutation_type} :: {sample} :: {transcript_id} :: position {aa_end} out of range :: {variant.aa} :: {variant.na} :: {variant.cosmic_id}"
                )

        else:
            err.fail(f"[ERROR] - {variant.mutation_type} - AA format unknwn:: {variant.aa} :: {variant.na} :: {variant.cosmic_id}")
            continue

        # Handle different patterns of in-frame insertions based on the NA change annotation, similar to the AA patterns.
        if m := re.match(r"c\.(\d+)_(\d+)dup$", variant.na):
            na_start = int(m.group(1)) - 1
            na_end = int(m.group(2)) - 1
            ins = "".join(combo.na[na_start:na_end + 1])
            if na_end < len(combo.na):
                combo.na[na_end] += ins
        elif m := re.match(r"c\.(\d+)_(\d+)ins(\D+)$", variant.na):
            na_end = int(m.group(2)) - 1
            ins = m.group(3)
            if na_end < len(combo.na):
                combo.na[na_end] += ins
        else:
            err.fail(f"[ERROR] - {variant.mutation_type} - NA format unknown:: {variant.na} :: {variant.aa} :: {variant.cosmic_id}")

    return False, "".join(combos[0].aa)

# Function to handle nonsense variants 
# by parsing the AA change annotation, verifying the reference amino acid, and writing the truncated sequence to the output FASTA file.
def handle_nonsense(combos: List[SequenceCombination], variant: VariantRecord, pos: int, out_handle, sample: str, transcript_id: str, gene: GeneInfo, tissue: str, err: ErrorWriter) -> Tuple[bool, str]:
    
    # Parse the amino acid annotation to extract the reference amino acid and position for the nonsense mutation, 
    # indicated by a pattern like "p.A123*".
    match = re.match(r"p\.(\D+)(\d+)\*$", variant.aa)
    if not match:
        err.fail(f"[ERROR] - {variant.mutation_type} - AA format unknown:: {variant.aa} :: {variant.na} :: {variant.cosmic_id}")
        return False, None
    
    ref_aa, aa_pos_s = match.groups()
    aa_pos = int(aa_pos_s) - 1 # Convert to 0-based index

    # Verify that the amino acid position is valid and matches the reference amino acid in the sequence combination. 
    if aa_pos >= len(combos[0].aa) or combos[0].aa[aa_pos] != ref_aa:
        err.fail(
            f"[ERROR] - NONSENSE AA :: {sample} :: {transcript_id} :: Positions: {pos} == {aa_pos} :: {variant.aa} :: {variant.na} :: {variant.cosmic_id} :: "
            f"does not match {ref_aa} to {combos[0].aa[aa_pos] if aa_pos < len(combos[0].aa) else 'OUT_OF_RANGE'}"
        )
        return False, None

    # Write the truncated sequence up to the nonsense mutation position to the output FASTA file, using the appropriate header format.
    outseq = ""
    for combo in combos:
        outseq = "".join(combo.aa[:aa_pos])
        write_record(
            out_handle,
            format_header(transcript_id, sample, f"{combo.var_label}_Nons:{aa_pos}", gene, tissue, combo.vlist + [variant], combo.poslist + [aa_pos]),
            outseq
        )

    return True, outseq

# Function to handle frameshift insertions
# by parsing the AA and NA change annotations, determining the affected positions, applying the frameshift insertion to the sequence combinations, and writing the resulting sequences to the output FASTA file.
def handle_frameshift_insertion( combos: List[SequenceCombination], variant: VariantRecord, pos: int, utr_seq: str, out_handle, sample: str, transcript_id: str, gene: GeneInfo, tissue: str, err: ErrorWriter) -> Tuple[bool, str]:
    
    base_na = combos[0].na.copy() + list(utr_seq) #Get the first ( unmutated) NA sequence plus utr extension to determine the inserted sequence based on the NA annotation)
    insert_pos = -1 #Initialize the position where the insertion will occur
    insert_seq = ""

    # Parse the NA annotation to determine the type of frameshift insertion and extract the relevant positions and inserted sequence.
    # The patterns include single nucleotide duplications (e.g., "c.123dup"), range duplications (e.g., "c.123_125dup"), and insertions of specific nucleotides (e.g., "c.123_125insATG").
    if m := re.match(r"c\.(\d+)dup$", variant.na):
        insert_pos = int(m.group(1)) - 1
        if insert_pos < len(base_na):
            insert_seq = base_na[insert_pos] + base_na[insert_pos]
        else:
            err.fail(
                f"[ERROR] - NA INSERT - {variant.mutation_type} :: {sample} :: {transcript_id} :: position {insert_pos} out of range :: {variant.aa} :: {variant.na} :: {variant.cosmic_id}"
            )
    elif m := re.match(r"c\.(\d+)_(\d+)dup$", variant.na):
        s = int(m.group(1)) - 1
        e = int(m.group(2)) - 1
        insert_pos = e
        insert_seq = base_na[e] + "".join(base_na[s:e + 1])
    elif m := re.match(r"c\.(\d+)_\d+ins(\D+)$", variant.na):
        insert_pos = int(m.group(1)) - 1
        ins = m.group(2)
        if insert_pos < len(base_na):
            insert_seq = base_na[insert_pos] + ins
        else:
            err.fail(
                f"[ERROR] - NA INSERT - {variant.mutation_type} :: {sample} :: {transcript_id} :: position {insert_pos} out of range :: {variant.aa} :: {variant.na} :: {variant.cosmic_id}"
            )
    elif m := re.match(r"c\.(\d+)del", variant.na):
        insert_pos = int(m.group(1)) - 1
        if insert_pos < len(base_na):
            insert_seq = ""
        else:
            err.fail(
                f"[ERROR] - NA DELETE - {variant.mutation_type} :: {sample} :: {transcript_id} :: position {insert_pos} out of range :: {variant.aa} :: {variant.na} :: {variant.cosmic_id}"
            )
    else:
        err.fail(f"[ERROR] - NA format unknown - {variant.mutation_type} :: {variant.na} :: {variant.aa} :: {variant.cosmic_id}")
        return True, None

    # Define a function to apply the frameshift insertion to a given sequence combination, 
    # update the variant label, and write the resulting sequence to the output FASTA file.
    def write_fs(combo: SequenceCombination) -> str:
        tmp = combo.na.copy() + list(utr_seq) #Get the NA sequence for the current combination plus utr extension to apply the insertion
        if insert_pos < len(tmp):
            tmp[insert_pos] = insert_seq
        else:
            err.fail(
                f"[ERROR] - NA INSERTION OUT OF RANGE - {variant.mutation_type} :: {sample} :: {transcript_id} :: position {insert_pos} out of range :: {variant.aa} :: {variant.na} :: {variant.cosmic_id}"
            )
            return None
          
        outseq = translate_na("".join(tmp))
        write_record(
            out_handle,
            format_header(transcript_id, sample, f"{combo.var_label}_InsFS:{pos}", gene, tissue, combo.vlist + [variant], combo.poslist + [pos]),
            outseq,
        )
        return outseq

    # Apply the frameshift insertion to the sequence combinations and write the results to the output FASTA file.
    oseq = ""
    for combo in combos:
        oseq = write_fs(combo)

    return True, oseq

# Function to handle frameshift deletions
# by parsing the NA change annotations, determining the affected positions, applying the frameshift deletion to the sequence combinations, and writing the resulting sequences to the output FASTA file.
def handle_frameshift_deletion( combos: List[SequenceCombination], variant: VariantRecord, pos: int, utr_seq: str, out_handle, sample: str, transcript_id: str, gene: GeneInfo, tissue: str, err: ErrorWriter ) -> Tuple[bool, str]:
    
    # Parse the annotation to determine the position of the frameshift deletion, 
    # indicated by a pattern like "p.A123fs".
    if m := re.match(r"c\.(\d+)del$", variant.na):
        na_start = na_end = int(m.group(1)) - 1
    elif m := re.match(r"c\.(\d+)_(\d+)del$", variant.na):
        na_start = int(m.group(1)) - 1
        na_end = int(m.group(2)) - 1
    else:
        err.fail(f"[ERROR] - NA format unknown - {variant.mutation_type} :: {variant.na} :: {variant.aa} :: {variant.cosmic_id}")
        return True, None

    # Define a function to apply the frameshift deletion to a given sequence combination,
    # update the variant label, and write the resulting sequence to the output FASTA file.
    def write_fs(combo: SequenceCombination) -> str:
        tmp = combo.na.copy() + list(utr_seq) #Get the NA sequence for the current combination plus utr extension to apply the deletion
        LOG.debug("Original NA sequence: %s", "".join(tmp))

        for j in range(na_start, min(na_end + 1, len(tmp))):
            LOG.debug("Deleting nucleotide at position %d: %s", j, tmp[j])
            tmp[j] = ""
        
        LOG.debug("NA sequence after deletion: %s", "".join(tmp))
        outseq = translate_na("".join(tmp))
        LOG.debug("Translated sequence after deletion: %s", outseq)

        LOG.debug("Check formmatting header with var_label: %s", combo.var_label)
        LOG.debug("...VLIST: %s", combo.vlist)
        LOG.debug("...variant: %s", variant)
        LOG.debug("...POSLIST: %s", combo.poslist)
        LOG.debug("...pos: %s", pos)

        write_record(
            out_handle,
            format_header(transcript_id, sample, f"{combo.var_label}_DelFS:{pos}", gene, tissue, combo.vlist + [variant], combo.poslist + [pos]),
            outseq,
        )
        return outseq

    # Apply the frameshift deletion to the sequence combinations and write the result to the output FASTA file.
    oseq = ""
    for combo in combos:
        oseq = write_fs(combo)

    return True, oseq

# Function to handle nonstop mutations
# by parsing the AA and NA change annotations, determining the affected positions, applying the nonstop extension via 3'UTR to the sequence combinations, and writing the resulting sequences to the output FASTA file.
def handle_nonstop( combos: List[SequenceCombination], variant: VariantRecord, pos: int, utr_seq: str, out_handle, sample: str, transcript_id: str,gene: GeneInfo, tissue: str, err: ErrorWriter) -> Tuple[bool, str]:
    
    #Function to build the extended sequence for a given sequence combination by applying the nonstop mutation and appending the 3'UTR sequence.
    def build_seq(combo: SequenceCombination) -> Optional[str]:
        tmp = combo.na.copy()

        # Parse the NA annotation to determine the type of nonstop mutation and extract the relevant position and new base.
        # The patterns include single nucleotide substitutions (e.g., "c.123A>T") and more complex changes (e.g., "c.123_125delinsATG").
        if m := re.match(r"c\.(\d+)\D+>(\D+)$", variant.na):
            na_pos = int(m.group(1)) - 1
            new_base = m.group(2)
            if na_pos < len(tmp):
                tmp[na_pos] = new_base

        elif m := re.match(r"c\.(\d+)_(\d+)delins(\D+)$", variant.na):
            start = int(m.group(1)) - 1
            end = int(m.group(2)) - 1
            ins = m.group(3)
            for j in range(start, min(end + 1, len(tmp))):
                tmp[j] = ""
            if end < len(tmp):
                tmp[end] = ins
        else:
            err.fail(f"[ERROR] - NA format unknown - {variant.mutation_type} :: {variant.na} :: {variant.aa} :: {variant.cosmic_id}")
            return None
        return "".join(tmp) + utr_seq

    # Build the extended sequence for the first sequence combination and write it to the output FASTA file. 
    # If there are additional sequence combinations (e.g., from branching due to close variants) and the current variant is within the distance threshold from the last processed variant, 
    # also build and write the extended sequence for the last combination.
    seq = ""
    for combo in combos:
        seq = build_seq(combo)
        if seq is not None:
            write_record(
                out_handle,
                format_header(transcript_id, sample, f"{combo.var_label}_NonstopEXT:{pos}", gene, tissue, combo.vlist + [variant], combo.poslist + [pos]),
                translate_na(seq),
            )
    
    return True, seq

# Function to handle complex variants that do not fit into the standard categories of missense, nonsense, in-frame indels, or frameshift mutations.
# This function includes handling for special cases such as stop-gain / nonsense insertions and complex frameshift mutations, as well as applying the appropriate changes to the sequence combinations and writing the results to the output FASTA file.
def handle_complex(combos: List[SequenceCombination], variant: VariantRecord, pos: int, utr_seq: str, out_handle, sample: str, transcript_id: str, gene: GeneInfo, tissue: str, err: ErrorWriter) -> Tuple[bool, str]:
    
    # First, check for a special case of a stop-gain / nonsense mutation indicated by an annotation like "p.A123*"
    if m := re.match(r"p\.\D+(\d+)\*$", variant.aa):
        stop_pos = int(m.group(1)) - 1
        outseq = ""
        for combo in combos:
            outseq = "".join(combo.aa[:stop_pos])
            write_record(
                out_handle,
                format_header(transcript_id, sample, f"{combo.var_label}_Comp:{stop_pos}*", gene, tissue, combo.vlist + [variant], combo.poslist + [stop_pos]),
                outseq,
            )
        return True, outseq

    # Next, check for complex frameshift mutations indicated by annotations containing "fs" in the AA change.
    if "fs*" in variant.aa:
        def build_fs(combo: SequenceCombination) -> Optional[str]:
            tmp = combo.na.copy() + list(utr_seq)
            if m := re.match(r"c\.(\d+)del$", variant.na):
                s = int(m.group(1)) - 1
                if s < len(tmp):
                    tmp[s] = ""
                for j in range(s, len(tmp)):
                    tmp[j] = tmp[j]
            elif m := re.match(r"c\.(\d+)_(\d+)del$", variant.na):
                s = int(m.group(1)) - 1
                e = int(m.group(2)) - 1
                for j in range(s, min(e + 1, len(tmp))):
                    tmp[j] = ""
                for j in range(e, len(tmp)):
                    if j >= 0:
                        tmp[j] = tmp[j]
            elif m := re.match(r"c\.(\d+)_(\d+)delins(\D+)$", variant.na):
                s = int(m.group(1)) - 1
                e = int(m.group(2)) - 1
                ins = m.group(3)
                for j in range(s, min(e + 1, len(tmp))):
                    tmp[j] = ""
                if e < len(tmp):
                    tmp[e] = ins
            else:
                err.fail(f"[ERROR] - NA format unknown - {variant.mutation_type} :: {variant.na} :: {variant.aa} :: {variant.cosmic_id}")
                return None
            return "".join(tmp)

        seq = ""
        for combo in combos:
            seq = build_fs(combo)
            if seq is not None:
                write_record(
                    out_handle,
                    format_header(transcript_id, sample, f"{combo.var_label}_CompFS:{pos}", gene, tissue, combo.vlist + [variant], combo.poslist + [pos]),
                    translate_na(seq),
                )
        
        return True, translate_na(seq)

    for combo in combos:

        # Handle different patterns of complex variants based on the AA change annotation, 
        # including single amino acid deletions, range deletions, and complex indels.
        if m := re.match(r"p\.(\D+)(\d+)del$", variant.aa):
            aa_pos = int(m.group(2)) - 1
            if aa_pos < len(combo.aa):
                combo.aa[aa_pos] = ""
            combo.var_label += f"_CompDel{aa_pos}"
            combo.vlist.append(variant)
            combo.poslist.append(aa_pos)
        elif m := re.match(r"p\.(\D+)(\d+)_(\D+)(\d+)del$", variant.aa):
            aa_start = int(m.group(2)) - 1
            aa_end = int(m.group(4)) - 1
            for j in range(aa_start, min(aa_end + 1, len(combo.aa))):
                combo.aa[j] = ""
            combo.var_label += f"_CompDel{pos}"
            combo.vlist.append(variant)
            combo.poslist.append(aa_start)
        elif m := re.match(r"p\.(\D+)(\d+)_(\D+)(\d+)delins(\D+)$", variant.aa):
            aa_start = int(m.group(2)) - 1
            aa_end = int(m.group(4)) - 1
            ins = m.group(5)
            for j in range(aa_start, min(aa_end + 1, len(combo.aa))):
                combo.aa[j] = ""
            if aa_end < len(combo.aa):
                combo.aa[aa_end] = ins
            combo.var_label += f"_CompDelIns{pos}"
            combo.vlist.append(variant)
            combo.poslist.append(aa_start)
        else:
            err.fail(f"[ERROR] - Complex AA format unmatched - {variant.mutation_type} :: {variant.aa} :: {variant.na} :: {variant.cosmic_id}")
            continue

        # Handle different patterns of complex variants based on the NA change annotation, similar to the AA patterns, 
        # including single nucleotide deletions, range deletions, and complex indels.
        if m := re.match(r"c\.(\d+)_(\d+)del$", variant.na):
            na_start = int(m.group(1)) - 1
            na_end = int(m.group(2)) - 1
            for j in range(na_start, min(na_end + 1, len(combo.na))):
                combo.na[j] = ""
        elif m := re.match(r"c\.(\d+)_(\d+)delins(\D+)$", variant.na):
            na_start = int(m.group(1)) - 1
            na_end = int(m.group(2)) - 1
            ins = m.group(3)
            for j in range(na_start, min(na_end + 1, len(combo.na))):
                combo.na[j] = ""
            if na_end < len(combo.na):
                combo.na[na_end] += ins
        else:
            err.fail(f"[ERROR] - Complex NA format unmatched - {variant.mutation_type} :: {variant.na} :: {variant.aa} :: {variant.cosmic_id}")

    return False, "".join(combos[0].aa)

# Function to process a single variant by determining its type, applying the appropriate changes to the sequence combinations, and writing the resulting sequences to the output FASTA file.
def process_variant( combos: List[SequenceCombination], variant: VariantRecord, pos: int, out_handle, err: ErrorWriter, sample: str, transcript_id: str, gene: GeneInfo, tissue: str, utr_seq: str) -> Tuple[bool, str]:
    
    #Check if the amino acid annotation can be parsed to extract the reference amino acid and position. 
    #If so, verify that the reference amino acid matches the corresponding position in the sequence combination.
    #Exclude special cases such as stop-gain mutations (e.g., "p.X123*") and complex indels (e.g., "c.123_125delinsATG") from this verification, as they may not follow the standard reference amino acid format.
    ref_aa, aa_index = None, None
    parsed = parse_variant_annotation(variant.aa)
    if parsed is not None:
        ref_aa, aa_index = parsed
        if aa_index >= len(combos[0].aa) or (combos[0].aa[aa_index] != ref_aa and not re.search(r"c.*del", variant.na, flags=re.IGNORECASE)):
            if not re.search(r"p\.\*", variant.aa, flags=re.IGNORECASE):
                lo = max(0, aa_index - 3)
                hi = min(len(combos[0].aa), aa_index + 4)
                window = "".join(combos[0].aa[lo:hi]) if aa_index < len(combos[0].aa) else "OUT_OF_RANGE"
                err.fail(
                    f"[ERROR] - Parsing AA annotation - {variant.mutation_type} :: {sample} :: {transcript_id} :: Positions: {pos} == {aa_index} :: {variant.aa} :: {variant.na} :: {variant.cosmic_id} :: "
                    f"does not match reference: {ref_aa} to position: {combos[0].aa[aa_index] if aa_index < len(combos[0].aa) else 'OUT_OF_RANGE'} :: {window}"
                )
                return False, None

    

    # Determine the type of mutation based on the mutation_type annotation and call the appropriate handler function to process the variant, 
    # apply the necessary changes to the sequence combinations, and write the results to the output FASTA file.
    mutation_type = variant.mutation_type
    emitted_special = False
    mutant_seq = None

    LOG.debug(f"Processing variant {variant.cosmic_id} with mutation type {mutation_type} for sample {sample} and transcript {transcript_id} :: {variant.aa} :: {variant.na}")

    #Handle Missense mutations
    if re.search(r"missense", mutation_type, flags=re.IGNORECASE):
        mutant_seq = handle_missense(combos, variant, pos, err, sample, transcript_id)

    #Handle Nonsense Stop-gains
    elif re.search(r"stop_gained|nonsense", mutation_type, flags=re.IGNORECASE):
        emitted_special, mutant_seq = handle_nonsense(
            combos, variant, pos, out_handle, sample, transcript_id, gene, tissue, err
        )

    #Handle Frameshift Mutations (both insertions and deletions) formatted differently for old and new cosmic annotations.
    elif re.search(r"frameshift", mutation_type, flags=re.IGNORECASE):

        if re.search(r"insertion", mutation_type, flags=re.IGNORECASE) or re.search(r"dup", variant.na, flags=re.IGNORECASE):
            emitted_special, mutant_seq = handle_frameshift_insertion(
                combos, variant, pos, utr_seq, out_handle, sample, transcript_id, gene, tissue, err
            )

        elif re.search(r"complex", mutation_type, flags=re.IGNORECASE) or re.search(r"delins", variant.na, flags=re.IGNORECASE):
            emitted_special, mutant_seq = handle_complex(
                combos, variant, pos, utr_seq, out_handle, sample, transcript_id, gene, tissue, err
            )

        elif re.search(r"deletion", mutation_type, flags=re.IGNORECASE) or re.search(r"del", variant.na, flags=re.IGNORECASE):
             emitted_special, mutant_seq = handle_frameshift_deletion(
                combos, variant, pos, utr_seq, out_handle, sample, transcript_id, gene, tissue, err
            )
        else:
            err.fail(f"[ERROR] - Frameshift type unrecognized - {variant.mutation_type} :: {variant.aa} :: {variant.na} :: {variant.cosmic_id}")
            return False, None

    #Handle In-frame Indels (both insertions and deletions) formatted differently for old and new cosmic annotations.
    elif re.search(r"inframe", mutation_type, flags=re.IGNORECASE):
        
        if re.search(r"insertion", mutation_type, flags=re.IGNORECASE):
            emitted_special, mutant_seq = handle_inframe_insertion(
                combos, variant, pos, out_handle, sample, transcript_id, gene, tissue, err
            )

        elif re.search(r"deletion", mutation_type, flags=re.IGNORECASE):

            #Handle Special Cases where although the delelion is in-frame the nucleotide change does not align with the reference amino acid sequence hence it may be more appropriate to treat it as a frameshift deletion. 
            #This is determined by checking if the reference amino acid from the annotation matches the corresponding amino acid in the sequence combination at the specified position.
            if ref_aa is not None and ref_aa == combos[0].aa[aa_index]:
                mutant_seq = handle_inframe_deletion(combos, variant, err)
            else:
                emitted_special, mutant_seq = handle_frameshift_deletion(
                    combos, variant, pos, utr_seq, out_handle, sample, transcript_id, gene, tissue, err
                )
        else:
            err.fail(f"[ERROR] - In-frame type unrecognized - {variant.mutation_type} :: {variant.aa} :: {variant.na} :: {variant.cosmic_id}")
            return False, None
        
    #Handle Nonstop / stop lost mutations
    elif re.search(r"stop_lost|nonstop", mutation_type, flags=re.IGNORECASE):
        emitted_special, mutant_seq = handle_nonstop(
            combos, variant, pos, utr_seq, out_handle, sample, transcript_id, gene, tissue, err
        )

    else:
        err.fail(f"[ERROR] - UNHANDLED VARIANT TYPE :: type = {mutation_type} :: AAannotation = {variant.aa} :: NAannotation =  {variant.na} :: Cosmic ID = {variant.cosmic_id} :: Sample = {sample} :: Transcript = {transcript_id}")
        return False, None

    modified = False
    # If the variant was not handled by a special case (e.g. nonsense, frameshift, nonstop), 
    # write the resulting sequences for all combinations to the output FASTA file using a standard header format.
    if not emitted_special:
        modified = True
        for combo in combos:
            if not combo.var_label:
                continue
            write_record(
                out_handle,
                format_header(transcript_id, sample, combo.var_label, gene, tissue, combo.vlist, combo.poslist),
                "".join(combo.aa),
            )

    return modified, mutant_seq

# Function to perform in silico tryptic digestion of a protein sequence
# Make sequence isobaric by replacing I with L, then split on K/R not followed by P, and return the resulting peptides.
# Only return peptides that are between 6 and 100 amino acids in length, which are typical for mass spectrometry analysis.
def digest_tryptic_peptides(sequence: str) -> List[str]:
    seq = sequence.replace("I", "L")
    seq = re.sub(r"([KR])", r"\1_", seq)
    seq = seq.replace("_P", "P")
    return [pep for pep in seq.split("_") if 6 <= len(pep) <= 100]

# Function to identify peptides that are affected by a mutation 
# by comparing the tryptic peptides generated from the reference and mutant protein sequences, 
# and returning the lists of peptides that are unique to each sequence.
def get_mutation_site_peptides(ref_sequence: str, mutant_sequence: str) -> Tuple[List[str], List[str]]:
    ref_peptides = set(digest_tryptic_peptides(ref_sequence))
    mut_peptides = set(digest_tryptic_peptides(mutant_sequence))
    ref_affected = sorted(ref_peptides - mut_peptides)
    mut_affected = sorted(mut_peptides - ref_peptides)
    return ref_affected, mut_affected
    
#Function to process a single sample by iterating through the associated protein variants, 
# applying the necessary changes to the sequences based on the variant annotations, 
# and writing the resulting mutant sequences to an output FASTA file.
def process_sample( sample: str, protein_variants: Dict[str, Dict[int, Dict[str, VariantRecord]]], proteins: Dict[str, str], transcripts: Dict[str, str], utr: Dict[str, str], gene_info: Dict[str, GeneInfo], tissue: str, output_prefix: Path, distance_threshold: int, mutation_tsv_writer: csv.writer, written_mutations: set[str], out_fasta: IO ) -> Tuple[int, int]:
    
    #Define the path for the error log file based on the sample name and output prefix.
    err_path = Path(f"{output_prefix}_{sample}_mutant_proteoforms.err")

    success_count = 0
    fail_count = 0
    err = ErrorWriter(err_path)

    #Iterate through the protein variants associated with the current sample
    for protein_count, transcript_id in enumerate(protein_variants, start=1):
        if protein_count % 1000 == 0:
            LOG.info("%s: processed %d transcripts", sample, protein_count)

        if transcript_id not in proteins or transcript_id not in transcripts:
            continue

        aa_seq = list(proteins[transcript_id])
        na_seq = list(transcripts[transcript_id])
        ref_aa_sequence = "".join(aa_seq)
        gene = gene_info.get(transcript_id, GeneInfo())
        positions = sorted(protein_variants[transcript_id])

        #Iterate through the variant positions for the current transcript and process each variant
        for index, pos in enumerate(positions):
            err.seen()
            success_count += 1
            
            LOG.debug(f"Processing sample {sample} :: transcript {transcript_id} :: position {pos} :: variants at this position {list(protein_variants[transcript_id][pos].values())}")
            
            #Iterate through the variants at the current position for the transcript 
            for variant in protein_variants[transcript_id][pos].values():

                dist = 0
                combo = [SequenceCombination(aa=aa_seq.copy(), na=na_seq.copy(), var_label="", vlist=[], poslist=[])]

                #Check if the current variant position is within the length of the amino acid sequence for the transcript.
                if pos > len(aa_seq) and not re.search(r"p\.\*", variant.aa, flags=re.IGNORECASE):
                    err.fail(
                        f"[ERROR] - AA OUT OF RANGE - {variant.mutation_type} :: Skipping Variant :: Current Variant Position {pos} is outside of protein {transcript_id} with length {len(aa_seq)} :: {sample} :: {variant.cosmic_id} :: {variant.aa} :: {variant.na} :: {aa_seq[-1]}"
                    )
                    fail_count += 1
                    continue

                before_failures = err.fail_count


                LOG.debug(f"Processing variant {variant.cosmic_id} at position {pos} for transcript {transcript_id} in sample {sample} with mutation type {variant.mutation_type} :: AA change {variant.aa} :: NA change {variant.na}")
                LOG.debug(f"Current sequence combination before processing variant: {''.join(combo[0].aa)} :: {''.join(combo[0].na)} :: var_label: {combo[0].var_label} :: vlist: {combo[0].vlist} :: poslist: {combo[0].poslist}")

                #Process the variant to apply the necessary changes to the sequence combinations 
                # write the resulting mutant sequences to the output FASTA file, and return the mutant combo
                m, mutant_sequence = process_variant(
                    combos=combo,
                    variant=variant,
                    pos=pos,
                    out_handle=out_fasta,
                    err=err,
                    sample=sample,
                    transcript_id=transcript_id,
                    gene=gene,
                    tissue=tissue,
                    utr_seq=utr.get(transcript_id, ""),
                )

                #Before processing the variant, build a unique identifier for the mutation sample 
                #Check if it has already been written to the combined TSV file.
                unique_mutation_sample_id = build_mutation_sample_id(sample, transcript_id, pos, variant)
                if unique_mutation_sample_id not in written_mutations:

                    written_mutations.add(unique_mutation_sample_id)

                    #Find the tryptic peptides that are affected by the mutation
                    ref_peptides: List[str] = []
                    mut_peptides: List[str] = []
                    LOG.debug(f"MUTANT SEQUENCE:: {mutant_sequence} :: for variant {variant.cosmic_id} :: {variant.aa} :: {variant.na}")
                    
                    if mutant_sequence is None:
                        #err.fail(f"[ERROR] - Failed to build mutant sequence for TSV output - {variant.mutation_type} :: {sample} :: {transcript_id} :: {variant.aa} :: {variant.na} :: {variant.cosmic_id}")
                        LOG.debug(f"[ERROR] - Failed to build mutant sequence for TSV output - {variant.mutation_type} :: {sample} :: {transcript_id} :: {variant.aa} :: {variant.na} :: {variant.cosmic_id}")

                    else:
                        ref_peptides, mut_peptides = get_mutation_site_peptides( ref_aa_sequence, mutant_sequence )

                        # Write the mutation information to the combined TSV file
                        mutation_tsv_writer.writerow(
                            [
                                unique_mutation_sample_id,
                                gene.ensg,
                                transcript_id,
                                gene.gene_name,
                                sample,
                                variant.mutation_type,
                                tissue,
                                variant.zygosity,
                                variant.cosmic_id,
                                pos,
                                variant.aa,
                                variant.na,
                                ";".join(ref_peptides),
                                ";".join(mut_peptides),
                            ]
                        )
                else:
                    err.warn(f"[WARNING] - Duplicate mutation sample ID {unique_mutation_sample_id} already written to TSV, skipping TSV write for this variant.")

                #After processing the variant, check if any new failures were recorded in the error log and update the failure count accordingly.
                if err.fail_count > before_failures:
                    fail_count += 1


                #If there are multiple variants in the transcript, check the distance to the next variant and if it is within the specified threshold, 
                # process the next variant as well to account for potential interactions between nearby variants.
                nindex = index + 1
                if nindex >= len(positions):
                    continue

                npos = positions[nindex]
                dist = npos - pos
                while dist < distance_threshold and nindex + 1 < len(positions):

                    LOG.debug(f"Next variant at position {npos} is within distance threshold {distance_threshold} of current position {pos} (distance {dist}), processing next variant for transcript {transcript_id} in sample {sample}")

                    if dist > 0:

                        #Create a store for modified sequence combinations after processing the next variant, which will be used to branch the sequence combinations if there are multiple variants at the same position.
                        current_combos = []

                        for vindex, nvariant in enumerate(protein_variants[transcript_id][npos].values()):
                                    
                            #Before processing the next variant, clone the current sequence combinations to preserve their state for potential branching if there are multiple variants at the same position.
                            tmp_combos = [c.clone() for c in combo]

                            #Process the next variant to apply the necessary changes to the current sequence comb 
                            #write the resulting mutant sequences to the output FASTA file
                            modifed, mseq = process_variant(
                                combos=tmp_combos,
                                variant=nvariant,
                                pos=npos,
                                out_handle=out_fasta,
                                err=err,
                                sample=sample,
                                transcript_id=transcript_id,
                                gene=gene,
                                tissue=tissue,
                                utr_seq=utr.get(transcript_id, ""),
                            )

                            #If the next variant modified the sequence combinations, update the current_combos with the modified combinations for further processing of subsequent variants.
                            if modifed:
                                LOG.debug(f"Variant at position {npos} modified the sequence combinations, updating current_combos for transcript {transcript_id} in sample {sample}")
                                for c in tmp_combos:
                                    current_combos.append(c.clone())

                        #If there were modifications from processing the next variant, update the main combo list with the current_combos to ensure that subsequent variants are processed with the most up-to-date sequence combinations.     
                        combo = current_combos if current_combos else combo

                    nindex += 1
                    if nindex >= len(positions):
                        break
                    dist = positions[nindex] - pos

    err.close()
    return success_count, fail_count

# Main function 
# #orchestrate the processing of COSMIC variants, including reading input data, processing each cell line and its associated variants, 
# and writing the results to output files.
def main() -> int:

    #Parse command-line arguments and configure logging based on the specified log level.
    args = parse_args()
    configure_logging(args.log_level)

    #Parse sample information from the provided sample info file and store it in a dictionary for easy access during processing.
    sample_info = parse_sample_info(args.sample_info)

    requested_samples = set()
    for experiment in sample_info:
        requested_samples.update(sample_info[experiment].keys())


    LOG.info("Reading COSMIC data")
    sampledata, tissue_map, selected_ids, required_utr = load_cosmic_variants( args.cosmic_tsv, requested_samples, args.lock_transcript_version, args.old_cosmic )

    LOG.info("Loading protein sequences")
    proteins, gene_info = load_protein_sequences(args.protein_fasta, args.lock_transcript_version, selected_ids)

    LOG.info("Loading transcript sequences")
    transcripts, utr = load_transcript_sequences(args.transcript_fasta, args.lock_transcript_version, selected_ids, required_utr)

    #Check all selected transcript IDs have corresponding protein and transcript sequences, and log any missing sequences.
    unmatched_ids = set()
    for tid in selected_ids:
        if tid not in proteins:
            LOG.debug(f"Selected transcript ID {tid} not found in protein FASTA")
            unmatched_ids.add(tid)
        if tid not in transcripts:
            LOG.debug(f"Selected transcript ID {tid} not found in transcript FASTA")
            unmatched_ids.add(tid)

    matched_id_count = len(selected_ids) - len(unmatched_ids)
    LOG.info(f"{matched_id_count} out of {len(selected_ids)} transcript IDs with COSMIC mutation have matching sequences in the FASTA files")

    combined_tsv_path = Path(f"{args.output_prefix}_mapped_mutations.tsv")

    with combined_tsv_path.open("w") as TSV:
        combined_tsv_writer = csv.writer(TSV, delimiter="\t")
        combined_tsv_writer.writerow(
            [
                "UNIQUE_MUTATION_ID",
                "GENE_ID",
                "TRANSCRIPT_ID",
                "GENE_NAME",
                "SAMPLE_NAME",
                "MUTATION_TYPE",
                "SAMPLE_TISSUE",
                "MUTATION_ZYGOSITY",
                "COSMIC_ID",
                "AA_POSITION",
                "AA_CHANGE",
                "NA_CHANGE",
                "OBSERVABLE_WT_PEPTIDES",
                "OBSERVABLE_MUTANT_PEPTIDES",
            ]
        )

        written_mutations: set[str] = set()

        total_seen = 0
        total_failed = 0

        for experiment in sorted(sample_info):
            LOG.info("Processing mutations for experiment %s", experiment)

            #Define the path for the experiment output FASTA file.
            fasta_path = Path(f"{args.output_prefix}_{experiment}_mutant_proteoforms.fasta")
            with fasta_path.open("w", encoding="utf-8") as OUTFASTA:

                for sample in sorted(sample_info[experiment]):
                    if sample not in sampledata:
                        LOG.warning("Sample %s not found in COSMIC data, skipping", sample)
                        continue
            
                    LOG.info("Current Sample == %s", sample)
                    seen, failed = process_sample(
                        sample=sample,
                        protein_variants=sampledata[sample],
                        proteins=proteins,
                        transcripts=transcripts,
                        utr=utr,
                        gene_info=gene_info,
                        tissue=tissue_map.get(sample, "NA"),
                        output_prefix=args.output_prefix,
                        distance_threshold=args.variant_distance,
                        mutation_tsv_writer=combined_tsv_writer,
                        written_mutations=written_mutations,
                        out_fasta=OUTFASTA,
                    )
                    total_seen += seen
                    total_failed += failed

    LOG.info("FAILED COSMIC MUTANTION MAPPINGS == %d out of %d", total_failed, total_seen)
    return 0


# Entry point of the script
if __name__ == "__main__":
    raise SystemExit(main())
