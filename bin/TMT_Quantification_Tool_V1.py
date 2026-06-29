#!/usr/bin/python3
# PIPPA TMT Quantifier
# Dr James Wright 2026

#used to get cmd line arguments 
import argparse, csv, re, glob, statistics
from typing import Dict, Set
import logging, sys
from dataclasses import dataclass

# =============================================================================
# Define data class
# =============================================================================

#Spectrum information data class
@dataclass
class SpectrumInfo:
    Peptide = str('')
    PrecursorScanID: str = ''
    PrecursorIntensity: float = None
    Charge: int = None
    SelectedMZ: float = None
    IsoWindow: float = None
    Title: str = ''
    BasePeakMZ: float = None
    BasePeakIntensity: float = None
    RawFile: str = ''
    RawScan: str = ''
    Info: str = ''
    MatchMZ: float = None
    RT: float = None
    qvalue: float = None
    score: float = None

class PeptideInfo:
    def __init__(self):
        self.Spectra: Dict[str, float] = {}
        self.rawTMT: Dict[int, float] = {}
        self.normTMT: Dict[int, float] = {}
        self.sumTMT: float = 0


# =============================================================================
# Get user parameters
# =============================================================================

# Function to load command line arguments
def parse_args() -> argparse.Namespace:
    #Read command line arguments and create help documentation using argparse
    parser = argparse.ArgumentParser(
        description='''TMT Quantification Tool - Quantifies TMT intensities for MS2 or MS3 experiments. 
        Maps ConsensusXML from OpenMS IsoBaric_Analyzer tool TMT intensities to PSMs, Peptides and Proteins in PIPPA. Outputs Spectral Tables with Normalised and Scaled TMT intensities.
        ---- James Wright (2025) - The Institute of Cancer Research, London ----
        ''')

    parser.add_argument('--psmtable', '-i', dest='psmtab', help='Txt File Containing PIPPA PSM Results', required=True)
    parser.add_argument('--peptidetable', '-p', dest='peptab', help='Txt File Containing PIPPA Peptide Results', required=True)
    parser.add_argument('--consensusxml', '-c', dest='conxmldir', help='Consensus XML directory Containing cxml files with OpenMS TMT intensities - Output for isobaric quant tool', required=True)
    parser.add_argument('--mzml', '-s', dest='mzmldir', help='mzML spectra directory', required=True)
    parser.add_argument('--proteintable', '-r', dest='protab', help='Txt File Containing PIPPA Protein Results (Optional)', required=False)
    parser.add_argument('--genetable', '-g', dest='genetab', help='Txt File Containing PIPPA Gene Level Results (Optional)', required=False)
    parser.add_argument('--fdr_threshold', '-f', type=float, dest='fdr_threshold', default=0.01, help='FDR threshold for PSM quantification (Default: 0.01)')
    parser.add_argument('--score_threshold', '-t', type=float, dest='score_threshold', default=0.01, help='PEP Score threshold for PSM quantification (Default: 0.01)')
    parser.add_argument('--tmt_threshold', '-x', type=float, dest='tmt_threshold', default=0.05, help='TMT channel intensity threshold for filtering low intensity channels from quantification (Default: 0.05)')
    parser.add_argument('--mslevel', '-m', type=int, dest='mslevel', default=2, choices=[2,3], help='Set MS TMT Quant Level (MS2=2, MS3=3). Default=2')
    parser.add_argument('--prefix', '-o', dest='pre', default='', help='Set output files prefix. Default=')
    parser.add_argument('--log_level', '-l', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], help='Logging level (Default: INFO)')
    
    return parser.parse_args()

# =============================================================================
# Logging
# =============================================================================

#Setup Logging Handle
LOG = logging.getLogger("tmt_quantification_tool")

# Function to configure logging based on the specified log level
def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="[%(levelname)s] %(message)s",
        stream=sys.stdout,
    )

# =============================================================================
# MZML PARSER and Spectral Table Builder
# =============================================================================

def parseMzML (spectrum_table: Dict[str, Dict[int, Dict[str, SpectrumInfo]]], sample: str, file: str) -> None:

    currentSpecID = ''
    currentSType = 0
    spectrum_count = 0

    LOG.info(f"Parsing mzML file: {file} for sample: {sample}")

    #Open Spectra 
    with open(file, 'r') as ML:

        #loop each line in the file
        for line in ML:
            
            #MATCH SPECTRUM ID HEADER AND EXTRACT CURRENT SPECTRUM ID
            if "<spectrum id" in line:
                m = re.search('<spectrum\s+id="spectrum=(\d+)"\s+index', line)
                if m:
                    currentSpecID = m.group(1)
                    currentSType = 0

                #<spectrum id="controllerType=0 controllerNumber=1 scan=25" index="24" defaultArrayLength="4">
                elif 'spectrum id="controllerType=' in line:
                    m = re.search('scan=(\d+)"', line)
                    if m:
                        currentSpecID = m.group(1)
                        currentSType = 0

            #MATCH MS LEVEl OF SPECTRUM
            #<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="3" />
            elif "ms level" in line:
                m = re.search('value="(\d+)"', line)
                if m:
                    currentSType = int(m.group(1))

            #Filter for MS2 and MS3 spectra
            if currentSType > 1:
                
                if sample not in spectrum_table:
                    spectrum_table[sample] = {}
                    
                if currentSType not in spectrum_table[sample]:
                    spectrum_table[sample][currentSType] = {}

                if currentSpecID not in spectrum_table[sample][currentSType]:
                    spectrum_table[sample][currentSType][currentSpecID]=SpectrumInfo()
                    spectrum_count += 1
                    
                if "base peak m/z" in line:
                    m = re.search('value="([^"]+)"', line)
                    if m:
                        spectrum_table[sample][currentSType][currentSpecID].BasePeakMZ = float(m.group(1))

                elif "base peak intensity" in line:
                    m = re.search('value="([^"]+)"', line)
                    if m:
                        spectrum_table[sample][currentSType][currentSpecID].BasePeakIntensity = float(m.group(1))

                elif "spectrum title" in line:
                    m = re.search('value="([^"]+)"', line)
                    if m:
                        spectrum_table[sample][currentSType][currentSpecID].Title = m.group(1)

                    m = re.search('File:&quot;(.*).raw&quot;.*scan=(.*)&quot', line)
                    if m:   
                        spectrum_table[sample][currentSType][currentSpecID].RawFile = m.group(1)
                        spectrum_table[sample][currentSType][currentSpecID].RawScan = m.group(2)   

                elif "filter string" in line:
                    m = re.search('value="([^"]+)"', line)
                    if m:
                        spectrum_table[sample][currentSType][currentSpecID].Info = m.group(1)   

                    m = re.search('ms\d+\s+(\S+)@[hcd|cid]', line)
                    #<userParam name="filter string" type="xsd:string" value="ITMS + c NSI r d Full ms2 484.8442@cid35.00 [128.0000-980.0000]"/>
                    #<userParam name="filter string" type="xsd:string" value="FTMS + p NSI sps d Full ms3 484.8442@cid35.00 407.0459@hcd45.00 [120.0000-140.0000]"/>
                    if m:
                        spectrum_table[sample][currentSType][currentSpecID].MatchMZ = str(m.group(1))

                elif "scan start time" in line:
                    m = re.search('value="([^"]+)"', line)
                    if m:
                        spectrum_table[sample][currentSType][currentSpecID].RT = round(float(m.group(1)), 7)  

                elif "precursor spectrumRef=" in line:
                    m = re.search('scan=([^"]+)"', line)
                    if m:
                        spectrum_table[sample][currentSType][currentSpecID].PrecursorScanID = m.group(1)  

                elif "isolation window target m/z" in line:
                    m = re.search('value="([^"]+)"', line)
                    if m:
                        spectrum_table[sample][currentSType][currentSpecID].IsoWindow = m.group(1)  

                elif "selected ion m/z" in line:
                    m = re.search('value="([^"]+)"', line)
                    if m:
                        spectrum_table[sample][currentSType][currentSpecID].SelectedMZ = float(m.group(1))   

                elif "charge state" in line:
                    m = re.search('value="([^"]+)"', line)
                    if m:
                        spectrum_table[sample][currentSType][currentSpecID].Charge = int(m.group(1))  

                elif "peak intensity" in line:
                    m = re.search('value="([^"]+)"', line)
                    if m:
                        spectrum_table[sample][currentSType][currentSpecID].PrecursorIntensity = float(m.group(1)) 

    LOG.info(f"Finished parsing mzML file: {file} for sample: {sample} - Total spectra parsed: {spectrum_count}")
            

# =============================================================================
# Spectral Table Indexing - Create indexes for fast mapping of spectra to TMT values
# =============================================================================
def indexSTable (spectrum_table: Dict[str, Dict[int, Dict[str, SpectrumInfo]]], sample: str) -> Dict[str, Dict[str, str]]:
    idx = {}
    for mstype in spectrum_table[sample]:
        for spectrumid, spectruminfo in spectrum_table[sample][mstype].items():
            if spectruminfo.RawFile and spectruminfo.RawScan:
                if spectruminfo.RawFile not in idx:
                    idx[spectruminfo.RawFile] = {}

                idx[spectruminfo.RawFile][spectruminfo.RawScan] = spectrumid

            else:
                if mstype == 2:
                    if spectruminfo.MatchMZ:
                        if spectruminfo.MatchMZ not in idx:
                            idx[str(spectruminfo.MatchMZ)] = {}
                        idx[str(spectruminfo.MatchMZ)][spectrumid] = 1
                        LOG.debug(f"DEBUG-indexSTable: MS2 Spectrum ID={spectrumid} has no RawFile/RawScan information for indexing. Using MatchMZ={spectruminfo.MatchMZ} for indexing.")
                    else:
                        LOG.warning(f"WARNING-indexSTable: MS2 Spectrum ID={spectrumid} has no RawFile/RawScan or MatchMZ information for indexing. This spectrum will not be matched to TMT values using spectral matching, only RT matching.")
                else:
                    if not spectruminfo.MatchMZ:
                        LOG.warning(f"WARNING-indexSTable: MS3 Spectrum ID={spectrumid} has no RawFile/RawScan or MatchMZ information for indexing. This spectrum will not be matched to TMT values using spectral matching, only RT matching.")

                       
    return idx


# =============================================================================
# RT Indexing - Create index for mapping RT to spectrum ID for MS spectra (This is used to map TMT values to spectra using RT matching)
# =============================================================================

def indexRT (spectrum_table: Dict[str, Dict[int, Dict[str, SpectrumInfo]]], sample: str) -> Dict[float, str]:
    idx = {}
    for mstype in spectrum_table[sample]:
        for spectrumid, spectrum in spectrum_table[sample][mstype].items():
            if spectrum.RT is None:
                LOG.warning(f"WARNING-indexRT: Spectrum ID={spectrumid} has no RT information for indexing. This spectrum will not be matched to TMT values using RT matching.")
                continue

            if spectrum.RT in idx:
                LOG.warning(f"WARNING-indexRT: RT={spectrum.RT} for Spectrum ID={spectrumid} is already indexed to Spectrum ID={idx[spectrum.RT]}. This may lead to incorrect mapping of TMT values to spectra using RT matching.")
            else:
                idx[spectrum.RT] = spectrumid
    return idx

# =============================================================================
# ConsensusXML Parser - Parse TMT intensities from consensusXML files and map to spectra in spectral table 
# =============================================================================

"""
EXAMPLE CONSENSUSXML TMT XML:

<consensusElement id="e_1603513335834568183" quality="0.0" charge="2">
    <centroid rt="795.8193359375" mz="405.830932617188012" it="9023.494141"/>
    <groupedElementList>
        <element map="0" id="2260" rt="795.8193359375" mz="126.127725999999996" it="781.687134"/>
        <element map="1" id="2260" rt="795.8193359375" mz="127.124761000000007" it="583.309021"/>
        <element map="2" id="2260" rt="795.8193359375" mz="127.131080999999995" it="0.0"/>
        <element map="3" id="2260" rt="795.8193359375" mz="128.128116000000006" it="1385.418457"/>
        <element map="4" id="2260" rt="795.8193359375" mz="128.134435999999994" it="0.0"/>
        <element map="5" id="2260" rt="795.8193359375" mz="129.131471000000005" it="829.05127"/>
        <element map="6" id="2260" rt="795.8193359375" mz="129.137789999999995" it="0.0"/>
        <element map="7" id="2260" rt="795.8193359375" mz="130.134825000000006" it="1431.373779"/>
        <element map="8" id="2260" rt="795.8193359375" mz="130.141144999999995" it="0.0"/>
        <element map="9" id="2260" rt="795.8193359375" mz="131.138180000000006" it="4012.653809"/>
        <element map="10" id="2260" rt="795.8193359375" mz="131.144499999999994" it="0.0"/>
    </groupedElementList>
    <UserParam type="float" name="precursor_purity" value="0.246882540297846"/>
    <UserParam type="string" name="scan_id" value="spectrum=5157"/>
    <UserParam type="float" name="precursor_intensity" value="0.0"/>
</consensusElement>
"""

def parseConsensusXML(spectrum_map: Dict[str, Dict[str, Dict[int, float]]], spectrum_table: Dict[str, Dict[int, Dict[str, SpectrumInfo]]], rt_index: Dict[float,str], spectrum_index: Dict[str, Dict[str, str]], sample: str, file: str, tmt_channels: Set[int], mslevel: int) -> None:

    LOG.info(f"Parsing ConsensusXML file: {file} for sample: {sample} and MS level: {mslevel}")
    
    if sample not in spectrum_table:
        LOG.warning(f"WARNING: Sample {sample} not found in spectral table. TMT values from ConsensusXML file {file} will not be mapped to spectra for this sample.")
        return
    
    if sample in spectrum_map:
        LOG.warning(f"WARNING: Sample {sample} already has TMT values mapped from a previous ConsensusXML file. TMT values from ConsensusXML file {file} will be added to existing mapped values for this sample, which may lead to incorrect mapping of TMT values to spectra if multiple ConsensusXML files are used for the same sample.")
    else:
        spectrum_map[sample] = {}
    

    #Open ConsensusXML 
    with open(file, 'r') as CX:

        element_map = {}
        element_map["elements"] = {}

        #loop each line in the file
        for line in CX:

            #Match element XML line with TMT event
            if "<element map=" in line:
                m = re.search('<element map="(\d+)" id=".*" rt="([^"]+)" mz=.*it="([^"]+)"', line)
                if m:

                    #Extract key values
                    TMTchannel = int(m.group(1)) + 1
                    rt = round(float(m.group(2)),7)
                    intensity = float(m.group(3))

                    tmt_channels.add(TMTchannel)
                    element_map["elements"][TMTchannel] = {"RT": rt, "Intensity": intensity}

            if "<UserParam" in line:
                m = re.search('<UserParam\s+type="([^"]+)"\s+name="([^"]+)"\s+value="([^"]+)"', line)
                if m:
                    param_name = m.group(2)
                    param_value = m.group(3)
                    element_map[param_name] = param_value

            if "</consensusElement>" in line:
                    
                matchRT = True
                if "scan_id" in element_map and mslevel == 2:
                    msID = element_map["scan_id"].replace("spectrum=","")
                    if msID in spectrum_table[sample][mslevel]:
                        matchRT = False
                        #Initialise spectral TMT map dictionary sMap if not previously done
                        if msID not in spectrum_map[sample]:
                            spectrum_map[sample][msID] = {}
                        #Map each TMT channel to MS2 spectrum ID  
                        for TMTchannel in element_map["elements"]:
                            LOG.debug(f"parseConsensusXML - SCAN Mapping TMT channel {TMTchannel} intensity {element_map['elements'][TMTchannel]['Intensity']} to spectrum ID {msID} for sample {sample} with RT {rt}")
                            spectrum_map[sample][msID][TMTchannel] = element_map["elements"][TMTchannel]["Intensity"]
                    else:
                        LOG.warning(f"WARNING: MS2 Spectrum ID={msID} not found in spectal table. Defaulting to RT matching.")

                if matchRT:
                    #Match Retention time of TMT peak to Retention time of MS2/MS3 Spectrum (Note: RT has be limited to 7 decimal places)
                    #Check if RT is matched rounding error can make RT fail to match, try fuzzy match if not found
                    if rt not in rt_index:
                        rtm = round( (rt / 60), 7)
                        if rtm in rt_index:
                            rt_index[rt] = rt_index[rtm]
                        else:
                            LOG.warning(f"WARNING_2: RT={rt} not found in RT Index. TMTchan:{TMTchannel} Int:{intensity} - trying fuzzy RT match with RT Index")
                            rtmatch = 0
                            for rtx in rt_index:
                                if abs(rtx-rt) <= 0.000001:
                                    LOG.info(f"RTx matched:{rtx} to {rt}:: xtimes:{rtmatch}")
                                    rtmatch = rtx
                            if rtmatch == 0:
                                bestdiff = 1000
                                for rtx in rt_index:
                                    diff = abs(rtx-rtm)
                                    if diff <= 0.000001:
                                        LOG.info(f"RTx matched:{rtx} to {rt}:: xtimes:{rtmatch}")
                                        rtmatch = rtx
                                    if diff < bestdiff:
                                        bestdiff = diff
                                        bestmatch = rtx
                            if rtmatch == 0:
                                LOG.warning(f"WARNING-21: unable to match {rt} or {rtm} with RT in RT Index")
                                LOG.info(f"BEST DIFF == {bestdiff} for RT:{bestmatch}")
                            else:
                                rt_index[rt] = rt_index[rtmatch]

                    if rt not in rt_index:
                        LOG.warning(f"WARNING: RT={rt} still not found in RT Index. TMTchan:{TMTchannel} Int:{intensity} - TMT values for this element will not be mapped to spectra using RT matching.")
                    else:
                        #Set spectrum scan ID of MS spectrum for this TMT element
                        msID = rt_index[rt]

                        if msID not in spectrum_table[sample][mslevel]:
                            LOG.warning(f"WARNING: RT={rt} matched to scan ID:{msID} but this scan ID not found in spectral table for MS level:{mslevel}")
                        else:
                            spectrum = spectrum_table[sample][mslevel][msID]
                            if msID not in spectrum_map[sample]:
                                spectrum_map[sample][msID] = {}
                            for TMTchannel in element_map["elements"]:
                                #Map each TMT channel to MS3 spectrum ID
                                LOG.debug(f"parseConsensusXML - RT Mapping TMT channel {TMTchannel} intensity {element_map['elements'][TMTchannel]['Intensity']} to spectrum ID {msID} for sample {sample} with RT {rt}")
                                spectrum_map[sample][msID][TMTchannel] = element_map["elements"][TMTchannel]["Intensity"]
                            if mslevel == 3:
                                ms3ID = msID    
                                ms2ID = "" 
                                #Use spectral index to extract MS2 scan ID  use MS3 Raw file name and Precusor Scan ID if available or use MzMatching
                                if spectrum.RawFile and spectrum.PrecursorScanID:
                                    if spectrum.RawFile in spectrum_index and spectrum.PrecursorScanID in spectrum_index[spectrum.RawFile]:
                                        ms2ID = spectrum_index[spectrum.RawFile][spectrum.PrecursorScanID]
                                if not ms2ID:
                                    LOG.debug(f"parseConsensusXML - MS3 ({ms3ID}) MatchMZ ({spectrum.MatchMZ}) not found in IDX using RawFile and PrecursorScanID. Trying to find MS2 scan using MatchMZ in IDX")
                                    if spectrum.MatchMZ in spectrum_index:
                                        diff = 1000
                                        for ms2scan in spectrum_index[spectrum.MatchMZ]:
                                            sd = int(ms3ID) - int(ms2scan)
                                            LOG.debug(f"parseConsensusXML - MS3 ({ms3ID}) MatchMZ ({spectrum.MatchMZ}) found in IDX. Checking MS2 scan {ms2scan} with difference {sd}")
                                            if sd > 0 and sd < diff:
                                                diff = sd
                                                ms2ID = ms2scan

                                        if ms2ID == '':
                                            LOG.warning(f"WARNING_81: MS3 ({ms3ID}) MatchMZ ({spectrum.MatchMZ}) found in IDX but no close matching MS2 for this spectrum! RT:{rt} TMTchan:{TMTchannel} Int:{intensity}")
                                    else:
                                        LOG.warning(f"WARNING_82: MS3 ({ms3ID}) MatchMZ ({spectrum.MatchMZ}) not found in IDX no matching MS2 for this spectrum! RT:{rt} TMTchan:{TMTchannel} Int:{intensity}")     
                                if not ms2ID:
                                    LOG.warning("WARNING-83: Unable to find MS2ID for MS3ID ("+ms3ID+") using RawFile and PrecursorScanID or MatchMZ")
                                else:
                                    #Initialise spectral TMT map dictionary sMap if not previously done
                                    if ms2ID not in spectrum_map[sample]:
                                        spectrum_map[sample][ms2ID] = {}
                                    #Map each TMT channel to MS2 spectrum ID
                                    for TMTchannel in element_map["elements"]:    
                                        spectrum_map[sample][ms2ID][TMTchannel] = element_map["elements"][TMTchannel]["Intensity"]

                    element_map = {}
                    element_map["elements"] = {}
                    


# =============================================================================
# Write Spectral Table with TMT values to file
# =============================================================================

def writeSTable(spectrum_table: Dict[str, Dict[int, Dict[str, SpectrumInfo]]], spectrum_map: Dict[str, Dict[str, Dict[int, float]]], prefix: str, tmt: Set[int]) -> None:

    outfile = prefix + "SpectralTable.tsv"
    LOG.info(f'Writing Spectral Table to file: {outfile}')

    OUT = open(outfile, 'w')
    headers = [
        "Experiment_Sample_ID",
        "Spectrum_ID",
        "MS_Level",
        "Base_Peak_MZ",
        "Base_Peak_Intensity",
        "Spectrum_Title",
        "Raw_File",
        "Raw_Scan",
        "Spectrum_Info",
        "Match_MZ",
        "RT",
        "Prescursor_Scan_ID",
        "Isolation_Window",
        "Selected_MZ",
        "Charge_State",
        "Precursor_Intensity",
        "Peptide",
        "q-value",
        "PEP"
    ]

    for t in sorted(tmt):
        headers.append(f"Raw_TMT_Intensity_{t}")
        
    OUT.write("\t".join(headers) + "\n")

    for sample in spectrum_table:
        for mstype in spectrum_table[sample]:
            for spectrumid, spectrum in spectrum_table[sample][mstype].items():

                line = [
                    sample,
                    spectrumid,
                    str(mstype),
                    str(spectrum.BasePeakMZ) if spectrum.BasePeakMZ else '-',
                    str(spectrum.BasePeakIntensity) if spectrum.BasePeakIntensity else '-',
                    spectrum.Title,
                    str(spectrum.RawFile) if spectrum.RawFile else '-',
                    str(spectrum.RawScan) if spectrum.RawScan else '-',
                    str(spectrum.Info) if spectrum.Info else '-',
                    str(spectrum.MatchMZ) if spectrum.MatchMZ else '-',
                    str(spectrum.RT) if spectrum.RT else '-',
                    str(spectrum.PrecursorScanID) if spectrum.PrecursorScanID else '-',
                    str(spectrum.IsoWindow) if spectrum.IsoWindow else '-',
                    str(spectrum.SelectedMZ) if spectrum.SelectedMZ else '-',
                    str(spectrum.Charge) if spectrum.Charge else '-',
                    str(spectrum.PrecursorIntensity) if spectrum.PrecursorIntensity else '-',
                    spectrum.Peptide if spectrum.Peptide else '-',
                    str(spectrum.qvalue) if spectrum.qvalue else '-',
                    str(spectrum.score) if spectrum.score else '-'
                ]
                for t in sorted(tmt):
                    if sample in spectrum_map and spectrumid in spectrum_map[sample] and t in spectrum_map[sample][spectrumid]:
                        line.append(str(spectrum_map[sample][spectrumid][t]))
                    else:
                        line.append('-')
                OUT.write("\t".join(line) + "\n")

    OUT.close()

# =============================================================================
# Parse PSMs and map to spectra in spectral table 
# =============================================================================
#Example Formats:
#0Sample  1Spectrum        2Peptide 3Type    4Engine  5PEPscore        6FDR     7Delta   8UnambiguousMatch        9Conflict        10Annotation      11Database        12Search  13Charge  14Proteins
#Prot_04_V1_Prot_04_     Prot_04_V1_Prot_04__470685      LLTASL  REF     mergedresults,comet,msgfplus    0.0338236       0.00413514      99      True            .(TMT6plex)LLTASL,n[+229.162932000000012]LLTASL Cosmic_Bridge_CellLines,gencode.v38.pc_translations_cRap,Cosmic_Prot_04_CellLines       V1      2       sp_Q9NXF7_DCA16_HUMAN,ENSP00000371682.1_ENST00000382247.6_ENSG00000163257.11_OTTHUMG00000128536.4_OTTHUMT00000250371.2_DCAF16-201_DCAF16_216

#Experiment	SpectrumID	PeptideSequence	Database	SearchEngine	Score(PEP)	FDR	DeltaScoreToBest	BestPeptide	ConflictingPeptidesAndScores	RetentionTime	MZ	Charge
#Lumos1	Lumos1_UPS1+UPS2+Ecoli_TMT10plex_20210604_R1_28153	PYTDYVVGSDQLLQESEDFFTLLESHEGKPLK	Human_Ecoli	mergedresults	0.1039	99.0000	-0.0839	True		4280.20	1458.7648	3

def parsePSMs(spectrum_table: Dict[str, Dict[int, Dict[str, SpectrumInfo]]], spectrum_map: Dict[str, Dict[str, Dict[int, float]]], file: str, tmt: Set[int], fdr_threshold: float =0.01, score_threshold: float =0.01) -> Dict[str, Dict[str, Dict[str, float]]]:

    LOG.info(f"Parsing PSM information file: {file} with FDR threshold: {fdr_threshold} and Score(PEP) threshold: {score_threshold}")

    peptide_map = {}

    with open(file, "r") as psminfo:
        reader = csv.DictReader(psminfo, delimiter="\t")
        if not {"SpectrumID", "PeptideSequence", "FDR", "Score(PEP)"}.issubset(reader.fieldnames) and not {"Spectrum", "Peptide", "FDR", "PEPscore"}.issubset(reader.fieldnames):
            LOG.warning(f"PSM information file {file} is missing required columns. Expected columns: (SpectrumID or Spectrum, PeptideSequence or Peptide, FDR, Score(PEP) or PEPscore). Found columns: {reader.fieldnames}. PSM information will not be loaded.")
            return
        for row in reader:

            peptide = row.get("PeptideSequence") or row.get("Peptide")
            psmFDR = row.get("FDR")
            score = row.get("Score(PEP)") or row.get("PEPscore")
            
            spectrumID = row.get("SpectrumID") or row.get("Spectrum")
    
            spectrum = ""
            sample = ""

            m = re.search('[^/_]+_([^/]+)_(\d+)$', spectrumID)
            if m:
                sample = m.group(1)
                spectrum = m.group(2)
            else:
                m = re.search('(.*)_\d+\s+controllerNumber=\d+\s+scan=(\d+)$', spectrumID)
                if m:
                    sample = m.group(1)
                    spectrum = m.group(2)

            if sample not in spectrum_table:
                m = re.search('(.*)_(\d+)$', spectrumID)
                if m:
                    sample = m.group(1)
                    spectrum = m.group(2)

            if sample not in spectrum_table:
                LOG.warning(f"Warning-61: Sample in PSMs not found in spectral table. Sample: {sample} Spectrum: {spectrum} Peptide: {peptide} \n {row}" )
                continue

            if peptide not in peptide_map:
                peptide_map[peptide] = PeptideInfo()
                
                for t in tmt:
                    peptide_map[peptide].rawTMT[t]=0

            peptide_map[peptide].Spectra[spectrumID] = score

            if spectrum in spectrum_table[sample].get(3, {}):
                LOG.warning(f"Warning-69: PSM spectrum is MS3 not MS2. Sample: {sample} Spectrum: {spectrum} Peptide: {peptide} \n {row}" )

            else:
                if spectrum in spectrum_map[sample]:
                    LOG.debug(f"Found match between PSM spectrum and spectrum map for sample {sample} spectrum {spectrum} peptide {peptide} with PSM FDR {psmFDR} and Score(PEP) {score}")
                    if float(psmFDR) <= fdr_threshold or float(score) <= score_threshold:
                        LOG.debug(f"PSM meets FDR or Score(PEP) threshold for quantification. Mapping TMT values to peptide {peptide} for sample {sample} spectrum {spectrum} with PSM FDR {psmFDR} and Score(PEP) {score}")
                        for t in tmt:
                            if t in spectrum_map[sample][spectrum]:
                                peptide_map[peptide].rawTMT[t] += spectrum_map[sample][spectrum][t]
                                LOG.debug(f"Mapping TMT channel {t} intensity {spectrum_map[sample][spectrum][t]} to peptide {peptide} for sample {sample} spectrum {spectrum} with PSM FDR {psmFDR} and Score(PEP) {score}")
                            else:
                                LOG.warning(f"Warning-63: No TMT channel match in MAP SID. Sample: {sample} Spectrum: {spectrum} Peptide: {peptide} TMTchan {t}" )
                else:
                    LOG.warning(f"Warning-62: No match between PSMs and spectrum TMT MAP SID. Sample: {sample} Spectrum: {spectrum} Peptide: {peptide}" )

                if spectrum in spectrum_table[sample][2]:
                    spectrum_table[sample][2][spectrum].qvalue = psmFDR
                    spectrum_table[sample][2][spectrum].score = score
                    spectrum_table[sample][2][spectrum].Peptide = peptide
                else:
                    LOG.warning(f"Warning-60: Not match between PSMs and spectrum table SID. Sample: {sample} Spectrum: {spectrum} Peptide: {peptide}" )

    return peptide_map

# =============================================================================
# Peptide Quantification - Quantify TMT intensities for peptides in peptide table using TMT values mapped to spectra and normalisation factors calculated from total channel intensities
# =============================================================================
def calculate_TMT_channels(peptide_map: Dict[str, PeptideInfo], tmt: Set[int], tmt_threshold: float = 0.05) -> Dict[str, Dict[int, float]]:
    
    LOG.info(f"Calculating total and median TMT channel intensities for TMT channels: {sorted(tmt)} with intensity threshold for filtering low intensity channels: {tmt_threshold}")
    
    TMTchannels = {}
    TMTchannels['Total'] = {}
    TMTchannels['Median'] = {}
    MaxINT = 0

    for tag in sorted(tmt):
        LOG.debug(f"Processing TMT channel: {tag}")
        for p in peptide_map:
            if tag in peptide_map[p].rawTMT:
                LOG.debug(f"Peptide: {p} - Raw TMT Intensity for channel {tag}: {peptide_map[p].rawTMT[tag]}")
        TMTchannels['Total'][tag] = sum(peptide_map[p].rawTMT[tag] for p in peptide_map if tag in peptide_map[p].rawTMT)
        LOG.debug(f"Total intensity for TMT channel {tag}: {TMTchannels['Total'][tag]}")
        
        if TMTchannels['Total'][tag] > 0:
            TMTchannels['Median'][tag] = statistics.median([peptide_map[p].rawTMT[tag] for p in peptide_map if tag in peptide_map[p].rawTMT and peptide_map[p].rawTMT[tag] > 0])
        else:
            TMTchannels['Median'][tag] = 0
        LOG.debug(f"Median intensity for TMT channel {tag}: {TMTchannels['Median'][tag]}")

        if TMTchannels['Total'][tag] > MaxINT: MaxINT = TMTchannels['Total'][tag]

    #Filter out TMT channels with very low total intensity as these channels are unlikely to provide reliable quantification and can skew normalisation factors. 
    #Threshold is set at 5% of max total intensity across channels but can be adjusted based on dataset and experimental design. 
    #Filtered channels will be logged and excluded from peptide quantification and output results.
    for tag in sorted(tmt):
        tmt_scaled = 0
        if TMTchannels['Total'][tag] > 0:
            tmt_scaled = TMTchannels['Total'][tag] / MaxINT

        if tmt_scaled < tmt_threshold:
            LOG.info("TMT Channel:" + str(tag) + " TOTAL INT == " + str(TMTchannels['Total'][tag]) + " MEDIAN INT ==" + str(TMTchannels['Median'][tag]) + " SCALE=" + str(tmt_scaled) + " -----THIS CHANNEL IS BELOW 5% OF MAX INTENSITY AND WILL BE FLAGGED FOR EXCLUSION FROM RESULTS-----")
            del TMTchannels['Total'][tag]
            del TMTchannels['Median'][tag]
        else:
            LOG.info("TMT Channel:" + str(tag) + " TOTAL INT == " + str(TMTchannels['Total'][tag]) + " MEDIAN INT ==" + str(TMTchannels['Median'][tag]) + " SCALE=" + str(tmt_scaled) )
        

    return TMTchannels


# =============================================================================
# Peptide Quantification - Quantify TMT intensities for peptides in peptide table using TMT values mapped to spectra and normalisation factors calculated from total channel intensities
# ============================================================================= 
def peptideQuant(peptide_map: Dict[str, PeptideInfo], prefix: str, file: str, tmt_channels: Dict[str, Dict[str, float]] ) -> None:
    
    outfile = prefix + "PeptideQuant.tsv"
    LOG.info(f"Quantifying peptides using PSM information from file: {file} and TMT channels: {sorted(tmt_channels['Total'].keys())} and writing results to file: {outfile}")

    OUT = open(outfile, 'w')
    with open(file, "r") as peptideinfo:
        reader = csv.DictReader(peptideinfo, delimiter="\t")
        if not {"Peptide"}.issubset(reader.fieldnames) and not {"PeptideSequence"}.issubset(reader.fieldnames):
            LOG.warning(f"Peptide information file {file} is missing required column 'Peptide' or 'PeptideSequence'. Found columns: {reader.fieldnames}. Peptide quantification will not be performed.")
            return
        
        headers = reader.fieldnames + ["QuantPSMs"] + [f"TMT_Abundance_Raw_{t}" for t in sorted(tmt_channels['Total'])] + [f"TMT_Abundance_Normalised_{t}" for t in sorted(tmt_channels['Total'])] + [f"TMT_Abundance_Scaled_{t}" for t in sorted(tmt_channels['Total'])]
        OUT.write("\t".join(headers) + "\n")

        for row in reader:
            peptide = row.get("Peptide") or row.get("PeptideSequence")
            OUT.write("\t".join(row[field] for field in reader.fieldnames))

            if peptide in peptide_map:

                OUT.write("\t" + str(len(peptide_map[peptide].Spectra)))

                for tag in sorted(tmt_channels['Total']):
                    if tag in peptide_map[peptide].rawTMT:
                        OUT.write("\t" + str(peptide_map[peptide].rawTMT[tag]))
                    else:
                        OUT.write("\t-")

                for tag in sorted(tmt_channels['Total']):
                    peptide_map[peptide].normTMT[tag] = 0
                    if tmt_channels['Median'][tag] > 0:
                        peptide_map[peptide].normTMT[tag] = (( peptide_map[peptide].rawTMT[tag] / tmt_channels['Median'][tag] ) * 100)

                    OUT.write("\t" + str(peptide_map[peptide].normTMT[tag]))
                    peptide_map[peptide].sumTMT += peptide_map[peptide].normTMT[tag]

                for tag in sorted(tmt_channels['Total']):
                    try: 
                        TMTscaled = (peptide_map[peptide].normTMT[tag] / (peptide_map[peptide].sumTMT / len(tmt_channels['Total']))) * 100
                    except ZeroDivisionError:
                        TMTscaled = '-'
                    OUT.write('\t' + str(TMTscaled))

            else:
                LOG.warning(f"Warning-PeptideQuant: Peptide not found in PSM map! Peptide: {peptide} \n {row}")

            OUT.write("\n")

    OUT.close()

# =============================================================================
# Protein Quantification - Quantify TMT intensities for proteins in protein table using peptide TMT values and normalisation factors calculated from total channel intensities
# =============================================================================

def proteinQuant(peptide_map: Dict[str, PeptideInfo], prefix: str, file: str, tmt_channels: Dict[str, Dict[str, float]], otype: str ) -> None:
    
    outfile = prefix + otype + "Quant.tsv"
    LOG.info(f"Quantifying {otype} using peptide information from file: {file} and TMT channels: {sorted(tmt_channels['Total'])} and writing results to file: {outfile}")

    OUT = open(outfile, 'w')
    with open(file, "r") as proteininfo:
        reader = csv.DictReader(proteininfo, delimiter="\t")
        if not {"ProteotypicPeptides"}.issubset(reader.fieldnames):
            LOG.warning(f"Proteotypic Peptide information file {file} is missing required column 'ProteotypicPeptides'. Found columns: {reader.fieldnames}. Protein quantification will not be performed.")
            return
        
        headers = reader.fieldnames + [f"TMT_Abundance_Raw_{t}" for t in sorted(tmt_channels['Total'])] + [f"TMT_Abundance_Normalised_{t}" for t in sorted(tmt_channels['Total'])] + [f"TMT_Abundance_Scaled_{t}" for t in sorted(tmt_channels['Total'])]
        OUT.write("\t".join(headers) + "\n")

        for row in reader:

            peptide_list = row.get("ProteotypicPeptides").split(',')
            tmt_quant = {t: 0 for t in sorted(tmt_channels['Total'])}
            for pep in peptide_list:
                if pep != '' and pep != '-':
                    if pep in peptide_map:
                        for tag in sorted(tmt_channels['Total']):
                            if tag in peptide_map[pep].rawTMT:
                                tmt_quant[tag] += peptide_map[pep].rawTMT[tag]
                            else:
                                LOG.warning(f"Warning-ProteinQuant: No TMT channel match in PSM map for peptide {pep} in protein: {row.get('ProteinID', 'Unknown Protein ID')} TMTchan {tag} \n {row}" )
                    else:
                        LOG.warning(f"Warning-ProteinQuant: Proteotypic peptide: {pep} not found in peptide map for protein: {row.get('ProteinID', 'Unknown Protein ID')} \n {row}")
                
            for tag in sorted(tmt_channels['Total']):
                row[f"TMT_Abundance_Raw_{tag}"] = tmt_quant[tag]
                if tmt_channels['Median'][tag] > 0:
                    normTMT = (( tmt_quant[tag] / tmt_channels['Median'][tag] ) * 100)
                else:
                    normTMT = 0
                row[f"TMT_Abundance_Normalised_{tag}"] = normTMT
            
            for tag in sorted(tmt_channels['Total']):
                try: 
                    TMTscaled = (row[f"TMT_Abundance_Normalised_{tag}"] / (sum(row[f"TMT_Abundance_Normalised_{t}"] for t in sorted(tmt_channels['Total'])) / len(tmt_channels['Total']))) * 100
                except ZeroDivisionError:
                    TMTscaled = '-'
                row[f"TMT_Abundance_Scaled_{tag}"] = TMTscaled

            OUT.write("\t".join(str(row[field]) for field in headers))
            OUT.write("\n")

    OUT.close()

# =============================================================================
# Main Script
# =============================================================================

def main() -> int:
    
    args = parse_args()
    configure_logging(args.log_level)

    LOG.info("Starting TMT Quantification Tool")
    LOG.info(f"Command line arguments: {args}")

    #initialize data structures for spectral table, spectral map and TMT channels
    spectralTable = {}
    spectralMap = {}
    tmtchans = set()

    #Get list of consensusXML files in directory
    conxmls = glob.glob(args.conxmldir + '/*.consensusXML')

    #Loop each mzML file in directory and parse spectra and TMT values from consensusXML files
    for mzml in glob.glob(args.mzmldir + '/*.mzML'):

        sample = mzml.split('/')[-1].replace('.mzML', '')
        cpat = sample + "_isoquant.consensusXML"
        conxml = ''
        for cxml in conxmls:
            if cpat in cxml:
                conxml = cxml

        if conxml == '':
            LOG.warning(f"WARNING: No consensusXML file found for mzML file {mzml} with expected pattern {cpat}. This mzml will be skipped.")
        else:
            LOG.info(f"Processing mzML file: {mzml} with sample name: {sample} and consensusXML: {conxml}")
            
            #1. Build Spectral Table from mzML - Scan Number, RT, MZ, Level, Spectrum Title
            parseMzML(spectralTable, sample, mzml)
            spectralIndex = indexSTable(spectralTable, sample)
            rtIndex = indexRT(spectralTable, sample)

            #2. Load TMT intensities from consensusXML into Spectral Table
            parseConsensusXML(spectralMap, spectralTable, rtIndex, spectralIndex, sample, conxml, tmtchans, args.mslevel)

    #check that any spectra were found in mzML files
    if len(spectralTable) == 0:
        LOG.error("No spectra found in mzML files! Check that your mzML files are correct and that the correct directory is specified! Exiting!")
        return 1
    else:
        LOG.info("Total Spectra Found in mzML files: " + str(len(spectralTable)))
        

    #check if any TMT channels were found
    if len(tmtchans) == 0:
        LOG.error("No TMT channels found in consensusXML files! Check that your consensusXML files contain TMT channels and that the correct directory is specified! Exiting!")
        return 1
    else:
        LOG.info("Total TMT channels found in consensusXML files: " + str(len(tmtchans)) + " Channels: " + str(sorted(tmtchans)))

    #3. Map Spectra in Table to PSM hits and peptide sequences 
    peptideMap = parsePSMs(spectralTable, spectralMap, args.psmtab, tmtchans, args.fdr_threshold, args.score_threshold)

    #4. Write Spectral Table with TMT data to file
    writeSTable(spectralTable, spectralMap, args.pre, tmtchans)

    #5. Calculate TMT channel normalisation factors from total channel intensities across all PSMs
    tmt_channels = calculate_TMT_channels(peptideMap, tmtchans, args.tmt_threshold)

    #6. Recalculate TMT peptide quantifications in PEPTIDE table (Rewrite Peptide Table)
    peptideQuant(peptideMap, args.pre, args.peptab, tmt_channels)

    if args.protab:
        #7. Protein quantification
        proteinQuant(peptideMap, args.pre, args.protab, tmt_channels, "Protein")

    if args.genetab:
        #8. Gene quantification
        proteinQuant(peptideMap, args.pre, args.genetab, tmt_channels, "Gene")

    LOG.info("TMT Quantification Tool finished successfully")

# Entry point of the script
if __name__ == "__main__":
    raise SystemExit(main())

