#!/usr/bin/env nextflow

/* 
----------------------------------------------------------------------------------------
 James Wright -  May 2025 - V1.0.1 
 PIPPA OpenMS Mutant Proteoform Analysis Pipeline
----------------------------------------------------------------------------------------
*/ 

nextflow.enable.dsl=2

// IMPORT MODULES AND PROCESSES
include {
    BUILD_SAINDEX;
    MERGE_MZMLS;
    MSGF;
    COMET;
    MSFRAGGER;
    FEATURE_EXTRACT;
    PERCOLATOR;
    PEP_CALC;
    MERGE_RESULTS;
    CONSENSUS_ID;
    RESULTS_FILTER;
    GENERATE_MZTAB;
    ISOBARIC_ANALYZER;
    EPIFANY;
    MAP_IDS;
    RESOLVE_CONFLICTS;
    PROTEIN_QUANT;
} from './modules/openms'

include {
    BUILDMUTANTDB;
    DECOYPYRAT;
    ANNOTATE_RESULTS;
    TMT_QUANT;
    MUTANT_PROTEOFORM_ANALYSIS;
} from './modules/icrtools'



//PIPPA MAIN WORKFLOW
workflow {

    //Workflow Header Print Summary of Settings
    log.info """\
         .                                                  .

         ===================================================
         I C R ---- PIPPA -> MUTANT PROTEOFORM <-  P I P E L I N E
         ===================================================

         mzML Path          : ${params.mzml_dir}
         Sample Grouping    : ${params.mzml_merge_pattern}
         mzML Merging       : ${params.mzml_merge ? "Yes" : "No"}

         Sample Table       : ${params.sample_info}
         Mutant Table       : ${params.mutation_table}

         Build Mutant Fasta : ${params.build_mutant_db ? "Yes" : "No"}

         Database Directory : ${params.database_dir}
         Reference FASTA    : ${params.reference_fasta}
         Decoy Strategy     : ${params.decoy_name ? params.decoy_name : "No decoy generation"}

         Output Directory   : ${params.out_dir}
         Log File Directory : ${params.log_dir}

         Search Engines     : ${params.use_msgf ? "MSGF+" : ""} ${params.use_comet ? "COMET" : ""} ${params.use_msfragger ? "MSFragger" : ""}
         PEP Calculation    : ${params.use_pep_calc ? "PEP-Calc" : "Percolator"}

         TMT Quant          : ${params.tmt ? "Yes" : "No"}
         MS3 Quant          : ${params.ms3_quant ? "Yes" : "No"}
         Labelling          : ${params.isotype}
         TMT Threshold      : ${params.tmt_channel_threshold}

         Enzyme             : ${params.enzyme}
         Cleavage           : ${params.cleavage_termini}
         Missed Cleavages   : ${params.missed_cleavages}
         Precursor Tol      : ${params.precursor_tolerance} ${params.precursor_tol_units}
         Fragment Tol       : ${params.fragment_tol} ${params.fragment_tol_units}
         Instrument Type    : ${params.instrument_type}

         Fixed Mods         : ${params.fixed_modifications}
         Variable Mods      : ${params.var_modifications}

         Peptide FDR        : ${params.peptide_fdr_threshold}
         PSM FDR            : ${params.psm_fdr_threshold}
         PSM PEP Threshold  : ${params.psm_pep_threshold}

         Mutant Analysis    : ${params.mutant_proteoforms ? "Yes" : "No"}
         Control Channel    : ${params.control_channel_name}
         
         ===================================================

         .                                                  .
        """
        .stripIndent()

    //Merge MZMLs if param is true otherwise pass mzMLs as they are
    ch_merged_mzml_files = Channel.empty()
    if (params.mzml_merge){    

        //Make channel with list of spectra file paths and sampleIDs based on pattern
        ch_mzml_files = Channel.fromPath("${params.mzml_dir}/*.mzML", checkIfExists: true)
            .map{ it -> [ ("${it.baseName}" =~ /($params.mzml_merge_pattern)(.*)/)[0][1], it ] }
            .view()

        ch_merged_mzml_files = MERGE_MZMLS(ch_mzml_files.groupTuple(sort:true)).merged_mzml
    } else {
        ch_merged_mzml_files = Channel.fromPath("${params.mzml_dir}/*.mzML", checkIfExists: true).view()
    }

    //Make channel for fasta databases
    //Regex explanation: match everything unless the _TDxxxx preceeds the .fasta
    ch_target_db = Channel.fromPath("${params.database_dir}/*.fasta", checkIfExists: true)
        .filter( file -> !file.name.matches(/.*_TD[^_]*\.fasta$/) ) 
       .view()

    //Create a channel with pairs of target and decoy database paths
    //check is decoy database exists, split into two branches  (decoy_exists and no_decoy)
    ch_target_db.map { it -> [ it, file("${it.Parent}/${it.baseName}_TD${params.decoy_name}.fasta") ] }
        .branch { 
            decoy_exists: it[1].exists() 
                return it[1]
            no_decoy: true 
                return it[0] 
        }
        .set{ ch_dbs }

    //If build_mutant_db is true, build mutant proteoform fasta database
    ch_mutant_db = Channel.empty()
    if (params.build_mutant_db) {
        ch_builddb_input = Channel.of(params.searchid)
            .map { it -> [ it, file(params.sample_info), file(params.mutation_table), file(params.cosmic_aa_reference), file(params.cosmic_na_reference) ] }
        ch_mutant_db = BUILDMUTANTDB(ch_builddb_input).fasta
        ch_mutant_db.view()
    }

    //Generate decoy sequences using decoyPYrat
    ch_new_decoys = DECOYPYRAT(ch_dbs.no_decoy.mix(ch_mutant_db)).fasta

    //Combine the new target+decoy databases with exisiting target_decoy databases into single channel
    ch_databases = ch_dbs.decoy_exists
        .mix(ch_new_decoys)
        .view()
    
    //If MSGF is selected, Index databases and run MSGF+ 
    ch_msgf_results = Channel.empty()
    if (params.use_msgf) {

        //Split databases into indexed and non-indexed
        ch_databases.map { it -> [ it, file("${it.Parent}/${it.baseName}.canno"), file("${it.Parent}/${it.baseName}.cnlcp"), file("${it.Parent}/${it.baseName}.csarr"), file("${it.Parent}/${it.baseName}.cseq") ] }
            .branch {  
                index_exists: it[1].exists() 
                    return it 
                no_index: true 
                    return it 
            }
            .set{ ch_msgf_dbs }

        //Index the non-indexed databases using MSGFs BuildSA tool
        ch_newly_indexed_msgf_dbs = BUILD_SAINDEX(ch_msgf_dbs.no_index.map{ item -> [ item[0] ] }).dbs

        //Combine the previously indexed databases with the newly indexed databases into a single channel
        ch_indexed_msgf_databases = ch_msgf_dbs.index_exists
            .mix(ch_newly_indexed_msgf_dbs)

        //Run MSGF+ on the indexed databases
        ch_msgf_results = MSGF(ch_merged_mzml_files.combine(ch_indexed_msgf_databases)).results   
    }

    //If COMET is selected run COMET
    ch_comet_results = Channel.empty()
    if (params.use_comet) {
        ch_comet_results = COMET(ch_merged_mzml_files.combine(ch_databases)).results
    }

    //If MSFragger is selected run MSFragger
    ch_msfragger_results = Channel.empty()
    if (params.use_msfragger) {
        ch_msfragger_results = MSFRAGGER(ch_merged_mzml_files.combine(ch_databases)).results
    }

    //Use PSM Feature Extractor to extract search engine specific features for percolator
    ch_psm_features = Channel.empty()
    ch_psm_features = FEATURE_EXTRACT(ch_msgf_results.mix(ch_comet_results).mix(ch_msfragger_results)).features

    //Run Percolator or PEP calc on the PSMs
    ch_pep_results = Channel.empty()
    if (params.use_pep_calc) {
        ch_pep_results = PEP_CALC(ch_psm_features).idxml
    } else {
        ch_pep_results = PERCOLATOR(ch_psm_features).idxml
    }

    //Merge the results from all search engines
    ch_merged_results = MERGE_RESULTS(ch_pep_results.groupTuple()).idxml

    //Run ConsensusID to combine the results from all search engines
    ch_consensus_results = CONSENSUS_ID(ch_merged_results).idxml

    //Run Results Filter to filter the results
    ch_filtered_results = RESULTS_FILTER(ch_consensus_results).idxml

    //Generate mztab file from the filtered results
    ch_mztab = GENERATE_MZTAB(ch_filtered_results).mztab

    //Run Epifany to generate a protein inference
    ch_epifany_results = EPIFANY(ch_filtered_results).idxml

    /////// QUANTIFICATION PROCESSES //////////
    ch_pc_mztab = Channel.empty()
    if (params.tmt) {
        //Run IsobaricAnalyzer to analyze the TMT peaks
        ch_isobaric_results = ISOBARIC_ANALYZER(ch_merged_mzml_files).consensusXML

        ch_iso_epifany_cross = ch_isobaric_results.cross(ch_epifany_results).map{ it -> [ it[0][0], it[1][1], it[0][1], it[1][2] ] } .view()

        //Map the IDs to the quant features
        //Cross the isobaric quant results and epifany results and map to input format
        ch_quant_mapped_ids = MAP_IDS(ch_iso_epifany_cross).consensusXML

        //Resolve conflicts in the quant features
        ch_resolved_quant = RESOLVE_CONFLICTS(ch_quant_mapped_ids).consensusXML

        //Run ProteinQuantifier to quantify the proteins
        ch_pc_mztab = PROTEIN_QUANT(ch_resolved_quant).mztab
    }

    //Run Results Summary Script
    ANNOTATE_RESULTS(ch_mztab.mix(ch_pc_mztab).collect())

    //run ICR TMT Quantification 
    if (params.tmt) {
        ch_summary_results = ANNOTATE_RESULTS.out.summary
        TMT_QUANT(ch_summary_results)
    }

    //run Mutant Proteoform Analysis
    if (params.mutant_proteoforms) {
        //tuple val(id), path(spectrum_table), path(peptide_table), path(mutations)
        if (params.build_mutant_db) {
            ch_mutation_info = BUILDMUTANTDB.out.mutations
        } else {
            ch_mutation_info = Channel.fromPath("${params.mutation_info}", checkIfExists: true)
        }
        ch_quant_tables = TMT_QUANT.out.quant_summary.combine(ch_mutation_info)
        MUTANT_PROTEOFORM_ANALYSIS(ch_quant_tables)
    }
    
}

//END OF WORKFLOW