// modules/openms.nf

//BUILD SAINDEX(); indexes the database using MSGF+ BuildSA tool
process BUILD_SAINDEX {

    clusterOptions '--job-name=nf-buildSA --output=nf_buildSA.txt'
    label 'Medium_Task'

    conda = "$baseDir/env/openms"

    input: 
    path database 

    output: 
    tuple path(database), path("${database.baseName}.canno"), path("${database.baseName}.cnlcp"), path("${database.baseName}.csarr"), path("${database.baseName}.cseq"), emit: dbs
    path "*.log", emit: hidden // only for logging into the log directory

    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.database_dir}", mode: 'copy', pattern: '*.canno'
    publishDir "${params.database_dir}", mode: 'copy', pattern: '*.cnlcp'
    publishDir "${params.database_dir}", mode: 'copy', pattern: '*.csarr'
    publishDir "${params.database_dir}", mode: 'copy', pattern: '*.cseq'

    script: 

    mem = task.cpus * 6000

    """
    ts=\$(date '+%Y%m%d_%H%M%S')

    java -Xmx${mem}M -cp ${params.msgf_path} \\
    edu.ucsd.msjava.msdbsearch.BuildSA \\
    -d $database \\
    -tda 0 \\
    > ${task.process}_\${ts}_${database.baseName}_BuildSA.log

    """
}

//MERGE_MZMLS(); takes a tuple of sample name and mzML files and merges mzMLs into a single file
process MERGE_MZMLS {

    clusterOptions  '--job-name=nf_mergeMzML --output=nf_mergeMzML.txt '
    label 'Medium_Task'

    conda = "$baseDir/env/openms"

    input:
    tuple val(sample), path(mzmls) 

    output:
    path "${sample}.mzML", emit: merged_mzml
    path "*.log", emit: hidden // only for logging into the log directory

    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.out_dir}/merged_mzml", mode: 'copy', pattern: '*.mzML'

    script:
     """
     ts=\$(date '+%Y%m%d_%H%M%S')

     FileMerger  -in ${mzmls.join(' ')} \\
                 -out ${sample}.mzML \\
                 -threads ${task.cpus} \\
                 -rt_concat:gap 0.1 \\
                 > ${task.process}_\${ts}_${sample}_mzML_FileMerge.log
     """
}

//MSFG(); Run MSGFplus Search of mzML files against an indexed database
process MSGF {

    clusterOptions  '--job-name=nf_msgfplus --output=nf_msgfplus.txt '
    label 'MSGF_Task'

    conda = "$baseDir/env/openms"

    input:
    tuple path(mzml), path(database), path(canno), path(cnlcp), path(csarr), path(cseq)

    output:
    tuple val("V${params.searchid}_${mzml.baseName}_DBx${database.baseName}"), val("V${params.searchid}_${mzml.baseName}"), path("V${params.searchid}_${mzml.baseName}_DBx${database.baseName}_msgfplus.idxml"), emit: results
    path "*.log", emit: hidden // only for logging into the log directory

    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.out_dir}/idxmls", mode: 'copy', pattern: '*.idxml'

    script:

    mem = task.cpus * 6000

     """
     ts=\$(date '+%Y%m%d_%H%M%S')

     MSGFPlusAdapter -in $mzml \\
                     -out V${params.searchid}_${mzml.baseName}_DBx${database.baseName}_msgfplus.idxml \\
		             -executable ${params.msgf_path} \\
                     -database $database \\
                     -instrument ${params.instrument_type} \\
                     -enzyme '${params.msgf_enzyme}' \\
                     -tryptic ${params.cleavage_termini} \\
                     -precursor_mass_tolerance ${params.precursor_tolerance} \\
                     -precursor_error_units ${params.precursor_tol_units} \\
                     -fixed_modifications ${params.fixed_modifications.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                     -variable_modifications ${params.var_modifications.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                     -max_missed_cleavages ${params.msgf_missed_cleavages} \\
                     -isotope_error_range ${params.msgf_isotope_error} \\
                     -fragment_method from_spectrum \\
                     -protocol automatic \\
                     -min_precursor_charge ${params.min_pre_charge} \\
                     -max_precursor_charge ${params.max_pre_charge} \\
                     -min_peptide_length ${params.min_pep_length} \\
                     -max_peptide_length ${params.max_pep_length} \\
                     -matches_per_spec ${params.hit_ranks} \\
                     -add_features true \\
                     -max_mods ${params.max_modifications} \\
                     -PeptideIndexing:unmatched_action ${params.unmatched_peptide_action} \\
                     -tasks ${params.msgf_tasks} \\
                     -threads ${task.cpus} \\
                     -java_memory $mem \\
                     > ${task.process}_\${ts}_${mzml.baseName}_DBx${database.baseName}_msgfplus.log
     """
}

    
//SAGE(); Run Sage Search of mzML files against database
process SAGE {

    clusterOptions  '--job-name=nf_sage --output=nf_sage.txt '
    label 'SAGE_Task'

    conda = "$baseDir/env/openms"

    input:
    tuple path(mzml), path(database) 

    output:
    tuple val("V${params.searchid}_${mzml.baseName}_DBx${database.baseName}"), val("V${params.searchid}_${mzml.baseName}"), path("V${params.searchid}_${mzml.baseName}_DBx${database.baseName}_sage.idxml"), emit: results
    path "*.log", emit: hidden // only for logging into the log directory
     
    //Copy output to log file folder
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.out_dir}/idxmls", mode: 'copy', pattern: '*.idxml'

    script:

    threads = task.cpus 
    neg_pre_tol = -1 * params.precursor_tolerance
    neg_frag_tol = -1 * params.sage_frag_tol

    """
    ts=\$(date '+%Y%m%d_%H%M%S')

    SageAdapter         -sage_executable ${params.sage_path} \\
                        -in $mzml \\
                        -out V${params.searchid}_${mzml.baseName}_DBx${database.baseName}_sage.idxml \\
                        -database $database \\
                        -decoy_prefix DECOY_ \\
                        -precursor_tol_left $neg_pre_tol \\
                        -precursor_tol_right ${params.precursor_tolerance} \\
                        -precursor_tol_unit ${params.precursor_tol_units} \\
                        -fragment_tol_left $neg_frag_tol \\
                        -fragment_tol_right ${params.sage_frag_tol} \\
                        -fragment_tol_unit ${params.sage_frag_tol_units} \\
                        -enzyme '${params.enzyme}' \\
                        -fixed_modifications ${params.fixed_modifications.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                        -variable_modifications ${params.var_modifications.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                        -reindex ${params.sage_reindex} \\
                        -threads $threads \\
                        -min_len ${params.min_pep_length} \\
                        -max_len ${params.max_pep_length} \\
                        -missed_cleavages ${params.missed_cleavages} \\
                        -fragment_min_mz ${params.sage_frag_min_mz} \\
                        -fragment_max_mz ${params.sage_frag_max_mz} \\
                        -peptide_min_mass ${params.sage_pep_min_mass} \\
                        -peptide_max_mass ${params.sage_pep_max_mass} \\
                        -max_variable_mods ${params.sage_max_var_mods} \\
                        -PeptideIndexing:unmatched_action ${params.unmatched_peptide_action} \\
                        -charges ${params.sage_charges} \\
                        > ${task.process}_\${ts}_${mzml.baseName}_DBx${database.baseName}_sage.log
     
    """
    //TODO: remove idxml file if it is empty - maybe add process to count peptide hits and report and only return idxmls with >n hits
}


//MSFRAGGER(); Run MSFragger Search of mzML files against a database
process MSFRAGGER {

    clusterOptions  '--job-name=nf_msfragger --output=nf_msfragger.txt '
    label 'MSFRAG_Task'

    conda = "$baseDir/env/openms"

    input:
    tuple path(mzml), path(database) 

    output:
    tuple val("V${params.searchid}_${mzml.baseName}_DBx${database.baseName}"), val("V${params.searchid}_${mzml.baseName}"), path("V${params.searchid}_${mzml.baseName}_DBx${database.baseName}_msfragger.idxml"), emit: results
    path "*.log", emit: hidden // only for logging into the log directory

    //Copy output to log file folder
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.out_dir}/idxmls", mode: 'copy', pattern: '*.idxml'

    script:

    mem = task.cpus * 7500

    flags = ""
    if (params.msfragger_clip_nterm_m) flags = flags + "-varmod:clip_nterm_m "
    if (params.msfragger_varmod_multi) flags = flags + "-varmod:not_allow_multiple_variable_mods_on_residue "
    if (params.msfragger_override_charge) flags = flags + "-spectrum:override_charge "

     """
    ts=\$(date '+%Y%m%d_%H%M%S')

    MSFraggerAdapter    -license yes \\
                        -java_executable ${params.java_path} \\
                        -java_heapmemory $mem \\
                        -executable ${params.msfragger_path} \\
                        -in \$PWD/$mzml \\
                        -out \$PWD/V${params.searchid}_${mzml.baseName}_DBx${database.baseName}_msfragger.idxml \\
                        -opt_out \$PWD/V${params.searchid}_${mzml.baseName}_DBx${database.baseName}_msfragger.pepXML \\
                        -database \$PWD/$database \\
                        -tolerance:precursor_mass_tolerance_lower ${params.precursor_tolerance} \\
                        -tolerance:precursor_mass_tolerance_upper ${params.precursor_tolerance} \\
                        -tolerance:precursor_mass_unit ${params.precursor_tol_units} \\
                        -tolerance:precursor_true_tolerance ${params.msfragger_pre_true_tol} \\
                        -tolerance:precursor_true_unit ${params.msfragger_pre_true_tol_units} \\
                        -tolerance:fragment_mass_tolerance ${params.msfragger_frag_tol} \\
                        -tolerance:fragment_mass_unit ${params.msfragger_frag_tol_units} \\
                        -tolerance:isotope_error ${params.msfragger_isotope_err} \\
                        -digest:search_enzyme_name ${params.msfragger_enzyme_name.tokenize(',').collect { "'${it}'" }.join(" ")} \\
                        -digest:search_enzyme_cutafter ${params.msfragger_enzyme_cutafter} \\
                        -digest:search_enzyme_nocutbefore ${params.msfragger_enzyme_nocutbefore} \\
                        -digest:num_enzyme_termini ${params.msfragger_enzyme_numtermini.tokenize(',').collect { "'${it}'" }.join(" ")} \\
                        -digest:allowed_missed_cleavage ${params.missed_cleavages} \\
                        -digest:min_length ${params.min_pep_length} \\
                        -digest:max_length ${params.max_pep_length} \\
                        -digest:mass_range_min ${params.msfragger_digest_min_mass} \\
                        -digest:mass_range_max ${params.msfragger_digest_max_mass} \\
                        -varmod:max_variable_mods_per_peptide ${params.msfragger_varmod_maxmod} \\
                        -varmod:max_variable_mods_combinations ${params.msfragger_varmod_maxcomb} \\
                        -spectrum:minimum_peaks ${params.msfragger_minpeaks} \\
                        -spectrum:use_topn_peaks ${params.msfragger_topnpeaks} \\
                        -spectrum:minimum_ratio ${params.msfragger_minratio} \\
                        -spectrum:clear_mz_range_min ${params.msfragger_clear_mz_min} \\
                        -spectrum:clear_mz_range_max ${params.msfragger_clear_mz_max} \\
                        -spectrum:max_fragment_charge ${params.msfragger_max_frag_charge} \\
                        -spectrum:precursor_charge_min ${params.msfragger_min_charge} \\
                        -spectrum:precursor_charge_max ${params.msfragger_max_charge} \\
                        -search:track_zero_topn ${params.msfragger_track_zero_topn} \\
                        -search:zero_bin_accept_expect ${params.msfragger_zero_bin_accept_expect} \\
                        -search:zero_bin_mult_expect ${params.msfragger_zero_bin_mult_expect} \\
                        -search:add_topn_complementary ${params.msfragger_add_topn_comp} \\
                        -search:min_fragments_modeling ${params.msfragger_min_frag_modeling} \\
                        -search:min_matched_fragments ${params.msfragger_min_frag_matched} \\
                        -search:output_report_topn ${params.msfragger_report_topn} \\
                        -search:output_max_expect ${params.msfragger_max_evalue} \\
                        -search:localize_delta_mass ${params.msfragger_local_deltamass} \\
                        -statmod:add_C_cysteine 0.0 \\
                        -varmod:unimod ${params.var_modifications.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                        -statmod:unimod ${params.fixed_modifications.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                        -reindex ${params.msfragger_reindex} \\
                        -PeptideIndexing:unmatched_action ${params.unmatched_peptide_action} \\
                        -threads ${task.cpus} \\
                        $flags \\
                        > \$PWD/${task.process}_\${ts}_${mzml.baseName}_DBx${database.baseName}_msfragger.log

     """
}

//COMET();  Run Comet Search of mzML files against a database
process COMET {

    clusterOptions  '--job-name=nf_comet --output=nf_comet.txt '
    label 'Comet_Task'

    conda = "$baseDir/env/openms"

    input:
    tuple path(mzml), path(database) 
    
    output:
    tuple val("V${params.searchid}_${mzml.baseName}_DBx${database.baseName}"), val("V${params.searchid}_${mzml.baseName}"), path("V${params.searchid}_${mzml.baseName}_DBx${database.baseName}_comet.idxml"), emit: results
    path "*.log", emit: hidden // only for logging into the log directory

    //Copy output to log file folder
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.out_dir}/idxmls", mode: 'copy', pattern: '*.idxml'

    script:
     """

        ts=\$(date '+%Y%m%d_%H%M%S')

         CometAdapter  -in $mzml \\
                       -out V${params.searchid}_${mzml.baseName}_DBx${database.baseName}_comet.idxml \\
                       -database $database \\
                       -instrument ${params.instrument_type} \\
                       -missed_cleavages ${params.missed_cleavages} \\
                       -num_enzyme_termini ${params.cleavage_termini} \\
                       -enzyme '${params.com_enzyme}' \\
		               -fixed_modifications ${params.fixed_modifications.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                       -variable_modifications ${params.var_modifications.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                       -precursor_mass_tolerance ${params.precursor_tolerance} \\
                       -precursor_error_units ${params.precursor_tol_units} \\
                       -fragment_mass_tolerance ${params.com_frag_tolerance} \\
                       -clip_nterm_methionine ${params.com_clip_met} \\
                       -fragment_bin_offset ${params.com_frag_bin_offset} \\
                       -isotope_error ${params.com_isotope_error} \\
                       -use_A_ions ${params.com_a_ion} \\
                       -use_B_ions ${params.com_b_ion} \\
                       -use_C_ions ${params.com_c_ion} \\
                       -use_X_ions ${params.com_x_ion} \\
                       -use_Y_ions ${params.com_y_ion} \\
                       -use_Z_ions ${params.com_z_ion} \\
                       -use_NL_ions ${params.com_nl_ion} \\
                       -second_enzyme ${params.com_second_enzyme} \\
                       -min_peptide_length ${params.min_pep_length} \\
                       -max_peptide_length ${params.max_pep_length} \\
                       -num_hits ${params.hit_ranks} \\
                       -override_charge '${params.com_override_charge}' \\
                       -ms_level ${params.com_mslevel} \\
                       -activation_method ALL \\
                       -digest_mass_range ${params.com_mass_range} \\
                       -max_fragment_charge ${params.com_max_frag_charge} \\
                       -max_precursor_charge ${params.max_pre_charge} \\
                       -spectrum_batch_size ${params.com_batch_size} \\
                       -minimum_peaks ${params.com_min_peaks} \\
                       -minimum_intensity ${params.com_min_intensity} \\
                       -remove_precursor_peak ${params.com_remove_pre} \\
                       -remove_precursor_tolerance ${params.com_remove_pre_tol} \\
                       -clear_mz_range ${params.com_clear_mz} \\
                       -binary_modifications ${params.com_binary_mod} \\
                       -max_variable_mods_in_peptide ${params.max_modifications} \\
                       -require_variable_mod ${params.com_require_mod} \\
                       -PeptideIndexing:unmatched_action ${params.unmatched_peptide_action} \\
                       -threads ${task.cpus} \\
                   > ${task.process}_\${ts}_${mzml.baseName}_DBx${database.baseName}_comet.log
     """
}



//FEATURE_EXTRACT(); extracts search engine specifc percolator features from the idXML files
process FEATURE_EXTRACT {

    clusterOptions  '--job-name=nf_psmextract --output=nf_psmextract.txt '

    conda = "$baseDir/env/openms"

    input:
    tuple val(id), val(mzml), path(input_idxml) 

    output:
    tuple val(id), val(mzml), path("${input_idxml.baseName}_psmfeatures.idXML"), emit: features
    path "*.log", emit: hidden // only for logging into the log directory

    //Copy output to log file folder
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    script:
     """
     ts=\$(date '+%Y%m%d_%H%M%S')

     PSMFeatureExtractor -in $input_idxml \\
                         -out ${input_idxml.baseName}_psmfeatures.idXML \\
                         -threads ${task.cpus} \\
                         > ${task.process}_\${ts}_${input_idxml.baseName}_psm_feature_extractor.log
     """
}

//PERCOLATOR(); Run Percolator on the target decoy search engine PSM results
process PERCOLATOR {

    clusterOptions  '--job-name=nf_perc --output=nf_perc.txt '

    conda = "$baseDir/env/openms"

    errorStrategy 'ignore'

    input:
    tuple val(id), val(mzml), path(input_idxml) 

    output:
    tuple val(id), val(mzml), path("${input_idxml.baseName}_percolator.idXML"), emit: idxml
    path "*.log", emit: hidden // only for logging into the log directory
    path "*.tab", emit: tab // only for copying into the mztab directory

    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.out_dir}/tabs", mode: 'copy', pattern: '*.tab'

    script:

      """

      ts=\$(date '+%Y%m%d_%H%M%S')

      PercolatorAdapter   -in ${input_idxml} \\
                          -out ${input_idxml.baseName}_percolator.idXML \\
                          -out_pout_target ${input_idxml.baseName}_percolator.tab \\
                          -decoy_pattern DECOY_ \\
                          -score_type pep \\
                          -threads ${task.cpus} \\
                          > ${task.process}_\${ts}_${input_idxml.baseName}_percolator.log
      """
}
   
//PEP_CALC(); Run Posterior Error Probability calculation on the target decoy search engine PSM results
process PEP_CALC {

    clusterOptions  '--job-name=nf_IDPEP --output=nf_IDPEP.txt '

    conda = "$baseDir/env/openms"

    input:
    tuple val(id), val(mzml), path(input_idxml) 

    output:
    tuple val(id), val(mzml), path("${input_idxml.baseName}_nopercolator.idXML") , emit: idxml
    path "*.log", emit: hidden // only for logging into the log directory

    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    script:

      """

      ts=\$(date '+%Y%m%d_%H%M%S')

      IDPosteriorErrorProbability   -in ${input_idxml} \\
                          -out ${input_idxml.baseName}_nopercolator.idXML \\
                          -threads ${task.cpus} \\
                          > ${task.process}_\${ts}_${input_idxml.baseName}_IDPosteriorProb.log
      """
}
   

//MERGE_RESULTS(); takes the idXML files from all search engines and merges them into a single idXML file
process MERGE_RESULTS {

    clusterOptions  '--job-name=nf_mergeIDs --output=nf_mergeIDs.txt '

    conda = "$baseDir/env/openms"

    input:
    tuple val(id), val(mzml), path(input_idxmls) 

    output:
    tuple val(id), val("${mzml[0]}"), path ("${id}_mergedresults.idXML") , emit: idxml
    path "*.log", emit: hidden // only for logging into the log directory
    
     //Copy output to log file folder
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    script:
    """

    ts=\$(date '+%Y%m%d_%H%M%S')

    IDMerger        -in $input_idxmls \\
                    -out ${id}_mergedresults.idXML \\
                    -threads ${task.cpus} \\
                    > ${task.process}_\${ts}_${id}_IDMerge.log
    """

}

//CONSENSUS_ID(); combines and merges the search engine specific PSMs into single set of PSMs
process CONSENSUS_ID {
    
    clusterOptions  '--job-name=nf_consensus --output=nf_consensus.txt '

    conda = "$baseDir/env/openms"

    input:
    tuple val(id), val(mzml), path(idxml)

    output:
    tuple val(id), val(mzml), path("${idxml.baseName}_consensus.idXML") , emit: idxml
    path "*.log", emit: hidden // only for logging into the log directory

    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.out_dir}/idxmls", mode: 'copy', pattern: '*.idXML'

    script:
     """
     ts=\$(date '+%Y%m%d_%H%M%S')

     ConsensusID    -in $idxml \\
                    -out ${idxml.baseName}_consensus.idXML \\
                    -per_spectrum \\
                    -algorithm best \\
                    -threads ${task.cpus} \\
                    > ${task.process}_\${ts}_${idxml.baseName}_ConsensusID.log
     """
}

// RESULTS_FILTER(); filters the consensus idXML PSMs based on the user defined parameters
//This is needed for FIdo and Epifany to work with large sets of spectra as if scores to extreme they crash. 
//This PEP threshold of 0.5 is not significant as is currently hard coded as setting it too low will also cause problems for inference.
process RESULTS_FILTER {

    clusterOptions  '--job-name=nf_IDfilter --output=nf_IDfilter.txt '

    conda = "$baseDir/env/openms"

    input:
    tuple val(id), val(mzml), path(idxml)

    output:
    tuple val(id), val(mzml), path("${idxml.baseName}_filtered.idXML") , emit: idxml
    path "*.log", emit: hidden // only for logging into the log directory

    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    script:
     """
     ts=\$(date '+%Y%m%d_%H%M%S')

     IDFilter   -in $idxml \\
                -out ${idxml.baseName}_filtered.idXML \\
                -score:psm 0.5 \\
                > ${task.process}_\${ts}_${idxml.baseName}_IDFilter.log
     """
}



//GENERATE_MZTAB(); takes the filtered idXML files and generates a mztab file
process GENERATE_MZTAB {

    clusterOptions  '--job-name=nf_writeMZTAB --output=nf_writeMZTAB.txt '

    conda = "$baseDir/env/openms"

    input:
    tuple val(id), val(mzml), path(idxml) 

    output:
    path ("${idxml.baseName}.mzTab") , emit : mztab
    path "*.log", emit: hidden // only for logging into the log directory
 
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.out_dir}/tabs", mode: 'copy', pattern: '*.mzTab'

    script:
    """
    ts=\$(date '+%Y%m%d_%H%M%S')

    MzTabExporter   -in $idxml \\
                    -out ${idxml.baseName}.mzTab \\
                    > ${task.process}_\${ts}_${idxml.baseName}_MzTabExport.log
    """
}

//ISOBARIC_ANALYZER(); Extracts TMT peaks from raw spectra and generates a consensusXML file
process ISOBARIC_ANALYZER {

    clusterOptions  '--job-name=nf_isoQuant --output=nf_isoQuant.txt '
    label 'Medium_mem'

    conda = "$baseDir/env/openms"

    input:
    path(mzml) 

    output:
    tuple val("V${params.searchid}_${mzml.baseName}"), path ("${mzml.baseName}_isoquant.consensusXML"), emit: consensusXML
    path "*.log", emit: hidden // only for logging into the log directory
 
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.out_dir}/consensusXMLs", mode: 'copy', pattern: '*.consensusXML'

    script:
     """
     ts=\$(date '+%Y%m%d_%H%M%S')

     IsobaricAnalyzer   -in $mzml \\
                        -out ${mzml.baseName}_isoquant.consensusXML \\
                        -type '${params.isotype}' \\
                        -threads ${task.cpus} \\
                        -extraction:select_activation '${params.qactivation}' \\
                 > ${task.process}_\${ts}_${mzml.baseName}_IsoANz.log
     """
}

//EPIFANY(); Run Epifany to calculate the protein inference and generate a idXML file
process EPIFANY {

    clusterOptions  '--job-name=nf_epifany --output=nf_epifany.txt '
    label 'Medium_Task'
    errorStrategy 'ignore' 

    conda = "$baseDir/env/openms"

    input:
    tuple val(rid),  val(mzml), path(idxml) 

    output:
    tuple val(mzml), val(rid), path("${idxml.baseName}_epifany.idXML") , emit: idxml
    path "*.log", emit: hidden // only for logging into the log directory

    //Copy output to log file folder
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    script:
     """
     ts=\$(date '+%Y%m%d_%H%M%S')

     Epifany   -in $idxml \\
                    -out ${idxml.baseName}_epifany.idXML \\
                    -threads ${task.cpus} \\
                 > ${task.process}_\${ts}_${idxml.baseName}_epifany.log
     """
}

//MAP_IDS(); takes the consensusXML file and the idXML file and maps the IDs to the quantification results
process MAP_IDS {

    clusterOptions  '--job-name=nf_IDmap --output=nf_IDmap.txt '
    label 'Medium_Task'

    conda = "$baseDir/env/openms"

    input:
    tuple val(mzml), val(pid), path(conxml), path(protxml) 

    output:
    tuple val(pid), path(protxml), path("${protxml.baseName}_isoquant_mapped.consensusXML"), emit: consensusXML
    path "*.log", emit: hidden // only for logging into the log directory

    //Copy output to log file folder
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    script:
     """
     ts=\$(date '+%Y%m%d_%H%M%S')

     IDMapper  	-id $protxml \\
		        -in $conxml \\
                -out ${protxml.baseName}_isoquant_mapped.consensusXML \\
                -threads ${task.cpus} \\
                 > ${task.process}_\${ts}_${protxml.baseName}_IQmapping.log
     """
}

//RESOLVE_CONFLICTS(); takes the mapped consensusXML file and the idXML file and resolves conflicts
//Required for protein quantification tool 
process RESOLVE_CONFLICTS {
    
    clusterOptions  '--job-name=nf_IDresolver --output=nf_IDResolve.txt '
    label 'Medium_Task'

    conda = "$baseDir/env/openms"

    input:
    tuple val(rid), path(protxml), path(conxml) 

    output:
    tuple val(rid), path(protxml), path("${conxml.baseName}_Resolved.consensusXML"), emit: consensusXML 
    path "*.log", emit: hidden // only for logging into the log directory

    //Copy output to log file folder
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    script:
     """
     ts=\$(date '+%Y%m%d_%H%M%S')

     IDConflictResolver  -in $conxml \\
			-out ${conxml.baseName}_Resolved.consensusXML \\
		> ${task.process}_\${ts}_${conxml.baseName}_IDResolver.log
     """
}

//PROTEIN_QUANT(); Run ProteinQuantifier on the resolved consensusXML file and the idXML file to get mzTab with Quant values
process PROTEIN_QUANT {

    clusterOptions  '--job-name=nf_proQuant --output=nf_proQuant.txt '
    label 'Medium_Task'

    conda = "$baseDir/env/openms"

    errorStrategy 'ignore'

    input:
    tuple val(rid), path(protxml), path(conxml) 

    output:
    path("${conxml.baseName}_ProteinQuant.mzTab"), emit: mztab
    path "*.log", emit: hidden // only for logging into the log directory

    //Copy output to log file folder
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.out_dir}/tabs", mode: 'copy', pattern: '*.mzTab'

    script:
     """
     ts=\$(date '+%Y%m%d_%H%M%S')

     ProteinQuantifier  -in $conxml \\
                        -protein_groups $protxml \\
                        -out ${conxml.baseName}_ProteinQuant.csv \\
                        -peptide_out ${conxml.baseName}_PeptideQuant.csv \\
                        -mztab ${conxml.baseName}_ProteinQuant.mzTab \\
                        -top:N  ${params.pq_top} \\
                        -top:aggregate ${params.pq_average} \\
                        -consensus:normalize \\
                        -ratios \\
                        -threads ${task.cpus} \\
                 > ${task.process}_\${ts}_${conxml.baseName}_ProteinQuant.log
     """
}

