// modules/icrtools.nf

// BUILDMUTANTDB(); takes a reference fasta file and a mutant table and generates a mutant proteoform database
process BUILDMUTANTDB {

    clusterOptions '--job-name=nf_builmutdb --output=nf_buildmutdb.txt '

    input:
    tuple val(id), path(sample_tsv), path(mutation_tsv), path(aa_ref), path(na_ref)

    output:
    path "*_mutant_proteoforms.fasta", emit: fasta
    path "*_mapped_mutations.tsv", emit: mutations
    path "*.log", emit: hidden // only for logging into the log directory

    publishDir "${params.database_dir}", mode: 'copy', pattern: '*_mutant_proteoforms.fasta'
    publishDir "${params.out_dir}/final_results", mode: 'copy', pattern: '*_mapped_mutations.tsv'
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    script:

    flags = ""
    if (params.cosmic_lockversion) flags = flags + "-v "
    if (params.cosmic_old) flags = flags + "-c "

    """
    ts=\$(date '+%Y%m%d_%H%M%S')

    python3 ${baseDir}/bin/Build_Cosmic_FASTA_Databases.py  -i ${mutation_tsv} \\
                                                            -p ${aa_ref} \\
                                                            -n ${na_ref} \\
                                                            -o ${id} \\
                                                            -s ${sample_tsv} \\
                                                            -d ${params.cosmic_dist} \\
                                                            $flags \\
                                                            > ${task.process}_\${ts}_${id}_BuildMutantDB.log


    """


}
 
// DECOYPYRAT(); takes a protein fasta file and generates a decoy database
process DECOYPYRAT {

    clusterOptions  '--job-name=nf_decoyPYrat --output=nf_decoyPYrat.txt '
    label 'Medium_mem'

    //Input channel is protein fasta file from params.database
    input:
    path target

    //Output channel contains path to target decoy database and the log file
    output:
    path "${target.baseName}_TD${params.decoy_name}.fasta" , emit: fasta
    path "*.log", emit: hidden // only for logging into the log directory

    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.database_dir}", mode: 'copy', pattern: '*.fasta'

    //Run decoyPYrate
    script:

    flags = "-ot "
    if (params.decoypyrate_pfasta) flags = flags + "-pf "
    if (params.decoypyrate_memsave) flags = flags + "-m "

     """
     ts=\$(date '+%Y%m%d_%H%M%S')

    python3 ${baseDir}/bin/decoyPYrateV2.py \\
                    -c ${params.decoypyrate_csites} \\
                    -a ${params.decoypyrate_non_csites} \\
                    -p ${params.decoypyrate_cpos} \\
                    -ml ${params.decoypyrate_min_length} \\
                    -mx ${params.decoypyrate_max_length} \\
                    -n ${params.decoypyrate_max_it} \\
                    -u ${params.decoypyrate_max_sub} \\
                    $flags \\
                    -o ${target.baseName}_TD${params.decoy_name}.fasta \\
                    $target \\
                    > ${task.process}_\${ts}_${target.baseName}_decoyPYrateV2.log
     """
}

process ANNOTATE_RESULTS {

    clusterOptions  '--job-name=nf_annotator --output=nf_annotator.txt '
    label 'Summary_Task'

    input: 
    path(mztabs)

    output: 
    path("${params.searchid}_Peptides.tsv"), emit: peptide_table
    path("${params.searchid}_PSMs.tsv"), emit: psm_table
    path("${params.searchid}_Proteins.tsv"), emit: protein_table
    path("${params.searchid}_Genes.tsv"), emit: gene_table
    path "${params.searchid}_*Annotated_Results.mzTab", emit: mztab optional true 

    tuple val("${params.searchid}"), path("${params.searchid}_PSMs.tsv"), path("${params.searchid}_Peptides.tsv"), path("${params.searchid}_Proteins.tsv"), path("${params.searchid}_Genes.tsv"), emit: summary
    path "*.log", emit: hidden // only for logging into the log directory
    path "*.tsv", emit: tables optional true 

    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.out_dir}/final_results", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.out_dir}/final_results", mode: 'copy', pattern: '*_Annotated_Results.mzTab'

    script: 

    flags = ""
    if (params.ms3_quant) flags = flags + "--ms3 "
    if (params.sample_info != "") flags = flags + "-s ${params.sample_info} "
    if (params.summary_write_mztab) flags = flags + "-mz "

    """
    ts=\$(date '+%Y%m%d_%H%M%S')

    python3 $projectDir/bin/Results_Annotator.py -d ${baseDir}/${params.out_dir}/tabs \\
                                        -r ${params.reference_fasta} \\
                                        -f ${params.peptide_fdr_threshold} \\
                                        -pm ${params.psm_pep_threshold} \\
                                        -fm ${params.psm_fdr_threshold} \\
                                        -o ${params.searchid}_ \\
                                        ${flags} \\
                                        > ${task.process}_\${ts}_${params.searchid}_Annotator.log

    """

}

//TMT_ReQuant(); runs the TMT MS2 or MS3 ReQuantification workflow on the peptide table and the original mzML files to attempt to rescue some of the missing values.
//This process is the same as TMT_REQUANT except it also takes the protein table as input and produces a quantified protein table as output. 
process TMT_QUANT {

    clusterOptions  '--job-name=nf_TMTquant --output=nf_TMTquant.txt '

    input:
    tuple val(id), path(psm_table), path(peptide_table), path(protein_table), path(gene_table)

    output:
    path("${params.searchid}_PeptideQuant.tsv"), emit: peptide_table
    path("${params.searchid}_SpectralTable.tsv"), emit: spectra_table
    path("${params.searchid}_ProteinQuant.tsv"), emit: protein_table optional true
    path("${params.searchid}_GeneQuant.tsv"), emit: gene_table optional true
    tuple val(id), path("${params.searchid}_PeptideQuant.tsv"), path("${params.searchid}_SpectralTable.tsv"), emit: quant_summary
    path "*.log", emit: hidden // only for logging into the log directory

    publishDir "${params.out_dir}/final_results", mode: 'copy', pattern: '*_PeptideQuant.tsv'
    publishDir "${params.out_dir}/final_results", mode: 'copy', pattern: '*_SpectralTable.tsv'
    publishDir "${params.out_dir}/final_results", mode: 'copy', pattern: '*_ProteinQuant.tsv'
    publishDir "${params.out_dir}/final_results", mode: 'copy', pattern: '*_GeneQuant.tsv'
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    script:

    mzmldir = "${params.mzml_dir}"
    if (params.mzml_merge) mzmldir = "${params.out_dir}/merged_mzml"

    ms_level = "2"
    if (params.ms3_quant) ms_level = "3"
    
    """
    ts=\$(date '+%Y%m%d_%H%M%S')

    python3 $projectDir/bin/TMT_Quantification_Tool_V1.py -i ${psm_table} \\
                                        -p ${peptide_table} \\
                                        -r ${protein_table} \\
                                        -g ${gene_table} \\
                                        -c $projectDir/${params.out_dir}/consensusXMLs \\
                                        -s $projectDir/${mzmldir} \\
                                        -f ${params.psm_fdr_threshold} \\
                                        -t ${params.psm_pep_threshold} \\
                                        -x ${params.tmt_channel_threshold} \\
                                        -m ${ms_level} \\
                                        -o ${id}_ \\
                                        > ${task.process}_\${ts}_${id}_TMT_Quant.log

    """

}

process MUTANT_PROTEOFORM_ANALYSIS {

    clusterOptions  '--job-name=nf_mutantProteoform --output=nf_mutantProteoform.txt '

    input: 
    tuple val(id), path(peptide_table), path(spectrum_table),  path(mutations)

    output: 
    path("${params.searchid}_Mutant_Proteoforms.tsv"), emit: mutant_proteoforms
    path "*.log", emit: hidden // only for logging into the log directory

    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.out_dir}/final_results", mode: 'copy', pattern: '*_Mutant_Proteoforms.tsv'

    script: 

    """
    ts=\$(date '+%Y%m%d_%H%M%S')

    python3 $projectDir/bin/Mutant_Proteoform_Analysis_V4.py --samplefile ${params.sample_info} \\
                                                             --pepfile ${peptide_table} \\
                                                             --mutfile ${mutations} \\
                                                             --fasta ${params.cosmic_aa_reference} \\
                                                             --spectrafile ${spectrum_table} \\
                                                             --control-channel-name ${params.control_channel_name} \\
                                                             --output ${params.searchid}_Mutant_Proteoforms.tsv \\
                                                            > ${task.process}_\${ts}_${params.searchid}_Mutant_Proteoform_Analysis.log

    """

}