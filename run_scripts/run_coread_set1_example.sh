#!/bin/bash
#SBATCH --job-name=COREAD_SET1 --partition=master-worker --ntasks=1 --mem=4021 --output=COREAD_SET1_TMT10_MS3_Log.txt --time=120:00:00

module load Nextflow

cd /data/scratch/DBI/DUDBI/FUNCPROT/jwright/Search_OpenMSV3/pippa-mutant-proteoforms-pipeline

nextflow run main.nf -profile slurm,TMT10MS3 \
                    --searchid='COREAD_SET1' \
                    --sample_info='/data/scratch/DBI/DUDBI/FUNCPROT/jwright/Search_OpenMSV3/pippa-mutant-proteoforms-pipeline/projects/COREAD_MS3/data/COREAD_SET1_sample_info.tsv' \
                    --mutation_table='/data/scratch/DBI/DUDBI/FUNCPROT/jwright/Search_OpenMSV3/pippa-mutant-proteoforms-pipeline/data/CellLinesProject_GenomeScreensMutant_v103_GRCh38.tsv' \
                    --cosmic_aa_reference='/data/scratch/DBI/DUDBI/FUNCPROT/jwright/Search_OpenMSV3/pippa-mutant-proteoforms-pipeline/data/gencode.v47.pc_translations.fa' \
                    --cosmic_na_reference='/data/scratch/DBI/DUDBI/FUNCPROT/jwright/Search_OpenMSV3/pippa-mutant-proteoforms-pipeline/data/gencode.v47.pc_transcripts.fa' \
                    --reference_fasta='/data/scratch/DBI/DUDBI/FUNCPROT/jwright/Search_OpenMSV3/pippa-mutant-proteoforms-pipeline/data/gencode.v47.pc_translations.fa' \
                    --control_channel_name='SW-48' \
					--mzml_merge_pattern='COREAD_MS3_SET\d+' \
					--mzml_dir='./projects/COREAD_MS3/mzml/' \
					--database_dir='./projects/COREAD_MS3/db/MPDB/' \
					--out_dir='./projects/COREAD_MS3/results_mp/' \
					-resume