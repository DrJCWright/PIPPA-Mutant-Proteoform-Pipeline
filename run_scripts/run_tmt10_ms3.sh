#!/bin/bash
#SBATCH --job-name=TMT10_MS3 --partition=master-worker --ntasks=1 --mem=4021 --output=TMT11_MS3_Log.txt --time=120:00:00

#module load is for ICR SLURM cluster. If you are running on your own machine, you can comment out this line.
#module load Nextflow

#Remove 'slurm' from the command below if you are running on your own machine.
nextflow run main.nf -profile slurm,TMT10_MS3 \
                    --searchid='ExampleTMT10MS3' \
					--mzml_merge_pattern='.*' \
					--mzml_dir='./mzml/' \
					--database_dir='./db/' \
					--out_dir='./results/' \
					-resume