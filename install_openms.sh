#!/bin/bash
#
#SBATCH --job-name=mamba
#SBATCH --output=env_creator.txt
#SBATCH --partition=compute
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=8000
#SBATCH --time=10:00:00

source ~/.bash_profile

mamba env create -f openms-conda-enviroment.yml -p ./env/openms
