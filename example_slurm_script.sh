#!/bin/bash

#SBATCH --time=1000
#SBATCH --ntasks=6
#SBATCH --mem=90000
#SBATCH --job-name=GP3
#SBATCH --output=GP3.log
#SBATCH --error=GP3.err


module load <python module name>
module load <R module name>

cd <full path to your virtual env> 
source bin/activate
cd <full path to GP3 directory/repo>

time chunky run run_GWAS_analysis_pipeline.py -inputPLINK <full/path/to/PLINK.bed> -phenoFile <full/path/to/sample_spreadsheet.xlsx> --outDir <full/path/to/existing/directory/to/output/results> --projectName <pick a project name(no spaces)>
