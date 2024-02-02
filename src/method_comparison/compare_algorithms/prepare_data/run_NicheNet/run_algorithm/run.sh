#!/bin/bash
#SBATCH -p slim18
#SBATCH -c 30
#SBATCH -o nichenet_analysis.log
#SBATCH -e nichenet_analysis.log

source /work/project/ladcol_011/conda_path/miniconda3/bin/activate NicheNet


Rscript run_NN.r

