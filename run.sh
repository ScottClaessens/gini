#!/usr/bin/env bash

#SBATCH --job-name=gini
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10-00:00:00
#SBATCH --mem=30G
#SBATCH --account=arch039044

module load languages/R/4.5.1
Rscript run.R
