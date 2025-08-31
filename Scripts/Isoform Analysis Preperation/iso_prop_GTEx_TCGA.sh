#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60g
#SBATCH --time=00:30:00
#SBATCH --job-name=iso_prop_fil
#SBATCH --output=/gpfs01/home/mbxjy4/logs/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/mbxjy4/logs/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbxjy4@exmail.nottingham.ac.uk

source $HOME/.bash_profile

conda activate R

Rscript iso_prop_GTEx_TCGA.R

conda deactivate
