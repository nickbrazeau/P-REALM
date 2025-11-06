#!/bin/bash
#SBATCH --job-name=prealm_bactopia
#SBATCH --output=logs/p_realm_%j.out
#SBATCH --error=logs/p_realm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nick.brazeau@duke.edu

# Load environment
source ~/.bashrc
conda init bash
conda env create -f env/environment.yml
conda activate prealm

# run snake
snakemake \
    --snakefile bactopia_run.snake \
    --cores 8 \
    --printshellcmds \
    --latency-wait 90 \
    --rerun-incomplete \
    --use-conda \
    --keep-going \
    --default-resources /work/nfb9/projects/P-REALM/bactopia_work/ \
    --workdir /work/nfb9/projects/P-REALM/bactopia_work/
