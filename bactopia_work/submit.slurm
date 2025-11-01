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
conda env create -f env/environment.yml
conda activate prealm

# run
snakemake \
    --snakefile bactopia_run.snake \
    --cores 16 \
    --jobs 20 \
    --latency-wait 90 \
    --rerun-incomplete \
    --use-conda \
    --cluster "sbatch --nodes=1 --ntasks=1 --cpus-per-task={threads} --mem={resources.mem_mb}M --time={resources.time}" \
    --cluster-config config/cluster.yaml \
    --keep-going \
    --directory /hpc/group/taylorlab/users/nfb/projects/P-REALM/
