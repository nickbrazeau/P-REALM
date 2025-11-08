#!/bin/bash
#SBATCH --job-name=prealm_bactopia
#SBATCH --output=logs/p_realm_%j.out
#SBATCH --error=logs/p_realm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=96:00:00
#SBATCH --export=ALL,APPTAINER_CACHEDIR=/hpc/group/taylorlab/containers/apptainer_cache,NXF_SINGULARITY_CACHEDIR=/hpc/group/taylorlab/containers/apptainer_cache,NXF_WORK=/work/nfb9/projects/P-REALM/nxf_work,TMPDIR=/work/nfb9/projects/P-REALM/tmp,NUMBA_CACHE_DIR=/work/nfb9/projects/P-REALM/tmp/numba_cache,XDG_CACHE_HOME=/work/nfb9/projects/P-REALM/tmp/.cache
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nick.brazeau@duke.edu

# Load environment
source ~/.bashrc
#micromamba env create -f env/environment.yml
micromamba activate prealm

# run snake
snakemake \
    --snakefile bactopia_run.snake \
    --cores 24 \
    --printshellcmds \
    --latency-wait 90 \
    --rerun-incomplete \
    --use-conda \
    --keep-going \
