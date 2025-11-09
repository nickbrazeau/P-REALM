#!/bin/bash
#SBATCH --job-name=prealm_bactopia
#SBATCH --output=logs/p_realm_%j.out
#SBATCH --error=logs/p_realm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nick.brazeau@duke.edu
#SBATCH --export=ALL,APPTAINER_CACHEDIR=/hpc/group/taylorlab/containers/apptainer_cache,NXF_APPTAINER_CACHEDIR=/hpc/group/taylorlab/containers/apptainer_cache
# (not: NXF_SINGULARITY_CACHEDIR; keep Singularity unset to avoid dual-engine)

# 🔒 Disable Numba JIT *globally* (container will inherit this)
export NUMBA_DISABLE_JIT=1
export APPTAINERENV_NUMBA_DISABLE_JIT=1
export NXF_APPTAINER_OPTS="--env NUMBA_DISABLE_JIT=1"

# Make sure Singularity vars are not fighting Apptainer
unset SINGULARITY_CACHEDIR NXF_SINGULARITY_CACHEDIR NXF_SINGULARITY_OPTS


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
