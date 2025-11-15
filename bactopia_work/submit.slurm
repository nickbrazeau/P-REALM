#!/bin/bash
#SBATCH --job-name=prealm_bactopia
#SBATCH --output=logs/p_realm_%j.out
#SBATCH --error=logs/p_realm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=512G
#SBATCH --time=96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nick.brazeau@duke.edu
#SBATCH --export=ALL,APPTAINER_CACHEDIR=/hpc/group/taylorlab/containers/apptainer_cache,NXF_APPTAINER_CACHEDIR=/hpc/group/taylorlab/containers/apptainer_cache

# env
source ~/.bashrc
micromamba activate prealm

# JIT off globally (Nextflow task-level config also sets this for Gubbins)
export NUMBA_DISABLE_JIT=1
export APPTAINERENV_NUMBA_DISABLE_JIT=1
export NXF_APPTAINER_OPTS="--env NUMBA_DISABLE_JIT=1"

# ensure no singularity vars shadow apptainer
unset SINGULARITY_CACHEDIR NXF_SINGULARITY_CACHEDIR NXF_SINGULARITY_OPTS

# run
snakemake \
  --snakefile bactopia_run.snake \
  --cores 48 \
  --printshellcmds \
  --latency-wait 90 \
  --rerun-incomplete \
  --keep-going
