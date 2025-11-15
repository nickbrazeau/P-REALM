#!/bin/bash
#SBATCH --job-name=prealm_bactopia
#SBATCH --output=logs/bactopia_%j.out
#SBATCH --error=logs/bactopia_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=72:00:00
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
cd /hpc/group/taylorlab/users/nfb/projects/P-REALM/bactopia_work
mkdir -p logs
bash bactopia_run.sh
