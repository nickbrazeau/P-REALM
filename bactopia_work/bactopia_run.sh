#!/usr/bin/env bash
set -euo pipefail

######################################################
#######  Nextflow Bactopia Pipeline          #########
######################################################

# ---------------------------
## project-specific configurations
# ---------------------------
indir="/hpc/group/taylorlab/users/nfb/projects/P-REALM/bactopia_work"
outdir="/work/nfb9/projects/P-REALM/bactopia_work"
REF="${indir}/ref/GCF_000013465.1_ASM1346v1_genomic.gbff"
DB="${indir}/baktadb/db-light"

# Configs and helper paths
MLST_CONFIG="${indir}/config/bactopia-mlst.config"
NF_CONFIG="${indir}/config/slurm-nextflow.config"
EXCLUDE_TSV="${indir}/bactopia-exclude.tsv"
LOGDIR="${indir}/logs"
PREALM_RET="${outdir}/prealm_ret"

# ---------------------------
# Environment
# ---------------------------
export TMPDIR="/work/nfb9/projects/P-REALM/tmp"
export PATH="/hpc/home/nfb9/miniconda3/envs/bactopia/bin:${PATH}"
export MPLCONFIGDIR="/work/nfb9/projects/P-REALM/tmp/matplotlib"
export NUMBA_DISABLE_JIT=1 # JIT off globally (Nextflow task-level config also sets this for Gubbins)
export APPTAINERENV_NUMBA_DISABLE_JIT=1
export NXF_APPTAINER_OPTS="--env NUMBA_DISABLE_JIT=1"
unset SINGULARITY_CACHEDIR NXF_SINGULARITY_CACHEDIR NXF_SINGULARITY_OPTS # ensure no singularity vars shadow apptainer

# dirs we need
mkdir -p "${LOGDIR}" "${outdir}" "${TMPDIR}" "${MPLCONFIGDIR}"


# -----------------------------------------------
# Pipeline Steps
# -----------------------------------------------
# ---------------------------
# Prep & Main
# ---------------------------
# prep
bactopia prepare \
--path "${indir}/fastq" \
--species "Staphylococcus aureus" \
--genome-size 2800000 \
> "${indir}/samples.fofn.txt"

# Main
# bactopia runs tools sequentially within its own run_pipeline, details below for flags I have incorporated:
# gather: increased coverage from default 100 to 300
# QC: adding fastqc and fastp (this is default)
#    adding scrubbr to remove human genome data
# Assembler: using Shovill-PE with default SKESA. Enable adapter trimming
# Annotator: using Bakta intead of Prokka default
# Sketcher: using default (really a sanity check)
# Sequence Typing: MLST module, default
# Antimicrobial Resistance: AMR module, default
# Merlin: uses sketches to run species-specific tools; letting it detect SA
bactopia \
-profile apptainer \
--samples "${indir}/samples.fofn.txt" \
--coverage 200 \
--shovill_assembler skesa \
--trim \
--use_bakta \
--bakta_db "${DB}" \
--nfconfig "${MLST_CONFIG}" \
--outdir "${PREALM_RET}" \
--max_cpus 32 \
--max_memory 400.GB \
2> "${LOGDIR}/bactopiamain.2out.txt"

# summary for exclude file
bactopia summary	--bactopia-path "${PREALM_RET}" \
--outdir "${indir}" \
2> "${LOGDIR}/bashreport.txt"


# ---------------------------
# Staph-specific tools
# ---------------------------
# staphtyper
bactopia \
  -profile apptainer \
  --wf staphtyper \
  --bactopia "${PREALM_RET}" \
  --nfconfig "${NF_CONFIG}" \
  --exclude "${EXCLUDE_TSV}" \
  2> "${LOGDIR}/staphtyper.2out.txt"

# staphopiasccmec
bactopia \
  -profile apptainer \
  --wf staphopiasccmec \
  --bactopia "${PREALM_RET}" \
  --nfconfig "${NF_CONFIG}" \
  --exclude "${EXCLUDE_TSV}" \
  2> "${LOGDIR}/staphopiasccmec.2out.txt"

# ---------------------------
# Genome tools
# ---------------------------
# Mobsuite
bactopia \
  -profile apptainer \
  --wf mobsuite \
  --bactopia "${PREALM_RET}" \
  --nfconfig "${NF_CONFIG}" \
  --exclude "${EXCLUDE_TSV}" \
  2> "${LOGDIR}/mobsuite.2out.txt"

# run_rgi
bactopia \
  -profile apptainer \
  --wf rgi \
  --bactopia "${PREALM_RET}" \
  --nfconfig "${NF_CONFIG}" \
  --exclude "${EXCLUDE_TSV}" \
  2> "${LOGDIR}/rgi.2out.txt"

# plasmid finder
 bactopia \
  -profile apptainer \
  --wf plasmidfinder \
  --bactopia "${PREALM_RET}" \
  --nfconfig "${NF_CONFIG}" \
  --exclude "${EXCLUDE_TSV}" \
  2> "${LOGDIR}/plasmidfinder.2out.txt"

# ---------------------------
# Modules
# ---------------------------
# snippy
bactopia \
  -profile apptainer \
  --wf snippy \
  --bactopia "${PREALM_RET}" \
  --nfconfig "${NF_CONFIG}" \
  --exclude "${EXCLUDE_TSV}" \
  --reference "${REF}" \
  2> "${LOGDIR}/snippy.2out.txt"

# mashtree
bactopia \
  -profile apptainer \
  --wf mashtree \
  --bactopia "${PREALM_RET}" \
  -qs 8 \
  --max_cpus 32 \
  --nfconfig "${NF_CONFIG}" \
  --exclude "${EXCLUDE_TSV}" \
  2> "${LOGDIR}/mashtree.2out.txt"

# pangenome panaroo
bactopia \
  -profile apptainer \
  --wf pangenome \
  --bactopia "${PREALM_RET}" \
  --exclude "${EXCLUDE_TSV}" \
  --nfconfig "${NF_CONFIG}" \
  --use_panaroo \
  2> "${LOGDIR}/pangenome_panaroo.2out.txt"

# pangenome pirate
bactopia \
  -profile apptainer \
  --wf pangenome \
  --bactopia "${PREALM_RET}" \
  --exclude "${EXCLUDE_TSV}" \
  --nfconfig "${NF_CONFIG}" \
  --use_pirate \
  2> "${LOGDIR}/pangenome_pirate.2out.txt"

# ---------------------------
# EOF
# ---------------------------
echo "[$(date)] All steps completed. Cheers 🍻!"
