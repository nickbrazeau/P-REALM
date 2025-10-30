project/
 ├─ data/
 │   ├─ raw/               

 │   ├─ derived/               

 │   ├─ meta/              # sample sheets, matching keys
 │   └─ refs/              # reference fasta, annotation, gene panels
 ├─ pipelines/
 │   ├─ bactopia/          # assemblies, pangenome, QC
 │   ├─ variants/          # vcfs, cohort matrices

├─ analysis/          # stats outputs (csv, rds, png)
 │   └─ models/            # trained models, metrics
 ├─ configs/               # yaml for refs, filters, thresholds
 ├─ scripts/               # R/Python/bash entrypoints
 ├─ prompts/               # Claude tool prompts & test cases
 ├─ Makefile               # or Snakemake/justfile
 └─ README.md