#!/bin/bash --login
#---------------
# gfastats2.sh: calculates assembly statistics for a single FASTA file
# Use this after scaffolding (e.g. yahs output) or on any final assembly FASTA
#---------------
#SBATCH --account=pawsey0964
#SBATCH --job-name=gfastats
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --export=ALL
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err


# ============================================================
# USER CONFIGURATION - edit these variables before submitting
# ============================================================

# Sample name / output prefix
SAMPLE=OG849

# Input FASTA file
FA=OG849.yahs_scaffolds_final.fa

# Output stats file (defaults to <SAMPLE>.gfastats.txt)
OUT="${SAMPLE}.gfastats.txt"

# gfastats container version
GFASTATS_SIF=gfastats:1.3.10.sif

# ============================================================
# END USER CONFIGURATION
# ============================================================

set -euo pipefail

echo "=== gfastats2: ${FA} -> ${OUT} | $(date) ==="

singularity run "$SING/${GFASTATS_SIF}" gfastats -f "$FA" > "$OUT"

echo "=== gfastats2 complete: $(date) ==="
