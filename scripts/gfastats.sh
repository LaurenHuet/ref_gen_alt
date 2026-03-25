#!/bin/bash --login
#---------------
# gfastats.sh: converts GFA to FASTA and calculates assembly summary statistics
# Processes all *ctg.gfa files in the working directory
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

# Estimated genome size in base pairs (use GenomeScope2 estimate)
GENOME_SIZE=525039910

# gfastats container version
GFASTATS_SIF=gfastats:1.3.6.sif

# ============================================================
# END USER CONFIGURATION
# ============================================================

set -euo pipefail

echo "=== gfastats: genome size ${GENOME_SIZE} bp | $(date) ==="

#---------------
# Convert GFA to FASTA
echo "--- Converting GFA to FASTA ---"
for gfa in *ctg.gfa; do
    [[ -f "$gfa" ]] || { echo "No *ctg.gfa files found"; exit 1; }
    fasta="${gfa%.*}.fasta"
    echo "  $gfa -> $fasta"
    singularity run "$SING/${GFASTATS_SIF}" gfastats --discover-paths "$gfa" -o fa > "$fasta"
done

#---------------
# Calculate summary statistics
echo "--- Calculating assembly statistics ---"
for gfa in *ctg.gfa; do
    [[ -f "$gfa" ]] || continue
    stats="${gfa%.*}.assembly.summary.txt"
    echo "  Stats: $stats"
    singularity run "$SING/${GFASTATS_SIF}" gfastats "$gfa" "${GENOME_SIZE}" \
        --discover-paths --tabular --nstar-report > "$stats"
done

echo "=== gfastats complete: $(date) ==="
