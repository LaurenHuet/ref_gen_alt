#!/bin/bash --login

#---------------
# purge_dups.sh: removes haplotypic duplications from genome assemblies
#
# Supports two modes:
#   primary    - purge a single primary assembly (p_ctg)
#   haplotypes - purge hap1 and hap2 independently
#
# GenomeScope2 cutoffs: run GenomeScope2 on your HiFi reads BEFORE purging.
# Set CALCUTS_LOW/MID/HIGH from the GenomeScope2 coverage histogram to ensure
# accurate purging thresholds rather than relying on automatic inference.
# See README.md for instructions on reading GenomeScope2 output.
#---------------
#SBATCH --account=pawsey0964
#SBATCH --job-name=purge_dups
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=10:00:00
#SBATCH --mem=150G
#SBATCH --export=ALL
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err


# ============================================================
# USER CONFIGURATION - edit these variables before submitting
# ============================================================

# Sample name / output prefix
SAMPLE=OG849

# Purging mode:
#   primary    - purge single primary assembly
#   haplotypes - purge hap1 and hap2 separately
MODE=primary

# PacBio HiFi reads (used for read depth calculation)
READS_HIFI=/scratch/pawsey0964/lhuet/deep-sea/OG849_OG849_m84154_241127_094241_s2.hifi_reads.bc2009.filt.fastq.gz

# --- Primary mode inputs ---
FASTA_PRIMARY=/scratch/pawsey0964/lhuet/deep-sea/primary/OG849.hic.p_ctg.fasta

# --- Haplotype mode inputs ---
FASTA_HAP1=/scratch/pawsey0964/lhuet/deep-sea/haplotypes/OG849.hic.hap1.p_ctg.fasta
FASTA_HAP2=/scratch/pawsey0964/lhuet/deep-sea/haplotypes/OG849.hic.hap2.p_ctg.fasta

# Threads (should match --cpus-per-task above)
THREADS=48

# --- GenomeScope2 coverage cutoffs ---
# Recommended: set these from your GenomeScope2 results for accurate purging.
# Leave ALL THREE blank to use automatic calcuts inference (less accurate).
#
# How to read GenomeScope2 output:
#   CALCUTS_LOW  = lower coverage cutoff; use ~1/3 of the heterozygous (1n) peak
#                  (filters out low-coverage noise/contamination)
#   CALCUTS_MID  = mid coverage cutoff; use the heterozygous (1n) peak value
#                  (separates haplotigs from true primary contigs)
#   CALCUTS_HIGH = upper coverage cutoff; use ~1.5-2x the homozygous (2n) peak
#                  (filters out collapsed repeats / high-copy regions)
#
# Example: if GenomeScope2 shows het peak at 30x and hom peak at 60x:
#   CALCUTS_LOW=10
#   CALCUTS_MID=30
#   CALCUTS_HIGH=120
CALCUTS_LOW=""
CALCUTS_MID=""
CALCUTS_HIGH=""

# ============================================================
# END USER CONFIGURATION
# ============================================================

set -euo pipefail

# Build calcuts options string if manual cutoffs provided
build_calcuts_opts() {
    local opts=""
    if [[ -n "${CALCUTS_LOW}" && -n "${CALCUTS_MID}" && -n "${CALCUTS_HIGH}" ]]; then
        opts="-l ${CALCUTS_LOW} -m ${CALCUTS_MID} -u ${CALCUTS_HIGH}"
        echo "  Using manual GenomeScope2 cutoffs: low=${CALCUTS_LOW} mid=${CALCUTS_MID} high=${CALCUTS_HIGH}"
    else
        echo "  Using automatic calcuts inference (set CALCUTS_LOW/MID/HIGH for better accuracy)"
    fi
    echo "$opts"
}

# Run purge_dups pipeline on a single assembly FASTA
# Usage: run_purge_dups <fasta> <output_dir>
run_purge_dups() {
    local FASTA="$1"
    local OUT="$2"

    echo "--- Purging: ${FASTA} -> ${OUT}/ ---"
    mkdir -p "${OUT}"
    pushd "${OUT}" > /dev/null

    # Map HiFi reads to assembly for depth estimation
    echo "  [1/6] Mapping HiFi reads (minimap2 map-hifi)..."
    singularity run "$SING/minimap2:2.26.sif" minimap2 \
        -x map-hifi -t "${THREADS}" "${FASTA}" "${READS_HIFI}" \
        | gzip -c > pb.paf.gz

    # Compute per-base and per-read depth statistics
    echo "  [2/6] Computing depth statistics (pbcstat)..."
    singularity run "$SING/purge_dups:1.2.6.sif" pbcstat pb.paf.gz -O pd

    # Infer coverage cutoffs (use GenomeScope2 values if provided)
    echo "  [3/6] Inferring coverage cutoffs (calcuts)..."
    local calcuts_opts
    calcuts_opts=$(build_calcuts_opts)
    singularity run "$SING/purge_dups:1.2.6.sif" calcuts \
        $calcuts_opts pd/PB.stat > cutoffs 2> calcuts.log
    echo "  Cutoffs:"
    cat cutoffs

    # Split contigs at N gaps for improved haplotig detection
    echo "  [4/6] Splitting contigs (split_fa)..."
    singularity run "$SING/purge_dups:1.2.6.sif" split_fa "${FASTA}" > contigs.split.fa

    # Self-alignment to identify duplicated regions
    echo "  [5/6] Self-alignment (minimap2 asm5)..."
    singularity run "$SING/minimap2:2.26.sif" minimap2 \
        -x asm5 -DP -t "${THREADS}" contigs.split.fa contigs.split.fa \
        | gzip -c > contigs.self.paf.gz

    # Identify duplicated/haplotig intervals
    echo "  [6/6] Calling duplicates (purge_dups)..."
    singularity run "$SING/purge_dups:1.2.6.sif" purge_dups \
        -2 -T cutoffs -c pd/PB.base.cov contigs.self.paf.gz > dups.bed

    # Extract purged primary assembly and haplotigs
    echo "  Extracting sequences (get_seqs)..."
    singularity run "$SING/purge_dups:1.2.6.sif" \
        get_seqs -s -e dups.bed "${FASTA}" > purged.fa

    singularity run "$SING/purge_dups:1.2.6.sif" \
        get_seqs -s dups.bed "${FASTA}" > purged_hap.fa

    # Assembly statistics before and after purging
    echo "  Assembly statistics:"
    singularity run "$SING/seqkit:2.8.2.sif" seqkit stats \
        "${FASTA}" purged.fa purged_hap.fa

    popd > /dev/null
    echo "  Done: ${OUT}/purged.fa"
}


# ============================================================
# Main
# ============================================================
echo "=== purge_dups: ${SAMPLE} | mode: ${MODE} | $(date) ==="

if [[ "${MODE}" == "primary" ]]; then
    run_purge_dups "${FASTA_PRIMARY}" "${SAMPLE}_pd_primary"

elif [[ "${MODE}" == "haplotypes" ]]; then
    run_purge_dups "${FASTA_HAP1}" "${SAMPLE}_pd_hap1"
    run_purge_dups "${FASTA_HAP2}" "${SAMPLE}_pd_hap2"

else
    echo "ERROR: MODE must be 'primary' or 'haplotypes'" >&2
    exit 1
fi

echo "=== purge_dups complete: $(date) ==="
