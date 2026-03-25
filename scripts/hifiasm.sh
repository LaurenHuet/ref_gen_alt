#!/bin/bash --login

#---------------
# hifiasm.sh: runs HiFiasm - a fast haplotype-resolved de novo assembler for PacBio HiFi reads
# Supports primary assembly mode or dual-haplotype mode (hap1 + hap2)
#---------------
#SBATCH --account=pawsey0964
#SBATCH --job-name=hifiasm_primary
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

# PacBio HiFi reads (fastq.gz)
READS_HIFI=/scratch/pawsey0964/lhuet/deep-sea/OG849_OG849_m84154_241127_094241_s2.hifi_reads.bc2009.filt.fastq.gz

# Hi-C reads (paired-end, fastq.gz)
HIC_R1=OG849_hic_R1.fastq.gz
HIC_R2=OG849_hic_R2.fastq.gz

# Assembly mode:
#   primary    - produce a single primary assembly (--primary flag); faster, less memory
#   haplotypes - produce phased hap1 + hap2 assemblies; requires Hi-C
MODE=primary

# Threads (should match --cpus-per-task above)
THREADS=48

# ============================================================
# END USER CONFIGURATION
# ============================================================

set -euo pipefail

echo "=== hifiasm: ${SAMPLE} | mode: ${MODE} | $(date) ==="

if [[ "${MODE}" == "primary" ]]; then
    singularity run $SING/hifiasm:0.25.0.sif hifiasm \
        -o "${SAMPLE}" -t "${THREADS}" --primary \
        --h1 "${HIC_R1}" --h2 "${HIC_R2}" \
        "${READS_HIFI}" 2>&1 | tee "${SAMPLE}.asm.log"

elif [[ "${MODE}" == "haplotypes" ]]; then
    singularity run $SING/hifiasm:0.25.0.sif hifiasm \
        -o "${SAMPLE}" -t "${THREADS}" \
        --h1 "${HIC_R1}" --h2 "${HIC_R2}" \
        "${READS_HIFI}" 2>&1 | tee "${SAMPLE}.asm.log"

else
    echo "ERROR: MODE must be 'primary' or 'haplotypes'" >&2
    exit 1
fi

echo "=== hifiasm complete: $(date) ==="
