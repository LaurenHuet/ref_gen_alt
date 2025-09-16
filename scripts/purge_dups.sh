#!/bin/bash --login

#---------------
#Requested resources:
#SBATCH --account=pawsey0812
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


# === INPUTS ===
READS_HIFI=/scratch/pawsey0964/lhuet/deep-sea/OG849_OG849_m84154_241127_094241_s2.hifi_reads.bc2009.filt.fastq.gz
FASTA=/scratch/pawsey0964/lhuet/deep-sea/primary/OG849.hic.p_ctg.fasta
OUT=OG849_pd                               

# === THREADS ===
T=48

mkdir -p "${OUT}"
cd "${OUT}"


# Long-read → contig mapping for depth (use map-hifi preset for PacBio HiFi)
singularity run $SING/minimap2:2.26.sif minimap2 -x map-hifi -t ${T} $FASTA ${READS_HIFI} | gzip -c > pb.paf.gz


# Compute depth stats
singularity run $SING/purge_dups:1.2.6.sif pbcstat pb.paf.gz -O pd

# Infer cutoffs from depth distribution
singularity run $SING/purge_dups:1.2.6.sif calcuts pd/PB.stat > cutoffs 2> calcuts.log

# Split contigs for sensitivity
singularity run $SING/purge_dups:1.2.6.sif split_fa $FASTA > contigs.split.fa

# Self PAF (asm5 is standard for assembly-vs-assembly)
singularity run $SING/minimap2:2.26.sif minimap2 -x asm5 -DP -t ${T} contigs.split.fa contigs.split.fa | gzip -c > contigs.self.paf.gz

# Call duplicated/haplotig segments
singularity run $SING/purge_dups:1.2.6.sif purge_dups -2 -T cutoffs -c pd/PB.base.cov contigs.self.paf.gz > dups.bed

# Extract purged primary assembly (-e removes dups; produces purged.fa and hap.fa)
singularity run $SING/purge_dups:1.2.6.sif get_seqs -e dups.bed $FASTA > purged.fa
singularity run $SING/purge_dups:1.2.6.sif get_seqs dups.bed $FASTA > purged_hap.fa   # optional: the removed bits

# # Size/N50 shift
singularity run $SING/seqkit:2.8.2.sif seqkit stats $FASTA purged.fa purged_hap.fa



# Primary (kept) assembly, with internal trims realized as breaks
singularity run $SING/purge_dups:1.2.6.sif \
  get_seqs -s -e dups.bed $FASTA > purged.fa

# Haplotig/removed fragments (now materialized)
singularity run $SING/purge_dups:1.2.6.sif \
  get_seqs -s dups.bed $FASTA > purged_hap.fa

# Sanity check
singularity run $SING/seqkit:2.8.2.sif seqkit stats $FASTA purged.fa purged_hap.fa
