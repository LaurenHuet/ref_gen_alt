#!/bin/bash --login

#---------------
#hifiasm.sh: runs HiFiasm - a fast haplotype-resolved de novo assembler for PacBio HiFi reads

#---------------
#Requested resources:
#SBATCH --account=pawsey0812
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


reads=/scratch/pawsey0964/lhuet/deep-sea/OG849_OG849_m84154_241127_094241_s2.hifi_reads.bc2009.filt.fastq.gz


THREADS=48
singularity run $SING/hifiasm:0.25.0.sif hifiasm \
  -o OG849 -t $THREADS --primary \
  --h1 OG849_hic_R1.fastq.gz --h2 OG849_hic_R2.fastq.gz \
  $reads 2>&1 | tee OG849.asm.log
