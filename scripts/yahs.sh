#!/bin/bash --login

#---------------
#Requested resources:
#SBATCH --account=pawsey0812
#SBATCH --job-name=purge_dups
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=05:00:00
#SBATCH --mem=120G
#SBATCH --export=ALL
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err


bam=/scratch/pawsey0964/lhuet/deep-sea/primary/OG849_pd/OG849_purged_omnic.mapped.PT.bam
fasta=/scratch/pawsey0964/lhuet/deep-sea/primary/OG849_pd/purged.fa

singularity run $SING/yahs:1.2a.2.sif yahs --no-contig-ec -o OG849.yahs $fasta $bam
