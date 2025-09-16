#!/bin/bash --login
#---------------
#gfastats.sh : converts gfa to fasta format, and calculates assembly summary statistics
#---------------
#Requested resources:
#SBATCH --account=pawsey0812
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


#---------------
# Convert gfa to fasta

FA=OG849.yahs_scaffolds_final.fa
LINEAGE=/scratch/references/busco_db/actinopterygii_odb10
T=48

singularity run "$SING/gfastats:1.3.10.sif" gfastats -f "$FA" > OG849.yahs.gfastats.txt
