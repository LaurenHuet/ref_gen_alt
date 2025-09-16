#!/bin/bash --login
 
#---------------
#tiara.sh : runs the program Tiara (v.1.0.3) to search for additional decontamination in the assemblies, particularity Mitochondrial 
 
#---------------
#Requested resources:
#SBATCH --account=pawsey0812
#SBATCH --job-name=tiara
#SBATCH --partition=work

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1:00:00
#SBATCH --mem=48G

#SBATCH --export=ALL
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

date=$(date +%y%m%d)

echo "========================================="
echo "SLURM_JOB_ID = $SLURM_JOB_ID"
echo "SLURM_NODELIST = $SLURM_NODELIST"
echo "DATE: $date"
echo "========================================="

#---------------
# Load modules
module load singularity/3.11.4-nompi


#---------------
# Define variables

# specify the output directory
out_dir="tiara"

#---------------
# Run tiara

fasta=/scratch/pawsey0964/lhuet/deep-sea/primary/OG849_pd/yahs/OG849.yahs_scaffolds_final.fa

sample=OG849.yahs.tiara

# Run tiara
singularity run $SING/tiara:1.0.3.sif tiara -i "$fasta" -o "$sample.tiara.txt" -m 1000 --tf mit pla pro -t 4 -p 0.65 0.60 --probabilities



