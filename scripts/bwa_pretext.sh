#!/bin/bash --login
 
#---------------
#omnic.sh : creates an index file for the reference and an associated genome file, aligns the Omni-C library to the reference, creates a sorted pairs file and removes PCR duplicates
 
#---------------
#Requested resources:
#SBATCH --account=pawsey0812
#SBATCH --job-name=bwa_pretext
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=15:00:00
#SBATCH --mem=64G
#SBATCH --export=ALL
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

date=$(date +%y%m%d)

echo "========================================="
echo "SLURM_JOB_ID = $SLURM_JOB_ID"
echo "SLURM_NODELIST = $SLURM_NODELIST"
echo "DATE: $date"
echo "========================================="

# Define path to temp directory
temp="/scratch/pawsey0964/lhuet/deep-sea/primary/OG849_pd"
[ -d "$temp" ] || { echo "Creating directory: $temp"; mkdir -p "$temp"; }
echo "Temp dir: $temp"

# Define Hi-C directory
hic_dir="/scratch/pawsey0964/lhuet/deep-sea/primary/OG849_pd/omnic_2"

# Define Hi-C files 
H1="${hic_dir}/OG849_hic_R1.fastq.gz"
echo "Hi-C forward: $H1"
H2="${hic_dir}/OG849_hic_R2.fastq.gz"
echo "Hi-C reverse: $H2"

# Define the assembly file 
assembly_file=OG849.yahs.tiara.scaffolds.fasta
output_file="OG849.yahs.tiara.omnic"

# Map Omni-C to the EXACT FASTA you plan to view
singularity run $SING/bwa:0.7.17.sif bwa index $assembly_file

# FIX: make a proper pipeline (mem → view -b -h → sort)
singularity run $SING/bwa:0.7.17.sif bwa mem -t 48 -5SP -T0 $assembly_file $H1 $H2 \
  | singularity run $SING/samtools_1.16.1.sif samtools view -b -h - \
  | singularity run $SING/samtools_1.16.1.sif samtools sort -@8 -o OG849.yahs.tiara.omnic.hic.bam -

singularity run $SING/samtools_1.16.1.sif samtools index OG849.yahs.tiara.omnic.hic.bam

# Build Pretext (light filtering)
# -F 0x904 removes unmapped(0x4) + secondary(0x100); keep supplementary
singularity run "$SING/samtools_1.16.1.sif" \
  samtools view -h -F 0x904 -q 1 OG849.yahs.tiara.omnic.hic.bam \
  | singularity run "$SING/pretextmap:0.1.9.sif" PretextMap -o OG849.yahs.tiara.omnic.hic.pretext
