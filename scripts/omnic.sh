#!/bin/bash --login
 
#---------------
#omnic.sh : creates an index file for the reference and an associated genome file, aligns the Omni-C library to the reference, creates a sorted pairs file and removes PCR duplicates
 
#---------------
#Requested resources:
#SBATCH --account=pawsey0812
#SBATCH --job-name=omnic
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

# Check if the directory exists, if not, create it
if [ ! -d "$temp" ]; then
    echo "Creating directory: $temp"
    mkdir -p "$temp"
fi

echo "Temp dir: $temp"

# Define Hi-C directory
hic_dir="/scratch/pawsey0964/lhuet/deep-sea/primary/fine_tuned_hifiasm"

# Define Hi-C files 
H1="${hic_dir}/OG849_hic_R1.fastq.gz"
for file in $H1; do
    echo "Hi-C forward: $file"
done

H2="${hic_dir}/OG849_hic_R2.fastq.gz"
for file in $H2; do
    echo "Hi-C reverse: $file"
done

# Define the assembly file 
assembly_file=purged.fa

# Remove the "_ctg.fasta" from the end of the string
output_file="OG849_purged_omnic"

# --------------
# Run the Omni-C pipeline



singularity run $SING/samtools_1.16.1.sif samtools faidx $assembly_file && \
cut -f1,2 *.fai > "${output_file}.genome" \
&& singularity run $SING/bwa:0.7.17.sif bwa index $assembly_file \
&& singularity run $SING/bwa:0.7.17.sif bwa mem -5SP -T0 -t64 $assembly_file $H1 $H2 -o "${output_file}.aligned.sam" \
&& singularity run $SING/pairtools:1.1.0.sif pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 32 --nproc-out 32 --chroms-path "${output_file}.genome" "${output_file}.aligned.sam" >  "${output_file}.parsed.pairsam" \
&& singularity run $SING/pairtools:1.1.0.sif pairtools sort --nproc 32 --tmpdir=$temp  "${output_file}.parsed.pairsam" > "${output_file}.sorted.pairsam" \
&& singularity run $SING/pairtools:1.1.0.sif pairtools dedup --nproc-in 32 --nproc-out 32 --mark-dups --output-stats "${output_file}.stats.txt" --output "${output_file}.dedup.pairsam" "${output_file}.sorted.pairsam" \
&& singularity run $SING/pairtools:1.1.0.sif pairtools split --nproc-in 32 --nproc-out 32 --output-pairs "${output_file}.mapped.pairs" --output-sam "${output_file}.unsorted.bam" "${output_file}.dedup.pairsam" \
&& singularity run $SING/samtools_1.16.1.sif samtools sort -@32 -T  "${temp}/${output_file}_temp.bam" -o "${output_file}.mapped.PT.bam" "${output_file}.unsorted.bam" \
&& singularity run $SING/samtools_1.16.1.sif samtools index "${output_file}.mapped.PT.bam"

## Each line in the above code performs the following:
## 1. Create an index for reference - samtools faidx
## 2. Use the index to generate a genome file - cut
## 3. Generate bwa index files - bwa index
## 4. Alignment - bwa mem
## 5. Recording valid ligation events - pairtools parse
## 6. Sorting the pairsam file - pairtools sort
## 7. Removig PCR duplicates - pairtools dedup
## 8. Generating .pairs and bam files
## 9. Generating the final bam file
## 10. Index the bam file


# --------------
# Cleanup: Delete unwanted files
rm -f "${output_file}.aligned.sam" "${output_file}.mapped.pairs" "${output_file}.genome" "${output_file}.pairsam" "${output_file}.parsed.pairsam" "${output_file}.sorted.pairsam" "${output_file}.unsorted.bam" "${assembly_file}.sa" "${assembly_file}.amb" "${assembly_file}.ann" "${assembly_file}.pac" "${assembly_file}.bwt"
echo "Processing complete for $output_file"


#---------------
#Successfully finished
echo "Done"
exit 0
