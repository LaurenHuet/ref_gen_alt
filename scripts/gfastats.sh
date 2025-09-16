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

# Iterate over the files in the directory
for filename in *ctg.gfa; do
    if [ -f "$filename" ]; then
        # Construct the output file path by replacing the extension with ".fasta"
        output_file="${filename%.*}.fasta"
        echo "Output file: $output_file"
        
        # Execute the gfastats command to convert GFA to FASTA
        singularity run $SING/gfastats:1.3.6.sif gfastats --discover-paths "$filename" -o fa > "$output_file"
        echo "Converted $filename to $output_file"

    fi
done


#---------------
# Calculate summary statistics
gsize=525039910
# Iterate over the files in the directory
for filename in *ctg.gfa; do
    # Construct the output file path by replacing the extension with ".assembly.summary.txt"
    output_file="${filename%.*}.assembly.summary.txt"

    # Execute the gfastats command to calculate assembly statistics
    singularity run $SING/gfastats:1.3.6.sif gfastats "$filename" $gsize --discover-paths --tabular --nstar-report > "$output_file"
done
