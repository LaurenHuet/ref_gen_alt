#!/bin/bash --login
 
#---------------
#fcsgx.sh : runs the NCBI Foreign Contamination Screen (fcs) tool for identifying and removing contaminant sequences in genome assemblies 


#---------------
#Requested resources:
#SBATCH --account=pawsey0812
#SBATCH --job-name=fcs-gx
#SBATCH --partition=highmem
 
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=512G
 
#SBATCH --export=ALL
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err


# Define variables

fasta=/scratch/pawsey0964/lhuet/deep-sea/primary/OG849_pd/yahs/OG849.yahs_scaffolds_final.fa


# Path to the database
gxdb="/scratch/references/Foreign_Contamination_Screening/gxdb"

# Specify the number of cores (GitHub recommends 48 cores which is approx 24 CPUs)
GX_NUM_CORES=48


    echo 'copying database files to /tmp/'
    mkdir -p /tmp/gxdb/
    cp -v ${gxdb}/all.gxi /tmp/gxdb/
    cp -v ${gxdb}/all.gxs /tmp/gxdb/
    cp -v ${gxdb}/all.meta.jsonl /tmp/gxdb/
    cp -v ${gxdb}/all.blast_div.tsv.gz /tmp/gxdb/
    cp -v ${gxdb}/all.taxa.tsv /tmp/gxdb/
    echo 'done copying database files'
    ls -l /tmp/gxdb/



#---------------
# Run fcs-gx
python3 $SING/fcs.py --image=$SING/fcs-gx.sif screen genome --fasta $fasta --out-dir "./gx_out/${fasta}" --gx-db "/tmp/gxdb/" --tax-id 3362464


#---------------
#Successfully finished
echo "Done"
exit 0
