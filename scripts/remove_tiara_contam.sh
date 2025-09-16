#!/bin/bash

# Make a file with the list of unwanted scaffolds
cat > to_remove.txt <<EOF
scaffold_585
scaffold_586
scaffold_771
EOF

# Exclude those scaffolds from your assembly
singularity run $SING/seqkit:2.8.2.sif seqkit grep -v -f to_remove.txt OG849.yahs_scaffolds_final.fa > OG849.yahs.tiara.scaffolds.fasta

