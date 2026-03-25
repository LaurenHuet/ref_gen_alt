# Reference Genome Assembly Pipeline

Primary assembly from PacBio HiFi + Hi-C reads using HiFiasm, with duplicate purging via purge_dups.

## Overview

```
HiFi reads + Hi-C reads
        |
   [1] hifiasm.sh          → GFA assembly files (primary or hap1/hap2)
        |
   [2] gfastats.sh         → Convert GFA → FASTA, assembly statistics
        |
   [3] GenomeScope2         → Coverage model, genome size estimate, cutoffs
   (manual step, no script)
        |
   [4] purge_dups.sh       → Purge haplotigs/duplications
        |
   [5] gfastats2.sh        → Final assembly statistics
```

---

## Prerequisites

- `$SING` environment variable must point to your Singularity/Apptainer image directory
- Required containers:
  - `hifiasm:0.25.0.sif`
  - `gfastats:1.3.6.sif` and/or `gfastats:1.3.10.sif`
  - `minimap2:2.26.sif`
  - `purge_dups:1.2.6.sif`
  - `seqkit:2.8.2.sif`

---

## Step 1: Assembly with HiFiasm

**Script:** `scripts/hifiasm.sh`

Edit the USER CONFIGURATION block at the top:

```bash
SAMPLE=MySpecies            # output prefix
READS_HIFI=/path/to/hifi.fastq.gz
HIC_R1=/path/to/hic_R1.fastq.gz
HIC_R2=/path/to/hic_R2.fastq.gz
MODE=primary                # 'primary' or 'haplotypes'
THREADS=48
```

**Mode options:**
- `primary` — produces a single collapsed primary assembly (`*.hic.p_ctg.gfa`). Faster and uses less memory. Recommended as a starting point.
- `haplotypes` — produces phased `hap1` and `hap2` assemblies (`*.hic.hap1.p_ctg.gfa`, `*.hic.hap2.p_ctg.gfa`). Requires Hi-C data for phasing.

Submit:
```bash
sbatch scripts/hifiasm.sh
```

**Outputs (primary mode):**
```
MySpecies.hic.p_ctg.gfa     # primary contig graph
MySpecies.hic.a_ctg.gfa     # alternate contigs
MySpecies.asm.log
```

**Outputs (haplotypes mode):**
```
MySpecies.hic.hap1.p_ctg.gfa
MySpecies.hic.hap2.p_ctg.gfa
MySpecies.asm.log
```

---

## Step 2: GFA to FASTA Conversion and Assembly Statistics

**Script:** `scripts/gfastats.sh`

Run from the directory containing the GFA files. Edit the USER CONFIGURATION block:

```bash
GENOME_SIZE=525039910       # estimated genome size in bp (use GenomeScope2 estimate)
```

Submit:
```bash
sbatch scripts/gfastats.sh
```

**Outputs** (for each `*ctg.gfa`):
```
*.fasta                     # FASTA sequences
*.assembly.summary.txt      # N50, L50, contig counts, etc.
```

---

## Step 3: GenomeScope2 — Coverage Model and Purge_Dups Cutoffs

> **This step is critical for accurate purge_dups results.**

GenomeScope2 models the k-mer frequency distribution of your HiFi reads to estimate genome size, heterozygosity, and ploidy. The coverage peaks it identifies are used to set accurate thresholds for purge_dups.

### 3a. Count k-mers with Jellyfish

```bash
jellyfish count -C -m 21 -s 1G -t 48 -o reads.jf <(zcat hifi_reads.fastq.gz)
jellyfish histo -t 48 reads.jf > reads.histo
```

### 3b. Run GenomeScope2

Upload `reads.histo` to http://genomescope.org/genomescope2/ (or run locally):

```bash
genomescope2 -i reads.histo -o genomescope_out -k 21 -p 2
```

### 3c. Reading the GenomeScope2 output for purge_dups cutoffs

Open `genomescope_out/model.txt` or inspect the plot. Identify:

| Value | Description | purge_dups parameter |
|-------|-------------|----------------------|
| Heterozygous peak (1n) | Coverage at the haplotig peak — roughly half the main peak | `CALCUTS_MID` |
| Homozygous peak (2n) | Coverage at the main diploid peak | — |
| Lower bound | ~1/3 of the het peak; excludes noise/contamination | `CALCUTS_LOW` |
| Upper bound | ~1.5–2× the hom peak; excludes collapsed repeats | `CALCUTS_HIGH` |

**Example:** GenomeScope2 reports het peak at 30×, hom peak at 60×:
```bash
CALCUTS_LOW=10
CALCUTS_MID=30
CALCUTS_HIGH=120
```

If you leave these blank, `calcuts` will try to infer thresholds automatically from the depth histogram, but this is less reliable — especially for heterozygous or repeat-rich genomes.

---

## Step 4: Purge Duplications

**Script:** `scripts/purge_dups.sh`

Edit the USER CONFIGURATION block. **Set the GenomeScope2 cutoffs from Step 3** for best results:

```bash
SAMPLE=MySpecies
MODE=primary                # 'primary' or 'haplotypes'

READS_HIFI=/path/to/hifi_reads.fastq.gz

# Primary mode
FASTA_PRIMARY=/path/to/MySpecies.hic.p_ctg.fasta

# Haplotype mode (used when MODE=haplotypes)
FASTA_HAP1=/path/to/MySpecies.hic.hap1.p_ctg.fasta
FASTA_HAP2=/path/to/MySpecies.hic.hap2.p_ctg.fasta

THREADS=48

# From GenomeScope2 (recommended)
CALCUTS_LOW=10
CALCUTS_MID=30
CALCUTS_HIGH=120
```

Submit:
```bash
sbatch scripts/purge_dups.sh
```

**Outputs (primary mode):**
```
MySpecies_pd_primary/
├── pb.paf.gz               # HiFi → assembly alignments
├── pd/PB.stat              # per-read depth
├── pd/PB.base.cov          # per-base coverage
├── cutoffs                 # coverage thresholds used
├── calcuts.log             # calcuts log
├── contigs.split.fa        # contigs split at N gaps
├── contigs.self.paf.gz     # self-alignment
├── dups.bed                # duplicated intervals
├── purged.fa               # purged primary assembly (use this)
└── purged_hap.fa           # removed haplotigs
```

**Outputs (haplotypes mode):** same structure under `MySpecies_pd_hap1/` and `MySpecies_pd_hap2/`.

### Checking purge_dups results

Inspect the seqkit stats output in the job log. A successful run typically shows:
- `purged.fa` total bases slightly smaller than input (haplotigs removed)
- `purged_hap.fa` contains the removed sequences
- N50 of `purged.fa` should be equal to or higher than the input

Also check `cutoffs` to confirm the thresholds are sensible relative to your GenomeScope2 peaks.

---

## Step 5: Final Assembly Statistics

**Script:** `scripts/gfastats2.sh`

Run on any final FASTA (e.g. after scaffolding with yahs). Edit the USER CONFIGURATION block:

```bash
SAMPLE=MySpecies
FA=MySpecies.yahs_scaffolds_final.fa
OUT=MySpecies.gfastats.txt
```

Submit:
```bash
sbatch scripts/gfastats2.sh
```

---

## Quick Reference: Pipeline Commands

```bash
# 1. Assembly
sbatch scripts/hifiasm.sh

# 2. GFA → FASTA + stats (run from hifiasm output directory)
cd /path/to/hifiasm_output
sbatch /path/to/scripts/gfastats.sh

# 3. GenomeScope2 (manual — see Step 3 above)

# 4. Purge duplications
sbatch scripts/purge_dups.sh

# 5. Final stats
sbatch scripts/gfastats2.sh
```

---

## Notes

- All scripts use `set -euo pipefail` — they will exit on the first error.
- The `$SING` environment variable must be set in your environment or loaded via a module.
- GenomeScope2 cutoffs are strongly recommended for heterozygous genomes; auto-inference works best for highly homozygous samples.
- For haplotype-mode purge_dups, each haplotype is purged independently. The two purged assemblies together should approximately equal one haploid genome size.
