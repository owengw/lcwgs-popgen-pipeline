# lcWGS Population Genomics Pipeline

A SLURM-based pipeline for population genomic analysis of low-coverage whole genome sequencing (lcWGS) data (1–5×). Designed for non-model organisms with fragmented genome assemblies. All analyses use genotype likelihood (GL) approaches via ANGSD to avoid biases introduced by hard genotype calling at low coverage.

Originally developed for *Plestiodon longirostris* (Bermuda skink), but applicable to any diploid organism with a reference genome.

---

## Overview

The pipeline is divided into numbered stages, each represented by a script. Scripts within the same stage can typically be run in parallel. Scripts in later stages depend on outputs from earlier stages.

```
p01                        Read QC and trimming
p02                        Genome indexing
p03                        Alignment
p04                        BAM merging, deduplication, overlap clipping
p05                        Depth assessment
p_sex                      Sex identification (optional)
p06                        GL generation, LD pruning, PCA, admixture, theta, FST
p07                        Individual heterozygosity, population theta, inbreeding, SFS
p08                        ROH, HWE, LD decay
p09                        lcMLkin relatedness (GL-aware, HWE-free)
plot_lcmlkin_relatedness.R Relatedness network plots and summary tables (local R)
Diff.R                     Selection analysis — Tajima's D and FST outliers, gene annotation, GO enrichment (local R)
```

---

## Requirements

### Software

| Tool | Version tested | Notes |
|---|---|---|
| ANGSD | 0.940 | Core GL engine |
| ngsLD | 1.2.1 | Use `/users/bi4og/ngsLD/ngsLD`, not conda version |
| PCAngsd | — | Via conda |
| NGSadmix | — | Via conda |

| ROHan | — | Compiled binary required |
| lcMLkin | v2.1 | Python script, cloned from GitHub |
| realSFS / thetaStat | (ANGSD suite) | |
| samtools | — | Via conda |
| bowtie2 | — | Via conda |
| Trimmomatic | 0.40 | Via conda |
| Picard | — | Via conda |
| bamUtil | — | Via conda (clipOverlap) |
| PLINK | 1.9 | Community genomics bin |
| bcftools | — | Via conda |
| R | 4.4.1 | System module |
| Python | 3.10.8 | System module (for lcMLkin) |

### R packages

- `HardyWeinberg` — installed automatically by `p08c_hwe.sh` to a user library
- Base R only required for all other R steps (no tidyverse dependency)

### Python packages (for lcMLkin)

```bash
pip install --user "numpy<2.0" pandas scipy
```

> **Note:** numpy ≥ 1.25 removes `np.warnings` — use numpy < 2.0.

### Cluster environment

Scripts are written for SLURM on a shared HPC with:
- A conda environment at `/mnt/parscratch/users/<user>/conda_envs/owenspopgen`
- Working directory: `/mnt/parscratch/users/<user>/`
- All scripts sourced from a `scripts/` subdirectory
- Logs written to a `logs/` subdirectory

**Before running:** update all hardcoded paths (search for `/mnt/parscratch/users/bi4og/` and replace with your own).

---

## Directory structure

```
project_root/
├── scripts/           # All pipeline scripts
├── logs/              # SLURM log output
├── raw_data/          # Raw FASTQ files
├── genome/            # Reference genome FASTA and indices
├── physalia/
│   ├── adapter_clipped/   # Trimmed FASTQs
│   ├── aligned/           # Per-lane BAMs
│   ├── merged/            # Per-individual merged BAMs
│   ├── prefix_list.txt    # List of FASTQ prefixes
│   └── angsd_results/     # All downstream outputs
│       ├── pop_map_lists/     # Per-population BAM lists
│       ├── theta_corrected/   # Population theta/pi
│       ├── theta_downsampled/ # Downsampled theta (equal n)
│       ├── heterozygosity_corrected/  # Individual het
│       ├── inbreeding_corrected/      # F_HET
│       ├── relatedness/              # ngsRelate output
│       ├── hwe/                      # HWE results
│       ├── ld_decay/                 # LD decay curves
│       ├── roh/                      # ROHan output
│       │   ├── individual/
│       │   └── population/
│       ├── ngsadmix/                 # NGSadmix results
│       └── lcmlkin/                  # lcMLkin results
└── metadata.tsv       # Sample metadata (see format below)
```

---

## Metadata file

Many scripts require a `metadata.tsv` at `physalia/angsd_results/metadata.tsv`. Required columns:

| Column | Description | Example |
|---|---|---|
| `sample_id` | Unique sample identifier matching BAM prefix | `285-CAI14` |
| `population` | Population group label | `pop_map_CAI` |
| `site` | Short site code | `CAI` |
| `year` | Sampling year | `2023` |

---

## Population BAM lists

Scripts expect per-population BAM lists at `physalia/angsd_results/pop_map_lists/pop_map_<POP>.list`, one absolute BAM path per line. Also required:

- `physalia/angsd_results/all_samples.list` — all individuals combined

These must be created manually before running p06 onwards.

---

## Stage-by-stage guide

### p01 — Read trimming

```bash
# Array job — one task per FASTQ prefix
sbatch --array=1-N scripts/p01_trim_array.sh \
    -p raw_data \
    -l physalia/prefix_list.txt \
    -f _R1.fastq.gz \
    -r _R2.fastq.gz \
    -k ILLUMINACLIP:scripts/NEBNext.fa:2:30:10 \
    -s SLIDINGWINDOW:4:20 \
    -L LEADING:3 \
    -t TRAILING:3 \
    -m MINLEN:36
```

**Adapt:** set `-p` to your raw data directory, `-l` to your prefix list, `-f`/`-r` to your FASTQ extensions, `-k` to your adapter file. Adjust `--array` to match the number of prefixes.

For SE-only libraries, only provide `-f`. The script auto-detects PE vs SE.

---

### p02 — Genome indexing

```bash
sbatch scripts/p02_index_genome.sh genome/your_reference.fasta
```

Runs samtools faidx, Picard CreateSequenceDictionary, and bowtie2-build. Only needs to run once.

---

### p03 — Alignment

```bash
sbatch --array=1-N scripts/p03_align_array.sh \
    -g genome/your_reference.fasta \
    -l physalia/prefix_list.txt \
    -p physalia/adapter_clipped
```

Aligns with bowtie2 `--very-sensitive`, filters MAPQ ≥ 20, outputs sorted BAMs.

---

### p04 — Merge, dedup, overlap clip

```bash
sbatch --array=1-N scripts/p04_merge_bam_array.sh \
    -q physalia/prefix_list.txt \
    -i physalia/aligned \
    -o physalia/merged \
    -r your_reference_name
```

Merges across lanes, marks and removes duplicates with Picard, clips overlapping read pairs with bamUtil. Output BAMs are used in all downstream analyses. Adjust `--array` to match the number of unique biological samples (not prefix count).

---

### p05 — Depth assessment

```bash
sbatch --array=1-N scripts/p05_depth.sh \
    -l physalia/sample_lists/bam_list_dedup_overlapclipped.list \
    -o physalia/depths
```

Calculates per-site depth for each individual. Use the output to set appropriate `setMinDepthInd` and `setMaxDepthInd` values in downstream ANGSD runs.

---

### p_sex — Sex identification (optional)

```bash
sbatch scripts/p_sex_identification.sh
```

Estimates sex by comparing mean depth on putative sex chromosome contigs vs autosomes. Requires a RagTag AGP file from scaffolding. Edit the script to set the correct chromosome and contig paths for your organism.

---

### p06 — GL generation, structure, diversity

This is the core ANGSD stage. Run scripts in the following order.

#### p06a — Genome-wide GL (for LD estimation)

```bash
sbatch --export=BAM_LIST=physalia/angsd_results/all_samples.list,\
REFERENCE=genome/your_reference.fasta,\
OUT_DIR=physalia/angsd_results \
scripts/p06a_all_samples.sh
```

Generates `all_samples.beagle.gz` across all individuals genome-wide. Used as input for LD pruning. Adjust `MIN_IND_RATIO` (default 0.20) and `MIN_MAPQ`/`MIN_Q` as needed.

#### p06b — LD estimation and pruning

```bash
sbatch --export=OUT_DIR=physalia/angsd_results,\
REFERENCE=genome/your_reference.fasta \
scripts/p06b_ld.sh
```

Runs ngsLD on a random subsample of sites (default 100,000), then graph-based LD pruning at r² ≥ 0.5 within 2 Mb. Produces:
- `LDpruned_snps.list` — LD-pruned site list for ANGSD `-sites`
- `LDpruned_contigs.txt` — contig restriction file for ANGSD `-rf`
- `all_samples_unlinked.id` — unlinked site IDs

> **Note:** LD decay results are unreliable with highly fragmented assemblies. Report as supplementary with caveat if contig N50 < 1 Mb.

#### p06c — Per-population SAF generation

```bash
for POP in pop_map_CAI pop_map_NS pop_map_SB pop_map_SI; do
    sbatch --export=BAM_LIST=physalia/angsd_results/pop_map_lists/${POP}.list,\
REFERENCE=genome/your_reference.fasta,\
OUT_DIR=physalia/angsd_results,\
USE_RANDOM_SITES=1 \
    scripts/p06c_angsd_subset.sh
done
```

Generates per-population SAF index files on random 10M sites (unbiased for theta/FST). Used by p06d and p06e. Set `USE_RANDOM_SITES=0` to use full genome instead.

#### p06c2 — PCA beagle file (LD-pruned sites)

```bash
sbatch --export=BAM_LIST=physalia/angsd_results/pop_map_lists/pop_map_ALL.list,\
REFERENCE=genome/your_reference.fasta,\
OUT_DIR=physalia/angsd_results \
scripts/p06c2_pca.sh
```

Generates beagle GL file on LD-pruned sites for PCA and admixture. Also produces IBS matrix.

#### p06c3 — PCAngsd PCA

```bash
sbatch --export=BAM_LIST=physalia/angsd_results/pop_map_lists/pop_map_ALL.list,\
REFERENCE=genome/your_reference.fasta,\
OUT_DIR=physalia/angsd_results,\
THREADS=40 \
scripts/p06c3_pcangsd.sh
```

Runs PCAngsd on the LD-pruned beagle file. Output: `pcangsd_ALL.cov` covariance matrix. Use `--iter 200` for convergence in structured populations (edit script if needed).

#### p06d — Theta and Tajima's D

```bash
for POP in pop_map_CAI pop_map_NS pop_map_SB pop_map_SI; do
    sbatch --export=BAM_LIST=physalia/angsd_results/pop_map_lists/${POP}.list,\
REFERENCE=genome/your_reference.fasta,\
OUT_DIR=physalia/angsd_results \
    scripts/p06d_merge_theta.sh
done
```

Runs realSFS and thetaStat per population. Produces windowed π, θw, and Tajima's D. **Critical:** realSFS stdout must go directly to the `.sfs` file — do not use `2>&1 | tee` as this corrupts the SFS.

#### p06e — Pairwise FST

```bash
sbatch --export=OUT_DIR=physalia/angsd_results scripts/p06e_fst.sh
```

Calculates pairwise weighted FST between all population pairs using 2D-SFS. Temporal comparisons (year × population) are only run if global temporal FST exceeds threshold (default 0.015).

#### p06f — Summary tables

```bash
sbatch --export=OUT_DIR=physalia/angsd_results scripts/p06f_summary.sh
```

Collates theta and FST results into summary TSV files.

#### p06g — NGSadmix

```bash
sbatch --array=1-200 --export=OUT_DIR=physalia/angsd_results scripts/p06g_ngsadmix.sh
```

Runs NGSadmix for K=1–10, 20 replicates each (200 array jobs). Best K determined by log-likelihood stability across replicates. Adjust array size if changing K range or replicate count.

---

### p07 — Diversity and inbreeding

#### p07a — Individual heterozygosity

```bash
sbatch --array=1-N \
    --export=OUT_DIR=physalia/angsd_results,\
BAM_LIST=physalia/angsd_results/all_samples.list \
    scripts/p07a_individual_het.sh
```

Estimates per-individual heterozygosity using 1-sample folded SFS via realSFS. Adjust `--array` to match the number of individuals. Excludes failed or duplicate samples by editing the `EXCLUDE` pattern in the script.

#### p07b — Population theta

```bash
# All populations
sbatch --export=OUT_DIR=physalia/angsd_results scripts/p07b_pop_theta.sh

# Single population (useful for reruns)
sbatch --export=OUT_DIR=physalia/angsd_results,POP_SINGLE=SI scripts/p07b_pop_theta.sh
```

Runs the full theta pipeline (SAF → SFS → saf2theta → thetaStat) per population. Collates results to `theta_corrected/theta_summary_corrected.tsv`. DH (n=1) is excluded by default — edit `POPULATIONS` array if your dataset differs.

#### p07c — Inbreeding coefficients (F_HET)

```bash
sbatch --export=OUT_DIR=physalia/angsd_results scripts/p07c_inbreeding.sh
```

Calculates F_HET = 1 − (H_obs / H_exp) per individual. Requires p07a and p07b to be complete. H_expected is taken from corrected population theta_pi.

#### p07d — Downsampled theta

```bash
sbatch --export=OUT_DIR=physalia/angsd_results scripts/p07d_theta_downsampled.sh
```

Repeats theta estimation with equal sample sizes (default n=17) to remove sample size bias from Tajima's D comparisons. Adjust `DOWNSAMPLE_N` and `POPULATIONS` for your dataset. Random seed = 42 (report in methods for reproducibility).

#### p07e — SFS summary

```bash
sbatch --export=OUT_DIR=physalia/angsd_results scripts/p07e_sfs.sh
```

Summarises singleton proportion and segregating site counts from existing SFS files.

---

### p08 — ROH, HWE, LD decay

#### p08a — ROHan (individual ROH)

```bash
sbatch --array=1-N scripts/p08a_rohan_individual.sh
```

Runs ROHan on each individual BAM using `--tstv 2.1` (Ts/Tv ratio mode, appropriate for non-model organisms). Adjust `--array` and `EXCLUDE` pattern. Requires HTSlib and GSL modules. ROHan binary path must be set in script.

> **Note:** ROH detection is conservative with fragmented assemblies. ROH spanning contig boundaries are undetectable. ROHan theta (estimated outside ROH) is more reliable than ROH metrics in this case.

#### p08b — ROHan collation

```bash
sbatch --export=OUT_DIR=physalia/angsd_results scripts/p08b_rohan_population.sh
```

Collates individual ROHan `.summary.txt` files into population-level summaries. Run after all p08a array tasks complete.

#### p08c — Hardy-Weinberg equilibrium

```bash
# All populations
sbatch --export=OUT_DIR=physalia/angsd_results scripts/p08c_hwe.sh

# Single population
sbatch --export=OUT_DIR=physalia/angsd_results,POP_SINGLE=SI scripts/p08c_hwe.sh
```

Generates per-population beagle GL files with ANGSD, then tests HWE per site using the `HardyWeinberg` R package (GL-based exact test). Beagle files are shared with p08d. Reports per-site F, uncorrected p-values, and Bonferroni-corrected counts.

> **Interpretation note:** Bonferroni correction is conservative when testing only a few thousand sites. The proportion of sites with p < 0.05 (expected 5% under HWE) is a more informative population-level summary than per-site significance.

#### p08d — LD decay

```bash
sbatch --export=OUT_DIR=physalia/angsd_results scripts/p08d_ld_decay.sh
```

Estimates within-population LD decay using ngsLD on the beagle files from p08c. Requires the ngsLD binary at `/users/bi4og/ngsLD/ngsLD` (v1.2.1) and a `libgsl.so.27` symlink (created automatically). Positions file must be gzipped and tab-separated — the script handles this conversion from the beagle marker column.

> **Note:** LD decay results are unreliable with highly fragmented genome assemblies where most contigs are shorter than the detection window. Report with caveat or omit if assembly N50 < 100 kb.

---

### p09 — lcMLkin relatedness

lcMLkin v2.1 is preferred over ngsRelate KING for highly inbred populations because it uses a likelihood approach that does not assume Hardy-Weinberg equilibrium.

#### p09a — Per-population VCF generation

```bash
sbatch --array=1-4 \
    --mem=128G \
    --export=OUT_DIR=physalia/angsd_results \
    scripts/p09a_angsd_vcf.sh
```

Generates per-population VCFs from ANGSD (`-dogeno 1 -doBcf 1`). Array tasks map to POPULATIONS=(CAI NS SB SI) — edit the array for your populations. Memory requirements scale with sample size; SI (n=84) required 128G.

> **Setup note:** lcMLkin v2.1 requires several patches to work with non-standard chromosome names and modern numpy. See the lcMLkin setup section below.

#### p09 — lcMLkin

```bash
sbatch --export=OUT_DIR=physalia/angsd_results \
    --mem=64G \
    --time=24:00:00 \
    scripts/p09_lcmlkin.sh
```

For each population: filters VCF to biallelic SNPs with MAF > 0.05, subsamples to ~50,000 SNPs (random seed 42), strips PL field (lcMLkin uses GL; PL causes errors on missing data), converts to plink with CHROM_POS IDs, pre-generates allele frequency file, then runs lcMLkin from the population output directory.

Output columns: `Z0ag` (K0), `Z1ag` (K1), `Z2ag` (K2), `PI_HATag` (r).

---

## lcMLkin setup

lcMLkin v2.1 requires the following one-time setup:

```bash
# Clone
git clone https://github.com/altinisik/lcMLkin-v2.1.git /path/to/lcMLkin-v2.1

# Install Python dependencies (numpy must be <2.0)
module load Python/3.10.8-GCCcore-12.2.0
pip install --user "numpy<2.0" pandas scipy

# Patch 1: np.warnings compatibility (removed in numpy >= 1.25)
sed -i 's/np\.warnings\.filterwarnings/import warnings; warnings.filterwarnings/' \
    lcMLkin-v2.1/lcmlkinv2.py

# Patch 2: Replace plink calls with full path + --allow-extra-chr
# (required for non-standard chromosome names)
sed -i "s|Popen('plink |Popen('/path/to/plink |g" lcmlkinv2.py
sed -i "s|--bfile '+filenameinplink+' --allow-no-sex|--bfile '+filenameinplink+' --allow-extra-chr --allow-no-sex|g" lcmlkinv2.py
# ... (see script comments for full patch list)

# Patch 3: Fix snpldtest file path (basename strips directory)
# Line ~389: change os.path.basename(filenameout) to filenameout

# Patch 4: Fix GQ handling when PL field is absent
# Line ~189: handle case where GQ_key == -9 and PL_key == -9
```

These patches are applied automatically by `p09_lcmlkin.sh` when it validates the environment.

---

## Key parameter decisions

| Parameter | Default | Notes |
|---|---|---|
| `MIN_IND_RATIO` | 0.20–0.30 | Adjust based on expected missing data rate |
| `minMapQ` | 10–20 | Use 20 for stricter alignment filtering |
| `minQ` | 20 | Base quality threshold |
| `SNP_pval` | 1e-6 | SNP calling threshold for beagle/SFS |
| `FOLD_SFS` | 1 | Folded SFS (no outgroup); set 0 if outgroup available |
| `WINDOW` | 50,000 bp | Theta/FST window size |
| `STEP` | 25,000 bp | Window step size |
| `DOWNSAMPLE_N` | 17 | Target n for downsampled theta (set to smallest population) |
| `--tstv` | 2.1 | ROHan Ts/Tv ratio; adjust for your taxon |
| LD pruning r² | 0.5 | Threshold for LD pruning |
| lcMLkin MAF | 0.05 | Minimum allele frequency for relatedness |
| lcMLkin SNPs | ~50,000 | After subsampling; sufficient for relatedness |

---

## Excluded samples

Edit the `EXCLUDE` pattern in each script to match your sample naming convention. The default pattern excludes:

```bash
EXCLUDE="o_merged\|T_merged\|SI43_merged\|SI45_merged\|SI83_merged"
```

- `o_merged` — original samples kept for QC comparison (duplicate pairs)
- `T_merged` — tissue samples kept for QC comparison
- `SI43`, `SI45` — failed samples (< 1,000 sites)
- `SI83` — confirmed duplicate of SI77

Replace these with your own failed/duplicate sample identifiers.

---

## Output files — key results

| File | Description |
|---|---|
| `angsd_results/LDpruned_snps.list` | LD-pruned SNP sites for PCA/admixture |
| `angsd_results/pcangsd_ALL.cov` | PCAngsd covariance matrix |
| `angsd_results/ngsadmix/loglikelihoods.txt` | NGSadmix log-likelihoods by K and rep |
| `angsd_results/theta_corrected/theta_summary_corrected.tsv` | π, θw, Tajima's D per population |
| `angsd_results/theta_downsampled/theta_downsampled_summary.tsv` | Equal-n theta comparison |
| `angsd_results/fst_summary.tsv` | Pairwise weighted FST |
| `angsd_results/heterozygosity_corrected/*.het` | Individual heterozygosity |
| `angsd_results/inbreeding_corrected/F_HET_individual_corrected.txt` | Per-individual F_HET |
| `angsd_results/hwe/hwe_summary_all_populations.tsv` | HWE summary per population |
| `angsd_results/roh/population/population_ROH_summary.tsv` | ROHan population summary |
| `angsd_results/roh/population/individual_ROH_all.tsv` | ROHan individual results |
| `angsd_results/ld_decay/ld_decay_all_populations.tsv` | LD decay by distance class |
| `angsd_results/lcmlkin/lcmlkin_all_populations.tsv` | lcMLkin pairwise relatedness (all populations) |
| `results/relatedness_lcmlkin/lcmlkin_relationship_summary.tsv` | Relationship class counts per population |
| `results/relatedness_lcmlkin/lcmlkin_close_pairs.tsv` | All non-unrelated pairs with K0/K1/K2/r |
| `results/relatedness_lcmlkin/lcmlkin_network_<POP>.pdf` | Per-population relatedness network |
| `results/relatedness_lcmlkin/lcmlkin_network_combined.pdf` | Combined four-population network panel |
| `results/selection difference/plestiodon_TajimaD_manhattan.pdf` | Tajima's D Manhattan plot |
| `results/selection difference/plestiodon_FST_upset.pdf` | FST outlier gene UpSet plot |
| `results/selection difference/plestiodon_FST_manhattan_pairs.pdf` | Per-pair FST Manhattan plots |
| `results/selection difference/plestiodon_FST_outlier_gene_products.csv` | FST outlier genes with product names |
| `results/selection difference/plestiodon_TajimaD_outlier_gene_products.csv` | Tajima's D outlier genes with product names |

---

### plot_lcmlkin_relatedness.R — Relatedness visualisation

Run locally in RStudio after downloading `lcmlkin_all_populations.tsv` from the cluster. Requires no cluster access — pure R.

**Before running**, update `base_dir` and `lcmlkin_dir` at the top of the script to point to your local copy of the lcMLkin results.

The script produces:

- A summary table of relationship class counts per population (Duplicate, 1st degree, Full sibling, 2nd degree, 3rd degree)
- A bar chart of close pair counts per population
- Per-population kinship network plots — nodes are individuals, edges are pairs with r ≥ 0.125 (3rd degree or closer), edge colour encodes relationship class, edge width encodes r
- A combined four-population network panel

**Why lcMLkin over ngsRelate KING for this pipeline:** KING assumes Hardy-Weinberg equilibrium and underestimates relatedness in inbred populations. lcMLkin uses a likelihood framework that directly models inbreeding, making it more appropriate for small isolated populations with elevated F values, even at low coverage.

**R packages required:**

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "igraph",
                   "stringr", "patchwork", "scales"))
```

**Outputs** (saved to `results/relatedness_lcmlkin/`):
- `lcmlkin_relationship_summary.tsv` — relationship class counts per population
- `lcmlkin_close_pairs.tsv` — all non-unrelated pairs with K0, K1, K2, r
- `lcmlkin_network_<POP>.pdf/.png` — per-population network
- `lcmlkin_network_combined.pdf` — all populations combined
- `lcmlkin_relationship_counts.pdf/.png` — bar chart summary

---

### Diff.R — Selection analysis

Run locally in RStudio after downloading results from the cluster. Requires:
- `pop_map_ALL.thetas.windowed.pestPG` — genome-wide windowed Tajima's D (from p07b with `POP_SINGLE=ALL`)
- `fst_*.fst.windowed.tsv` files — per-pair windowed FST (from p06e)
- `ragtag.scaffold.agp` — AGP scaffolding file mapping contigs to chromosomes (from RagTag)
- A GFF3 annotation file for your species

**Before running**, update the file paths at the top of the script under `# --- 2. FILE PATHS ---` to point to your local copies of these files.

The script performs five analyses in sequence:

**1. Tajima's D outlier detection** — identifies windows in the top and bottom 1% of the genome-wide Tajima's D distribution as candidates for balancing and positive selection respectively. Windows are translated from contig coordinates to chromosome coordinates via the AGP file.

**2. FST outlier detection** — identifies windows in the top 1% of FST per population pair as candidates for divergent selection. Runs across all pairwise comparisons simultaneously.

**3. Gene annotation overlap** — uses GenomicRanges to find annotated genes overlapping outlier windows (±25 kb), extracts gene product names from the GFF3, and identifies genes flagged in multiple population comparisons.

**4. GO enrichment** — tests for enrichment of Gene Ontology terms among outlier gene sets using clusterProfiler. Requires GO annotations in the GFF3 (Dbxref field).

**5. Visualisation** — produces three figures:
- Tajima's D Manhattan plot across all chromosomes, with outlier windows highlighted
- FST UpSet plot showing overlap of divergent genes across all population pairs
- Per-pair FST Manhattan plots (one panel per comparison, combined with patchwork)

**R packages required:**

```r
install.packages(c("ggplot2", "dplyr", "ggrepel", "data.table",
                   "tidyr", "stringr", "conflicted", "patchwork", "UpSetR"))
BiocManager::install(c("rtracklayer", "GenomicRanges", "GO.db", "clusterProfiler"))
```

**Outputs** (saved to the path set in `OUT_DIR`):
- `plestiodon_TajimaD_manhattan.pdf/.png`
- `plestiodon_FST_upset.pdf/.png`
- `plestiodon_FST_manhattan_pairs.pdf/.png`
- `plestiodon_FST_outlier_gene_products.csv`
- `plestiodon_TajimaD_outlier_gene_products.csv`

> **Note:** The script uses chromosome-level coordinates from RagTag scaffolding for Manhattan plots. If your assembly is not scaffolded to chromosome level, plots will show per-contig results instead — remove the AGP translation step (sections 4–5) and plot directly from contig coordinates.

---



## Known limitations

- **Fragmented assembly:** ROH detection, LD decay, and within-population IBD analysis are unreliable when assembly N50 is much shorter than the expected ROH length or LD window. Report with caveats or exclude from main results.
- **Bonferroni HWE:** Conservative with few tested sites per population. Report uncorrected proportion significant alongside Bonferroni results.
- **KING bias in inbred populations:** KING-based relatedness estimators assume Hardy-Weinberg equilibrium and are biased downward in inbred populations. lcMLkin is preferred for small isolated populations with significant HWE deviation as it uses a likelihood approach that directly accounts for inbreeding.
- **lcMLkin SNP requirement:** Fewer than ~10,000 SNPs per population pair will produce unreliable K0/K1/K2 estimates. Ensure the per-population ANGSD VCF has sufficient SNP density (typically 5–8 million SNPs before filtering, yielding ~45,000–55,000 after MAF filtering and subsampling).
- **realSFS stdout corruption:** Never pipe realSFS through `tee` — log messages are written to stdout and corrupt the `.sfs` file. Always redirect stderr separately (`2> file.log`).

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **ANGSD:** Korneliussen et al. (2014) BMC Bioinformatics
- **ngsLD:** Fox et al. (2019) Bioinformatics
- **PCAngsd:** Meisner & Albrechtsen (2018) Genetics
- **NGSadmix:** Skotte et al. (2013) Genetics
- **ngsRelate:** Hanghøj et al. (2019) GigaScience — retained in pipeline for cross-validation but lcMLkin preferred as primary relatedness estimator for inbred populations
- **ROHan:** Renaud et al. (2019) PLOS Genetics
- **lcMLkin:** Lipatov et al. (2015) Bioinformatics; Žegarac et al. (2021) Bioinformatics
- **HardyWeinberg R package:** Graffelman & Morales-Camarena (2008) Human Heredity
- **Trimmomatic:** Bolger et al. (2014) Bioinformatics
- **bowtie2:** Langmead & Salzberg (2012) Nature Methods
- **clusterProfiler:** Wu et al. (2021) The Innovation
- **ggplot2:** Wickham (2016) Springer

---

## Contact

Pipeline developed for PhD research at Manchester Metropolitan University. For questions about adapting the pipeline to new datasets, open a GitHub issue.
