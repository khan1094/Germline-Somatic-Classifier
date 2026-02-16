# Somatic Variant Classifier

A lightweight, research-oriented command-line tool for classifying genomic variants as **LIKELY_SOMATIC**, **LIKELY_GERMLINE**, **CONFLICTING**, or **UNKNOWN** by integrating evidence from **Mutect2**, **gnomAD**, and **COSMIC**.

This project was implemented as part of a technical task and emphasizes **correctness**, **interpretability**, and **testability** over production-scale optimization.

---

## Overview

In cancer genomics, distinguishing somatic mutations from germline variants is critical for downstream interpretation.  
This tool applies simple, explainable scoring rules based on:

- **Population frequency** (gnomAD)
- **Cancer recurrence evidence** (COSMIC)

Each variant in a Mutect2 VCF is annotated with external evidence, assigned germline and somatic scores, and classified using a configurable decision matrix.

---

## Features

- Command-line interface with sensible defaults
- Deterministic, interpretable scoring (no machine learning)
- Allele-level processing of multi-allelic variants
- Unit-tested core classification logic
- Designed for small-to-medium research datasets

---

## Input & Output

### Input

- **Sample VCF**: Mutect2 output (tumor or control)
- **gnomAD VCF**: Population allele frequencies
- **COSMIC TSV (gzipped)**: Cancer recurrence database

### Output

A tab-separated (TSV) file with one row per variant allele.

| Column | Description |
|------|------------|
| chrom | Chromosome |
| pos | Genomic position |
| ref | Reference allele |
| alt | Alternate allele |
| filter | Mutect2 FILTER field |
| sample_dp | Sample read depth |
| sample_af | Sample allele frequency |
| gnomad_af | gnomAD allele frequency |
| gnomad_an | gnomAD allele number |
| cosmic_count | Total COSMIC sample count |
| cosmic_tissues | COSMIC tissue types |
| germline_score | Calculated germline evidence score |
| somatic_score | Calculated somatic evidence score |
| classification | Final classification |

---

## Classification Logic

### Germline Evidence Score

```
germline_score = AF * min(AN / AN_CONFIDENCE_THRESHOLD, 1.0)
```

- High allele frequency suggests germline origin
- Allele number (AN) is used to weight confidence in AF

### Somatic Evidence Score

```
somatic_score = log10(cosmic_count + 1)
```

- Log scaling handles the wide dynamic range of COSMIC recurrence counts

### Decision Matrix

| Germline Score | Somatic Score | Classification |
|---------------|---------------|----------------|
| High | Low | LIKELY_GERMLINE |
| Low | High | LIKELY_SOMATIC |
| High | High | CONFLICTING |
| Low | Low | UNKNOWN |

Thresholds are configurable via command-line arguments.

---

## Installation

### Requirements

- Python ≥ 3.9
- `cyvcf2==0.31.4`

Install dependencies:

```bash
pip install -r requirements.txt
```

### Platform Notes

This tool was developed and tested on **Linux (Linux Mint 22.2 / Ubuntu-based)**.

`cyvcf2` depends on **htslib**:

- cyvcf2 ≥ 0.20.0 requires **htslib ≥ 1.10**

#### Linux
On most Linux systems, htslib is available via system packages or conda and works out-of-the-box.

#### Windows
On Windows, `cyvcf2` may fail to install due to missing htslib binaries.

Recommended options:
- Use **WSL (Windows Subsystem for Linux)**
- Use **conda-forge** in a conda environment
- Run the tool inside a **Linux container or VM**

---

## Usage

## Directory Structure and Data Organization

The project expects input data to be organized in a predefined directory structure.
All required reference datasets and sample VCF files must be placed in the corresponding folders
before running the classification pipeline.

### Recommended Directory Layout

```
project_root/
├── data/
│   ├── gnomad4_1pct.vcf.gz
│   └── hg38_cosmic91.txt.gz
│
├── samples/
│   ├── CO8-MA-27_mutect2_vlod.vcf.gz
│   ├── CO8-MA-27_mutect2_vlod.vcf.gz.tbi
│   ├── CO8-MA-35_mutect2_vlod.vcf.gz
│   ├── CO8-MA-35_mutect2_vlod.vcf.gz.tbi
│   ├── CO8-PA-26_mutect2_vlod.vcf.gz
│   ├── CO8-PA-26_mutect2_vlod.vcf.gz.tbi
│   └── ...
│
├── results/
│   └── sample_name.tsv
│
├── tests/
│   └── test_classification.py
│
├── somatic_variant_classifier.py
└── README.md
```

### Data Directory (`data/`)

This directory contains reference datasets required for variant annotation:

- **gnomad4_1pct.vcf.gz**  
  Population allele frequency reference used for germline filtering and scoring.

- **hg38_cosmic91.txt.gz**  
  COSMIC v91 mutation database (hg38), used to estimate somatic recurrence across tissues.

These files are expected to be indexed and accessible locally. Their paths are resolved relative
to the project root.

### Sample Directory (`samples/`)

The `samples/` directory contains input VCF files to be analyzed.
Each sample must be provided as a compressed VCF file (`.vcf.gz`) together with its index file (`.tbi`).

The pipeline processes one sample at a time and generates a corresponding TSV output file in the
`results/` directory.

Example:
- Input: `samples/CO8-MA-27_mutect2_vlod.vcf.gz`
- Output: `results/CO8-MA-27_mutect2_vlod.tsv`

> Note: Reference datasets are not distributed with this repository and must be provided by the user
> due to licensing and size constraints.

### Basic Example

```bash
python somatic_variant_classifier.py   --sample-vcf samples/CO8-MA-27_mutect2_vlod.vcf.gz   --gnomad-vcf data/gnomad4_1pct.vcf.gz   --cosmic-tsv data/hg38_cosmic91.txt.gz   --output results/CO8-MA-27.tsv
```

This command:
1. Loads population data from gnomAD
2. Loads cancer recurrence data from COSMIC
3. Processes a Mutect2 tumor VCF
4. Writes a classified TSV output

### Optional Parameters

- `--germline-threshold` (default: 0.1)
- `--somatic-threshold` (default: 1.0)
- `--an-confidence-threshold` (default: 100000)

Use `--help` to see all available options.

---

## Example Outputs

### Tumor Sample

A variant example with low population frequency and high COSMIC recurrence are predominantly classified as **LIKELY_SOMATIC**.

chrom | pos | ref | alt | filter | sample_dp | sample_af | gnomad_af | gnomad_an | cosmic_count | cosmic_tissues | germline_score | somatic_score | classification | 
|---------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------| 
19 | 49033133 | C | A | haplotype;map_qual	| 1 | 0.684 |  |  | 21 | central_nervous_system,large_intestine,soft_tissue,stomach,thyroid,upper_aerodigestive_tract | 0 | 1.3424 | LIKELY_SOMATIC


### Control Sample (NC / NL)

Control samples are expected to yield mostly **LIKELY_GERMLINE** or **UNKNOWN** variants and serve as a validation of the classification logic.

chrom | pos | ref | alt | filter | sample_dp | sample_af | gnomad_af | gnomad_an | cosmic_count | cosmic_tissues | germline_score | somatic_score | classification | 
|---------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------| 
1 | 6124032 | A | G |	| 21 | 0.949 | 0.5968179703 | 1611694 | 1 | thyroid | 0.5968 | 0.301 | LIKELY_GERMLINE


---

## Testing

Core scoring and classification logic is unit-tested using `pytest`.

Run tests with the bash command of **(in 'germline_somatic_classifier' folder)**:

```bash
PYTHONPATH=$(pwd) pytest
```

All tests validate score calculations, threshold behavior, and edge cases such as missing annotations.


Example output:
```
============================= test session starts ==============================

collected 11 items

tests/test_classification.py ...........                                 [100%]

============================== 11 passed in 0.15s ==============================
```

---

## Design Decisions & Assumptions

- Reference datasets are loaded into in-memory lookup tables for simplicity
- Chromosome names are normalized (e.g. `chr1` → `1`)
- Multi-allelic variants are processed allele-by-allele
- For multi-allelic gnomAD records, only the first ALT allele frequency is used
- COSMIC indel representations are partially normalized; complex indels requiring left-alignment or anchor-base normalization may not match perfectly
- The tool prioritizes clarity and correctness over performance

---

## Intended Use

The output TSV is intended for **downstream filtering, review, or visualization** in research or clinical genomics pipelines.

---

## Status

✔ Functional requirements implemented  
✔ Technical requirements satisfied  
✔ Unit-tested core logic  
✔ Example usage and outputs provided  


# Somatic vs Germline Variant Classifier

## 🚀 Version 2 -- Performance & Architecture Update

Version 2 introduces a significant architectural improvement to the
pipeline.\
The scoring model and output format remain unchanged, but reference data
access has been redesigned for performance, stability, and scalability.

------------------------------------------------------------------------

## 🔧 What Changed in v2?

In Version 1:

-   gnomAD and COSMIC files were opened inside each worker process
-   Reference indexes were repeatedly initialized
-   High I/O overhead significantly increased runtime

In Version 2:

-   gnomAD and COSMIC files are initialized **once**
-   Indexed region-based queries are reused
-   Redundant file initialization is eliminated
-   Memory usage is reduced
-   Runtime is dramatically improved

Example benchmark (\~2000 variants):

  Version   Runtime
  --------- ------------------
  v1        \~1.5--2 minutes
  v2        \< 1 second

------------------------------------------------------------------------

## 📦 Reference Data Requirements

This pipeline requires **indexed reference files** to enable fast
region-based queries.

Required files inside the `data/` directory:

    data/
    ├── gnomad4_1pct.vcf.gz
    ├── gnomad4_1pct.vcf.gz.tbi
    ├── hg38_cosmic91.txt.gz
    └── hg38_cosmic91.txt.gz.tbi

If `.tbi` index files are missing, they must be created before running
the classifier.

------------------------------------------------------------------------

## 🛠 Preparing Reference Files

### 1️⃣ gnomAD (VCF)

If you already have:

    gnomad4_1pct.vcf.gz

Create the index:

``` bash
tabix -p vcf data/gnomad4_1pct.vcf.gz
```

This will generate:

    gnomad4_1pct.vcf.gz.tbi

------------------------------------------------------------------------

### 2️⃣ COSMIC

If you have:

    hg38_cosmic91.txt

Convert and index:

``` bash
bgzip data/hg38_cosmic91.txt
tabix -s 1 -b 2 -e 2 data/hg38_cosmic91.txt.gz
```

If you already have:

    hg38_cosmic91.txt.gz

But no index, run:

``` bash
tabix -s 1 -b 2 -e 2 data/hg38_cosmic91.txt.gz
```

This will generate:

    hg38_cosmic91.txt.gz.tbi

------------------------------------------------------------------------

## ▶ Running the Classifier

``` bash
python updated_somatic_variant_classifier.py   --sample-vcf samples/CO8-MA-27_mutect2_vlod.vcf.gz   --gnomad-vcf data/gnomad4_1pct.vcf.gz   --cosmic-tsv data/hg38_cosmic91.txt.gz   --output results/output.tsv
```

⚠️ Note: You do NOT need to explicitly pass `.tbi` files.\
They are automatically detected when present in the same directory.

------------------------------------------------------------------------

## 🧠 Design Principles

-   Memory-safe streaming access
-   Indexed region-based querying (Tabix / cyvcf2)
-   No full VCF loading into RAM
-   Production-ready architecture
-   Identical output format to Version 1

------------------------------------------------------------------------

## 📌 Output Format

The TSV output contains:

-   chrom
-   pos
-   ref
-   alt
-   filter
-   sample_dp
-   sample_af
-   gnomad_af
-   gnomad_an
-   cosmic_count
-   cosmic_tissues
-   germline_score
-   somatic_score
-   classification

Output schema is unchanged from Version 1.
