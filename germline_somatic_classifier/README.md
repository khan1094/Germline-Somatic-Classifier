# Somatic vs Germline Variant Classifier

A lightweight, research-oriented command-line tool for classifying
genomic variants as **LIKELY_SOMATIC**, **LIKELY_GERMLINE**,
**CONFLICTING**, or **UNKNOWN** by integrating evidence from
**Mutect2**, **gnomAD**, and **COSMIC**.

This project was implemented as part of a technical task and emphasizes
**correctness, interpretability, testability, and architectural
robustness**.

------------------------------------------------------------------------

# Version 2 -- Performance & Architecture Update

Version 2 introduces a major architectural improvement while keeping the
**scoring logic and output format fully identical to Version 1**.

## What Changed in v2?

### Version 1

-   gnomAD and COSMIC files were opened inside each worker process\
-   Reference indexes were repeatedly initialized\
-   High I/O overhead significantly increased runtime\
-   \~1.5--2 minutes runtime for \~2000 variants

### Version 2

-   gnomAD and COSMIC are initialized **once**
-   Indexed region-based queries are reused
-   Redundant file initialization eliminated
-   Significantly reduced memory footprint
-   Runtime improved to **\< 1 second** for \~2000 variants
-   Fully identical TSV output schema

This redesign makes the tool scalable and safe to run on older hardware
without OOM issues.

------------------------------------------------------------------------

# Overview

In cancer genomics, distinguishing somatic mutations from germline
variants is critical for downstream interpretation.

This tool applies simple, explainable scoring rules based on:

-   **Population frequency (gnomAD)**
-   **Cancer recurrence evidence (COSMIC)**

Each variant in a Mutect2 VCF is annotated with external evidence,
assigned germline and somatic scores, and classified using a
configurable decision matrix.

------------------------------------------------------------------------

# Features

-   Deterministic, interpretable scoring (no ML)
-   Allele-level processing of multi-allelic variants
-   Indexed region-based reference lookups
-   Memory-safe streaming access (no full VCF loading)
-   Identical output format between versions
-   Unit-tested classification logic
-   Suitable for research and evaluation pipelines

------------------------------------------------------------------------

# Input & Output

## Input

-   **Sample VCF** -- Mutect2 output (.vcf.gz + .tbi)
-   **gnomAD VCF** -- Population allele frequencies (.vcf.gz + .tbi)
-   **COSMIC TSV (bgzipped + indexed)** -- Cancer recurrence database

## Output

A tab-separated (TSV) file with one row per variant allele.

Columns:

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

------------------------------------------------------------------------

# Classification Logic

## Germline Evidence Score

germline_score = AF × min(AN / AN_CONFIDENCE_THRESHOLD, 1.0)

-   High population frequency suggests germline origin
-   Allele number (AN) weights confidence

## Somatic Evidence Score

somatic_score = log10(cosmic_count + 1)

-   Log scaling handles wide COSMIC recurrence ranges

## Decision Matrix

  Germline   Somatic   Classification
  ---------- --------- -----------------
  High       Low       LIKELY_GERMLINE
  Low        High      LIKELY_SOMATIC
  High       High      CONFLICTING
  Low        Low       UNKNOWN

Thresholds are configurable via CLI.

------------------------------------------------------------------------

# Installation

## Requirements

-   Python ≥ 3.9
-   cyvcf2 == 0.31.4
-   pysam
-   tabix / bgzip (htslib)

Install:

``` bash
pip install -r requirements.txt
```

------------------------------------------------------------------------

# Reference File Preparation

The pipeline requires indexed reference files inside the `data/`
directory:

    data/
    ├── gnomad4_1pct.vcf.gz
    ├── gnomad4_1pct.vcf.gz.tbi
    ├── hg38_cosmic91.txt.gz
    └── hg38_cosmic91.txt.gz.tbi

If index files are missing:

## gnomAD

``` bash
tabix -p vcf data/gnomad4_1pct.vcf.gz
```

## COSMIC

If starting from .txt:

``` bash
bgzip data/hg38_cosmic91.txt
tabix -s 1 -b 2 -e 2 data/hg38_cosmic91.txt.gz
```

If already gzipped but not indexed:

``` bash
tabix -s 1 -b 2 -e 2 data/hg38_cosmic91.txt.gz
```

You do NOT need to pass `.tbi` files explicitly --- they are
automatically detected.

------------------------------------------------------------------------

# Usage

``` bash
python updated_somatic_variant_classifier.py   --sample-vcf samples/sample.vcf.gz   --gnomad-vcf data/gnomad4_1pct.vcf.gz   --cosmic-tsv data/hg38_cosmic91.txt.gz   --output results/output.tsv
```

Optional parameters:

-   --germline-threshold (default: 0.1)
-   --somatic-threshold (default: 1.0)
-   --an-confidence-threshold (default: 100000)
-   --threads (default: CPU count)

------------------------------------------------------------------------

# Testing

Run unit tests:

``` bash
PYTHONPATH=$(pwd) pytest
```

All tests validate scoring logic, threshold behavior, and edge cases.

------------------------------------------------------------------------

# Design Philosophy

-   Transparent, rule-based classification
-   Reproducible deterministic results
-   Efficient indexed data access
-   Architecture designed for scalability
-   Research-first implementation with production-aware design

------------------------------------------------------------------------

# Status

✔ Functional requirements implemented\
✔ Architecture optimized (v2)\
✔ Unit-tested classification logic\
✔ Stable output schema\
✔ Ready for integration into research pipelines
