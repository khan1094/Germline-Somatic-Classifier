#!/usr/bin/env python3

import argparse
import gzip
import math
import os
import re
import sys
import csv
from collections import defaultdict
from cyvcf2 import VCF

####################################
# Utility functions
####################################

def normalize_chrom(chrom: str) -> str:
    return chrom.replace("chr", "")

def compute_germline_score(af, an, an_threshold):
    if af is None or an is None:
        return 0.0
    confidence = min(an / an_threshold, 1.0)
    return af * confidence

def compute_somatic_score(cosmic_count):
    if cosmic_count is None:
        return 0.0
    return math.log10(cosmic_count + 1)

####################################
# Load gnomAD
####################################

def load_gnomad(gnomad_vcf_path):
    if not os.path.exists(gnomad_vcf_path):
        sys.exit(f"ERROR: gnomAD VCF not found: {gnomad_vcf_path}")

    gnomad = {}
    vcf = VCF(gnomad_vcf_path)

    for var in vcf:
        chrom = normalize_chrom(var.CHROM)
        pos = var.POS
        ref = var.REF

        af = var.INFO.get("AF")
        an = var.INFO.get("AN")

        if isinstance(af, (list, tuple)):
            af = af[0]

        for alt in var.ALT:
            key = (chrom, pos, ref, alt)
            gnomad[key] = {"af": af, "an": an}

    return gnomad

####################################
# Load COSMIC
####################################

def load_cosmic(cosmic_path):
    if not os.path.exists(cosmic_path):
        sys.exit(f"ERROR: COSMIC file not found: {cosmic_path}")

    cosmic = defaultdict(lambda: {"count": 0, "tissues": set()})

    with gzip.open(cosmic_path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue

            chrom = normalize_chrom(parts[0])
            start = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            info = parts[5]

            matches = re.findall(r"(\d+)\(([^)]+)\)", info)
            total = 0
            tissues = set()

            for c, t in matches:
                total += int(c)
                tissues.add(t)

            key = (chrom, start, ref, alt)
            cosmic[key]["count"] += total
            cosmic[key]["tissues"].update(tissues)

    return cosmic


####################################
# Classify Variant
####################################

def classify_variant(
    germline_score,
    somatic_score,
    germline_threshold,
    somatic_threshold
):
    """
    Classification matrix based on task specification.
    """
    g_high = germline_score >= germline_threshold
    s_high = somatic_score >= somatic_threshold

    if g_high and not s_high:
        return "LIKELY_GERMLINE"
    if not g_high and s_high:
        return "LIKELY_SOMATIC"
    if g_high and s_high:
        return "CONFLICTING"
    return "UNKNOWN"

def main():
    parser = argparse.ArgumentParser(
        description="Somatic Variant Classifier using gnomAD and COSMIC"
    )

    parser.add_argument("--sample-vcf", required=True, help="Mutect2 VCF file")
    parser.add_argument("--gnomad-vcf", required=True, help="gnomAD reference VCF")
    parser.add_argument("--cosmic-tsv", required=True, help="COSMIC TSV (gzipped)")
    parser.add_argument("--output", required=True, help="Output TSV file")

    parser.add_argument("--germline-threshold", type=float, default=0.1)
    parser.add_argument("--somatic-threshold", type=float, default=1.0)
    parser.add_argument("--an-confidence-threshold", type=int, default=100000)

    args = parser.parse_args()

    if not os.path.exists(args.sample_vcf):
        sys.exit(f"ERROR: Sample VCF not found: {args.sample_vcf}")

    print("Loading gnomAD...")
    gnomad = load_gnomad(args.gnomad_vcf)

    print("Loading COSMIC...")
    cosmic = load_cosmic(args.cosmic_tsv)

    print("Processing sample VCF...")
    sample_vcf = VCF(args.sample_vcf)

    output_fields = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "filter",
        "sample_dp",
        "sample_af",
        "gnomad_af",
        "gnomad_an",
        "cosmic_count",
        "cosmic_tissues",
        "germline_score",
        "somatic_score",
        "classification",
    ]

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    with open(args.output, "w", newline="") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()

        for var in sample_vcf:
            dp = var.format("DP")[0][0] if var.format("DP") is not None else None
            af = var.format("AF")[0][0] if var.format("AF") is not None else None

            for alt in var.ALT:
                chrom = normalize_chrom(var.CHROM)
                pos = var.POS
                ref = var.REF

                key = (chrom, pos, ref, alt)

                g = gnomad.get(key)
                c = cosmic.get(key)

                g_af = g["af"] if g else None
                g_an = g["an"] if g else None
                c_count = c["count"] if c else None
                c_tissues = ",".join(sorted(c["tissues"])) if c else ""

                germline_score = compute_germline_score(
                    g_af, g_an, args.an_confidence_threshold
                )
                somatic_score = compute_somatic_score(c_count)

                classification = classify_variant(
                    germline_score,
                    somatic_score,
                    args.germline_threshold,
                    args.somatic_threshold)

                writer.writerow({
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "filter": var.FILTER,
                    "sample_dp": dp if dp is not None else "",
                    "sample_af": af if af is not None else "",
                    "gnomad_af": g_af if g_af is not None else "",
                    "gnomad_an": g_an if g_an is not None else "",
                    "cosmic_count": c_count if c_count is not None else "",
                    "cosmic_tissues": c_tissues,
                    "germline_score": round(germline_score, 4),
                    "somatic_score": round(somatic_score, 4),
                    "classification": classification,
                })

    print(f"Done. Output written to: {args.output}")

if __name__ == "__main__":
    main()
