#!/usr/bin/env python3

import os
import gzip
import pandas as pd

# === READ TARGET GENES FROM BED FILE ===
# This list is used to filter variants called by MuTect2.
# Assumes BED file is located in the sibling 'resources' folder.
BED_FILE = "../resources/genes_interes_2.bed" 
try:
    genes_df = pd.read_csv(BED_FILE, sep="\t", header=None, names=["chr", "start", "end", "gene"])
    GENES_INTERES = set(genes_df["gene"])
    print(f"✔️ {len(GENES_INTERES)} genes loaded from {BED_FILE}")
except FileNotFoundError:
    print(f"❌ Error: BED file not found at {BED_FILE}")
    GENES_INTERES = set()

# === CORE FUNCTIONS ===
def parse_ann_field(ann_str):
    """Parses the VCF's ANN field (from SnpEff annotation) to extract gene, effect, and impact."""
    fields = ann_str.split(",")
    genes = []
    for entry in fields:
        parts = entry.split("|")
        if len(parts) > 3:
            gene = parts[3]
            effect = parts[1]
            impact = parts[2]
            genes.append((gene, effect, impact))
    return genes

def process_vcf(file_path):
    """Processes an annotated VCF and extracts variants for genes of interest."""
    sample = os.path.basename(file_path).replace(".filtered.vcf.gz", "").replace(".ann.vcf", "")
    records = []

    # Handle both compressed and uncompressed VCFs
    if file_path.endswith(".gz"):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'rt'

    with opener(file_path, mode) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            # CHROM, POS, ID (skip), REF, ALT, QUAL (skip), FILTER (skip), INFO
            chrom, pos, ref, alt, info = cols[0], cols[1], cols[3], cols[4], cols[7]
            
            # Extract ANN field (assumes VCF has been annotated, e.g., with SnpEff)
            ann_data = [x for x in info.split(";") if x.startswith("ANN=")]
            if not ann_data:
                continue
            
            ann_str = ann_data[0][4:]
            gene_effects = parse_ann_field(ann_str)
            
            for gene, effect, impact in gene_effects:
                if gene in GENES_INTERES:
                    records.append({
                        "Sample": sample,
                        "Chrom": chrom,
                        "Pos": pos,
                        "Ref": ref,
                        "Alt": alt,
                        "Gene": gene,
                        "Effect": effect,
                        "Impact": impact
                    })

    return pd.DataFrame(records)

# === PROCESS ALL ANNOTATED VCF FILES ===
# Adjust this path if VCF files are not directly in the 'results' folder
INPUT_DIR = "../results"
output_file = "../data/processed/variant_summary_genes_2.xlsx"
all_records = []

# Ensure the results folder exists before listing files
if not os.path.exists(INPUT_DIR):
    print(f"❌ Error: Input directory {INPUT_DIR} not found. Ensure 05_mutect.sh ran successfully.")
    exit()

for fname in os.listdir(INPUT_DIR):
    # Assuming the final VCFs from GATK were annotated to produce .ann.vcf or similar
    if fname.endswith(".ann.vcf") or fname.endswith(".filtered.vcf.gz"): # Check both annotated VCFs and the raw filtered GATK output
        print(f"🔍 Processing {fname}...")
        df = process_vcf(os.path.join(INPUT_DIR, fname))
        if not df.empty:
            all_records.append(df)

# === SAVE RESULTS ===
if all_records:
    final_df = pd.concat(all_records, ignore_index=True)
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    final_df.to_excel(output_file, index=False)
    print(f"✅ Summary file generated: {output_file}")
else:
    print("⚠️ No variants found for the genes of interest.")
