#!/usr/bin/env python3

import os
import pandas as pd

# --- CNV Thresholds ---
# These numeric thresholds define an initial, permissive CNV call.
# The final biological filtering will be done in the R script (07_oncoprint_visualization.R).
AMP_THRESHOLD = 0.3
DEL_THRESHOLD = -0.3

# --- Path Variables ---
# Directory containing GATK CNV segment files (.cr.igv.seg)
CNV_DIR = "../cnv_gatk" 
# Target gene BED file
GENE_BED = "../resources/genes_interes_2.bed"
# Final output file for the R script
OUTFILE = "../data/processed/CNV_genes_overlap_3.csv" 

def load_bed(bedfile):
    """Loads target genes from the BED file."""
    genes = pd.read_csv(bedfile, sep='\t', header=None, names=["chr", "start", "end", "gene"])
    return genes

def load_segments(segfile):
    """Loads and cleans the GATK CNV segment file (.cr.igv.seg)."""
    df = pd.read_csv(segfile, sep='\t')
    df = df.rename(columns=lambda x: x.strip()) # Clean column names
    df["Segment_Mean"] = pd.to_numeric(df["Segment_Mean"], errors="coerce")
    df = df[df["Segment_Mean"].notnull()]
    return df

def overlaps(seg, gene):
    """Checks if a CNV segment overlaps with a gene region."""
    # Check if chromosomes match AND they are NOT completely disjoint
    return seg["chr"] == gene["chr"] and not (seg["end"] < gene["start"] or seg["start"] > gene["end"])

def main():
    try:
        gene_df = load_bed(GENE_BED)
    except FileNotFoundError:
        print(f"❌ Error: Gene BED file not found at {GENE_BED}.")
        return

    # Ensure chromosome format consistency (e.g., remove "chr" prefix if present)
    gene_df["chr"] = gene_df["chr"].astype(str).str.replace("^chr", "", regex=True)

    results = []
    
    # Ensure the CNV directory exists before listing files
    if not os.path.exists(CNV_DIR):
        print(f"❌ Error: CNV input directory {CNV_DIR} not found. Ensure 06_CNV.sh ran successfully.")
        return

    for fname in os.listdir(CNV_DIR):
        if fname.endswith(".cr.igv.seg"): # Expected GATK CNV segment file format
            sample = fname.replace(".cr.igv.seg", "")
            print(f"🔍 Processing CNV segments for: {sample}")
            
            seg_df = load_segments(os.path.join(CNV_DIR, fname))
            seg_df = seg_df.rename(columns={"Chromosome": "chr", "Start": "start", "End": "end"})
            seg_df["chr"] = seg_df["chr"].astype(str).str.replace("^chr", "", regex=True)

            for _, seg in seg_df.iterrows():
                cnv_type = None
                
                # Check against initial numeric thresholds
                if seg["Segment_Mean"] >= AMP_THRESHOLD:
                    cnv_type = "Amplification"
                elif seg["Segment_Mean"] <= DEL_THRESHOLD:
                    cnv_type = "Deletion"
                    
                if cnv_type:
                    # Check overlap with genes of interest
                    for _, gene in gene_df.iterrows():
                        if overlaps(seg, gene):
                            results.append({
                                "Sample": sample,
                                "Gene": gene["gene"],
                                "CNV_Type": cnv_type,
                                "Segment_Mean": seg["Segment_Mean"]
                            })

    if results:
        df = pd.DataFrame(results)
        # Ensure the output directory exists
        os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)
        df.to_csv(OUTFILE, index=False)
        print(f"✅ Summary file generated: {OUTFILE}")
    else:
        print("⚠️ No relevant CNVs found for the genes of interest using the specified thresholds.")

if __name__ == "__main__":
    main()
