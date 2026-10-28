#!/usr/bin/env python3
"""
Fetch AlphaMissense scores using the GRCh38 VEP API.
"""
import requests
import pandas as pd
import time
import os

# We will use our master table from Project 4 as the input
MASTER = os.path.expanduser("~/projects/egfr-variant-structural-analysis/data/variants/luad_egfr_variants_master.tsv")
OUTPUT = "data/alphamissense_scores.tsv"
URL = "https://rest.ensembl.org/vep/human/hgvs?AlphaMissense=1"

def main():
    print("Loading variants...")
    df = pd.read_csv(MASTER, sep="\t")
    
    # Get unique missense substitutions
    variants = df[df["variant_type"] == "substitution"]["protein_change"].dropna().unique().tolist()
    print(f"Fetching AlphaMissense for {len(variants)} variants...")
    
    results = []
    
    for i in range(0, len(variants), 10):
        batch = variants[i:i+10]
        # Use ENSP00000275493 (no version number)
        payload = {"hgvs_notations": [f"ENSP00000275493:p.{v}" for v in batch]}
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        
        try:
            r = requests.post(URL, headers=headers, json=payload, timeout=30)
            r.raise_for_status()
            data = r.json()
            
            for item in data:
                pc = item.get("input", "").split(":p.")[-1]
                am_score = None
                am_class = None
                
                for tc in item.get("transcript_consequences", []):
                    if tc.get("gene_symbol") == "EGFR":
                        if "alphamissense" in tc and tc["alphamissense"]:
                            am_score = tc["alphamissense"][0].get("am_pathogenicity")
                            am_class = tc["alphamissense"][0].get("am_class")
                
                results.append({
                    "protein_change": pc,
                    "alphamissense_score": am_score,
                    "alphamissense_class": am_class
                })
            
            print(f"  Processed {min(i+10, len(variants))}/{len(variants)}")
            time.sleep(1)
            
        except Exception as e:
            print(f"  Error on batch {i}: {e}")
            time.sleep(2)

    out_df = pd.DataFrame(results)
    out_df.to_csv(OUTPUT, sep="\t", index=False)
    print(f"\nWrote {OUTPUT}")
    print(out_df.head(10).to_string(index=False))

if __name__ == "__main__":
    main()

