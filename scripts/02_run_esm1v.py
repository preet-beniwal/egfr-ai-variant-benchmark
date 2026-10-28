#!/usr/bin/env python3
"""
Fetch the canonical UniProt sequence for EGFR and run ESM-1v.
"""
import torch
import esm
import pandas as pd
import requests
import os

UNIPROT_ID = "P00533"

def fetch_sequence():
    print(f"Fetching sequence for {UNIPROT_ID} from UniProt...")
    url = f"https://rest.uniprot.org/uniprotkb/{UNIPROT_ID}.fasta"
    r = requests.get(url)
    r.raise_for_status()
    lines = r.text.splitlines()
    seq = "".join([line.strip() for line in lines if not line.startswith(">")])
    print(f"  Sequence length: {len(seq)} amino acids")
    return seq

def main():
    seq = fetch_sequence()
    
    print("Loading ESM-1v model (this may take a minute)...")
    # Load the first of the 5 ESM-1v models
    model, alphabet = esm.pretrained.esm1v_t33_650M_UR90S_1()
    batch_converter = alphabet.get_batch_converter()
    model.eval()
    
    print("Loading variants from Project 4...")
    master = pd.read_csv(os.path.expanduser("~/projects/egfr-variant-structural-analysis/data/variants/luad_egfr_variants_master.tsv"), sep="\t")
    variants = master[master["variant_type"] == "substitution"]["protein_change"].dropna().unique().tolist()
    
    print(f"Scoring {len(variants)} variants using ESM-1v...")
    data = [("EGFR_wt", seq)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    
    results = []
    
    with torch.no_grad():
        results_dict = model(batch_tokens, repr_layers=[33], return_contacts=False)
        logits = results_dict["logits"].squeeze(0)
        
    for pc in variants:
        wt_aa = pc[0]
        pos = int(pc[1:-1])
        mut_aa = pc[-1]
        
        # UniProt is 1-based. The ESM token index corresponds directly to the position.
        if pos < 1 or pos > len(seq):
            print(f"  {pc} out of bounds.")
            continue
            
        seq_wt_aa = seq[pos - 1]
        if seq_wt_aa != wt_aa:
            print(f"  Warning: WT mismatch for {pc}. Sequence has {seq_wt_aa}, expected {wt_aa}.")
            continue
            
        wt_idx = alphabet.get_idx(wt_aa)
        mut_idx = alphabet.get_idx(mut_aa)
        
        # ESM-1v score is the log-odds ratio: log P(mut) - log P(wt)
        wt_logit = logits[pos, wt_idx].item()
        mut_logit = logits[pos, mut_idx].item()
        score = mut_logit - wt_logit
        
        results.append({
            "protein_change": pc,
            "esm1v_score": score
        })
        
    out_df = pd.DataFrame(results)
    out_df.to_csv("data/esm1v_scores.tsv", sep="\t", index=False)
    print(f"\nWrote data/esm1v_scores.tsv ({len(out_df)} variants scored).")
    print("\nTop 10 most pathogenic by ESM-1v:")
    print(out_df.sort_values("esm1v_score", ascending=False).head(10).to_string(index=False))

if __name__ == "__main__":
    main()
