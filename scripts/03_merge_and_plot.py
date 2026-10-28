#!/usr/bin/env python3
"""
Merge AlphaMissense (AI) with Project 4 (FoldX + Conservation) and plot.
"""
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os

# Paths to the data
PROJ4_MASTER = os.path.expanduser("~/projects/egfr-variant-structural-analysis/data/variants/luad_egfr_variants_master.tsv")
ALPHAMISSENSE = "data/alphamissense_egfr.tsv"
FIG_DIR = "figures"
os.makedirs(FIG_DIR, exist_ok=True)

def main():
    print("Loading Project 4 data (FoldX + Conservation)...")
    df_proj4 = pd.read_csv(PROJ4_MASTER, sep="\t")
    
    # We only want the 40 substitutions
    subs = df_proj4[df_proj4["variant_type"] == "substitution"].drop_duplicates("protein_change").copy()
    
    print("Loading AlphaMissense data...")
    # Read TSV, skipping any weird headers if necessary
    df_am = pd.read_csv(ALPHAMISSENSE, sep="\t")
    
    # Clean up column names just in case
    df_am.columns = [c.strip() for c in df_am.columns]
    
    # Filter for canonical EGFR transcript if needed
    if "Ensemble id" in df_am.columns:
        df_am = df_am[df_am["Ensemble id"].str.startswith("ENST00000275493", na=False)]
    
    # Create a 'protein_change' string (e.g., L858R) to match Project 4
    # AlphaMissense columns: a.a.1 (WT), position (num), a.a.2 (Mut)
    df_am["protein_change"] = df_am["a.a.1"].astype(str) + df_am["position"].astype(str) + df_am["a.a.2"].astype(str)
    
    # Rename the AI score column
    df_am = df_am.rename(columns={"pathogenicity score": "alphamissense_score"})
    
    # Merge
    merged = subs.merge(df_am[["protein_change", "alphamissense_score"]], on="protein_change", how="inner")
    print(f"\nSuccessfully merged {len(merged)} variants with AI scores.")
    
    # Drop NAs
    merged = merged.dropna(subset=["alphamissense_score", "ddg_kcal_per_mol", "conservation_score"])
    
    # Calculate correlations
    rho_ai_phys, p_ai_phys = stats.spearmanr(merged["alphamissense_score"], merged["ddg_kcal_per_mol"])
    rho_ai_evo, p_ai_evo = stats.spearmanr(merged["alphamissense_score"], merged["conservation_score"])
    
    print("\n=== CORRELATION RESULTS ===")
    print(f"AlphaMissense (AI) vs FoldX (Physics):  rho = {rho_ai_phys:.3f}, p = {p_ai_phys:.4f}")
    print(f"AlphaMissense (AI) vs Conservation (Evo): rho = {rho_ai_evo:.3f}, p = {p_ai_evo:.4f}")
    
    # --- Plot 1: AI vs Physics ---
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(merged["ddg_kcal_per_mol"], merged["alphamissense_score"], alpha=0.7, color="#1f77b4", edgecolor="black")
    
    # Highlight key variants
    for var, color in {"L858R": "red", "T790M": "green", "G719S": "purple"}.items():
        row = merged[merged["protein_change"] == var]
        if not row.empty:
            ax.scatter(row["ddg_kcal_per_mol"], row["alphamissense_score"], color=color, s=100, label=var, zorder=5)
            
    ax.set_xlabel("FoldX ΔΔG (kcal/mol) - Physics", fontsize=11)
    ax.set_ylabel("AlphaMissense Score - AI", fontsize=11)
    ax.set_title(f"AI vs Physics\nSpearman ρ = {rho_ai_phys:.3f}, p = {p_ai_phys:.4f}", fontsize=12)
    ax.axhline(0.56, color="grey", linestyle="--", alpha=0.5, label="Pathogenic threshold")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/ai_vs_physics.png", dpi=300)
    
    # --- Plot 2: AI vs Evolution ---
    fig, ax = plt.subplots(figsize=(6, 5))
    # Add jitter to conservation score since it's discrete (1, 2, 3)
    import numpy as np
    jitter = np.random.RandomState(42).uniform(-0.1, 0.1, len(merged))
    ax.scatter(merged["conservation_score"] + jitter, merged["alphamissense_score"], alpha=0.7, color="#d62728", edgecolor="black")
    
    for var, color in {"L858R": "blue", "T790M": "green", "G719S": "purple"}.items():
        row = merged[merged["protein_change"] == var]
        if not row.empty:
            ax.scatter(row["conservation_score"], row["alphamissense_score"], color=color, s=100, label=var, zorder=5)
            
    ax.set_xticks([1, 2, 3])
    ax.set_xlabel("Conservation Score (0-3) - Evolution", fontsize=11)
    ax.set_ylabel("AlphaMissense Score - AI", fontsize=11)
    ax.set_title(f"AI vs Evolution\nSpearman ρ = {rho_ai_evo:.3f}, p = {p_ai_evo:.4f}", fontsize=12)
    ax.axhline(0.56, color="grey", linestyle="--", alpha=0.5, label="Pathogenic threshold")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}/ai_vs_evolution.png", dpi=300)

    print(f"\nWrote {FIG_DIR}/ai_vs_physics.png and {FIG_DIR}/ai_vs_evolution.png")
    
    # Save the merged table for the paper
    merged.to_csv("data/merged_ai_physics_evolution.tsv", sep="\t", index=False)
    print("Wrote data/merged_ai_physics_evolution.tsv")

if __name__ == "__main__":
    main()
