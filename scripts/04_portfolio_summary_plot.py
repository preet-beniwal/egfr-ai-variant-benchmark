import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import os

DATA = "data/merged_ai_physics_evolution.tsv"
os.makedirs("figures", exist_ok=True)

df = pd.read_csv(DATA, sep="\t").dropna(subset=["conservation_score", "ddg_kcal_per_mol", "alphamissense_score"])
np.random.seed(42)
df["cons_jitter"] = df["conservation_score"] + np.random.uniform(-0.1, 0.1, len(df))

highlights = {
    "L858R": {"color": "red", "marker": "*"},
    "T790M": {"color": "green", "marker": "D"},
    "G719S": {"color": "purple", "marker": "^"}
}

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

# Panel 1
ax = axes[0]
rho1, p1 = stats.spearmanr(df["conservation_score"], df["ddg_kcal_per_mol"])
ax.scatter(df["cons_jitter"], df["ddg_kcal_per_mol"], alpha=0.5, color="grey", edgecolor="black")
for var, style in highlights.items():
    row = df[df["protein_change"] == var]
    if not row.empty:
        ax.scatter(row["cons_jitter"], row["ddg_kcal_per_mol"], color=style["color"], marker=style["marker"], s=150, label=var, zorder=5)
ax.set_xticks([1, 2, 3])
ax.set_xlabel("Conservation Score (Evolution)")
ax.set_ylabel("FoldX ΔΔG (Physics)")
ax.set_title(f"Evolution vs Physics\nSpearman ρ = {rho1:.2f}, p = {p1:.3f}")
ax.axhline(0, color="black", linestyle="--", alpha=0.3)
ax.grid(True, alpha=0.3)
ax.legend()

# Panel 2
ax = axes[1]
rho2, p2 = stats.spearmanr(df["ddg_kcal_per_mol"], df["alphamissense_score"])
ax.scatter(df["ddg_kcal_per_mol"], df["alphamissense_score"], alpha=0.5, color="grey", edgecolor="black")
for var, style in highlights.items():
    row = df[df["protein_change"] == var]
    if not row.empty:
        ax.scatter(row["ddg_kcal_per_mol"], row["alphamissense_score"], color=style["color"], marker=style["marker"], s=150, label=var, zorder=5)
ax.set_xlabel("FoldX ΔΔG (Physics)")
ax.set_ylabel("AlphaMissense Score (AI)")
ax.set_title(f"Physics vs AI\nSpearman ρ = {rho2:.2f}, p = {p2:.3f}")
ax.axhline(0.56, color="red", linestyle=":", alpha=0.5, label="Pathogenic Threshold")
ax.grid(True, alpha=0.3)
ax.legend()

# Panel 3
ax = axes[2]
rho3, p3 = stats.spearmanr(df["conservation_score"], df["alphamissense_score"])
ax.scatter(df["cons_jitter"], df["alphamissense_score"], alpha=0.5, color="grey", edgecolor="black")
for var, style in highlights.items():
    row = df[df["protein_change"] == var]
    if not row.empty:
        ax.scatter(row["cons_jitter"], row["alphamissense_score"], color=style["color"], marker=style["marker"], s=150, label=var, zorder=5)
ax.set_xticks([1, 2, 3])
ax.set_xlabel("Conservation Score (Evolution)")
ax.set_ylabel("AlphaMissense Score (AI)")
ax.set_title(f"Evolution vs AI\nSpearman ρ = {rho3:.2f}, p = {p3:.3f}")
ax.axhline(0.56, color="red", linestyle=":", alpha=0.5, label="Pathogenic Threshold")
ax.grid(True, alpha=0.3)
ax.legend()

fig.suptitle("The Limits of Single-Method Variant Interpretation\n(EGFR LUAD Variants: Evolution vs Physics vs AI)", fontsize=16, fontweight="bold", y=1.05)
plt.tight_layout()
plt.savefig("figures/portfolio_summary_figure.png", dpi=300, bbox_inches='tight')
print("Successfully generated figures/portfolio_summary_figure.png")
