import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load both files
predictions_path = "" # path to predicted peptides
strong_binders_path = "" # path to strong binders

predictions_df = pd.read_csv(predictions_path)
strong_df = pd.read_csv(strong_binders_path)

# Create output directory for plots
output_dir = "" # Set output directory
os.makedirs(output_dir, exist_ok=True)

# Plot 1: Distribution of binding affinities
plt.figure(figsize=(10, 6))
sns.histplot(predictions_df["affinity"], bins=50, kde=True)
plt.axvline(500, color="red", linestyle="--", label="Strong Binder Threshold (500 nM)")
plt.title("Distribution of Binding Affinities (All Predictions)")
plt.xlabel("Predicted Affinity (nM)")
plt.ylabel("Peptide Count")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "affinity_distribution_all.png"))
plt.close()

# Plot 2: Sample count distribution for strong binders
plt.figure(figsize=(10, 6))
sns.histplot(strong_df["sample_count"], bins=30, kde=False)
plt.title("Distribution of Sample Counts (Strong Binders)")
plt.xlabel("Number of Samples Sharing Peptide")
plt.ylabel("Peptide Count")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "sample_count_strong_binders.png"))
plt.close()

# Plot 3: Top 20 most common strong binder peptides
top_common = strong_df.sort_values("sample_count", ascending=False).head(20)
plt.figure(figsize=(12, 8))
sns.barplot(data=top_common, x="sample_count", y="peptide", palette="viridis")
plt.title("Top 20 Most Common Strong Binder Peptides")
plt.xlabel("Number of Samples")
plt.ylabel("Peptide")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "top20_strong_binders.png"))
plt.close()
