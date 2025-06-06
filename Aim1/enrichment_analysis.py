#Author: Yugendran Rajaendran
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
from collections import defaultdict, Counter
from matplotlib_venn import venn2

try:
    import tensorflow as tf
    gpus = tf.config.list_physical_devices("GPU")
    if gpus:
        print(f"GPU detected: {len(gpus)} device(s) available.")
    else:
        print("No GPU detected.")
except ImportError:
    print("TensorFlow not installed. GPU check skipped.")

# SETUP
BASE_DIR = ""  # UPDATE this to your local folder
OUTPUT_DIR = "new_neoantigen_analysis_outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# HELPERS
def extract_metadata(filename):
    name = os.path.basename(filename).replace(".xlsx", "")
    parts = name.split("_")
    tool = parts[2].lower()
    cancer_type = parts[3].upper()
    threshold = parts[4]
    hla_type = parts[-1] if parts[-1] in ['A', 'B'] else 'NA'
    return tool, cancer_type, threshold, hla_type

# DATA CONTAINERS
peptide_sample_map = defaultdict(lambda: defaultdict(set))
sample_peptide_burden = defaultdict(Counter)
tool_peptide_set = defaultdict(set)

# INGEST DATA
file_list = glob.glob(os.path.join(BASE_DIR, "*.xlsx"))
for filepath in file_list:
    try:
        df = pd.read_excel(filepath)
        tool, cancer_type, threshold, hla_type = extract_metadata(filepath)
        key = (tool, cancer_type, threshold, hla_type)

        for _, row in df.iterrows():
            peptide = str(row['peptide_8to14']).strip().upper()
            sample = str(row['sample_id'])

            peptide_sample_map[key][peptide].add(sample)
            sample_peptide_burden[key][sample] += 1
            tool_peptide_set[(tool, cancer_type, threshold)].add(peptide)
    except Exception as e:
        print(f"Error processing {filepath}: {e}")

# CANCER-TOOL-TABLE
df_cancer_stats = []
for (tool, cancer_type, threshold, hla_type), peptides in peptide_sample_map.items():
    total_peptides = len(peptides)
    shared = sum(1 for p in peptides.values() if len(p) > 1)
    total_samples = len(set(s for p in peptides.values() for s in p))
    df_cancer_stats.append({
        "Tool": tool, "CancerType": cancer_type,
        "Threshold": threshold, "HLAType": hla_type,
        "TotalPeptides": total_peptides,
        "SharedPeptides": shared,
        "TotalSamples": total_samples
    })
df_cancer_stats = pd.DataFrame(df_cancer_stats)
df_cancer_stats.to_csv(os.path.join(OUTPUT_DIR, "cancer_type_summary.csv"), index=False)

# PEPTIDE BURDEN TABLE
sample_burden = []
for (tool, cancer_type, threshold, hla_type), counts in sample_peptide_burden.items():
    for sample, burden in counts.items():
        sample_burden.append({
            "Tool": tool,
            "CancerType": cancer_type,
            "Threshold": threshold,
            "SampleID": sample,
            "PeptideBurden": burden
        })
df_burden = pd.DataFrame(sample_burden)
df_burden.to_csv(os.path.join(OUTPUT_DIR, "sample_burden_summary.csv"), index=False)

# ------------------------------------------
# PLOT 1: Peptide Burden by Tool (SCLC only)
# ------------------------------------------
df_plot1 = df_burden[df_burden["CancerType"] == "SCLC"]
plt.figure(figsize=(10, 6))
sns.boxplot(data=df_plot1, x="Tool", y="PeptideBurden")
plt.yscale("log")
plt.xlabel("Prediction Tool")
plt.ylabel("Number of Predicted Neoantigenic Peptides (Log Scale)")
plt.title("Neoantigen Burden per Sample in SCLC")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "burden_boxplot_SCLC.png"))
plt.close()


# ------------------------------------------
# PLOT 2: Threshold Analysis (DeepMHCi only)
# ------------------------------------------
df_plot2 = df_cancer_stats[(df_cancer_stats["Tool"] == "deepmhci") & (df_cancer_stats["CancerType"] == "SCLC")]
plt.figure(figsize=(10, 6))
sns.barplot(data=df_plot2, x="Threshold", y="TotalPeptides", order=["0.1", "0.2", "0.3"])
plt.xlabel("Binding Score Threshold")
plt.ylabel("Total Predicted Neoantigenic Peptides")
plt.title("Effect of Threshold on Peptide Predictions (DeepMHCi, SCLC)")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "threshold_deepmhci_SCLC.png"))
plt.close()

# ------------------------------------------
# PLOT 3: Cancer Comparison (Threshold = 0.3)
# ------------------------------------------
df_plot3 = df_cancer_stats[df_cancer_stats["Threshold"] == "0.3"]
plt.figure(figsize=(10, 6))
sns.barplot(data=df_plot3, x="CancerType", y="TotalPeptides", hue="Tool")
plt.xlabel("Cancer Type")
plt.ylabel("Total Predicted Neoantigenic Peptides")
plt.title("Tool Comparison at Threshold 0.3")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "tool_comparison_threshold_0.3.png"))
plt.close()

# ------------------------------------------
# PLOT 4: Venn Overlap for SCLC @ 0.1/500
# ------------------------------------------
try:
    peptides1 = tool_peptide_set[("deepmhci", "SCLC", "0.1")]
    peptides2 = tool_peptide_set[("mhcflurry", "SCLC", "500")]
    plt.figure(figsize=(6, 6))
    venn2([peptides1, peptides2], set_labels=["DeepMHCi 0.1", "MHCflurry 500"])
    plt.title("Peptide Overlap in SCLC: DeepMHCi vs MHCflurry")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "venn_sclc_overlap.png"))
    plt.close()
except KeyError:
    print("One or both peptide sets missing for SCLC Venn comparison.")

# ------------------------------------------
# PLOT 5: Number of Samples per Cancer Type
# ------------------------------------------
df_sample_counts = df_burden.groupby("CancerType")["SampleID"].nunique().reset_index()
df_sample_counts = df_sample_counts.rename(columns={"SampleID": "NumSamples"})

plt.figure(figsize=(8, 6))
sns.barplot(data=df_sample_counts, x="CancerType", y="NumSamples", order=["SCLC", "NSCLC", "PAN"])
plt.title("Number of Unique Samples per Cancer Type")
plt.xlabel("Cancer Type")
plt.ylabel("Number of Samples")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "samples_per_cancer_type.png"))
plt.close()

# ------------------------------------------
# PLOT 6: Peptide length distribution
# ------------------------------------------
length_data = []
for (tool, cancer_type, threshold, hla_type), peptides in peptide_sample_map.items():
    for peptide in peptides:
        length_data.append({
            "Tool": tool,
            "CancerType": cancer_type,
            "Threshold": threshold,
            "PeptideLength": len(peptide)
        })
df_lengths = pd.DataFrame(length_data)

plt.figure(figsize=(10, 6))
sns.histplot(data=df_lengths, x="PeptideLength", hue="CancerType", multiple="stack", binwidth=1)
plt.title("Peptide Length Distribution by Cancer Type")
plt.xlabel("Peptide Length")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "peptide_length_distribution.png"))
plt.close()

# ------------------------------------------
# PLOT 7: Top Peptides and Genes per cancer type
# ------------------------------------------
# Plot 1: Top Peptides
peptide_counter = defaultdict(Counter)

for (tool, cancer_type, threshold, hla_type), peptides in peptide_sample_map.items():
    for peptide, samples in peptides.items():
        peptide_counter[cancer_type][peptide] += len(samples)

top_peptides_df = pd.DataFrame([
    {"CancerType": ct, "Peptide": p, "Count": c}
    for ct, counter in peptide_counter.items()
    for p, c in counter.most_common(10)
])

plt.figure(figsize=(12, 8))
sns.barplot(data=top_peptides_df, x="Count", y="Peptide", hue="CancerType")
plt.title("Top Predicted Neoantigenic Peptides per Cancer Type")
plt.xlabel("Peptide Occurrence Count")
plt.ylabel("Peptide")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "top_neoantigen_peptides.png"))
plt.close()

# Plot 2: Top Genes
gene_counter = defaultdict(Counter)

for filepath in file_list:
    try:
        df = pd.read_excel(filepath)
        tool, cancer_type, threshold, hla_type = extract_metadata(filepath)
        for _, row in df.iterrows():
            gene = str(row['gene_symbol']).strip().upper()
            gene_counter[cancer_type][gene] += 1
    except Exception as e:
        print(f"Error processing {filepath}: {e}")

top_genes_df = pd.DataFrame([
    {"CancerType": ct, "Gene": g, "Count": c}
    for ct, counter in gene_counter.items()
    for g, c in counter.most_common(10)
])

plt.figure(figsize=(12, 8))
sns.barplot(data=top_genes_df, x="Count", y="Gene", hue="CancerType")
plt.title("Top Genes Contributing to Neoantigen Predictions per Cancer Type")
plt.xlabel("Gene Occurrence Count")
plt.ylabel("Gene")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "top_neoantigen_genes.png"))
plt.close()

# ------------------------------------------
# PLOT 8: Affinity per cancer type
# ------------------------------------------
# Step 1: Load scores and retain per-cancer threshold info
affinity_data = []

for filepath in file_list:
    tool, cancer_type, threshold, hla_type = extract_metadata(filepath)
    if tool == "deepmhci" and (cancer_type in ["SCLC", "NSCLC", "PAN"] and threshold == "0.1"):
        df = pd.read_excel(filepath)
        if "score" in df.columns:
            for val in df["score"]:
                affinity_data.append({
                    "CancerType": cancer_type,
                    "Score": val
                })

df_affinity = pd.DataFrame(affinity_data)
print(df_affinity["CancerType"].value_counts())

# Step 2: Individual plots
for cancer in df_affinity["CancerType"].unique():
    plt.figure(figsize=(10, 6))
    sub = df_affinity[df_affinity["CancerType"] == cancer]
    sns.histplot(data=sub, x="Score", bins=50, kde=True)

    plt.axvline(0.1, color="red", linestyle="--", label="Threshold = 0.1")
    plt.axvline(0.2, color="green", linestyle="--", label="Threshold = 0.2")
    plt.axvline(0.3, color="blue", linestyle="--", label="Threshold = 0.3")

    plt.title(f"Binding Score Distribution (DeepMHCi) - {cancer}")
    plt.xlabel("Predicted Binding Score (Affinity)")
    plt.ylabel("Neoantigen Count")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"deepmhci_score_distribution_{cancer}.png"))
    plt.close()


# Step 3: Combined Plots
plt.figure(figsize=(10, 6))
# Plot each cancer type manually to ensure proper labeling
cancer_types = df_affinity["CancerType"].unique()

for i, cancer in enumerate(cancer_types):
    subset = df_affinity[df_affinity["CancerType"] == cancer]
    sns.histplot(
        data=subset,
        x="Score",
        bins=50,
        kde=True,
        element="step",
        label=cancer
    )

# Threshold lines (no label here; legend will be added manually)
plt.axvline(0.1, color="red", linestyle="--")
plt.axvline(0.2, color="green", linestyle="--")
plt.axvline(0.3, color="blue", linestyle="--")

# Create custom threshold handles
threshold_handles = [
    Line2D([0], [0], color="red", linestyle="--", label="Threshold = 0.1"),
    Line2D([0], [0], color="green", linestyle="--", label="Threshold = 0.2"),
    Line2D([0], [0], color="blue", linestyle="--", label="Threshold = 0.3")
]

# First legend: Cancer types
legend1 = plt.legend(title="Cancer Type", loc="upper right")

# Second legend: Threshold lines
plt.gca().add_artist(legend1)  # Keep first legend
plt.legend(handles=threshold_handles, title="Thresholds", loc="upper center")

# Final polish
plt.title("Binding Score Distribution by Cancer Type (DeepMHCi)")
plt.xlabel("Predicted Binding Score (Affinity)")
plt.ylabel("Neoantigen Count")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "deepmhci_score_distribution_combined.png"))
plt.close()


# Step 4: Cumulative Plot
# Define cumulative thresholds to evaluate
import numpy as np
thresholds = np.arange(0.0, 1.01, 0.05)  # You can change the step size if needed

# Prepare data structure
cancer_types = df_affinity["CancerType"].unique()
threshold_counts = {cancer: [] for cancer in cancer_types}

# Count scores above each threshold per cancer type
for t in thresholds:
    for cancer in cancer_types:
        count = df_affinity[
            (df_affinity["CancerType"] == cancer) & (df_affinity["Score"] > t)
        ].shape[0]
        threshold_counts[cancer].append(count)

# Plot cumulative counts
plt.figure(figsize=(10, 6))

for cancer in cancer_types:
    plt.plot(thresholds, threshold_counts[cancer], label=cancer)

# Threshold markers
plt.axvline(0.1, color="red", linestyle="--")
plt.axvline(0.2, color="green", linestyle="--")
plt.axvline(0.3, color="blue", linestyle="--")

# Legends
plt.title("Cumulative Neoantigen Counts Above Thresholds (DeepMHCi)")
plt.xlabel("Predicted Binding Score (Affinity)")
plt.ylabel("Cumulative Neoantigen Count")
plt.legend(title="Cancer Type", loc="upper right")

# Save the cumulative plot
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "deepmhci_cumulative_thresholds.png"))
plt.close()

# ------------------------------------------
# PLOT 9: Shared Peptides by Samples
# ------------------------------------------
cancer_types = ["SCLC", "NSCLC", "PAN"]
shared_counts = {}
# Threshold settings
thresholds = {
    "SCLC": "0.1",
    "NSCLC": "0.1",
    "PAN": "0.1"
}
for cancer in cancer_types:
    peptide_to_samples = defaultdict(set)

    for hla in ["A", "B"]:
        filename = f"neoantigen_candidates_deepmhci_{cancer}_{thresholds[cancer]}_{hla}.xlsx"
        filepath = os.path.join(BASE_DIR, filename)
        if not os.path.exists(filepath):
            print(f"File not found: {filepath}")
            continue

        df = pd.read_excel(filepath)
        for _, row in df.iterrows():
            peptide = row["peptide_8to14"]
            sample = row["sample_id"]
            peptide_to_samples[peptide].add(sample)

    shared_distribution = defaultdict(int)
    for peptide, samples in peptide_to_samples.items():
        shared_distribution[len(samples)] += 1

    # Store results in DataFrame
    shared_df = pd.DataFrame({
        "NumSamples": list(shared_distribution.keys()),
        "NumPeptides": list(shared_distribution.values()),
        "CancerType": cancer
    })
    shared_counts[cancer] = shared_df

sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))

df_sorted = []
for cancer, df in shared_counts.items():
    df_sorted = df.sort_values("NumSamples")
    sns.lineplot(data=df_sorted, x="NumSamples", y="NumPeptides", marker="o", label=cancer)

plt.title("Distribution of Shared Neoantigens Across Samples")
plt.xlabel("Number of Samples Sharing a Peptide")
plt.ylabel("Number of Peptides")
plt.yscale("log")  # Optional: useful for long-tail distributions
plt.xticks(ticks=range(0, max(df_sorted["NumSamples"]) + 2, 2))  # set x-axis ticks every 2
plt.legend(title="Cancer Type")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "shared_peptides_distribution.png"))
plt.close()
