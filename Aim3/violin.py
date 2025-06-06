import matplotlib.pyplot as plt
import seaborn as sns
import re
import numpy as np
import pandas as pd

# File paths (modify these if your files are in different directories)
gene_counts_file = 'results/SCLC/gene_counts_per_sample_deepmhci_SCLC_limit=2_combined.csv'
tpms_file = '../data/E-MTAB-2770-query-results.tpms.tsv'  # Assuming it's a TSV file

# 1. Read the gene counts Excel file
print("Reading gene counts data...")
gene_counts_df = pd.read_csv(gene_counts_file, index_col=0).transpose()
print(f"Gene counts data shape: {gene_counts_df.shape}")

# 2. Read the TPMS file
print("Reading TPMS data...")
# Determine if it's a TSV or CSV based on the file extension or content
# Here, assuming it's a TSV. Modify 'sep' if it's different.
if tpms_file.endswith('.tsv') or tpms_file.endswith('.tpms'):
    tpms_df = pd.read_csv(tpms_file, sep='\t', index_col=0)
else:
    tpms_df = pd.read_csv(tpms_file, sep=',', index_col=0)
print(f"TPMS data shape: {tpms_df.shape}")

# 3. Process sample names in tpms_df
print("Processing sample names in TPMS data...")

def process_tpms_sample_name(name):
    # Remove hyphen between 'NCI' and sample identifier
    name_no_hyphen = re.sub(r'^(NCI)-', r'\1', name)
    # Alternatively, remove all hyphens if needed
    name_no_hyphen = name_no_hyphen.replace('-', '')
    # Remove everything after the comma
    name_clean = name_no_hyphen.split(',')[0]
    return name_clean

processed_tpms_samples = [process_tpms_sample_name(sample) for sample in tpms_df.columns]
tpms_df.columns = processed_tpms_samples

# 4. Find common samples
print("Finding common samples...")
common_samples = list(set(gene_counts_df.columns) & set(tpms_df.columns))
common_samples.sort()
print(f"Number of common samples: {len(common_samples)}")

if len(common_samples) == 0:
    raise ValueError("No common samples found after processing sample names.")

# 5. Subset both DataFrames to common samples
gene_counts_common = gene_counts_df[common_samples]
tpms_common = tpms_df[common_samples]

# 6. Find common genes
print("Finding common genes...")
common_genes = list(set(gene_counts_common.index) & set(tpms_common.index))
common_genes.sort()
print(f"Number of common genes: {len(common_genes)}")

if len(common_genes) == 0:
    raise ValueError("No common genes found between the two datasets.")

# 7. Subset both DataFrames to common genes
gene_counts_common = gene_counts_common.loc[common_genes]
tpms_common = tpms_common.loc[common_genes]

# **NEW: Set the index name to 'Gene'**
gene_counts_common.index.name = 'Gene'
tpms_common.index.name = 'Gene'

# **NEW: Classify genes into 'With Neoantigen' and 'Without Neoantigen'**

print("Classifying genes based on neoantigen presence...")

# A gene is considered to have neoantigens if it has at least one neoantigen in any sample
genes_with_neo = gene_counts_common[(gene_counts_common >= 1).any(axis=1)].index.tolist()
genes_without_neo = gene_counts_common[(gene_counts_common == 0).all(axis=1)].index.tolist()

print(f"Number of genes with neoantigens: {len(genes_with_neo)}")
print(f"Number of genes without neoantigens: {len(genes_without_neo)}")

# 8. Calculate average TPMs for each group

print("Calculating average TPMs for each gene...")

# Calculate average TPM across all samples for each gene
tpms_common['Avg_TPM'] = tpms_common.mean(axis=1)

# Create a DataFrame with gene and its average TPM
avg_tpm_df = tpms_common[['Avg_TPM']].copy()

# Check for duplicate genes
duplicate_genes = avg_tpm_df.index.duplicated().sum()
if duplicate_genes > 0:
    print(f"Found {duplicate_genes} duplicate gene names. Aggregating TPMs by gene.")
    # Perform grouping and averaging only on numeric columns to avoid issues with non-numeric data
    avg_tpm_df = avg_tpm_df.groupby(avg_tpm_df.index).mean()

# Reset index to ensure 'Gene' is a column
avg_tpm_df = avg_tpm_df.reset_index()

# **Move the Group assignment here after aggregation**

# Add a new column indicating the group
avg_tpm_df['Group'] = avg_tpm_df['Gene'].isin(genes_with_neo)
avg_tpm_df['Group'] = avg_tpm_df['Group'].map({True: 'With Neoantigen', False: 'Without Neoantigen'})

# Ensure there are genes in both groups
if avg_tpm_df['Group'].nunique() != 2:
    raise ValueError("One of the groups ('With Neoantigen' or 'Without Neoantigen') has no genes.")

# 9. Plotting the Violin Plot

print("Creating the violin plot...")

plt.figure(figsize=(8, 6))  # Adjust the size as needed

# Set Seaborn style
sns.set(style="whitegrid")

# Optionally, define a palette
palette = {'With Neoantigen': 'skyblue', 'Without Neoantigen': 'lightgreen'}

# Create the violin plot
sns.violinplot(x='Group', y='Avg_TPM', data=avg_tpm_df, palette=palette, inner="box")

# Set labels and title
plt.xlabel('Gene Group')
plt.ylabel('Average TPM')
plt.title('Comparison of Average TPMs: Genes With vs. Without Neoantigens')

# Improve layout
plt.tight_layout()

# Save the plot (optional)
plt.savefig('violin_plot_avg_tpm_with_vs_without_neoantigens.png', dpi=300)

# Show the plot
plt.show()

# Optional: Save the average TPM dataframe for reference
avg_tpm_df.to_csv('average_tpm_with_vs_without_neoantigens.csv', index=False)
print("Average TPM dataframe saved to 'average_tpm_with_vs_without_neoantigens.csv'.")
print("Violin plot saved to 'violin_plot_avg_tpm_with_vs_without_neoantigens.png'.")









