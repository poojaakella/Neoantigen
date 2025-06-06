#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 21:19:59 2025
GPT was used for debugging/organization/library suggestions
@author: junwkim
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from scipy.stats import pearsonr  # Import pearsonr for correlation
import numpy as np

# File paths (modify these if your files are in different directories)
gene_counts_file = 'results/pan/gene_counts_per_sample_deepmhci_pan_limit=2_combined.csv'
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

# 3. Process sample names in gene_counts_df
print("Processing sample names in gene counts data...")
gene_counts_samples = gene_counts_df.columns.tolist()
processed_gene_counts_samples = [sample.replace('_LUNG', '') for sample in gene_counts_samples]
gene_counts_df.columns = processed_gene_counts_samples

# 4. Process sample names in tpms_df
print("Processing sample names in TPMS data...")
tpms_samples = tpms_df.columns.tolist()

def process_tpms_sample_name(name):
    # Remove hyphen between 'NCI' and sample identifier
    name_no_hyphen = re.sub(r'^(NCI)-', r'\1', name)
    # Alternatively, remove all hyphens if needed
    name_no_hyphen = name_no_hyphen.replace('-', '')
    # Remove everything after the comma
    name_clean = name_no_hyphen.split(',')[0]
    return name_clean

processed_tpms_samples = [process_tpms_sample_name(sample) for sample in tpms_samples]
tpms_df.columns = processed_tpms_samples

# 5. Find common samples
print("Finding common samples...")
common_samples = list(set(gene_counts_df.columns) & set(tpms_df.columns))
common_samples.sort()
print(f"Number of common samples: {len(common_samples)}")

if len(common_samples) == 0:
    raise ValueError("No common samples found after processing sample names.")

# 6. Subset both DataFrames to common samples
gene_counts_common = gene_counts_df[common_samples]
tpms_common = tpms_df[common_samples]

# 7. Find common genes
print("Finding common genes...")
common_genes = list(set(gene_counts_common.index) & set(tpms_common.index))
common_genes.sort()
print(f"Number of common genes: {len(common_genes)}")

if len(common_genes) == 0:
    raise ValueError("No common genes found between the two datasets.")

# 8. Subset both DataFrames to common genes
gene_counts_common = gene_counts_common.loc[common_genes]
tpms_common = tpms_common.loc[common_genes]

# **NEW: Set the index name to 'Gene'**
gene_counts_common.index.name = 'Gene'
tpms_common.index.name = 'Gene'

# 9. Melt the DataFrames to long format
print("Transforming data to long format...")
gene_counts_long = gene_counts_common.reset_index().melt(id_vars='Gene', var_name='Sample', value_name='Gene_Count')
tpms_long = tpms_common.reset_index().melt(id_vars='Gene', var_name='Sample', value_name='TPM')

# 10. Merge the two DataFrames on Gene and Sample
print("Merging datasets...")
merged_df = pd.merge(gene_counts_long, tpms_long, on=['Gene', 'Sample'])

# 11. Drop any rows with missing values
merged_df.dropna(subset=['Gene_Count', 'TPM'], inplace=True)

# **NEW: Filter out rows where Gene_Count is zero**
merged_df = merged_df[merged_df['Gene_Count'] > 2 ]

# Optional: Limit the number of points if the dataset is too large
# For example, randomly sample 10000 points
# merged_df = merged_df.sample(n=10000, random_state=42)

# 12. Plotting the dot plot
print("Creating the dot plot...")
plt.figure(figsize=(10, 8))
sns.set(style="whitegrid")
merged_df['log_TPM'] = np.log10(merged_df['TPM'] + 1)
merged_df['log_Gene_Count'] = np.log10(merged_df['Gene_Count'] + 1)
# Create scatter plot
ax = sns.scatterplot(
    data=merged_df, 
    x='log_TPM', 
    y='log_Gene_Count', 
    hue='Gene', 
    alpha=0.6, 
    edgecolor=None, 
    palette='viridis'
)
#ax.set_xscale('log', base=2)
plt.xlabel('TPM (E-MTAB-2770)')
plt.ylabel('Gene Counts (gene_counts_per_sample_limit=0)')
plt.title('Dot Plot: TPM vs Gene Counts per Sample and Gene')

# Optional: Improve legend (if too many genes, it may be cluttered)
# You might want to remove hue or use different visualization
# Here, we'll remove hue if too many unique genes
if merged_df['Gene'].nunique() > 20:
    plt.legend([], [], frameon=False)
else:
    plt.legend(title='Gene', bbox_to_anchor=(1.05, 1), loc='upper left')

# **NEW: Calculate and overlay correlation coefficient and p-value**
# Calculate Pearson correlation
print("Calculating Pearson correlation...")
corr, pval = pearsonr(merged_df['log_TPM'], merged_df['log_Gene_Count'])
print(f"Pearson correlation coefficient: {corr:.4f}, p-value: {pval:.4e}")

# Prepare the text to display
correlation_text = f'Pearson r = {corr:.2f}\nP-value = {pval:.2e}'

# Add the text box to the plot
# Positioning the text in the top-left corner inside the plot
plt.text(
    0.05, 0.95, correlation_text, 
    transform=ax.transAxes,
    fontsize=12, 
    verticalalignment='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.5)
)

plt.tight_layout()
plt.show()



