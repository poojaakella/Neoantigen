#GPT was used for debugging/organization/library suggestions
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import numpy as np

# File paths (modify these if your files are in different directories)
gene_counts_file = 'results/SCLC/gene_counts_per_sample_limit=9_combined.csv'
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

# 8. Initialize list to store per-gene results
print("Calculating per-gene expression ratios...")
results = []

epsilon = 1e-9  # To avoid division by zero

for gene in common_genes:
    # Get neoantigen counts for the current gene across all samples
    neo_counts = gene_counts_common.loc[gene]
    
    # Define samples with at least one neoantigen for this gene
    samples_with_neo = neo_counts[neo_counts >= 1].index.tolist()
    
    # Define samples without any neoantigen for this gene
    samples_without_neo = neo_counts[neo_counts == 0].index.tolist()
    
    # Ensure there are samples in both groups
    if len(samples_with_neo) == 0 or len(samples_without_neo) == 0:
        # You can choose to skip the gene or assign NaN
        """ratio = np.nan
        log2_ratio = np.nan
        avg_tpm_without = tpms_common.loc[gene].mean() if len(samples_without_neo) == 0 else tpms_common.loc[gene, samples_without_neo].mean()
        avg_tpm_with = tpms_common.loc[gene].mean() if len(samples_with_neo) == 0 else tpms_common.loc[gene, samples_with_neo].mean()"""
        continue
    else:
        # Calculate average TPM in samples without neoantigens
        avg_tpm_without = tpms_common.loc[gene, samples_without_neo].mean()
        
        # Calculate average TPM in samples with at least one neoantigen
        avg_tpm_with = tpms_common.loc[gene, samples_with_neo].mean()
        
        # Calculate ratio
        ratio = avg_tpm_with / (avg_tpm_without + epsilon)
        
        # Calculate Log2 ratio
        log2_ratio = np.log2(ratio + epsilon)
    
    # Append the results
    results.append({
        'Gene': gene,
        'Avg_TPM_Without_Neo': avg_tpm_without,
        'Avg_TPM_With_Neo': avg_tpm_with,
        'Ratio_Without_With': ratio,
        'Log2_Ratio': log2_ratio
    })

# Create DataFrame from results
ratio_df = pd.DataFrame(results)

# Optional: Drop genes where ratio could not be calculated
ratio_df_clean = ratio_df.dropna(subset=['Ratio_Without_With'])
ratio_df_clean["Ratio_Without_With"] = pd.to_numeric(ratio_df_clean["Ratio_Without_With"], errors="coerce")
ratio_df_clean["Log2_Ratio"] = pd.to_numeric(ratio_df_clean["Log2_Ratio"], errors="coerce")
ratio_df_clean = ratio_df_clean.dropna(subset=['Log2_Ratio'])

print(f"Number of genes with valid ratio calculations: {len(ratio_df_clean)}")

# 9. Plotting the Ratios
print("Creating the ratio plot...")

#plt.figure(figsize=(max(12, len(ratio_df_clean)*0.4), max(6, len(ratio_df_clean)*0.3)))  # Adjust size based on number of genes


# Sort genes by ratio for better visualization
ratio_df_sorted = ratio_df_clean.sort_values(by='Ratio_Without_With', ascending=True)

ratio_df_sorted.to_csv('gene_expression_ratios_without_vs_with_neoantigens.csv', index=False)
print("Ratio dataframe saved to 'gene_expression_ratios_without_vs_with_neoantigens.csv'.")

"""genes_ordered = ratio_df_sorted['Gene'].tolist()

# Convert 'Gene' column to a categorical type with the specified order
ratio_df_sorted['Gene'] = pd.Categorical(ratio_df_sorted['Gene'], categories=genes_ordered, ordered=True)

# Initialize the matplotlib figure
max_pixels = 65500
dpi = 100
# Calculate maximum height in inches based on dpi
max_height_inch = max_pixels / dpi  # 65500 / 100 = 655 inches

# Calculate required height based on number of genes
# Each gene needs a minimal amount of space, e.g., 0.01 inches
required_height = len(genes_ordered) * 0.01
# Set the figure height to the required height or the maximum allowed, whichever is smaller
fig_height = min(max(required_height, 6), 20)  # Setting an upper limit of 20 inches

plt.figure(figsize=(10, fig_height))

# Set Seaborn style
sns.set(style="whitegrid")

# Create a scatter plot (dot plot) with Genes on the Y-axis and Ratio on the X-axis
ax = sns.scatterplot(
    x='Ratio_Without_With', 
    y='Gene', 
    data=ratio_df_sorted, 
    color='steelblue',
    s=10  # Reduced size of the dots to accommodate more genes
)

# Add a vertical line at x=1 for reference
ax.axvline(x=1, color='red', linestyle='--', linewidth=1)

# Set labels and title
ax.set_xlabel('Expression Ratio (With Neoantigens / Without Neoantigens)')
ax.set_ylabel('Gene')
ax.set_title('Gene Expression Ratios Ordered in Ascending Order')

# Improve layout without using tight_layout to prevent automatic resizing issues
plt.tight_layout()"""


"""
plt.figure(figsize=(12, 8))
sns.set(style="whitegrid")
# Using Seaborn's barplot
ax = sns.barplot(
    x='Gene', 
    y='Log2_Ratio', 
    data=ratio_df_sorted, 
    palette='viridis'
)
ax.set_xlabel('Expression Ratio (With Neoantigens / Without Neoantigens)')
ax.set_yticklabels([])
ax.set_xlim(0, 1)
plt.tight_layout()
plt.show()"""

"""# 9. Choose a Visualization Method

# Option A: Histogram of Log2 Ratios
plt.figure(figsize=(10, 6))
sns.set(style="whitegrid")
sns.histplot(ratio_df_clean['Ratio_Without_With'], bins=5, kde=True, color='skyblue')
plt.xlabel('Log2 Ratio (With / Without Neoantigens)')
plt.ylabel('Number of Genes')
plt.title('Distribution of Gene Expression Ratios')
plt.axvline(0, color='red', linestyle='--', linewidth=1)
plt.tight_layout()
plt.show()
"""

"""# Create a bar plot
ax = sns.barplot(
    data=ratio_df_sorted,
    x='Gene',
    y='Ratio_Without_With',
    palette='viridis'
)

# Rotate x-axis labels for better readability if there are many genes
if len(ratio_df_sorted) > 20:
    plt.xticks(rotation=90)
else:
    plt.xticks(rotation=45)

plt.xlabel('Gene')
plt.ylabel('Average TPM Ratio (Without / With Neoantigens)')
plt.title('Gene-Specific Expression Ratios: Without vs With Neoantigens')

# Add a horizontal line at ratio=1 for reference
plt.axhline(1, color='red', linestyle='--', linewidth=1)

# Optional: Annotate bars with ratio values
for container in ax.containers:
    ax.bar_label(container, fmt='%.2f', rotation=90, fontsize=8)

plt.tight_layout()
plt.show()

# Optional: Save the ratio dataframe to a CSV file
ratio_df_sorted.to_csv('gene_expression_ratios_without_vs_with_neoantigens.csv', index=False)
print("Ratio dataframe saved to 'gene_expression_ratios_without_vs_with_neoantigens.csv'.")"""
