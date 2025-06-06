#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 14:02:41 2025

@author: junwkim
"""

import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from mpl_toolkits.mplot3d import Axes3D  # Import for 3D plotting

def sanitize_sample_id(sample_id):
    """
    Sanitize the sample ID by removing hyphens and the '_LUNG' suffix.
    """
    return sample_id.replace('-', '').replace('_LUNG', '')

def load_data(file_path, sheet_name=0, sample_id_col=None):
    """
    Load data from an Excel file.
    :param file_path: Path to the Excel file.
    :param sheet_name: Sheet name or index to load.
    :param sample_id_col: Column name for sample identifiers. If None, sample indices are used.
    :return: Tuple of (DataFrame with gene data, list of sample IDs)
    """
    try:
        df = pd.read_excel(file_path, sheet_name=sheet_name)
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        sys.exit(1)

    if sample_id_col and sample_id_col in df.columns:
        sample_ids = df[sample_id_col].tolist()
        gene_data = df.drop(columns=[sample_id_col])
    else:
        gene_data = df.copy()
        gene_data = gene_data.set_index(gene_data.columns[0])  # Assuming first column is gene names
        sample_ids = gene_data.columns.tolist()
        gene_data = gene_data.transpose().reset_index(drop=True)

    return gene_data, sample_ids

def load_subtypes(subtype_file):
    """
    Load subtype data and assign subtypes based on gene expression.
    :param subtype_file: Path to the subtype Excel file.
    :return: Dictionary mapping sanitized sample IDs to subtypes.
    """
    try:
        df = pd.read_excel(subtype_file, sheet_name=0)
    except Exception as e:
        print(f"Error reading subtype Excel file: {e}")
        sys.exit(1)

    # Assume the first column contains gene names
    gene_column = df.columns[0]
    df_genes = df.set_index(gene_column)

    # Transpose to have samples as rows
    df_samples = df_genes.transpose()

    # Sanitize sample IDs
    df_samples.index = [sanitize_sample_id(sid) for sid in df_samples.index]

    # Define target genes
    target_genes = ['YAP1', 'POU2F3', 'ASCL1', 'NEUROD1']

    # Check if all target genes are present
    missing_genes = [gene for gene in target_genes if gene not in df_samples.columns]
    if missing_genes:
        print(f"Missing target genes in subtype data: {missing_genes}")
        sys.exit(1)

    # Assign subtype based on highest expression among target genes
    def assign_subtype(row):
        expressions = row[target_genes]
        max_gene = expressions.idxmax()
        max_value = expressions.max()
        if expressions.tolist().count(max_value) == 1:
            if max_gene == 'YAP1':
                return 'y'
            elif max_gene == 'POU2F3':
                return 'p'
            elif max_gene == 'ASCL1':
                return 'a'
            elif max_gene == 'NEUROD1':
                return 'n'
        return 'unknown'

    df_samples['subtype'] = df_samples.apply(assign_subtype, axis=1)

    subtype_mapping = df_samples['subtype'].to_dict()
    return subtype_mapping

def preprocess_data(gene_data):
    """
    Preprocess the gene expression data: handle missing values and standardize.
    :param gene_data: DataFrame with gene expression data.
    :return: Standardized data as a NumPy array.
    """
    # Handle missing values by dropping them
    gene_data = gene_data.dropna()

    # Standardize the data (mean=0, variance=1)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(gene_data)

    return scaled_data

def perform_pca(scaled_data, n_components=3):
    """
    Perform PCA on the scaled data.
    :param scaled_data: Numpy array of scaled gene expression data.
    :param n_components: Number of principal components to generate.
    :return: PCA object and transformed data.
    """
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(scaled_data)
    return pca, principal_components

def plot_pca_2d(principal_components, sample_ids, subtypes, pca, title='PCA 2D Plot'):
    """
    Plot the first two principal components with subtype coloring.
    :param principal_components: Numpy array with principal components.
    :param sample_ids: List of sample identifiers.
    :param subtypes: List of subtype labels.
    :param pca: PCA object (to extract explained variance).
    :param title: Title of the plot.
    """
    pc_df = pd.DataFrame(data=principal_components[:, :2], columns=['PC1', 'PC2'])
    pc_df['Sample'] = sample_ids
    pc_df['Subtype'] = subtypes

    # Define color palette
    palette = {'y': 'red', 'p': 'blue', 'a': 'green', 'n': 'orange', 'unknown': 'gray'}

    plt.figure(figsize=(10, 8))
    sns.scatterplot(x='PC1', y='PC2', data=pc_df, hue='Subtype', palette=palette, s=100)

    for i in range(pc_df.shape[0]):
        plt.text(x=pc_df.PC1[i]+0.02, y=pc_df.PC2[i]+0.02, 
                 s=pc_df.Sample[i], fontdict=dict(color='black', size=10),
                 horizontalalignment='left', verticalalignment='bottom')

    plt.title(f"{title}\n"
              f"PC1: {pca.explained_variance_ratio_[0]*100:.2f}% "
              f"PC2: {pca.explained_variance_ratio_[1]*100:.2f}% Variance")
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.legend(title='Subtype', loc='best')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_pca_3d(principal_components, sample_ids, subtypes, pca, title='PCA 3D Plot'):
    """
    Plot the first three principal components in 3D with subtype coloring.
    :param principal_components: Numpy array with principal components.
    :param sample_ids: List of sample identifiers.
    :param subtypes: List of subtype labels.
    :param pca: PCA object (to extract explained variance).
    :param title: Title of the plot.
    """
    pc_df = pd.DataFrame(data=principal_components[:, :3], columns=['PC1', 'PC2', 'PC3'])
    pc_df['Sample'] = sample_ids
    pc_df['Subtype'] = subtypes

    # Define color palette
    palette = {'y': 'red', 'p': 'blue', 'a': 'green', 'n': 'orange', 'unknown': 'gray'}
    colors = pc_df['Subtype'].map(palette)

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(pc_df.PC1, pc_df.PC2, pc_df.PC3, 
                         c=colors, s=100, depthshade=True)

    for i in range(pc_df.shape[0]):
        ax.text(pc_df.PC1[i]+0.02, pc_df.PC2[i]+0.02, pc_df.PC3[i]+0.02, 
                s=pc_df.Sample[i], fontdict=dict(color='black', size=10))

    ax.set_title(f"{title}\n"
                 f"PC1: {pca.explained_variance_ratio_[0]*100:.2f}% "
                 f"PC2: {pca.explained_variance_ratio_[1]*100:.2f}% "
                 f"PC3: {pca.explained_variance_ratio_[2]*100:.2f}% Variance")
    ax.set_xlabel('Principal Component 1')
    ax.set_ylabel('Principal Component 2')
    ax.set_zlabel('Principal Component 3')

    # Create custom legend
    handles = [plt.Line2D([0], [0], marker='o', color='w', label=key,
                          markerfacecolor=value, markersize=10) 
               for key, value in palette.items()]
    ax.legend(handles=handles, title='Subtype', loc='best')

    plt.tight_layout()
    plt.show()

def perform_umap(scaled_data, n_components=2, n_neighbors=15, min_dist=0.1, random_state=42):
    """
    Perform UMAP dimensionality reduction.
    :param scaled_data: Numpy array of scaled gene expression data.
    :param n_components: Number of dimensions for UMAP.
    :param n_neighbors: The size of local neighborhood used for UMAP.
    :param min_dist: The effective minimum distance between embedded points.
    :param random_state: Seed used by the random number generator.
    :return: UMAP object and transformed data.
    """
    umap_reducer = umap.UMAP(n_components=n_components, n_neighbors=n_neighbors, 
                             min_dist=min_dist, random_state=random_state)
    umap_embedding = umap_reducer.fit_transform(scaled_data)
    return umap_reducer, umap_embedding

def plot_umap(umap_embedding, sample_ids, subtypes, title='UMAP Plot'):
    """
    Plot the UMAP embedding with subtype coloring.
    :param umap_embedding: Numpy array with UMAP dimensions.
    :param sample_ids: List of sample identifiers.
    :param subtypes: List of subtype labels.
    :param title: Title of the plot.
    """
    umap_df = pd.DataFrame(data=umap_embedding, columns=['UMAP1', 'UMAP2'])
    umap_df['Sample'] = sample_ids
    umap_df['Subtype'] = subtypes

    # Define color palette
    palette = {'y': 'red', 'p': 'blue', 'a': 'green', 'n': 'orange', 'unknown': 'gray'}

    plt.figure(figsize=(10, 8))
    sns.scatterplot(x='UMAP1', y='UMAP2', data=umap_df, hue='Subtype', palette=palette, s=100)

    for i in range(umap_df.shape[0]):
        plt.text(x=umap_df.UMAP1[i]+0.02, y=umap_df.UMAP2[i]+0.02, 
                 s=umap_df.Sample[i], fontdict=dict(color='black', size=10),
                 horizontalalignment='left', verticalalignment='bottom')

    plt.title(title)
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.legend(title='Subtype', loc='best')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def perform_tsne(scaled_data, n_components=2, perplexity=30, learning_rate=200, n_iter=1000, random_state=42):
    """
    Perform t-SNE dimensionality reduction.
    :param scaled_data: Numpy array of scaled gene expression data.
    :param n_components: Number of dimensions for t-SNE.
    :param perplexity: The perplexity parameter for t-SNE.
    :param learning_rate: The learning rate for t-SNE optimization.
    :param n_iter: Number of iterations for t-SNE.
    :param random_state: Seed used by the random number generator.
    :return: t-SNE object and transformed data.
    """
    tsne = TSNE(n_components=n_components, perplexity=perplexity, learning_rate=learning_rate, 
                n_iter=n_iter, random_state=random_state)
    tsne_embedding = tsne.fit_transform(scaled_data)
    return tsne, tsne_embedding

def plot_tsne(tsne_embedding, sample_ids, subtypes, title='t-SNE Plot'):
    """
    Plot the t-SNE embedding with subtype coloring.
    :param tsne_embedding: Numpy array with t-SNE dimensions.
    :param sample_ids: List of sample identifiers.
    :param subtypes: List of subtype labels.
    :param title: Title of the plot.
    """
    tsne_df = pd.DataFrame(data=tsne_embedding, columns=['tSNE1', 'tSNE2'])
    tsne_df['Sample'] = sample_ids
    tsne_df['Subtype'] = subtypes

    # Define color palette
    palette = {'y': 'red', 'p': 'blue', 'a': 'green', 'n': 'orange', 'unknown': 'gray'}

    plt.figure(figsize=(10, 8))
    sns.scatterplot(x='tSNE1', y='tSNE2', data=tsne_df, hue='Subtype', palette=palette, s=100)

    for i in range(tsne_df.shape[0]):
        plt.text(x=tsne_df.tSNE1[i]+0.5, y=tsne_df.tSNE2[i]+0.5, 
                 s=tsne_df.Sample[i], fontdict=dict(color='black', size=10),
                 horizontalalignment='left', verticalalignment='bottom')

    plt.title(title)
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.legend(title='Subtype', loc='best')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def main():
    # -------------------------------
    # User Input Parameters
    # -------------------------------
    gene_counts_file = 'gene_counts_per_sample_limit=1_combined.xlsx'  # Replace with your file path
    subtype_file = 'sclc_subtype.xlsx'  # Subtype definition file
    sheet_name = 0  # You can specify sheet name or index
    sample_id_col = "sample_id" # Replace with your sample ID column name if exists, else None
    title_pca = 'PCA of Gene Expression Data'
    title_umap = 'UMAP of Gene Expression Data'
    title_tsne = 't-SNE of Gene Expression Data'
    # -------------------------------

    # Load main gene expression data
    gene_data, sample_ids = load_data(gene_counts_file, sheet_name, sample_id_col)

    # Sanitize sample IDs
    sanitized_ids = [sanitize_sample_id(sid) for sid in sample_ids]

    # Load subtype mappings
    subtype_mapping = load_subtypes(subtype_file)

    # Assign subtypes to samples
    subtypes = [subtype_mapping.get(sid, 'unknown') for sid in sanitized_ids]
    subtype_df = pd.DataFrame({
        'Original_Sample_ID': sample_ids,
        'Sanitized_Sample_ID': sanitized_ids,
        'Subtype': subtypes
    })
    
    # Save the subtype mapping to an Excel file
    subtype_output_file = 'subtype_mapping.xlsx'  # You can change the file name/path as needed
    try:
        subtype_df.to_excel(subtype_output_file, index=False)
        print(f"Subtype mapping successfully saved to {subtype_output_file}")
    except Exception as e:
        print(f"Error saving subtype mapping to Excel: {e}")
    # Preprocess data
    scaled_data = preprocess_data(gene_data)

    # Perform PCA
    n_components_pca = 3  # Number of components for 3D plot
    pca, principal_components = perform_pca(scaled_data, n_components=n_components_pca)

    # Plot PCA (2D)
    plot_pca_2d(principal_components, sample_ids, subtypes, pca, title_pca)

    # Plot PCA (3D)
    plot_pca_3d(principal_components, sample_ids, subtypes, pca, title_pca)

    # Perform UMAP
    umap_reducer, umap_embedding = perform_umap(scaled_data, n_components=2)

    # Plot UMAP
    plot_umap(umap_embedding, sample_ids, subtypes, title_umap)

    # Perform t-SNE
    tsne_obj, tsne_embedding = perform_tsne(scaled_data, n_components=2, perplexity=30, 
                                           learning_rate=200, n_iter=1000)

    # Plot t-SNE
    plot_tsne(tsne_embedding, sample_ids, subtypes, title_tsne)

if __name__ == "__main__":
    main()

