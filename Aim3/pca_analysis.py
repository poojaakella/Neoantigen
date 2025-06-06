#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 10:37:26 2025

@author: junwkim
"""
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from mpl_toolkits.mplot3d import Axes3D  # Import for 3D plotting

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
        sample_ids = df.index.astype(str).tolist()
        gene_data = df.copy()
    
    return gene_data, sample_ids

def preprocess_data(gene_data):
    """
    Preprocess the gene expression data: handle missing values and standardize.

    :param gene_data: DataFrame with gene expression data.
    :return: Standardized data.
    """
    # Handle missing values if any (optional: here we drop them)
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

def plot_pca_2d(principal_components, sample_ids, pca, title='PCA 2D Plot'):
    """
    Plot the first two principal components.

    :param principal_components: Numpy array with principal components.
    :param sample_ids: List of sample identifiers.
    :param pca: PCA object (to extract explained variance).
    :param title: Title of the plot.
    """
    pc_df = pd.DataFrame(data=principal_components[:, :2], columns=['PC1', 'PC2'])
    pc_df['Sample'] = sample_ids
    
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x='PC1', y='PC2', data=pc_df, s=100)
    
    for i in range(pc_df.shape[0]):
        plt.text(x=pc_df.PC1[i]+0.02, y=pc_df.PC2[i]+0.02, s=pc_df.Sample[i], 
                 fontdict=dict(color='black', size=10),
                 horizontalalignment='left', verticalalignment='bottom')
    
    plt.title(f"{title}\n"
              f"PC1: {pca.explained_variance_ratio_[0]*100:.2f}% "
              f"PC2: {pca.explained_variance_ratio_[1]*100:.2f}% Variance")
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_pca_3d(principal_components, sample_ids, pca, title='PCA 3D Plot'):
    """
    Plot the first three principal components in 3D.

    :param principal_components: Numpy array with principal components.
    :param sample_ids: List of sample identifiers.
    :param pca: PCA object (to extract explained variance).
    :param title: Title of the plot.
    """
    pc_df = pd.DataFrame(data=principal_components[:, :3], columns=['PC1', 'PC2', 'PC3'])
    pc_df['Sample'] = sample_ids
    
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    scatter = ax.scatter(pc_df.PC1, pc_df.PC2, pc_df.PC3, s=100, c='b', depthshade=True)
    
    for i in range(pc_df.shape[0]):
        ax.text(pc_df.PC1[i]+0.02, pc_df.PC2[i]+0.02, pc_df.PC3[i]+0.02, 
                s=pc_df.Sample[i], 
                fontdict=dict(color='black', size=10))
    
    ax.set_title(f"{title}\n"
                 f"PC1: {pca.explained_variance_ratio_[0]*100:.2f}% "
                 f"PC2: {pca.explained_variance_ratio_[1]*100:.2f}% "
                 f"PC3: {pca.explained_variance_ratio_[2]*100:.2f}% Variance")
    ax.set_xlabel('Principal Component 1')
    ax.set_ylabel('Principal Component 2')
    ax.set_zlabel('Principal Component 3')
    plt.tight_layout()
    plt.show()

def main():
    # ------------------------------- #
    # User Input Parameters
    # ------------------------------- #
    excel_file_path = 'gene_counts_per_sample_limit=0.3.xlsx'  # Replace with your file path
    sheet_name = 0  # You can specify sheet name or index
    sample_id_col = "sample_id"  # Replace with your sample ID column name if exists, else None
    title = 'PCA of Gene Expression Data'
    # ------------------------------- #

    # Load data
    gene_data, sample_ids = load_data(excel_file_path, sheet_name, sample_id_col)
    
    # Preprocess data
    scaled_data = preprocess_data(gene_data)
    
    # Perform PCA
    n_components = 3  # Number of components for 3D plot
    pca, principal_components = perform_pca(scaled_data, n_components=n_components)
    
    # Plot PCA (2D)
    plot_pca_2d(principal_components, sample_ids, pca, title)
    
    # Plot PCA (3D)
    plot_pca_3d(principal_components, sample_ids, pca, title)

if __name__ == "__main__":
    main()



