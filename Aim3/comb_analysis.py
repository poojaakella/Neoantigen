"""
Created on Wed May 21 11:24:29 2025
@author: junwkim
GPT was used for debugging/organization/library suggestions
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

def load_data(file_path, sheet_name=0, sample_id_col=None):
    """
    Load data from an Excel file.

    :param file_path: Path to the Excel file.
    :param sheet_name: Sheet name or index to load.
    :param sample_id_col: Column name for sample identifiers. If None, sample indices are used.
    :return: Tuple of (DataFrame with gene data, list of sample IDs)
    """
    try:
        df = pd.read_csv(file_path)
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

def perform_umap(scaled_data, n_components=2, n_neighbors=15, min_dist=0.1, random_state=42):
    """
    Perform UMAP dimensionality reduction.

    :param scaled_data: Numpy array of scaled gene expression data.
    :param n_components: Number of dimensions for UMAP.
    :param n_neighbors: The size of local neighborhood (in terms of number of neighboring sample points) used for UMAP.
    :param min_dist: The effective minimum distance between embedded points.
    :param random_state: Seed used by the random number generator.
    :return: UMAP object and transformed data.
    """
    umap_reducer = umap.UMAP(n_components=n_components, n_neighbors=n_neighbors, 
                             min_dist=min_dist, random_state=random_state)
    umap_embedding = umap_reducer.fit_transform(scaled_data)
    return umap_reducer, umap_embedding

def plot_umap(umap_embedding, sample_ids, title='UMAP Plot'):
    """
    Plot the UMAP embedding.

    :param umap_embedding: Numpy array with UMAP dimensions.
    :param sample_ids: List of sample identifiers.
    :param title: Title of the plot.
    """
    umap_df = pd.DataFrame(data=umap_embedding, columns=['UMAP1', 'UMAP2'])
    umap_df['Sample'] = sample_ids
    
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x='UMAP1', y='UMAP2', data=umap_df, s=100, palette='viridis')
    
    for i in range(umap_df.shape[0]):
        plt.text(x=umap_df.UMAP1[i]+0.02, y=umap_df.UMAP2[i]+0.02, s=umap_df.Sample[i], 
                 fontdict=dict(color='black', size=10),
                 horizontalalignment='left', verticalalignment='bottom')
    
    plt.title(title)
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
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
    tsne = TSNE(n_components=n_components, perplexity=perplexity, 
                learning_rate=learning_rate, n_iter=n_iter, random_state=random_state)
    tsne_embedding = tsne.fit_transform(scaled_data)
    return tsne, tsne_embedding

def plot_tsne(tsne_embedding, sample_ids, title='t-SNE Plot'):
    """
    Plot the t-SNE embedding.

    :param tsne_embedding: Numpy array with t-SNE dimensions.
    :param sample_ids: List of sample identifiers.
    :param title: Title of the plot.
    """
    tsne_df = pd.DataFrame(data=tsne_embedding, columns=['tSNE1', 'tSNE2'])
    tsne_df['Sample'] = sample_ids
    
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x='tSNE1', y='tSNE2', data=tsne_df, s=100, palette='plasma')
    
    for i in range(tsne_df.shape[0]):
        plt.text(x=tsne_df.tSNE1[i]+0.5, y=tsne_df.tSNE2[i]+0.5, s=tsne_df.Sample[i], 
                 fontdict=dict(color='black', size=10),
                 horizontalalignment='left', verticalalignment='bottom')
    
    plt.title(title)
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def main():
    # ------------------------------- #
    # User Input Parameters
    # ------------------------------- #
    excel_file_path = 'gene_counts_per_sample_deepmhci_pan_0.2_B.csv'  # Replace with your file path
    sheet_name = 0  # You can specify sheet name or index
    sample_id_col = "sample_id"  # Replace with your sample ID column name if exists, else None
    title_pca = 'PCA of Gene Expression Data'
    title_umap = 'UMAP of Gene Expression Data'
    title_tsne = 't-SNE of Gene Expression Data'
    # ------------------------------- #

    # Load data
    gene_data, sample_ids = load_data(excel_file_path, sheet_name, sample_id_col)
    
    # Preprocess data
    scaled_data = preprocess_data(gene_data)
    
    # Perform PCA
    n_components_pca = 3  # Number of components for 3D plot
    pca, principal_components = perform_pca(scaled_data, n_components=n_components_pca)
    
    # Plot PCA (2D)
    plot_pca_2d(principal_components, sample_ids, pca, title_pca)
    
    # Plot PCA (3D)
    plot_pca_3d(principal_components, sample_ids, pca, title_pca)
    
    # Perform UMAP
    umap_reducer, umap_embedding = perform_umap(scaled_data, n_components=3)
    
    # Plot UMAP
    plot_umap(umap_embedding, sample_ids, title_umap)
    
    # Perform t-SNE
    tsne_obj, tsne_embedding = perform_tsne(scaled_data, n_components=3, perplexity=30, learning_rate=200, n_iter=1000)
    
    # Plot t-SNE
    plot_tsne(tsne_embedding, sample_ids, title_tsne)

if __name__ == "__main__":
    main()
