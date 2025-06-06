#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 18:27:49 2025
GPT was used for debugging/organization/library suggestions
@author: junwkim
"""
import pandas as pd
import glob
import os

def combine_gene_count_files(input_directory, output_directory):
    """
    Combines Excel files with the same 'limit' in their filenames by adding overlapping cells
    and taking the union of unique cells.

    Parameters:
    - input_directory: Path to the directory containing the input Excel files.
    - output_directory: Path to the directory where combined Excel files will be saved.
    """

    # Ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Define the pattern to match all relevant Excel files
    pattern = os.path.join(input_directory, "gene_counts_per_sample_deepmhci_SCLC_limit=*")

    # Find all matching files
    all_files = glob.glob(pattern)

    if not all_files:
        print("No files found matching the pattern.")
        return

    # Extract unique limits by removing the '_A' suffix if present
    limits = set()
    for file in all_files:
        basename = os.path.basename(file)
        # Remove prefix and possible '_A' suffix
        if basename.startswith("gene_counts_per_sample_deepmhci_SCLC_limit="):
            limit_part = basename.replace("gene_counts_per_sample_deepmhci_SCLC_limit=", "")
            limit = limit_part.split('_A')[0]
            limits.add(limit)
    
    print(f"Found limits: {sorted(limits)}")

    # Process each limit
    for limit in sorted(limits):
        # Find all files corresponding to this limit (including those with '_A')
        limit_pattern = os.path.join(input_directory, f"gene_counts_per_sample_deepmhci_SCLC_limit={limit}*")
        limit_files = glob.glob(limit_pattern)

        if not limit_files:
            print(f"No files found for limit={limit}. Skipping.")
            continue

        print(f"\nProcessing limit={limit} with files:")
        for lf in limit_files:
            print(f"  - {os.path.basename(lf)}")

        # Initialize an empty DataFrame for accumulation
        combined_df = None

        for file in limit_files:
            # Read the Excel file
            try:
                df = pd.read_csv(file, index_col=0)
                print(f"Loaded {os.path.basename(file)} with shape {df.shape}")
            except Exception as e:
                print(f"Error reading {file}: {e}")
                continue

            # Ensure that all data is numeric; non-numeric cells are treated as zero
            df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

            if combined_df is None:
                combined_df = df
            else:
                # Combine by adding, aligning on both index and columns
                combined_df = combined_df.add(df, fill_value=0).astype(int)

        if combined_df is not None:
            # Optional: Sort the DataFrame by sample and gene names
            combined_df = combined_df.sort_index().sort_index(axis=1)

            # Define the output file path
            output_file = os.path.join(output_directory, f"gene_counts_per_sample_limit={limit}_combined.csv")

            # Save to Excel
            try:
                combined_df.to_csv(output_file)
                print(f"Saved combined data to {output_file} with shape {combined_df.shape}")
            except Exception as e:
                print(f"Error saving combined file for limit={limit}: {e}")
        else:
            print(f"No valid data to combine for limit={limit}.")

if __name__ == "__main__":
    # Define the input directory containing the Excel files
    input_dir = "/Users/junwkim/Documents/BIOMEDIN212/analysis/results/SCLC"  # <-- Replace with your input directory path

    # Define the output directory where combined files will be saved
    output_dir = "/Users/junwkim/Documents/BIOMEDIN212/analysis/results/SCLC"  # <-- Replace with your desired output directory path

    combine_gene_count_files(input_dir, output_dir)
