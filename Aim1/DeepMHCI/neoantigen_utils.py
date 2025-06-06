#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 20 22:48:16 2025
GPT was used for debugging/organization/library suggestions
@author: junwkim
"""
import re
import pandas as pd
import logging
from collections import defaultdict

def read_variants_from_excel(excel_file: str, sheet_name=0) -> list:
    """Reads an Excel file with MAF-like columns and returns a list of dictionaries that match the structure of MAF_VARIANTS."""
    df = pd.read_excel(excel_file, sheet_name=sheet_name, dtype=str)
    variants_list = []

    for _, row in df.iterrows():
        variant_dict = {
            "Hugo_Symbol": row.get("Hugo_Symbol", ""),
            "Chromosome": row.get("Chromosome", ""),
            "Start_Position": convert_to_int(row.get("Start_Position", "")),
            "End_Position": convert_to_int(row.get("End_Position", "")),
            "Strand": row.get("Strand", ""),
            "Variant_Classification": row.get("Variant_Classification", ""),
            "Variant_Type": row.get("Variant_Type", ""),
            "Reference_Allele": row.get("Reference_Allele", ""),
            "Tumor_Seq_Allele2": row.get("Tumor_Seq_Allele2", ""),
            "Transcript_ID": row.get("Transcript_ID", ""),
            "Protein_position": convert_to_int(row.get("Protein_position", "")),
            "HGVSp_Short": row.get("HGVSp_Short", ""),
            "Tumor_Sample_Barcode": row.get("Tumor_Sample_Barcode", "")  # <-- ADDED
        }
        variants_list.append(variant_dict)

    return variants_list

def convert_to_int(value: str):
    """Helper function to convert string to integer if possible."""
    try:
        return int(float(value))
    except (ValueError, TypeError):
        return None

def parse_hgvsp_short(hgvsp_short):
    """Parse an HGVS protein string like 'p.E90G', returning (ref_aa, position, mut_aa)."""
    if hgvsp_short.startswith("p."):
        hgvsp_short = hgvsp_short[2:]
    match = re.match(r"([A-Za-z*])(\d+)([A-Za-z*?=])", hgvsp_short)
    if not match:
        return None, None, None
    ref_aa, pos_str, mut_aa = match.groups()
    pos = int(pos_str)
    return ref_aa, pos, mut_aa

def generate_peptides(ref_aa, mut_aa, prot_seq, prot_pos, flank_min=8, flank_max=14):
    """Generate 8–14mers around the mutated residue (prot_pos, 1-based)."""
    results = []
    seq_length = len(prot_seq)
    zero_idx = prot_pos - 1

    for length in range(flank_min, flank_max + 1):
        half_len = length // 2
        start = zero_idx - half_len
        end = start + length

        # Clamp boundaries
        if start < 0:
            start = 0
        if end > seq_length:
            end = seq_length
        start = end - length if start < 0 else start

        snippet = list(prot_seq[start:end])
        loc_in_snippet = zero_idx - start

        if 0 <= loc_in_snippet < len(snippet):
            snippet[loc_in_snippet] = mut_aa
            mutant_peptide = "".join(snippet)
            results.append({
                "peptide": mutant_peptide,
                "start_in_protein": start + 1,
                "end_in_protein": end,
                "center_index": loc_in_snippet
            })

    return results

def convert_allele_to_HLA_B(allele_str):
    """Given a string like "30:01:00", return "HLA-B30:01"."""
    allele_str = allele_str.strip("'*")
    parts = allele_str.split(":")
    if len(parts) >= 2:
        trimmed_allele = parts[0] + ":" + parts[1]
    else:
        trimmed_allele = parts[0]
    return f"HLA-B{trimmed_allele}"

def convert_allele_to_HLA_A(allele_str):
    """Given a string like "30:01:00", return "HLA-B30:01"."""
    allele_str = allele_str.strip("'*")
    parts = allele_str.split(":")
    if len(parts) >= 2:
        trimmed_allele = parts[0] + ":" + parts[1]
    else:
        trimmed_allele = parts[0]
    return f"HLA-A{trimmed_allele}"


def read_variants_from_txt(path, sep="\t", skip_rows=2):
    """Reads mutations from a text file and returns a list of dictionaries."""
    try:
        df = pd.read_csv(path, sep=sep, dtype=str, skiprows=skip_rows).fillna("")
        logging.info(f"Successfully read variants from {path}")
        return df.to_dict(orient="records")
    except Exception as e:
        logging.error(f"Error reading variants from {path}: {e}")
        return []

def save_neoantigen_candidates(candidates, output_filepath):
    """Prints the number of neoantigen candidates and saves them to Excel."""
    print("Neoantigen candidates with DeepMHCI predictions:")
    print(f"Total candidates: {len(candidates)}\n")
    df_candidates = pd.DataFrame(candidates)
    print("Sample of candidates:")
    print(df_candidates.head)
    try:
        df_candidates.to_excel(output_filepath, index=False)
        print(f"\nResults successfully saved to '{output_filepath}'.")
    except Exception as e:
        print(f"\nAn error occurred while saving to Excel: {e}")
