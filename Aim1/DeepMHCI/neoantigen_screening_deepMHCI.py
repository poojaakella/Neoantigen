#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 20 22:48:58 2025

@author: junwkim
"""
import logging
import pyensembl
import pandas as pd
import os
from collections import defaultdict
import argparse
from neoantigen_utils import (
    read_variants_from_txt,
    parse_hgvsp_short,
    convert_allele_to_HLA_B,
    convert_allele_to_HLA_A,
    generate_peptides,
    save_neoantigen_candidates
)
from predict_re import deepmhci_predict

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def main(cancer, threshold, hla):
    # Initialize PyEnsembl for GRCh37 (Ensembl release 75 recommended).
    ensembl = pyensembl.EnsemblRelease(release=75, species='homo_sapiens')
    ensembl.download()
    ensembl.index()

    # ------------------------------------------------------------------------
    # Filtering Steps:
    # 1. Get the SCLC sample IDs from ccle_broad_2019_clinical_data.tsv
    # 2. Filter the mutation for SCLC samples based on the sample IDs
    # 3. Get the HLA-B corresponding to the SCLC samples
    # ------------------------------------------------------------------------

    # 1. Read SCLC sample IDs
    if cancer == "SCLC":
        ccle_clinical_path = "../data/SCLC_ccle_broad_2019_clinical_data.tsv"
    elif cancer == "NSCLC":
        ccle_clinical_path = "../data/NSCLC_ccle_broad_2019_clinical_data.tsv"
    elif cancer == "pan":
        ccle_clinical_path = "../data/Pan_ccle_broad_2019_clinical_data.tsv"
    ccle_df = pd.read_csv(ccle_clinical_path, sep="\t", dtype=str)
    ccle_sample_ids = set(ccle_df["Sample ID"].unique())
    #ccle_sample_ids = set(["NCIH1048_LUNG"])  # For testing purposes

    # 2. Read mutations from data_mutations.txt
    mutations_path = "data_mutations.txt"
    maf_variants = read_variants_from_txt(mutations_path, sep="\t")

    # Filter variants for SCLC samples
    maf_variants_filtered = [
        v for v in maf_variants if v.get("Tumor_Sample_Barcode") in ccle_sample_ids
    ]

    # 3. Read HLA data
    hla_data_path = "20240226_TCLP_HLA_data.csv"
    hla_df = pd.read_csv(hla_data_path, dtype=str)
    hla_df_B = hla_df[hla_df["type"] == hla]

    # Build a map: sample_name -> [allele1, allele2]
    hla_map = {}
    for _, row in hla_df_B.iterrows():
        sample_name = row["name"]
        allele1 = row["allele1"]
        allele2 = row["allele2"]
        if hla == "A":
            allele1_formatted = convert_allele_to_HLA_A(allele1)
            allele2_formatted = convert_allele_to_HLA_A(allele2)
        elif hla == "B":
            allele1_formatted = convert_allele_to_HLA_B(allele1)
            allele2_formatted = convert_allele_to_HLA_B(allele2)

        hla_map[sample_name] = [allele1_formatted, allele2_formatted]

    # 3.2) Filter for nonsynonymous variants
    nonsynonymous = []
    for row in maf_variants_filtered:
        if row.get("Variant_Classification") != "Missense_Mutation":
            continue

        transcript_id = row["Transcript_ID"]
        pos_str = row.get("Protein_position", "")
        try:
            protein_pos = int(pos_str)
        except ValueError:
            protein_pos = None

        hgvsp_short = row.get("HGVSp_Short", "")
        ref_aa, pos, mut_aa = parse_hgvsp_short(hgvsp_short)
        sample_id = row.get("Tumor_Sample_Barcode", "")

        # Skip incomplete or synonymous
        if not ref_aa or not mut_aa or not sample_id:
            continue

        nonsynonymous.append({
            "gene_symbol": row["Hugo_Symbol"],
            "transcript_id": transcript_id,
            "protein_pos": protein_pos,
            "ref_aa": ref_aa,
            "mut_aa": mut_aa,
            "sample_id": sample_id
        })

    # 3.3) Retrieve protein sequences and generate peptides
    final_neoantigen_candidates = []
    peptide_to_samples = defaultdict(set)
    peptide_variant_sample_map = defaultdict(list)  # Key: peptide, Value: list of (variant, sample_id, allele)
#    peptide_allele_pairs = []
    peptide_allele_pairs = set()
    for var in nonsynonymous:
        t_id = var["transcript_id"]
        sample_id = var["sample_id"].split('_')[0]

        # Skip if there's no HLA type B data for this sample
        if sample_id not in hla_map:
            continue
        patient_hla = hla_map[sample_id]
        try:
            transcript_obj = ensembl.transcript_by_id(t_id)
            protein_seq = transcript_obj.protein_sequence
        except KeyError:
            logging.warning(f"Transcript ID {t_id} not found in Ensembl.")
            continue

        if not protein_seq:
            logging.warning(f"No protein sequence found for transcript ID {t_id}.")
            continue

        # Check reference amino acid for consistency
        idx = var["protein_pos"]
        if not isinstance(idx, int) or idx < 1 or idx > len(protein_seq):
            logging.warning(f"Invalid protein position {idx} for transcript ID {t_id}.")
            continue

        #logging.info(f"Processing sample id: {sample_id}")
        zero_idx = idx - 1
        ref_aa_ensembl = protein_seq[zero_idx]
        if ref_aa_ensembl != var["ref_aa"]:
            logging.warning(
                f"Reference amino acid mismatch at position {idx} for transcript ID {t_id}. "
                f"MAF: {var['ref_aa']}, Ensembl: {ref_aa_ensembl}"
            )
            continue

        # 3.4) Generate 8–14mers
        short_peptides = generate_peptides(
            ref_aa=var["ref_aa"],
            mut_aa=var["mut_aa"],
            prot_seq=protein_seq,
            prot_pos=var["protein_pos"],
            flank_min=8,
            flank_max=14
        )

        for sp in short_peptides:
            peptide = sp["peptide"]
            peptide_to_samples[peptide].add(sample_id)
            # Add all allele pairs for this peptide along with variant and sample info
            for allele in patient_hla:
               # peptide_allele_pairs.append((peptide, allele))
                peptide_allele_pairs.add((peptide, allele))
                peptide_variant_sample_map[peptide].append((var, sample_id, allele))

    # Filter peptides present in at least two samples
    filtered_peptides = {pep for pep, samples in peptide_to_samples.items() if len(samples) >= 0}
    filtered_peptide_allele_pairs = [
        (pep, allele) for (pep, allele) in peptide_allele_pairs if pep in filtered_peptides
    ]

    if not filtered_peptide_allele_pairs:
        logging.info("No peptides found present in at least two samples.")
        return

    # 3.5) Run DeepMHCI predictions
    data_yaml_path = "config/data.yaml"
    model_yaml_path = "config/model.yaml"

    logging.info("Running DeepMHCI predictions...")
    predictions = deepmhci_predict(
        peptide_allele_pairs=filtered_peptide_allele_pairs,
        data_cnf_path=data_yaml_path,
        model_cnf_path=model_yaml_path,
        start_id=0,
        num_models=1
    )

    # 3.6) Filter predicted binders with score thresholds
    strong_binders = [pred for pred in predictions if pred["score"] >= threshold]
    logging.info(f"Strong binders found: {len(strong_binders)}")

    # 3.7) Collect results
    
    for sb in strong_binders:
#here
        peptide = sb["peptide"]
        allele = sb["allele"]
        score = sb["score"] 
        # Retrieve all (variant, sample_id) pairs associated with this peptide and allele
        associated_variants_samples = [
            (var, sample_id)
            for var, sample_id, a in peptide_variant_sample_map[peptide]
            if a == allele
        ]

        for var, sample_id in associated_variants_samples:
            final_neoantigen_candidates.append({
                "gene_symbol": var["gene_symbol"],
                "transcript_id": var["transcript_id"],
                "peptide_8to14": peptide,
                "allele": allele,
                "sample_id": sample_id,
                "score": score
            })

    # 3.8) Print/save final results
    output_excel_path = f"../results/neoantigen_candidates_deepmhci_{cancer}_{threshold}_{hla}.xlsx"
    save_neoantigen_candidates(final_neoantigen_candidates, output_excel_path)

    df_candidates = pd.DataFrame(final_neoantigen_candidates)
    if not df_candidates.empty:
        # Get all unique genes from the nonsynonymous variants
        all_genes = sorted(set(v['gene_symbol'] for v in nonsynonymous))
        # Get all unique samples
        all_samples = sorted(set(df_candidates['sample_id']))
        # Group by sample_id and gene_symbol, then count the number of neoantigens
        gene_counts = df_candidates.groupby(['sample_id', 'gene_symbol']).size().unstack(fill_value=0)
        # Ensure all genes are represented as columns
        gene_counts = gene_counts.reindex(columns=all_genes, fill_value=0)
        # Optionally, ensure all samples are represented as rows
        gene_counts = gene_counts.reindex(index=all_samples, fill_value=0)
        # Save the gene count matrix to an Excel file
        gene_counts_path = f"../results/gene_counts_per_sample_deepmhci_{cancer}_{threshold}_{hla}.csv"
        gene_counts.to_csv(gene_counts_path)
        logging.info(f"Gene counts per sample successfully saved to '{gene_counts_path}'.")
    else:
        logging.info("No neoantigen candidates found; gene count vectors were not generated.")

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process neoantigen candidates with a specified model."
    )
    parser.add_argument(
        "--cancer",
        type=str,
        help="cancer type, e.g. NSCLC, SCLC, etc"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.2,
        help="threhold affinity"
    )
    parser.add_argument(
        "--hla",
        type=str,
        default="A",
        help="HLA type"
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    main(
        cancer = args.cancer,
        threshold = args.threshold,
        hla = args.hla
    )

