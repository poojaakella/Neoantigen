#Author: Yugendran Rajaendran
import os
import re
import pandas as pd
import pyensembl
import subprocess
from tqdm import tqdm
from collections import defaultdict
import logging

try:
    import tensorflow as tf
    gpus = tf.config.list_physical_devices("GPU")
    if gpus:
        print(f"GPU detected: {len(gpus)} device(s) available.")
    else:
        print("No GPU detected.")
except ImportError:
    print("TensorFlow not installed. GPU check skipped.")

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

BASE_DIR = ""  # Update 
OUTPUT_DIR = "" # Update

def convert_to_int(value):
    try:
        return int(float(value))
    except (ValueError, TypeError):
        return None

def parse_hgvsp_short(hgvsp_short):
    if hgvsp_short.startswith("p."):
        hgvsp_short = hgvsp_short[2:]
    match = re.match(r"([A-Za-z*])(\d+)([A-Za-z*?=])", hgvsp_short)
    if not match:
        return None, None, None
    ref_aa, pos_str, mut_aa = match.groups()
    return ref_aa, int(pos_str), mut_aa

def convert_allele_to_hla_B(allele):
    allele = str(allele).replace("'", "").strip()
    parts = allele.split(":")
    if len(parts) >= 2:
        return f"HLA-B{parts[0]}:{parts[1]}"
    return f"HLA-B{allele}"

def get_protein_sequence_cached(ensembl, transcript_id, cache):
    if transcript_id not in cache:
        try:
            cache[transcript_id] = ensembl.transcript_by_id(transcript_id).protein_sequence
        except Exception:
            cache[transcript_id] = None
    return cache[transcript_id]

def generate_27mer(ref_aa, mut_aa, prot_seq, prot_pos):
    idx = prot_pos - 1
    flank = 13
    start = max(0, idx - flank)
    end = min(len(prot_seq), idx + flank + 1)
    seq = list(prot_seq[start:end])
    if 0 <= idx - start < len(seq):
        seq[idx - start] = mut_aa
    return "".join(seq)

def generate_15to27mer(ref_aa, mut_aa, prot_seq, prot_pos, flank_min=15, flank_max=27):
    results = []
    seq_length = len(prot_seq)
    zero_idx = prot_pos - 1
    for length in range(flank_min, flank_max + 1):
        half_len = length // 2
        start = zero_idx - half_len
        end = start + length
        if start < 0:
            start = 0
        if end > seq_length:
            end = seq_length
        start = end - length if start < 0 else start
        snippet = list(prot_seq[start:end])
        loc_in_snippet = zero_idx - start
        if 0 <= loc_in_snippet < len(snippet):
            snippet[loc_in_snippet] = mut_aa
            results.append("".join(snippet))
    return results

def generate_subpeptides(peptide, min_len=8, max_len=14):
    seen = set()
    subs = []
    for length in range(min_len, max_len + 1):
        for start in range(0, len(peptide) - length + 1):
            sub = peptide[start:start + length]
            if sub not in seen:
                seen.add(sub)
                subs.append(sub)
    return subs


def save_neoantigen_candidates(candidates, output_file):
    df = pd.DataFrame(candidates)
    df.to_csv(output_file, index=False)
    print(f"Saved {len(df)} neoantigen candidates to: {output_file}")

def main():
    ensembl = pyensembl.EnsemblRelease(release=75, species='homo_sapiens')
    ensembl.download()
    ensembl.index()
    print("Ensembl loaded.")

    print("Loading CCLE sample data...")
    ccle_df = pd.read_csv(f"{BASE_DIR}/ccle_broad_2019_clinical_data.tsv", sep="\t", dtype=str)
    ccle_sample_ids = set(ccle_df["Sample ID"].unique())

    print("Loading HLA data...")
    hla_df = pd.read_csv(f"{BASE_DIR}/20240226_TCLP_HLA_data.csv", dtype=str).fillna("")
    #hla_df_B = hla_df[hla_df["type"] == "B"]
    hla_df_class1 = hla_df[hla_df["type"].isin(["A", "B", "C"])]

    valid_hla_rows = hla_df_class1[hla_df_class1["name"].apply(lambda x: f"{x}_LUNG" in ccle_sample_ids)].copy()
    valid_hla_rows[["name", "type", "allele1", "allele2"]].to_csv(os.path.join(BASE_DIR, "samples_with_hla_class1.csv"), index=False)
    print(f"Found {len(valid_hla_rows)} samples and HLA Class 1 data and in CCLE.")

    hla_df_class1["sample_id"] = hla_df_class1["name"] + "_LUNG"
    hla_df_class1["in_ccle"] = hla_df_class1["sample_id"].isin(ccle_sample_ids)

    summary = hla_df_class1.groupby("type")["in_ccle"].sum()
    print("HLA Class I matches by type:\n", summary)

    hla_map = defaultdict(list)
    for _, row in valid_hla_rows.iterrows():
        sid = row["name"] + "_LUNG"
        hla_map[sid].extend([
            convert_allele_to_hla_B(row["allele1"]),
            convert_allele_to_hla_B(row["allele2"])
        ])

    print("Loading mutation data...")
    maf = pd.read_csv(f"{BASE_DIR}/data_mutations.txt", sep="\t", skiprows=2, dtype=str,
                      engine="python").fillna("")
    mutations_by_sample = defaultdict(list)
    for _, row in maf.iterrows():
        sid = row.get("Tumor_Sample_Barcode")
        if sid in ccle_sample_ids and sid in hla_map:
            mutations_by_sample[sid].append(row)

    print(f"Found {len(mutations_by_sample)} samples with mutations and HLA data.")

    cache = {}

    for sample_id in tqdm(sorted(mutations_by_sample.keys()), desc="Processing samples"):
        output_path = os.path.join(OUTPUT_DIR, f"{sample_id}_neoantigens.csv")
        if os.path.exists(output_path):
            continue

        seen_sites = set()
        peptides = set()
        meta_map = {}

        for row in mutations_by_sample[sample_id]:
            if row.get("Variant_Classification") != "Missense_Mutation":
                continue

            ref, pos, mut = parse_hgvsp_short(row.get("HGVSp_Short", ""))
            if not ref or not mut:
                continue

            tid = row["Transcript_ID"]
            pos_int = convert_to_int(row.get("Protein_position"))
            if not tid or not pos_int:
                continue

            if (tid, pos_int) in seen_sites:
                continue
            seen_sites.add((tid, pos_int))

            prot_seq = get_protein_sequence_cached(ensembl, tid, cache)
            if not prot_seq or pos_int < 1 or pos_int > len(prot_seq):
                continue

            if prot_seq[pos_int - 1] != ref:
                continue

            long_peptides = generate_15to27mer(ref, mut, prot_seq, pos_int)
            for peptide in long_peptides:
                subpeps = generate_subpeptides(peptide)
                for sp in subpeps:
                    if sp in peptides:
                        continue
                    peptides.add(sp)
                    meta_map[sp] = {
                        "gene_symbol": row["Hugo_Symbol"],
                        "transcript_id": tid,
                        "peptide_15to27": peptide,
                        "sample_id": sample_id
                    }

        if not peptides:
            continue

        neo_output_path = os.path.join(OUTPUT_DIR, f"{sample_id}_neoantigens.pep")
        with open(neo_output_path, "w") as f:
            for pep in peptides:
                f.write(f"{pep}\n")

if __name__ == "__main__":
    main()
