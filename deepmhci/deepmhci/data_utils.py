import numpy as np
import torch


ACIDS = '0-ACDEFGHIKLMNPQRSTVWY'

        
def get_mhc_name_seq(mhc_name_seq_file):
    mhc_name_seq = {}
    with open(mhc_name_seq_file) as fp:
        for line in fp:
            mhc_name, mhc_seq = line.split()
            mhc_name_seq[mhc_name] = mhc_seq
    return mhc_name_seq


def get_data(data_file, mhc_name_seq_file):
    mhc_name_seq = get_mhc_name_seq(mhc_name_seq_file)
    data_list = []
    with open(data_file) as fp:
        for line in fp:
            peptide_seq, score, mhc_name = line.split()
            if len(peptide_seq) >= 8:
                data_list.append((mhc_name, peptide_seq, mhc_name_seq[mhc_name], float(score)))
    return data_list


def get_epitopes_data(data_file, mhc_name_seq_file):
    mhc_name_seq = get_mhc_name_seq(mhc_name_seq_file)
    epitopes_data = []
    with open(data_file) as fp:
        for line in fp:
            epitope, target, mhc_name = line.split()
            data_list = [(mhc_name, cur_:=epitope[i: i + len(target)], mhc_name_seq[mhc_name], float(cur_ == target))
                         for i in range(len(epitope) - len(target) + 1)]
            epitopes_data.append((mhc_name, data_list, len(target)))
    return epitopes_data


def get_raw_data(data_file, mhc_name, mhc_seq, peptide_pad=3):
    data_list, peptide_list = [] , []
    with open(data_file) as fp:
        for k, line in enumerate(fp):
            peptide_seq = line.strip()
            # pad_peptide_seq = ACIDS[0] * peptide_pad + peptide_seq + ACIDS[0]*(15-len(peptide_seq)) + ACIDS[0] * peptide_pad
            pad_peptide_seq = peptide_seq + ACIDS[0]*(15-len(peptide_seq))
            data_list += [(mhc_name, pad_peptide_seq, mhc_seq, 0.0)]
            peptide_list.append(peptide_seq)
    return peptide_list, data_list

def get_data_inference(txt_file, mhc_name_seq):
    """
    Reads a text file with lines of the form:
    PEPTIDE_SEQUENCE ALLELE

    Returns a list of [allele, peptide_seq, mhc_seq, dummy_score=0.0].
    """
    data_list = []
    with open(txt_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                # Expecting at least: PEPTIDE ALLELE
                raise ValueError(f"Invalid line (need at least 2 items): {line}")
            peptide_seq, allele = parts[0], parts[1]
            if allele not in mhc_name_seq:
                raise ValueError(f"Unknown allele: {allele}")
            # Assigning a dummy score of 0.0 as a placeholder
            data_list.append([allele, peptide_seq, mhc_name_seq[allele], 0.0])
    return data_list
