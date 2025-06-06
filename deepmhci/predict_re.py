#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 20 17:46:28 2025
GPT was used for debugging/organization/library suggestions
@author: junwkim
"""
import click
from pathlib import Path
from ruamel.yaml import YAML
import logging
import numpy as np
from torch.utils.data import DataLoader
import csv  # For writing CSV/TSV files

from deepmhci.datasets import MHCIDataset
from deepmhci.data_utils import get_mhc_name_seq, get_data_inference
from deepmhci.models import Model
from deepmhci.networks import DeepMHCI 

# Configure logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

def deepmhci_predict(peptide_allele_pairs, data_cnf_path, model_cnf_path, start_id=0, num_models=20):
    """
    Perform prediction using DeepMHCI on a list of (peptide, allele) tuples.

    Parameters:
    - peptide_allele_pairs: List of tuples [(peptide1, allele1), (peptide2, allele2), ...]
    - data_cnf_path: Path to data configuration YAML.
    - model_cnf_path: Path to model configuration YAML.
    - start_id: Starting ID for model ensemble.
    - num_models: Number of models in ensemble.

    Returns:
    - List of dictionaries containing 'peptide', 'allele', and 'score'.
    """
    yaml = YAML(typ='safe')

    # Load configurations
    data_config = yaml.load(Path(data_cnf_path))
    model_config = yaml.load(Path(model_cnf_path))

    model_name = model_config['name']
    logger.info(f'Model Name: {model_name}')

    model_path = Path(model_config['path']) / f'{model_name}'

    # Load MHC allele → sequence map
    mhc_name_seq = get_mhc_name_seq(data_config['mhc_seq'])

    # Prepare inference data
    inference_data_list = []
    valid_pairs = []
    for peptide, allele in peptide_allele_pairs:
    # If mhc allele sequence is not found, just skip
        try: 
            inference_data_list.append((allele, peptide, mhc_name_seq[allele], 0.0))  # Adjust fields if necessary
            valid_pairs.append((peptide, allele))
        except KeyError:
            print(f"Allele '{allele}' not found in mhc_name_seq. Skipping.")
            continue
    print(inference_data_list[0:5])
    # Create PyTorch dataset
    inference_dataset = MHCIDataset(
        inference_data_list,
        **model_config['padding']
    )

    # Initialize ensemble
    ensemble_scores = []
    for model_id in range(start_id, start_id + num_models):
        model_file = model_path.parent / f'{model_path.stem}-{model_id}{model_path.suffix}'
        logger.info(f'Loading model {model_id} from {model_file}')
        model_instance = Model(DeepMHCI, model_path=model_file, **model_config['model'])

        # Create DataLoader
        data_loader = DataLoader(
            inference_dataset,
            batch_size=model_config['test']['batch_size']
        )

        # Predict scores
        pred_scores = model_instance.predict(data_loader)
        ensemble_scores.append(pred_scores)

    # Average ensemble scores
    final_scores = np.mean(ensemble_scores, axis=0)

    # Prepare results
    results = []
#here
    for i, (peptide, allele) in enumerate(valid_pairs):
        score = final_scores[i]
        results.append({
            "peptide": peptide,
            "allele": allele,
            "score": float(score)
        })

    return results

# If you still want to keep the CLI functionality, retain the Click command
@click.command()
@click.option(
    '-d', '--data-cnf',
    type=click.Path(exists=True),
    required=True,
    help='Path to data configuration YAML.'
)
@click.option(
    '-m', '--model-cnf',
    type=click.Path(exists=True),
    required=True,
    help='Path to model configuration YAML.'
)
@click.option(
    '-s', '--start-id',
    default=0,
    type=int,
    show_default=True,
    help='Starting ID for model ensemble.'
)
@click.option(
    '-n', '--num_models',
    default=20,
    type=int,
    show_default=True,
    help='Number of models in ensemble.'
)
@click.option(
    '-o', '--output-file',
    type=click.Path(writable=True),
    default='predictions.tsv',
    show_default=True,
    help='Path to the output file where predictions will be saved.'
)
def predict_cli(data_cnf, model_cnf, start_id, num_models, output_file):
    """
    Perform prediction using a model ensemble and save the results to a file.
    """
    # Example usage: Load peptides and alleles from a file or other sources
    # Here, it's assumed that 'data/my_peptides_for_inference.txt' exists
    # You can adjust as needed
    peptide_allele_pairs = []
    with open(data_cnf, 'r') as f:
        for line in f:
            peptide, allele = line.strip().split()
            peptide_allele_pairs.append((peptide, allele))

    predictions = deepmhci_predict(peptide_allele_pairs, data_cnf, model_cnf, start_id, num_models)

    # Save predictions to the specified file
    try:
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')  # Using TSV format
            writer.writerow(['Peptide_Sequence', 'Allele', 'Score'])  # Header
            for pred in predictions:
                writer.writerow([pred['peptide'], pred['allele'], f"{pred['score']:.4f}"])
        logger.info(f"Predictions have been saved to {output_file}")
    except IOError as e:
        logger.error(f"Failed to write predictions to {output_file}: {e}")

if __name__ == '__main__':
    predict_cli()
