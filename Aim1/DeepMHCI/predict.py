#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 15 11:03:21 2025

@author: junwkim
"""

import click
from pathlib import Path
from ruamel.yaml import YAML
import logging
import numpy as np
from torch.utils.data import DataLoader
import csv  # Added for writing CSV/TSV files

from deepmhci.datasets import MHCIDataset
from deepmhci.data_utils import get_mhc_name_seq, get_data_inference
from deepmhci.models import Model
from deepmhci.networks import DeepMHCI  # Adjust import as needed

# Configure logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


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
def predict(data_cnf, model_cnf, start_id, num_models, output_file):
    """
    Perform prediction using a model ensemble and save the results to a file.
    """
    yaml = YAML(typ='safe')

    # Load configurations
    data_config = yaml.load(Path(data_cnf))
    model_config = yaml.load(Path(model_cnf))

    model_name = model_config['name']
    logger.info(f'Model Name: {model_name}')

    model_path = Path(model_config['path']) / f'{model_name}'

    # Load MHC allele → sequence map
    mhc_name_seq = get_mhc_name_seq(data_config['mhc_seq'])

    # Read inference data
    inference_data_list = get_data_inference(data_config['pred_file'], mhc_name_seq)

    # Create PyTorch dataset
    inference_dataset = MHCIDataset(
        inference_data_list,
        **model_config['padding']
    )

    # Initialize ensemble
    ensemble_scores = []
    for model_id in range(start_id, start_id + num_models):
        new_stem = f'{model_path.stem}-{model_id}'
        new_name = f'{new_stem}{model_path.suffix}' 
        model_file = model_path.with_name(new_name)
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

    # Output predictions to logger
    logger.info("Predictions for inference dataset:")
    for i, row in enumerate(inference_data_list):
        allele, peptide_seq, _, _ = row
        logger.info(f"{peptide_seq}\t{allele}\tScore={final_scores[i]:.4f}")

    # Save predictions to the specified file
    try:
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')  # Using TSV format
            writer.writerow(['Peptide_Sequence', 'Allele', 'Score'])  # Header
            for i, row in enumerate(inference_data_list):
                allele, peptide_seq, _, _ = row
                writer.writerow([peptide_seq, allele, f"{final_scores[i]:.4f}"])
        logger.info(f"Predictions have been saved to {output_file}")
    except IOError as e:
        logger.error(f"Failed to write predictions to {output_file}: {e}")


if __name__ == '__main__':
    predict()
