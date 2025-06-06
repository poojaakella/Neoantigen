#Author: Yugendran Rajaendran
import os
from collections import defaultdict
from tqdm import tqdm
import pandas as pd

try:
    import tensorflow as tf
    gpus = tf.config.list_physical_devices("GPU")
    if gpus:
        print(f"GPU detected: {len(gpus)} device(s) available.")
    else:
        print("⚠No GPU detected. MHCflurry will run on CPU.")
except ImportError:
    print("ℹTensorFlow not installed. GPU check skipped.")

PEP_FOLDER = ""  # Folder containing peptides
OUTPUT_CSV = ""  # Set output filename
MIN_SAMPLES = 3  # Minimum number of samples a peptide must appear in


def find_common_subpeptides():
    subpeptide_to_samples = defaultdict(set)
    pep_files = [f for f in os.listdir(PEP_FOLDER) if f.endswith(".pep")]
    print(f"Found {len(pep_files)} .pep files")

    for fname in tqdm(pep_files, desc=" Parsing peptide files"):
        #sample_id = os.path.splitext(fname)[0]
        sample_id = fname.replace("_neoantigens.pep", "")
        with open(os.path.join(PEP_FOLDER, fname)) as f:
            for line in f:
                peptide = line.strip()
                if peptide:
                    subpeptide_to_samples[peptide].add(sample_id)

    # Filter and collect common subpeptides
    common_peptides = [
        {
            "subpeptide": pep,
            "sample_count": len(samples),
            "samples": ";".join(sorted(samples))
        }
        for pep, samples in subpeptide_to_samples.items()
        if len(samples) >= MIN_SAMPLES
    ]

    # Save to CSV
    df = pd.DataFrame(common_peptides)
    df.sort_values(by="sample_count", ascending=False, inplace=True)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Saved {len(df)} common subpeptides to {OUTPUT_CSV}")

if __name__ == "__main__":
    find_common_subpeptides()
