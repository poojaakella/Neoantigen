from mhcflurry.class1_affinity_predictor import Class1AffinityPredictor
import pandas as pd
from tqdm import tqdm
from itertools import product

try:
    import tensorflow as tf
    gpus = tf.config.list_physical_devices("GPU")
    if gpus:
        print(f"GPU detected: {len(gpus)} device(s) available.")
    else:
        print("No GPU detected.")
except ImportError:
    print("TensorFlow not installed. GPU check skipped.")

# Load peptides
df = pd.read_csv("") # Update: common neopeptide
peptides = df["subpeptide"].tolist()
alleles = [
    "HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-A*11:01",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01",
    "HLA-C*07:01", "HLA-C*04:01"
] # Set alleles

# Prepare input combinations (cartesian product)
input_df = pd.DataFrame(
    list(product(peptides, alleles)),
    columns=["peptide", "allele"]
)

print("Loading MHCflurry model...")
predictor = Class1AffinityPredictor.load()

print("Predicting binding affinities...")
predictions = predictor.predict_to_dataframe(
    peptides=input_df["peptide"].tolist(),
    alleles=input_df["allele"].tolist()
)

# Rename column for clarity
predictions = predictions.rename(columns={"prediction": "affinity"})

# Merge metadata
predictions["sample_count"] = predictions["peptide"].map(df.set_index("subpeptide")["sample_count"])
predictions["samples"] = predictions["peptide"].map(df.set_index("subpeptide")["samples"])

# Save all predictions
predictions.to_csv("", index=False) # Set file name

strong = predictions[predictions["affinity"] <= 500].copy()
strong.to_csv("common_strong_binders_allClass1HLA.csv", index=False)

print("Done. Results saved to:")
