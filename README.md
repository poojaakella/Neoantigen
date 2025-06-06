# Neoantigen


Neoantigen_screening_deepMHCI allows users to screen for potential neoantigens in the Cancer Cell Line Encyclopedia (Broad, 2019) using DeepMHCI. This repository supplies an environment.yml file for Conda-based installation, and a Python script (`neoantigen_screening_deepMHCI.py`) to run the analysis.


## Installation


1. Clone or download this repository.  

2. Create a Conda environment using the provided YAML file:

   ```bash
   conda env create -f environment.yml
   ```

3. Activate the environment:

    ```bash
    conda activate deepmhci
    ```
4. Prepare the dataset (could not be uploaded due to the large size)
   4.1 Download the dataset (https://cbioportal-datahub.s3.amazonaws.com/ccle_broad_2019.tar.gz)
   4.2 Move data_mutations.txt file to the data folder in Aim1 directory

5. Usage:

Once installed and activated, run the screening script in a terminal:

  ```bash
  cd deepmhci \
 
  python neoantigen_screening_deepMHCI.py \

  --cancer <SCLC|NSCLC|pan> \

  --threshold <affinity score between 0.0 and 1.0> \

  --hla <A or B>
  ```
Example usage
  ```bash
  python neoantigen_screening_deepMHCI.py \

  --cancer SCLC \

  --threshold 0.5 \

  --hla A
 ```

### Contributions

Aim 1: Ady, Jun, Yugendran
Aim 2: Pooja, Gloria
Aim 3: Jun
   
