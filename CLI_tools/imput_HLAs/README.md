# HLA Imputation Command Line Tool

This tool is a command-line utility designed to impute Human Leukocyte Antigen (HLA) alleles from T-cell Receptor (TCR) repertoire sequencing data. It processes raw repertoire data, unifies gene naming conventions, and applies a pre-trained model to predict HLA alleles, outputting the results in a TSV format.

---

## How it Works

The tool follows a structured workflow to perform HLA imputation:

### Argument Parsing and Validation

- Uses `argparse` to define and parse command-line arguments, ensuring all mandatory inputs are provided.
- **Input validation checks include:**
  - `gene_naming` must be either `'Adaptive'` or `'IMGT'`.
  - `model` must be either `'TRA'` (TCR Alpha) or `'TRB'` (TCR Beta).
  - Paths to the input repertoire CSV file and the models directory are verified for existence.
  - The output directory for results is created if it does not already exist.

### Data Loading and Basic Statistics

- Loads the TCR repertoire data from the specified CSV file into a pandas DataFrame.
- Mandatory columns: `v_gene`, `CDR3`, `j_genes`, and `sample_name`.
- Prints basic statistics:
  - Total number of entries
  - Unique samples
  - Average clonotypes per sample
  - Min/Max clonotypes per sample
- Issues a warning if any sample has fewer than 25,000 clonotypes (may impact imputation recall).

### Gene Naming Unification

- Unifies the V and J gene naming conventions for compatibility with imputation models.
- Conversion logic:
  - If `model` is `'TRA'` and `gene_naming` is `'Adaptive'`, converts from Adaptive to IMGT using `adaptive_to_imgt_builder`.
  - If `model` is `'TRB'` and `gene_naming` is `'IMGT'`, converts from IMGT to Adaptive using `imgt_to_adaptive_builder`.
- Gene conversion mapping files must exist in the `path2models` directory.

### HLA Allele Prediction

- Initializes an `HLAPredictor` instance using `path2models` and the selected `model` type.
- Calls `predict_HLA_for_a_sample_collection` on the processed DataFrame to perform HLA imputation using pre-trained models.

### Results Translation and Writing

- Translates results (likely in dict/object form) to a DataFrame using `translate_to_tsv`.
- Saves the DataFrame to a TSV file specified by `--path2write_results`.

---

## Installation

Ensure the following Python packages are installed:

```bash
pip install pandas numpy
```

The tool also requires custom modules:
- `CLI_tools/query_database/gene_formating_converter.py`
- `CLI_tools/imput_HLAs/impute_HLA_alleles_API.py`

Make sure these are correctly structured relative to your `hla_imputer.py` script or available in your Python path.

---

## Usage

Run the tool from your terminal:

```bash
python hla_imputer.py \
    --input_repertoire /path/to/your/input_data.csv \
    --gene_naming IMGT \
    --model TRB \
    --path2models /path/to/your/models_directory \
    --path2write_results /path/to/your/output_results.tsv
```

### Arguments

- `--input_repertoire` (str, Mandatory):  
  Path to the CSV file containing TCR repertoire data. Must include: `v_gene`, `CDR3`, `j_genes`, and `sample_name`.

- `--gene_naming` (str, Mandatory):  
  Specifies the V and J gene naming format. Accepted values: `'Adaptive'` or `'IMGT'`.

- `--model` (str, Mandatory):  
  The model to use for imputation. Accepted values: `'TRA'` or `'TRB'`.

- `--path2models` (str, Mandatory):  
  Path to the directory with HLA models and gene conversion mappings.

- `--path2write_results` (str, Mandatory):  
  Path (including filename) to write the output TSV with imputed HLA alleles.

---

## Input File Format (`input_repertoire.csv`)

The CSV input must include these columns:

| Column Name  | Description                                           | Example Value         |
|--------------|-------------------------------------------------------|------------------------|
| `v_gene`     | The V gene segment of the TCR clonotype              | `TRBV1-1` or `TRBV1`   |
| `CDR3`       | The CDR3 amino acid sequence                         | `CASSLGQYEQYF`         |
| `j_genes`    | The J gene segment of the TCR clonotype              | `TRBJ1-1` or `TRBJ01`  |
| `sample_name`| Unique identifier for the sample                     | `patient_A_day_0`      |

---

## Output File Format (`output_results.tsv`)

A TSV file containing the imputed HLA alleles for each sample. Columns may vary based on the `translate_to_tsv` function.

### Example:

```
sample_name	imputed_hla_alleles
sample_A	HLA-A*01:01,HLA-A*02:01,HLA-B*07:02
sample_B	HLA-A*03:01,HLA-B*08:01,HLA-C*07:01
```
---
