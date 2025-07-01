# HLA Allele Imputation from TCR Repertoire Data

A command-line tool designed to impute HLA (Human Leukocyte Antigen) alleles directly from T-cell receptor (TCR) repertoire sequencing data. This tool leverages specific models (TRA or TRB clonotypes) to predict HLA allele carriership, providing insights into an individual's immune profile.

## Features

- **Flexible Input**: Accepts TCR repertoire data in CSV format with specified mandatory columns.
- **Gene Naming Support**: Accommodates both 'Adaptive' and 'IMGT' gene naming conventions.
- **Model Selection**: Allows users to choose between 'TRA' and 'TRB' clonotype models for imputation.
- **Configurable Output**:
  - Option to return raw probabilities for allele carriership.
  - Option to normalize probabilities among different alleles of a given HLA locus.
- **Results Export**: Saves imputation results to a specified TSV file.

## Installation

(Please replace this section with actual installation instructions for your project, e.g., pip install, cloning the repo, setting up a virtual environment, dependencies, etc.)

### Example:
```bash
git clone https://github.com/ikmb/TCR2HLA
cd TCR2HLA
pip install -r requirements.txt
cd TCR2HLA/CLI_tools/impute_HLA
```

## Usage

To run the HLA allele imputation tool, execute the script from your command line, providing the necessary arguments.

### Basic Usage Example
```bash
python imput_HLA_CLI.py \
    --input_repertoire /path/to/your/tcr_data.csv \
    --gene_naming Adaptive \
    --model TRB \
    --path2models ../../models/ \
    --path2write_results /path/to/output_results.tsv
```

### Returning Probabilities
To get probabilities instead of hard (0/1) predictions:
```bash
python imput_HLA_CLI.py \
    --input_repertoire /path/to/your/tcr_data.csv \
    --gene_naming Adaptive \
    --model TRB \
    --path2models ../../models/ \
    --return_prob \
    --path2write_results /path/to/output_probabilities.tsv
```

### Normalizing Probabilities
To normalize probabilities (this flag requires --return_prob to be set):
```bash
python imput_HLA_CLI.py \
    --input_repertoire /path/to/your/tcr_data.csv \
    --gene_naming Adaptive \
    --model TRB \
    --path2models ../../models/ \
    --return_prob \
    --normalize_prob \
    --path2write_results /path/to/output_normalized_probabilities.tsv
```

## Arguments

| Argument                | Type   | Required | Choices         | Description                                                                 |
|-------------------------|--------|----------|------------------|-----------------------------------------------------------------------------|
| `--input_repertoire`    | str    | Yes      | -                | Path to the CSV file containing TCR repertoire data (mandatory columns: `v_gene`, `CDR3`, `j_gene`, `sample_name`, `count`). |
| `--gene_naming`         | str    | Yes      | Adaptive, IMGT   | The formatting of V and J gene naming. Must be 'Adaptive' or 'IMGT'.        |
| `--model`               | str    | Yes      | TRA, TRB         | Which model shall be used for imputing HLA alleles. Either 'TRA' or 'TRB'.  |
| `--path2models`         | str    | Yes      | -                | Path to the directory where HLA imputation models are defined.              |
| `--return_prob`         | flag   | No       | -                | If set, the program will return probabilities instead of hard (0/1) predictions. |
| `--normalize_prob`      | flag   | No       | -                | If set, normalizes the probabilities among the different alleles of a given locus (e.g., among HLA-A, HLA-B, etc.). This flag requires --return_prob to be set. |
| `--path2write_results`  | str    | Yes      | -                | Base path to write results a TSV file.                                      |

## Output

The program will generate a TSV file at the specified `--path2write_results`. The content of this file will vary based on the `--return_prob` and `--normalize_prob` flags:

- **Default (no flags)**: Hard (0/1) predictions for allele carriership. 
- **--return_prob**: Raw probability scores for allele carriership.
- **--return_prob --normalize_prob**: Normalized probability scores for allele carriership, where normalization is applied per HLA locus.

## Contributing

Contributions are welcome! Please feel free to open issues or submit pull requests.