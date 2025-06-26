# TCR-HLA Query Tool

This document provides comprehensive information on the `tcr_query_tool.py` command-line utility, designed for querying TCR-HLA (T-cell Receptor - Human Leukocyte Antigen) databases. The tool supports various query types, allowing users to search for TCRs associated with specific HLA alleles, query individual TCR sequences, or process a table of TCRs.

---

## Table of Contents
- [Purpose](#1-purpose)
- [Installation](#2-installation)
- [Usage](#3-usage)
  - [General Syntax](#general-syntax)
  - [Querying an HLA Allele](#querying-an-hla-allele)
  - [Querying a Single TCR](#querying-a-single-tcr)
  - [Querying a Table of TCRs](#querying-a-table-of-tcrs)
- [Arguments](#4-arguments)
- [Input and Output Formats](#5-input-and-output-formats)
- [Modes of Operation](#6-modes-of-operation)

---

## 1. Purpose

The `tcr_query_tool.py` facilitates the interaction with TCR-HLA association databases. It streamlines the process of extracting relevant TCR data based on user-defined criteria, supporting three primary querying methods: by HLA allele name, by a single TCR sequence, or by a batch of TCR sequences provided in a table. The results are outputted to a specified TSV (Tab Separated Values) file.

---

## 2. Installation

This tool is a standalone Python script and requires Python 3.x. No external libraries beyond standard Python modules (`argparse`, `os`, `sys`) are strictly necessary for the command-line parsing itself. However, your internal "wiring" for database interaction may require additional libraries (e.g., `pandas` for table handling, database connectors).

### To use:
1. Ensure you have Python 3 installed.
2. Save the provided script as `tcr_query_tool.py`.
3. *(Optional)* Install any additional Python packages required by your internal database interaction logic.

---

## 3. Usage

The tool operates based on **mutually exclusive query types**. You must specify exactly one primary query argument: `--input_allele`, `--input_tcr`, or `--query_tcr_tables`.

### General Syntax
```bash
python tcr_query_tool.py --path2database <path_to_db> --database <db_type> --path2write_results <output_path_prefix> [query_type_arguments]
```

---

### Querying an HLA Allele

This mode allows you to retrieve all TCRs associated with a specific HLA allele from your database.

**Mandatory Arguments**: `--path2database`, `--database`, `--input_allele`, `--path2write_results`

**Example**:
```bash
python tcr_query_tool.py \
    --path2database /data/tcr_hla_db \
    --database TRA \
    --input_allele HLA-A*02:01 \
    --path2write_results /results/hla_a0201_matches
```

---

### Querying a Single TCR

This mode is used to find HLA alleles or other associated data for a single, specific TCR sequence.

**Mandatory Arguments**: `--path2database`, `--database`, `--input_tcr`, `--path2write_results`

**--input_tcr Format**:  
The TCR string must follow the format: `v_gene+CDR3+j_gene:format`

- `v_gene`: The name of the V gene (e.g., `TRBV12-3*01`)
- `CDR3`: The amino acid sequence of the CDR3 region (e.g., `CASSPGASGYTF`)
- `j_gene`: The name of the J gene (e.g., `TRBJ2-7*01`)
- `format`: The naming convention, either `'Adaptive'` or `'IMGT'`

**Example**:
```bash
python tcr_query_tool.py \
    --path2database /data/tcr_hla_db \
    --database TRB \
    --input_tcr "TRBV12-3*01+CASSPGASGYTF+TRBJ2-7*01:IMGT" \
    --path2write_results /results/single_tcr_match
```

---

### Querying a Table of TCRs

This mode processes a batch of TCRs provided in a TSV file, querying the database for each entry.

**Mandatory Arguments**: `--path2database`, `--database`, `--query_tcr_tables`, `--query_table_format`, `--path2write_results`

**--query_tcr_tables Requirements**:  
The input TSV file must contain at least three columns: `v_gene`, `CDR3`, and `j_gene`.

**--query_table_format Values**:  
Must be either `'Adaptive'` or `'IMGT'`.

**Example**:
```bash
python tcr_query_tool.py \
    --path2database /data/tcr_hla_db \
    --database TRA \
    --query_tcr_tables /inputs/my_tcr_list.tsv \
    --query_table_format Adaptive \
    --path2write_results /results/batch_query_output
```

---

## 4. Arguments

| Argument               | Type | Mandatory | Description |
|------------------------|------|-----------|-------------|
| `--path2database`      | str  | Yes       | Path to the TCR-HLA databases directory. |
| `--database`           | str  | Yes       | Which database to query: `'TRA'` or `'TRB'`. |
| `--input_allele`       | str  | No        | HLA allele to query. Required for HLA mode. |
| `--input_tcr`          | str  | No        | TCR string in format `"v+CDR3+j:format"`. Required for TCR mode. |
| `--query_tcr_tables`   | str  | No        | Path to TSV table of TCRs. Required for batch mode. |
| `--query_table_format` | str  | No        | Either `'Adaptive'` or `'IMGT'`. Required with `--query_tcr_tables`. |
| `--path2write_results` | str  | Yes       | Directory + prefix where output TSV(s) will be written. |

---

## 5. Input and Output Formats

### Input Database
The database should be structured and located according to the `--path2database` path. The script will use the selected `--database` type (`TRA` or `TRB`) to query appropriately.

### Input TCR Table (`--query_tcr_tables`)
Must be a TSV file with at least the following headers:

```
v_gene	CDR3	j_gene	other_data
TRBV20-1*01	CASSLAGATTYEQYF	TRBJ2-7*01	sample1
TRBV29-1*01	CASSQEGVGYTF	TRBJ2-1*01	sample2
```

Other columns are allowed and will be passed through.

### Output Results (`--path2write_results`)
The results will be saved as one or more TSV files in the location defined by the prefix. For example:
```
/path/to/output/results_prefix_matches.tsv
```

Output content may include:
- Matching TCRs and their HLA associations
- Additional metadata (if implemented)
- Echoed input information (in table mode)

---

## 6. Modes of Operation

The script detects its query mode based on the arguments you provide:

### I. Querying an HLA allele
- Triggered by: `--input_allele`
- Requires: `--path2database`, `--database`, `--input_allele`, `--path2write_results`
- Purpose: Return all TCRs associated with a specific HLA allele.

### II. Querying a single TCR
- Triggered by: `--input_tcr`
- Requires: `--path2database`, `--database`, `--input_tcr`, `--path2write_results`
- Purpose: Return all HLA alleles associated with a specific TCR.

### III. Querying a table of TCRs
- Triggered by: `--query_tcr_tables`
- Requires: `--path2database`, `--database`, `--query_tcr_tables`, `--query_table_format`, `--path2write_results`
- Purpose: Process a batch of TCRs and return matching results for each.