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

The `tcr_query_tool.py` facilitates the interaction with TCR-HLA association databases. It streamlines the process of extracting relevant TCR data based on user-defined criteria, supporting three primary querying methods:

- By HLA allele name  
- By a single TCR sequence  
- By a batch of TCR sequences provided in a table  

The results are outputted to a specified TSV (Tab Separated Values) file.

---

## 2. Installation

This tool is a standalone Python script and requires Python 3.x.  

### To use:
1. Ensure you have Python 3 installed.  
2. Save the provided script as `tcr_query_tool.py`.  
3. *(Optional)* Install any additional Python packages required by your internal database interaction logic (e.g., `pandas`, database connectors).

---

## 3. Usage

The tool operates based on mutually exclusive query types. You must specify **exactly one** primary query argument (`--input_allele`, `--input_tcr`, or `--query_tcr_tables`).

### General Syntax

```bash
python tcr_query_tool.py --path2database <path_to_db> --path2write_results <output_path_prefix> [query_type_arguments]
```

---

### Querying an HLA Allele

This mode retrieves all TCRs associated with a specific HLA allele.

**Mandatory Arguments:**  
`--path2database`, `--input_allele`, `--path2write_results`

**Example:**

```bash
python tcr_query_tool.py \
    --path2database /data/tcr_hla_db \
    --input_allele HLA-A*02:01 \
    --path2write_results /results/hla_a0201_matches
```

---

### Querying a Single TCR

This mode finds HLA alleles or other associated data for a specific TCR sequence.

**Mandatory Arguments:**  
`--path2database`, `--input_tcr`, `--path2write_results`

**`--input_tcr` Format:**

```text
"v_gene+CDR3+j_gene:format"
```

- `v_gene`: V gene name (e.g., `TRBV12-3*01`)  
- `CDR3`: Amino acid sequence (e.g., `CASSPGASGYTF`)  
- `j_gene`: J gene name (e.g., `TRBJ2-7*01`)  
- `format`: Either `'Adaptive'` or `'IMGT'`

**Example:**

```bash
python tcr_query_tool.py \
    --path2database /data/tcr_hla_db \
    --input_tcr "TRBV12-3*01+CASSPGASGYTF+TRBJ2-7*01:IMGT" \
    --path2write_results /results/single_tcr_match
```

---

### Querying a Table of TCRs

This mode processes a batch of TCRs from a TSV file.

**Mandatory Arguments:**  
`--path2database`, `--query_tcr_tables`, `--query_table_format`, `--path2write_results`

**TSV File Requirements:**

- Must contain `v_gene`, `CDR3`, and `j_gene` columns.  
- `--query_table_format` must be `'Adaptive'` or `'IMGT'`.

**Example:**

```bash
python tcr_query_tool.py \
    --path2database /data/tcr_hla_db \
    --query_tcr_tables /inputs/my_tcr_list.tsv \
    --query_table_format Adaptive \
    --path2write_results /results/batch_query_output
```

---

## 4. Arguments

| Argument | Type | Mandatory | Description |
|----------|------|-----------|-------------|
| `--path2database` | str | Yes | Path to your TCR-HLA database directory. |
| `--input_allele` | str | No | HLA allele to query (e.g., `HLA-A*02:01`). Mutually exclusive. |
| `--input_tcr` | str | No | TCR string in `"v_gene+CDR3+j_gene:format"` format. Mutually exclusive. |
| `--query_tcr_tables` | str | No | Path to a TSV file with TCRs. Mutually exclusive. |
| `--query_table_format` | str | No | Format for V and J genes: `'Adaptive'` or `'IMGT'`. Required with `--query_tcr_tables`. |
| `--path2write_results` | str | Yes | Path to write output files (with filename prefix). |

---

## 5. Input and Output Formats

### Input Database

Must be accessible at the `--path2database` path. Internal logic (to be implemented) should parse the database contents.

### Input TCR Table (`--query_tcr_tables`)

**Example `my_tcr_list.tsv`:**

```tsv
v_gene	CDR3	j_gene	other_data
TRBV20-1*01	CASSLAGATTYEQYF	TRBJ2-7*01	sample1
TRBV29-1*01	CASSQEGVGYTF	TRBJ2-1*01	sample2
```

---

### Output Results (`--path2write_results`)

Output will be one or more TSV files named like:

```
results_prefix_matches.tsv
```

Content will include:

- Matched TCRs and their associated HLA alleles  
- Other relevant data, depending on your database query implementation  
- For table queries, original TCR info with results

---

## 6. Modes of Operation

The script determines the query mode from the arguments:

### I. Querying an HLA Allele

**Trigger:** `--input_allele`  
**Required:** `--path2database`, `--input_allele`, `--path2write_results`  
**Purpose:** Find TCRs associated with an HLA allele

---

### II. Querying a Single TCR

**Trigger:** `--input_tcr`  
**Required:** `--path2database`, `--input_tcr`, `--path2write_results`  
**Purpose:** Find HLA alleles or records associated with the specific TCR

---

### III. Querying a Table of TCRs

**Trigger:** `--query_tcr_tables`  
**Required:** `--path2database`, `--query_tcr_tables`, `--query_table_format`, `--path2write_results`  
**Purpose:** Query a batch of TCRs from a TSV file

---