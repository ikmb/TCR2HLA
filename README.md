# TCR2HLA Tools and Databases

This repository provides tools and resources for analyzing T-Cell Receptor (TCR) repertoires in relation to Human Leukocyte Antigen (HLA) proteins. It includes functionalities for HLA imputation from TCR repertoire data and for querying a database of TCR-associated alleles.

## Preprint

For a detailed description of the methods and data presented in this repository, please refer to our preprint:

[Decoding the restriction of T cell receptors to human leukocyte antigen alleles using statistical learning](https://doi.org/10.1101/2025.02.06.636910)

## Overview

This repository is structured as follows:

* **`/CLI_tools/impute_HLA`**: Contains the `impute_HLAs` tool for imputing HLA alleles from TCR repertoire data.
* **`/CLI_tools/query_database`**: Contains the `query_database` tool for interfacing with a database of HLA-associated TCR clonotypes.
* **`/databases`**: Stores various databases, including TRA and TRB clonotypes strongly associated with specific HLA proteins.
* **`/models`**: Stores the models and the weights for imputing HLA alleles from TCR repertoire informations.
* **`test_dataset.tsv`**: A test TRB repertoire derived from the [immuneCODE database](https://www.researchsquare.com/article/rs-51964/v1). This dataset can be used to check the installation and test the functionality of the provided code.
* **`requirements.txt`**: Lists the necessary Python packages for running the tools.


## Tools

We provide two primary tools, each available as both a command-line interface (CLI) and a programmable API:

### 1. `impute_HLAs` (HLA Imputation from TCR Repertoire)

This tool allows users to predict HLA alleles based on their TCR repertoire data.

* **Functionality:** Imputes HLA from TCR repertoire.
* **Location:** `/CLI_tools/impute_HLA/`
* **Usage:**
    * **Command Line:** Detailed usage instructions can be found in `/CLI_tools/impute_HLA/README.md`.
    * **Programmable Interface (API):** For scripting and integration into other Python projects, refer to the documentation within the `/CLI_tools/impute_HLA/` directory.

### 2. `query_database` (TCR-HLA Database Interface)

This tool enables users to query our curated database of HLA-associated alleles.

* **Functionality:** Queries a database of HLA-associated alleles.
* **Location:** `/CLI_tools/query_database/`
* **Usage:**
    * **Command Line:** Detailed usage instructions can be found in `/CLI_tools/query_database/README.md`.
    * **Programmable Interface (API):** For scripting and integration into other Python projects, refer to the documentation within the `/CLI_tools/query_database/` directory.

## Databases

The `/databases` directory contains valuable resources, including:

* **`TRA` and `TRB` clonotypes:** These files contain TCR Alpha (TRA) and TCR Beta (TRB) clonotypes that have been identified as strongly associated with specific HLA proteins.
* **`adaptive_imgt_mapping.tsv`:** This file provides a mapping between Adaptive Biotechnologies' nomenclature and IMGT (International Immunogenetics Information System) standard nomenclature for TCR gene names. The file is obtained from the `tcrdist3` tool.

## Installation

To set up the project and install the required dependencies, follow these steps:

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/ikmb/TCR2HLA.git
    cd TCR2HLA
    ```
    
2.  **Install dependencies:**
    Ensure you have `pip` installed. Then, install the packages listed in `requirements.txt`:
    ```bash
    pip install -r requirements.txt
    ```

## Contributing

We welcome contributions! Please see `CONTRIBUTING.md` for guidelines on how to contribute to this project.


## Contact

For any questions or issues, please open an issue on the GitHub repository or contact h.elabd@ikmb.uni-kiel.de.