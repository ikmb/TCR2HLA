"""
@author: Hesham ElAbd
@brief: The code contain gene formatting converter, where genes can be convereted from the IMGT to Adaptive naming convention and vice versa. 
The converter functions are based on the tables released by tcrdist3 at: https://github.com/kmayerb/tcrdist3/blob/master/tcrdist/db/adaptive_imgt_mapping.csv
"""
# Load the modules
import os
import pandas as pd
from typing import Dict

def imgt_to_adaptive_builder(path2database:str)->Dict[str,str]:
    """
    Builds a mapping dictionary from IMGT gene names to Adaptive gene names for human species.

    This function reads a TSV file named 'adaptive_imgt_mapping.csv' from the
    specified database path, filters it for human species, and creates a dictionary
    where IMGT gene names are keys and Adaptive gene names are values.

    Args:
        path2database (str): The path to the directory where the 'adaptive_imgt_mapping.csv'
                             file is located.

    Returns:
        Dict[str, str]: A dictionary mapping IMGT gene names (keys) to Adaptive gene names (values).
                        Returns an empty dictionary if the file is not found or no human data is present.

    Raises:
        FileNotFoundError: If 'adaptive_imgt_mapping.csv' is not found at the specified path.
        KeyError: If expected columns ('species', 'imgt', 'adaptive') are missing from the CSV.
    """
    # Sanity check for file existence
    mapping_file_path = f'{path2database}/adaptive_imgt_mapping.csv'
    if not os.path.exists(mapping_file_path):
        raise FileNotFoundError(f"Mapping file not found at: {mapping_file_path}")

    input_map = pd.read_csv(mapping_file_path, sep=',')

    # Sanity check for required columns
    required_cols = ['species', 'imgt', 'adaptive']
    if not all(col in input_map.columns for col in required_cols):
        raise KeyError(f"Mapping file is missing one or more required columns: {required_cols}")

    input_map = input_map.loc[input_map['species']=='human',]

    # Handle case where no human data is found after filtering
    if input_map.empty:
        print("Warning: No 'human' species data found in the adaptive_imgt_mapping.csv file.")
        return {}

    return {imgt:adaptive for imgt, adaptive in zip(input_map['imgt'], input_map['adaptive'])}

def adaptive_to_imgt_builder(path2database:str)->Dict[str,str]:
    """
    Builds a mapping dictionary from Adaptive gene names to IMGT gene names for human species.

    This function reads a TSV file named 'adaptive_imgt_mapping.csv' from the
    specified database path, filters it for human species, and creates a dictionary
    where Adaptive gene names are keys and IMGT gene names are values.

    Args:
        path2database (str): The path to the directory where the 'adaptive_imgt_mapping.csv'
                             file is located.

    Returns:
        Dict[str, str]: A dictionary mapping Adaptive gene names (keys) to IMGT gene names (values).
                        Returns an empty dictionary if the file is not found or no human data is present.

    Raises:
        FileNotFoundError: If 'adaptive_imgt_mapping.csv' is not found at the specified path.
        KeyError: If expected columns ('species', 'imgt', 'adaptive') are missing from the CSV.
    """
    # Sanity check for file existence
    mapping_file_path = f'{path2database}/adaptive_imgt_mapping.csv'
    if not os.path.exists(mapping_file_path):
        raise FileNotFoundError(f"Mapping file not found at: {mapping_file_path}")
    # Loading 
    input_map = pd.read_csv(mapping_file_path, sep=',')

    # Sanity check for required columns
    required_cols = ['species', 'imgt', 'adaptive']
    if not all(col in input_map.columns for col in required_cols):
        raise KeyError(f"Mapping file is missing one or more required columns: {required_cols}")

    input_map = input_map.loc[input_map['species']=='human',]

    # Handle case where no human data is found after filtering
    if input_map.empty:
        print("Warning: No 'human' species data found in the adaptive_imgt_mapping.csv file.")
        return {}

    return {adaptive:imgt for adaptive, imgt in zip(input_map['adaptive'], input_map['imgt'])}
