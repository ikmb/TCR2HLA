#!/usr/bin/env python3
"""
@author: Hesham ElAbd
@brief: a tool that can be used for imputing HLA from TCR repertoire
"""
## Load the modules
#------------------
import argparse
import pandas as pd
import numpy as np
import json
import os
import sys
from gene_formating_converter import *
from impute_HLA_alleles_API import *

# define the function that is doing the work
#-------------------------------------------
def validate_arguments(args):
    """
    Validates the input arguments for correctness.
    """
    if args.gene_naming not in ['Adaptive', 'IMGT']:
        print(f"Error: Invalid value for --gene_naming. Must be 'Adaptive' or 'IMGT'. Found: {args.gene_naming}", file=sys.stderr)
        sys.exit(1)

    if args.model not in ['TRA', 'TRB']:
        print(f"Error: Invalid value for --model. Must be 'TRA' or 'TRB'. Found: {args.model}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.input_repertoire):
        print(f"Error: Input repertoire file not found at '{args.input_repertoire}'", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(args.path2models):
        print(f"Error: Models directory not found at '{args.path2models}'", file=sys.stderr)
        sys.exit(1)

    if (args.normalize_prob) and (not args.return_prob):
        print("Error: to return normalised probabilities, probabilities needs to be first set, please use the following flag: --return_prob", file=sys.stderr)
        sys.exit(1)

    # Ensure the output directory exists or can be created
    output_dir = os.path.dirname(args.path2write_results)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
        except OSError as e:
            print(f"Error: Could not create output directory '{output_dir}': {e}", file=sys.stderr)
            sys.exit(1)

def load_repertoire_data(file_path):
    """
    Loads the repertoire data from the given CSV file.
    Checks for mandatory columns: v_gene, CDR3, j_gene, sample_name.
    """
    try:
        df = pd.read_csv(file_path,sep='\t')
    except FileNotFoundError:
        print(f"Error: Input file not found at {file_path}", file=sys.stderr)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: Input file {file_path} is empty.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading CSV file {file_path}: {e}", file=sys.stderr)
        sys.exit(1)

    mandatory_columns = ['v_gene', 'CDR3', 'j_gene', 'sample_name', 'count']
    missing_columns = [col for col in mandatory_columns if col not in df.columns]
    if missing_columns:
        print(f"Error: Missing mandatory columns in '{file_path}': {', '.join(missing_columns)}", file=sys.stderr)
        sys.exit(1)
    return df

def main():
    parser = argparse.ArgumentParser(
        description="A command-line tool for HLA allele imputation from TCR repertoire data."
    )

    parser.add_argument(
        '--input_repertoire',
        type=str,
        required=True,
        help="Path to the CSV file containing TCR repertoire data (mandatory columns: v_gene, CDR3, j_gene, sample_name, count)."
    )
    parser.add_argument(
        '--gene_naming',
        type=str,
        required=True,
        choices=['Adaptive', 'IMGT'],
        help="The formatting of V and J gene naming. Must be 'Adaptive' or 'IMGT'."
    )
    parser.add_argument(
        '--model',
        type=str,
        required=True,
        choices=['TRA', 'TRB'],
        help="Which model shall be used for imputing HLA alleles. Either 'TRA' or 'TRB' clonotypes."
    )
    parser.add_argument(
        '--path2models',
        type=str,
        required=True,
        help="Path to the directory where HLA imputation models are defined."
    )
    parser.add_argument(
        '--return_prob',
        action='store_true',  # This is the key change!
        required=False,       # Can be omitted as False is the default for optional args
        help="If set, the program will return probabilities instead of hard (0/1) predictions for allele carriership."
    )
    
    parser.add_argument(
        '--normalize_prob',
        action='store_true', # This is the key for a boolean flag
        help="If set, normalizes the probabilities among the different alleles of a given locus (e.g., among HLA-A, HLA-B, etc.)."
    )
    
    parser.add_argument(
        '--path2write_results',
        type=str,
        required=True,
        help="Base path to write results a TSV file"
    )

    args = parser.parse_args()

    # 1. Validate arguments
    validate_arguments(args)

    # 2. Load input repertoire data and perfom some basic statstics
    #--------------------------------------------------------------
    print(f"Loading repertoire data from: {args.input_repertoire}")
    repertoire_df = load_repertoire_data(args.input_repertoire)
    print(f"Successfully loaded: {len(repertoire_df)} entries derived from: {len(repertoire_df.sample_name.unique())}")
    num_clonotypes_per_samples = repertoire_df.groupby('sample_name').size().reset_index().rename(columns={0:'num_clonotypes'})
    print(f"Average number of clonotypes per sample: {np.mean(num_clonotypes_per_samples.num_clonotypes)}\n the lowest number of clonotypes per sample is : {min(num_clonotypes_per_samples.num_clonotypes)} and the maximum is : {max(num_clonotypes_per_samples.num_clonotypes)}")

    if min(num_clonotypes_per_samples.num_clonotypes)<25_000:
        print(f"Warning, found {num_clonotypes_per_samples.loc[num_clonotypes_per_samples.num_clonotypes<25_000].sample_name.unique().shape[0]} samples with less than 25,000 clonotypes, shallow sequencing might impact the recall of the generated models.")

    ## 3. Unify gene names
    #---------------------
    if args.model == 'TRA' and args.gene_naming == 'Adaptive':
        Adaptive2IMGT = adaptive_to_imgt_builder(path2database=args.path2models)
        repertoire_df['v_gene'] = repertoire_df['v_gene'].replace(Adaptive2IMGT)
        repertoire_df['j_gene'] = repertoire_df['j_gene'].replace(Adaptive2IMGT)
    
    if args.model == 'TRB' and args.gene_naming == 'IMGT':
        IMGT2Adaptive = imgt_to_adaptive_builder(path2database=args.path2models)
        repertoire_df['v_gene'] = repertoire_df['v_gene'].replace(IMGT2Adaptive)
        repertoire_df['j_gene'] = repertoire_df['j_gene'].replace(IMGT2Adaptive)

    ## 4. Run the predictions
    #------------------------
    imputer = HLAPredictor(path2models=args.path2models,model=args.model)
    if args.return_prob:
        imputed_alleles = imputer.predict_HLA_for_a_sample_collection_prob(repertoire_df)
        if args.normalize_prob:
            imputed_alleles_tsv = normalise_by_softmax(translate_probabilities_to_tsv(imputed_alleles))
        else:
            imputed_alleles_tsv = translate_probabilities_to_tsv(imputed_alleles)
    else:
        imputed_alleles = imputer.predict_HLA_for_a_sample_collection(repertoire_df)
        imputed_alleles_tsv = translate_prediction_to_tsv(imputed_alleles)
    
    # 5. write the results 
    #---------------------
    imputed_alleles_tsv.to_csv(args.path2write_results,sep='\t',index=False)
    
    # 6. report the end of execusion
    #------------------------------
    print("\nTool execution finished successfully!")

if __name__ == "__main__":
    main()