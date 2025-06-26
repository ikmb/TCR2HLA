#!/usr/bin/env python3

"""
@author: Hesham ElAbd
@brief: the script is a command line tool that can be used for getting clonotypes from the dataset
"""

## import needed modules
import argparse
import os
from query_API import *

def main():
    """
    Main function to parse command-line arguments and validate input.
    """
    parser = argparse.ArgumentParser(
        description="A command-line tool for querying TCR-HLA databases based on HLA alleles, single TCRs, or TCR tables.",
        formatter_class=argparse.RawTextHelpFormatter # Preserve newlines in help
    )

    # Define arguments
    parser.add_argument(
        '--path2database',
        type=str,
        required=True,
        help="[MANDATORY] The path to TCR-HLA databases."
    )
    parser.add_argument(
        '--database',
        type=str,
        required=True,
        choices=['TRA', 'TRB'],
        help="[MANDATORY] The specific database to query: 'TRA' for TRA-HLA dataset or 'TRB' for TRB-HLA dataset."
    )
    parser.add_argument(
        '--input_allele',
        type=str,
        help="[OPTIONAL] The HLA allele name (e.g., 'HLA-A*02:01').\n"
             "  Mandatory for 'Querying an HLA allele' mode."
    )
    parser.add_argument(
        '--input_tcr',
        type=str,
        help="[OPTIONAL] A string representing a single TCR. Structured as:\n"
             "  'v_gene+CDR3+j_gene:format'\n"
             "  Where 'v_gene' is the V gene, 'CDR3' is the amino acid sequence,\n"
             "  'j_gene' is the J gene, and 'format' is the naming convention\n"
             "  ('Adaptive' or 'IMGT'). '+' and ':' are mandatory separators.\n"
             "  Example: 'TRBV12-3*01+CASSPGASGYTF+TRBJ2-7*01:IMGT'\n"
             "  Mandatory for 'Querying a single TCR' mode."
    )
    parser.add_argument(
        '--query_tcr_tables',
        type=str,
        help="[OPTIONAL] Path to a TSV table containing TCRs to query.\n"
             "  Must contain three mandatory columns: 'v_gene', 'CDR3', and 'j_gene'.\n"
             "  Mandatory for 'Querying a table of TCRs' mode."
    )
    parser.add_argument(
        '--query_table_format',
        type=str,
        help="[OPTIONAL] A string specifying the naming convention format of the\n"
             "  V and J genes in the '--query_tcr_tables' file. Can be 'Adaptive' or 'IMGT'.\n"
             "  Mandatory when '--query_tcr_tables' is provided."
    )
    parser.add_argument(
        '--path2write_results',
        type=str,
        required=True,
        help="[MANDATORY] The path (directory and filename prefix) to write results\n"
             "  as TSV tables. E.g., '/path/to/output/results_prefix'."
    )

    args = parser.parse_args()

    # --- Input Validation Logic ---

    # Determine which mode the user intends based on provided arguments
    is_allele_query = args.input_allele is not None
    is_single_tcr_query = args.input_tcr is not None
    is_table_tcr_query = args.query_tcr_tables is not None

    # Count how many primary query arguments are provided
    query_modes_active = sum([is_allele_query, is_single_tcr_query, is_table_tcr_query])

    if query_modes_active == 0:
        parser.error(
            "Error: You must specify one of the following query types:\n"
            "  --input_allele (for HLA allele query)\n"
            "  --input_tcr (for single TCR query)\n"
            "  --query_tcr_tables (for TCR table query)"
        )
    elif query_modes_active > 1:
        parser.error(
            "Error: You can only perform one type of query at a time.\n"
            "  Please provide only one of --input_allele, --input_tcr, or --query_tcr_tables."
        )

    # Mode I: Querying an HLA allele
    if is_allele_query:
        print("Mode: Querying an HLA allele")
        if args.input_tcr or args.query_tcr_tables or args.query_table_format:
            parser.error(
                "Error: For HLA allele query, only --path2database, --database, --input_allele, "
                "and --path2write_results are expected."
            )
        # All mandatory args for this mode are handled by required=True or existence check
        # No specific extra checks for 'database' here as it's universally required.

    # Mode II: Querying a single TCR
    elif is_single_tcr_query:
        print("Mode: Querying a single TCR")
        if args.input_allele or args.query_tcr_tables or args.query_table_format:
            parser.error(
                "Error: For single TCR query, only --path2database, --database, --input_tcr, "
                "and --path2write_results are expected."
            )
        # Basic format check for input_tcr (more robust parsing would be needed in the actual logic)
        if '+' not in args.input_tcr or ':' not in args.input_tcr:
            parser.error(
                "Error: --input_tcr must be in 'v_gene+CDR3+j_gene:format' format."
            )

    # Mode III: Querying a table of TCRs
    elif is_table_tcr_query:
        print("Mode: Querying a table of TCRs")
        if args.input_allele or args.input_tcr:
            parser.error(
                "Error: For TCR table query, only --path2database, --database, --query_tcr_tables, "
                "--query_table_format, and --path2write_results are expected."
            )
        if not args.query_table_format:
            parser.error(
                "Error: --query_table_format is mandatory when --query_tcr_tables is provided."
            )
        if args.query_table_format not in ['Adaptive', 'IMGT']:
            parser.error(
                "Error: --query_table_format must be either 'Adaptive' or 'IMGT'."
            )
        if not os.path.exists(args.query_tcr_tables):
            parser.error(f"Error: The file specified by --query_tcr_tables does not exist: {args.query_tcr_tables}")

    # Validate existence of database path (required for all modes)
    if not os.path.exists(args.path2database):
        parser.error(f"Error: The database path does not exist: {args.path2database}")

    # Validate output path parent directory
    output_dir = os.path.dirname(args.path2write_results)
    if output_dir and not os.path.exists(output_dir): # output_dir might be empty if path2write_results is just a filename
        try:
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
        except OSError as e:
            parser.error(f"Error: Could not create output directory {output_dir}: {e}")


    # If all validations pass, you can proceed with the core logic
    print("\nAll arguments validated successfully!")
    print(f"Database Path: {args.path2database}")
    print(f"Database Type: {args.database}")
    print(f"Output Path Prefix: {args.path2write_results}")

    if is_allele_query:
        print(f"Querying HLA Allele: {args.input_allele}")
        # Your wiring for HLA allele query goes here
        # For example: process_hla_allele_query(args.path2database, args.database, args.input_allele, args.path2write_results)
    elif is_single_tcr_query:
        print(f"Querying Single TCR: {args.input_tcr}")
        # Your wiring for single TCR query goes here
        # For example: process_single_tcr_query(args.path2database, args.database, args.input_tcr, args.path2write_results)
    elif is_table_tcr_query:
        print(f"Querying TCR Table: {args.query_tcr_tables}")
        print(f"TCR Table Format: {args.query_table_format}")
        # Your wiring for TCR table query goes here
        # For example: process_tcr_table_query(args.path2database, args.database, args.query_tcr_tables, args.query_table_format, args.path2write_results)

if __name__ == "__main__":
    main()
