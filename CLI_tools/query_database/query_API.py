"""
@author: Hesham ElAnd
@brief: The file contain an interface class that can be used to query the HLA-TCR database,
either it can query the HLA restriction associated with a given HLA or the TCR associated with a given HLA.
"""
## Load the modules
import os
import time
import pandas as pd
import numpy as np
import Levenshtein
from gene_formating_converter import imgt_to_adaptive_builder, adaptive_to_imgt_builder

class HLATCRDBInterface:
    """
    A class to interface with TCR-HLA association databases (TRA or TRB).

    This class loads a specified TCR-HLA database and provides methods to query
    clonotypes based on HLA alleles or matching TCR sequences using Levenshtein distance.

    Attributes:
        _database (pd.DataFrame): The loaded TCR-HLA database.
        _DATA_TYPE (str): The type of database loaded ('TRA' or 'TRB').
    """
    def __init__(self, database_dir: str, database_name: str):
        """
        Initializes the HLATCRDBInterface by loading the specified database.

        Args:
            database_dir (str): The path to the directory where the database parquet files are located.
            database_name (str): The name of the database to load. Must be 'TRB_DB' for TRB-HLA
                                 database or 'TRA_DB' for TRA-HLA database.

        Raises:
            ValueError: If an unsupported `database_name` is provided.
            FileNotFoundError: If the specified database file does not exist.
        """
        self._database_path=database_dir
        match database_name:
            case 'TRB':
                database_path = f'{database_dir}/databases/TRB_database.parquet'
                self._DATA_TYPE = 'TRB'
            case 'TRA':
                database_path = f'{database_dir}/databases/TRA_database.parquet'
                self._DATA_TYPE = 'TRA'
            case _:
                raise ValueError(
                    "Unsupported database_name. Only 'TRB_DB' (for TRB-HLA database) "
                    "and 'TRA_DB' (for TRA-HLA database) are supported."
                )
        
        if not os.path.exists(database_path):
            raise FileNotFoundError(f"Database file not found at: {database_path}")
        
        self._database = pd.read_parquet(database_path)
        
        # Ensure essential columns exist in the loaded database
        required_cols = ['loci', 'allele_name', 'v_gene', 'j_gene', 'CDR3']
        if not all(col in self._database.columns for col in required_cols):
            raise ValueError(
                f"Loaded database is missing one or more required columns: {required_cols}"
            )
        
        # Verify that 'loci' and 'allele_name' columns are present for allele queries
        if 'loci' not in self._database.columns or 'allele_name' not in self._database.columns:
            raise ValueError("Database must contain 'loci' and 'allele_name' columns for allele queries.")

    def get_clonotypes_associated_to_alleles(self, allele_name: str) -> pd.DataFrame:
        """
        Retrieves restricted clonotypes associated with a given HLA allele.

        Args:
            allele_name (str): The HLA allele name. Expected format for single alleles
                               is 'Locus*AlleleName', e.g., 'A*02:01'.
                               For DQ and DP alleles, format can be 'DQ*01:02+04:01' or 'DP*01:01+04:01'.

        Returns:
            pd.DataFrame: A table containing all clonotypes associated with the given allele.
                          The table is composed of columns such as 'V genes', 'J genes', and 'CDR3'.

        Raises:
            ValueError: If the provided locus or allele name is not supported/found in the database.
        """
        if '-' not in allele_name:
            raise ValueError("Invalid allele_name format. Expected 'Locus*AlleleName', e.g., 'A*02:01'.")

        loci_name, allele = allele_name.split('-')

        print('Checking your input ...')
        
        unique_loci_in_db = self._database.loci.unique()
        if loci_name not in unique_loci_in_db:
            raise ValueError(
                f"Locus: {loci_name} is not supported. Only these loci are supported: {list(unique_loci_in_db)}"
            )
        
        # Filter database for the specific locus first for efficiency
        locus_data = self._database.loc[self._database.loci == loci_name]
        if locus_data.empty:
            raise ValueError(f"No data found for locus: {loci_name} in the database.")

        unique_alleles_in_locus = locus_data.allele_name.unique()
        
      
        if allele_name not in unique_alleles_in_locus:
            raise ValueError(
                f"Allele: {allele} from Locus: {loci_name} is not supported. "
                f"Only these alleles from this locus are supported: {list(unique_alleles_in_locus)}"
            )
        
        print('Checks finished ... supported allele(s).')
        
        return self._database.loc[
            (self._database.loci == loci_name) & (self._database.allele_name==allele_name)
        ]
    
    def get_restricted_clonotypes_matching_a_clonotype(self,
        v_gene: str, j_gene: str, format: str, CDR3_amino_acids: str, num_mismatches: int = 1
    ) -> pd.DataFrame:
        """
        Finds restricted clonotypes in the database that match a single input clonotype.

        Matches are determined by identical V and J genes and a Levenshtein distance
        within `num_mismatches` for the CDR3 amino acid sequence. This function
        handles gene format conversion if needed based on the `_DATA_TYPE` of the
        loaded database and the provided `format`.

        Args:
            v_gene (str): The V gene name of the query clonotype.
            j_gene (str): The J gene name of the query clonotype.
            format (str): The naming convention of the V and J genes provided.
                          Must be 'Adaptive' or 'IMGT'.
            CDR3_amino_acids (str): The amino acid sequence of the CDR3 region of the query clonotype.
            num_mismatches (int, optional): The maximum allowed Levenshtein distance
                                            for CDR3 sequences to be considered a match. Defaults to 1.

        Returns:
            pd.DataFrame: A table containing matching clonotypes from the database,
                          including an additional 'distance' column showing the Levenshtein distance
                          to the query CDR3.

        Raises:
            ValueError: If an unsupported format is provided, or if the V/J gene is not found
                        in the (potentially converted) database.
        """
        if format not in ['Adaptive', 'IMGT']:
            raise ValueError(f"Unsupported format: '{format}'. The function supports V and J gene "
                             "name formatting in either 'IMGT' or 'Adaptive' gene formatting.")
        
        work_database = self._database.copy(deep=True)

        # Apply gene format conversion if necessary
        if self._DATA_TYPE == 'TRB' and format == 'IMGT':
            # Assuming adaptive_to_imgt_builder() returns a dict for replacement
            if 'adaptive_to_imgt_builder' not in globals() and 'adaptive_to_imgt_builder' not in locals():
                 raise NotImplementedError("adaptive_to_imgt_builder() is not defined or imported. "
                                           "Please implement or import this function.")
            Adaptive2IMGT = adaptive_to_imgt_builder(path2database=self._database_path)
            work_database['v_gene'] = work_database['v_gene'].replace(Adaptive2IMGT)
            work_database['j_gene'] = work_database['j_gene'].replace(Adaptive2IMGT)
        
        if self._DATA_TYPE == 'TRA' and format == 'Adaptive':
            # Assuming imgt_to_adaptive_builder() returns a dict for replacement
            if 'imgt_to_adaptive_builder' not in globals() and 'imgt_to_adaptive_builder' not in locals():
                raise NotImplementedError("imgt_to_adaptive_builder() is not defined or imported. "
                                          "Please implement or import this function.")
            IMGT2Adaptive = imgt_to_adaptive_builder(path2database=self._database_path)
            work_database['v_gene'] = work_database['v_gene'].replace(IMGT2Adaptive)
            work_database['j_gene'] = work_database['j_gene'].replace(IMGT2Adaptive)
        
        if v_gene not in work_database.v_gene.unique():
            raise ValueError(
                f"The provided V-gene: '{v_gene}' is not defined in the database (after potential format conversion). "
                f"Alternatively, it might not be in the correct format: '{format}'."
            )
        
        if j_gene not in work_database.j_gene.unique():
            raise ValueError(
                f"The provided J-gene: '{j_gene}' is not defined in the database (after potential format conversion). "
                f"Alternatively, it might not be in the correct format: '{format}'."
            )
        
        candiate_table = work_database.loc[
            (work_database.v_gene == v_gene) & (work_database.j_gene == j_gene)
        ].copy(deep=True)

        # Ensure CDR3 column exists before trying to calculate distance
        if 'CDR3' not in candiate_table.columns:
            raise KeyError("The 'CDR3' column is missing from the database subset.")

        # Using Levenshtein.distance requires the library to be imported
        if not hasattr(Levenshtein, 'distance'):
            raise ImportError("Levenshtein library is required but 'distance' method not found.")

        candiate_table['distance'] = [Levenshtein.distance(str(c), str(CDR3_amino_acids)) for c in candiate_table.CDR3]
        target_clonotypes = candiate_table.loc[candiate_table['distance'] <= num_mismatches].copy(deep=True)
        
        # Explicitly delete large DataFrames to free memory, especially important for large datasets
        del work_database
        del candiate_table 
        
        return target_clonotypes

    def get_restricted_clonotypes_matching_clonotypes(self,
        query_clonotypes: pd.DataFrame, format: str, num_mismatches: int = 1
    ) -> pd.DataFrame:
        """
        Finds restricted clonotypes in the database that match clonotypes from an input DataFrame.

        This function performs a cross-merge between the loaded database and the
        `query_clonotypes` DataFrame, then filters for V/J gene matches and
        calculates Levenshtein distance for CDR3 sequences within `num_mismatches`.
        It handles gene format conversion as per `_DATA_TYPE` and `format`.

        Args:
            query_clonotypes (pd.DataFrame): A DataFrame containing clonotypes to query.
                                             Must have 'v_gene', 'CDR3', and 'j_gene' columns.
            format (str): The naming convention of the V and J genes in `query_clonotypes`.
                          Must be 'Adaptive' or 'IMGT'.
            num_mismatches (int, optional): The maximum allowed Levenshtein distance
                                            for CDR3 sequences to be considered a match. Defaults to 1.

        Returns:
            pd.DataFrame: A table containing matching clonotypes from the database,
                          merged with the original query clonotype data, and including
                          a 'distance' column showing the Levenshtein distance between
                          the matched CDR3s.

        Raises:
            ValueError: If an unsupported format is provided, or if `query_clonotypes`
                        is missing required columns.
            ImportError: If the Levenshtein library or builder functions are not available.
        """
        if format not in ['Adaptive', 'IMGT']:
            raise ValueError(f"Unsupported format: '{format}'. The function supports V and J gene "
                             "name formatting in either 'IMGT' or 'Adaptive' gene formatting.")
        
        # Ensure query_clonotypes has the mandatory columns
        required_query_cols = ['v_gene', 'CDR3', 'j_gene']
        if not all(col in query_clonotypes.columns for col in required_query_cols):
            raise ValueError(
                f"The 'query_clonotypes' DataFrame must contain the following columns: "
                f"{required_query_cols}"
            )

        work_database = self._database.copy(deep=True)
        
        # Apply gene format conversion if necessary
        if self._DATA_TYPE == 'TRB' and format == 'IMGT':
            if 'adaptive_to_imgt_builder' not in globals() and 'adaptive_to_imgt_builder' not in locals():
                 raise NotImplementedError("adaptive_to_imgt_builder() is not defined or imported. "
                                           "Please implement or import this function.")
            Adaptive2IMGT = adaptive_to_imgt_builder(path2database=self._database_path)
            work_database['v_gene'] = work_database['v_gene'].replace(Adaptive2IMGT)
            work_database['j_gene'] = work_database['j_gene'].replace(Adaptive2IMGT)
        
        if self._DATA_TYPE == 'TRA' and format == 'Adaptive':
            if 'imgt_to_adaptive_builder' not in globals() and 'imgt_to_adaptive_builder' not in locals():
                raise NotImplementedError("imgt_to_adaptive_builder() is not defined or imported. "
                                          "Please implement or import this function.")
            IMGT2Adaptive = imgt_to_adaptive_builder(path2database=self._database_path)
            work_database['v_gene'] = work_database['v_gene'].replace(IMGT2Adaptive)
            work_database['j_gene'] = work_database['j_gene'].replace(IMGT2Adaptive)
        
        # Rename columns to avoid ambiguity after cross-merge
        work_database_renamed = work_database.rename(
            columns={'v_gene':'DB_v_gene', 'j_gene':'DB_j_gene', 'CDR3': 'DB_CDR3'}
        )
        # Create a copy of query_clonotypes to avoid modifying the original in-place
        query_clonotypes_renamed = query_clonotypes.copy(deep=True).rename(
            columns={'v_gene':'QW_v_gene', 'j_gene':'QW_j_gene', 'CDR3': 'QW_CDR3'}
        )

        # Performing a cross merge can be very memory intensive for large datasets.
        # Consider optimizing this if performance becomes an issue, e.g., by iterating
        # through query_clonotypes and calling get_restricted_clonotypes_matching_a_clonotype
        # or using a more targeted merge strategy.
        print("Processing your request... We are now cross-referencing your input with our databases."+
              f" This step can take more time if you have a large number of sequences. Current number of input sequences: {query_clonotypes_renamed.shape[0]}")
        print(f"Started crossing at : {time.ctime()}")
        crossed_dataset = work_database_renamed.merge(query_clonotypes_renamed, how='cross')
        print(f"Finished crossing at : {time.ctime()}")
        # Filter for matching V and J genes
        crossed_dataset = crossed_dataset.loc[
            (crossed_dataset.DB_v_gene == crossed_dataset.QW_v_gene) &
            (crossed_dataset.DB_j_gene == crossed_dataset.QW_j_gene)
        ]

        if crossed_dataset.shape[0] == 0:
            print(f"No match found in the database. This is because there are no matching V and/or J gene segments between"+
                  " your query sequences and the database entries. Please ensure you have used the correct naming convention,"+
                  " as the current version supports either IMGT or Adaptive conventions")

        # Ensure CDR3 columns exist before trying to calculate distance
        if 'DB_CDR3' not in crossed_dataset.columns or 'QW_CDR3' not in crossed_dataset.columns:
            raise KeyError("Missing 'DB_CDR3' or 'QW_CDR3' columns after merge.")

        # Using Levenshtein.distance requires the library to be imported
        if not hasattr(Levenshtein, 'distance'):
            raise ImportError("Levenshtein library is required but 'distance' method not found.")

        # Calculate Levenshtein distance
        crossed_dataset['distance'] = [
            Levenshtein.distance(str(db_cdr3), str(qw_cdr3))
            for db_cdr3, qw_cdr3 in zip(crossed_dataset.DB_CDR3, crossed_dataset.QW_CDR3)
        ]
        
        target_clonotypes = crossed_dataset.loc[crossed_dataset['distance'] <= num_mismatches].copy(deep=True)
        
        # Explicitly delete large DataFrames to free memory, especially important for large datasets
        del work_database
        del work_database_renamed
        del query_clonotypes_renamed
        del crossed_dataset

        return target_clonotypes