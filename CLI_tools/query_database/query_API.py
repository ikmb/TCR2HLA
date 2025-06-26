"""
@author: Hesham ElAnd
@brief: The file contain an interface class that can be used to query the HLA-TCR database,
either it can query the HLA restriction associated with a given HLA or the TCR associated with a given HLA.
"""
## Load the modules
import pandas as pd
import numpy as np
import Levenshtein
from CLI_tools.query_database.gene_formating_converter import imgt_to_adaptive_builder, adaptive_to_imgt_builder


## Define the class that acts as an API
#--------------------------------------
class HLATCRDBInterface: 
    def __init__(self,database_dir: str, database_name:str):
        """
        database_dir (str, mandatory): the path to where the database is locatted
        database_name (str, madatory): the name of the database, can be either TRB or TRA depending on database to load
        """
        match database_name:
            case 'TRB_DB':
                self._database = pd.read_parquet(f'{database_dir}/databases/TRB_database.parquet')
                self._DATA_TYPE = 'TRB'
            case 'TRA_DB':
                self._database = pd.read_parquet(f'{database_dir}/databases/TRA_database.parquet')
                self._DATA_TYPE = 'TRA'
            case _:
                raise ValueError("unknown database, the model only supports TRB_DB which is the TRB-HLA database and TRA_DB which is the TRA-HLA database")                
        return

    def get_clonotypes_associated_to_alleles(self,allele_name:str)->pd.DataFrame: 
        """
        The function recieves an allele name and return restricted clonotypes associated with this alleles

        allele_name (str): the allele name shall be provided as follow, e.g. 'A*02:01' where A represent the loci name, '*' is a sperator and '02:01' is the allele name. 
        For DQ and DP alleles, the following format can be used, DQ*01:02+04:01 or DP*01:01+04:01. 

        returns (pd.DataFrame): the fuction returns a table that contains all the clonotypes associated with a given allele, the table is composite of three column
        1. V genes --> which is either provided in the IMGT or Adaptive format
        2. J genes --> which is either provided in the IMGT or Adaptive format
        3. CDR3 --> the amino acid sequence of the CDR3
        """
        loci_name,allele = allele_name.split('*')
        print('Checking your input ...')
        if loci_name not in self._database.loci.unique():
            raise ValueError(f"Locus: {loci_name} is not support, only these loci are supported: {self._database.loci.unique()}")
        if allele not in self.loc[self.loci==loci_name].allele_name.unique(): 
            raise ValueError(f"Allele: {allele} from the Loci: {loci_name} is not support, only these alleles from this loci are supported: {self.loc[self.loci==loci_name].allele_name.unique()}")
        print('checks finished ...supported allele')
        return self._database.loc[
            (self._database.loci==loci_name) & (self._database.allele_name==allele)
        ]
    
    def get_restricted_clonotypes_matching_a_clonotype(self,
        v_gene:str, j_gene:str, format:str, CDR3_amino_acids:str, num_mismatches:int=1)->pd.DataFrame:
        """
        
        
        
        """
        if format not in ['Adaptive', 'IMGT']:
            raise ValueError(f"Unsupported format, the current function support V and J name formatting in either the 'IMGT' or 'Adaptive' gene formatting")
        
        if self._DATA_TYPE=='TRB' and format=='IMGT': 
            Adaptive2IMGT = adaptive_to_imgt_builder()
            work_database = self._database.copy(deep=True)
            work_database['v_gene'].replace(Adaptive2IMGT,inplace=True)
            work_database['j_gene'].replace(Adaptive2IMGT,inplace=True)
        else:
            work_database = self._database.copy(deep=True)
        
        if self._DATA_TYPE=='TRA' and format=='Adaptive':
            IMGT2Adaptive = imgt_to_adaptive_builder()
            work_database = self._database.copy(deep=True)
            work_database['v_gene'].replace(IMGT2Adaptive,inplace=True)
            work_database['j_gene'].replace(IMGT2Adaptive,inplace=True)
        else: 
            work_database = self._database.copy(deep=True)
        
        if v_gene not in work_database.v_gene.unique():
            raise ValueError(f"The provide V-gene: {v_gene} is not defined in the database, alternativly is might not be in the correct format. Current format is {format}")
        
        if j_gene not in work_database.j_gene.unique():
            raise ValueError(f"The provide V-gene: {j_gene} is not defined in the database, alternativly is might not be in the correct format. Current format is {format}")
        
        candiate_table = work_database.loc[
            (work_database.v_gene==v_gene) & (work_database.j_gene==j_gene)
        ].copy(deep=True)

        del work_database
        candiate_table['distance']=[Levenshtein.distance(c, CDR3_amino_acids) for c in candiate_table.CDR3]
        target_clonotypes = candiate_table.loc[candiate_table['distance']<=num_mismatches].copy(deep=True)
        del candiate_table
        return target_clonotypes

    def get_restricted_clonotypes_matching_clonotypes(self,query_clonotypes:pd.DataFrame, format:str,  num_mismatches:int=1)->pd.DataFrame:
        """
        """
        if format not in ['Adaptive', 'IMGT']:
            raise ValueError(f"Unsupported format, the current function support V and J name formatting in either the 'IMGT' or 'Adaptive' gene formatting")
        
        if self._DATA_TYPE=='TRB' and format=='IMGT': 
            Adaptive2IMGT = adaptive_to_imgt_builder()
            work_database = self._database.copy(deep=True)
            work_database['v_gene'].replace(Adaptive2IMGT,inplace=True)
            work_database['j_gene'].replace(Adaptive2IMGT,inplace=True)
        else:
            work_database = self._database.copy(deep=True)
        
        if self._DATA_TYPE=='TRA' and format=='Adaptive':
            IMGT2Adaptive = imgt_to_adaptive_builder()
            work_database = self._database.copy(deep=True)
            work_database['v_gene'].replace(IMGT2Adaptive,inplace=True)
            work_database['j_gene'].replace(IMGT2Adaptive,inplace=True)
        else: 
            work_database = self._database.copy(deep=True)
        
        ## merging the two datasets
        work_database.rename(columns={'v_gene':'DB_v_gene', 'j_gene':'DB_j_gene', 'CDR3': 'DB_CDR3'},inplace=True)
        query_clonotypes.rename(columns={'v_gene':'QW_v_gene', 'j_gene':'QW_j_gene', 'CDR3': 'QW_CDR3'},inplace=True)

        # doing the crossing 
        crossed_dataset = work_database.merge(query_clonotypes,how='cross')
        crossed_dataset = crossed_dataset.loc[crossed_dataset.DB_v_gene==crossed_dataset.QW_v_gene]
        crossed_dataset = crossed_dataset.loc[crossed_dataset.DB_j_gene==crossed_dataset.QW_j_gene]
        crossed_dataset['distance'] = [Levenshtein.distance(e, j) for e,j in zip(crossed_dataset.DB_CDR3,crossed_dataset.QW_CDR3)]
        target_clonotypes = crossed_dataset.loc[crossed_dataset['distance']<=num_mismatches].copy(deep=True)
        del crossed_dataset
        return target_clonotypes