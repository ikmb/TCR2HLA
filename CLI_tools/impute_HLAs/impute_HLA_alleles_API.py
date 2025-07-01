"""
@author: Hesham ElAbd
@brief: performing HLA imputation on a given TCR repertoire or a collection of repertoires
"""
## Load the modules 
#-------------------
import gzip
import time
import pickle
import pandas as pd
import numpy as np
from typing import Dict
from tqdm import tqdm
from scipy.special import softmax

class HLAPredictor:
    def __init__(self,path2models:str,model:str) -> None:
        """The class is used for predicting the HLA of a sample

        Args:
            path2trained_models (str, mandatory): the path where the models and weights are definestored
            loci (str; mandatory): which model to use for imputation, can be either TRA or TRB
            
        """
        match model:
            
            case 'TRA': 
                with open(f'{path2models}/TRA_models.pickle','rb') as reader:
                    self.hla_models = pickle.load(reader)
                with gzip.open(f'{path2models}/TRA_weights.pickle.gz','rb') as reader:
                    self.weighted_HLA_associated_clonotypes = pickle.load(reader)
            
            case 'TRB':
                with open(f'{path2models}/TRB_models.pickle','rb') as reader:
                    self.hla_models = pickle.load(reader)
                with gzip.open(f'{path2models}/TRB_weights.pickle.gz','rb') as reader:
                    self.weighted_HLA_associated_clonotypes = pickle.load(reader)
            
            case _:
                raise ValueError(f"The current version support imputation based on two models, either 'TRA' or 'TRB' models, however, your model is: {model}")
        return
    
    def predict_samples_HLA(self,samples_TCR:pd.DataFrame)->Dict[str,float]:
        """Predict the HLA of a sample

        Args:
            samples_TCR (pd.DataFrame): A pandas data frame that contains a list of clonotypes
            and their expansion

        Returns:
            Dict[str,float]: a dictionary contain the probability score that the samples has said HLA allele
        """
        samples_TCR=samples_TCR.copy(deep=True)
        samples_TCR['bioId'] = [f'{c}+{v}+{j}' for c,v,j in zip(samples_TCR['CDR3'],
                                                                samples_TCR['v_gene'],
                                                                samples_TCR['j_gene'])]
        sequencing_depth=np.log10(1+samples_TCR.shape[0])
        predicted_alleles = dict()
        for allele, clonotype_weights in self.weighted_HLA_associated_clonotypes.items():
            # 1. Calculate the weighted expansion of the allele associated clonotypes |
            #-------------------------------------------------------------------------|
            clonotypes = [e[0] for e in clonotype_weights]
            clonotype_to_weight_dc = {c:abs(w) for (c,w) in clonotype_weights}
            target_tcr_table = samples_TCR.loc[samples_TCR.bioId.isin(clonotypes)].copy(deep=True)
            target_tcr_table['weight']=[clonotype_to_weight_dc[b] for b in target_tcr_table.bioId]
            target_tcr_table['trust_weighted_expansion']=target_tcr_table['weight']*target_tcr_table['count']
            weighted_expansion = np.log10( 1+ np.sum(target_tcr_table['trust_weighted_expansion']))
            
            # 2. Prepare the input array to the model and predict the HLA |
            #-------------------------------------------------------------|
            input_array=np.array([weighted_expansion,sequencing_depth]).reshape(1,2)
            allele_score=self.hla_models[allele].predict(input_array)

            # 3. store the scores into the results dict
            #-------------------------------------------------------------|
            predicted_alleles[allele]=allele_score
        # Return the results 
        #-------------------
        return predicted_alleles
    
    def predict_samples_HLA_prob(self,samples_TCR:pd.DataFrame)->Dict[str,float]:
        """Predict the HLA of a sample

        Args:
            samples_TCR (pd.DataFrame): A pandas data frame that contains a list of clonotypes
            and their expansion

        Returns:
            Dict[str,float]: a dictionary contain the probability score that the samples has said HLA allele
        """
        samples_TCR=samples_TCR.copy(deep=True)
        samples_TCR['bioId'] = [f'{c}+{v}+{j}' for c,v,j in zip(samples_TCR['CDR3'],
                                                                samples_TCR['v_gene'],
                                                                samples_TCR['j_gene'])]
        sequencing_depth=np.log10(1+samples_TCR.shape[0])
        predicted_alleles = dict()
        for allele, clonotype_weights in self.weighted_HLA_associated_clonotypes.items():
            # 1. Calculate the weighted expansion of the allele associated clonotypes |
            #-------------------------------------------------------------------------|
            clonotypes = [e[0] for e in clonotype_weights]
            clonotype_to_weight_dc = {c:abs(w) for (c,w) in clonotype_weights}
            target_tcr_table = samples_TCR.loc[samples_TCR.bioId.isin(clonotypes)].copy(deep=True)
            target_tcr_table['weight']=[clonotype_to_weight_dc[b] for b in target_tcr_table.bioId]
            target_tcr_table['trust_weighted_expansion']=target_tcr_table['weight']*target_tcr_table['count']
            weighted_expansion = np.log10( 1+ np.sum(target_tcr_table['trust_weighted_expansion']))
            
            # 2. Prepare the input array to the model and predict the HLA |
            #-------------------------------------------------------------|
            input_array=np.array([weighted_expansion,sequencing_depth]).reshape(1,2)
            
            allele_score={c:s for c,s in zip (
               list(self.hla_models[allele].classes_),  list(self.hla_models[allele].predict_proba(input_array).reshape(-1,)), )
                }

            # 3. store the scores into the results dict
            #-------------------------------------------------------------|
            predicted_alleles[allele]=allele_score
        # Return the results 
        #-------------------
        return predicted_alleles
    
    def predict_HLA_for_a_sample_collection(self,collections_TCR:pd.DataFrame)->Dict[str,float]:
        """
        Predict HLA for a sample collection 

        Args:
            samples_TCR (pd.DataFrame): Predict HLA for a sample collection

        Returns:
            Dict[str,float]: Returns HLA for a sample collection
        """
        predicted_HLAs=dict()
        unique_sample_names=list(collections_TCR.sample_name.unique())
        for sample_name in tqdm(unique_sample_names): 
            print(f'{time.ctime()}:: Predicting the HLA for sample: {sample_name}')
            samples_TCR = collections_TCR.loc[collections_TCR.sample_name==sample_name]
            HLA_predictions_scores = self.predict_samples_HLA(samples_TCR)
            predicted_HLAs[sample_name]=HLA_predictions_scores
        return predicted_HLAs
    
    def predict_HLA_for_a_sample_collection_prob(self,collections_TCR:pd.DataFrame)->Dict[str,float]:
        """
        Predict HLA for a sample collection 

        Args:
            samples_TCR (pd.DataFrame): Predict HLA for a sample collection

        Returns:
            Dict[str,float]: Returns HLA for a sample collection
        """
        predicted_HLAs=dict()
        unique_sample_names=list(collections_TCR.sample_name.unique())
        for sample_name in tqdm(unique_sample_names): 
            print(f'{time.ctime()}:: Predicting the HLA for sample: {sample_name}')
            samples_TCR = collections_TCR.loc[collections_TCR.sample_name==sample_name]
            HLA_predictions_scores = self.predict_samples_HLA_prob(samples_TCR)
            predicted_HLAs[sample_name]=HLA_predictions_scores
        return predicted_HLAs

def translate_prediction_to_tsv(input_dict:Dict[str,Dict[str,np.ndarray]])->pd.DataFrame:
    sample_names = sorted(list(input_dict.keys()))
    column_names = sorted(list(input_dict[sample_names[0]].keys()))

    # add a sanity check to make sure the same column names
    for sname, res in input_dict.items(): 
        if sorted(list(res.keys())) != column_names:
            raise ValueError(f"Something went wron with decoding the model's output, please contact the developer at: h.elabd@ikmb.uni-kiel.de")
    
    columns = [list() for _ in range(len(column_names)+1)]
    for sname in sample_names:
        columns[0].append(sname)
        for idx, col_name in enumerate(column_names):
            carrier_ship_results =  input_dict[sname][col_name]
            if carrier_ship_results:
                columns[idx+1].append(1)
            else:
                columns[idx+1].append(0)
    database_as_table = pd.DataFrame({cn:c for (cn, c) in zip(['sample_name']+column_names, columns) })
    return database_as_table

def normalise_by_softmax(input_df:pd.DataFrame)->pd.DataFrame:
    """
    Normalize the predictions using a softmax
    """
    # 1. normalise the A alleles 
    A_alleles = [e for e in list(input_df.columns)[1:] if 'A' in e]
    normalised_A_loci = pd.DataFrame(softmax(input_df[A_alleles],axis=1))
    normalised_A_loci.columns=A_alleles
    # 2. normalise the B alleles
    B_alleles = [e for e in list(input_df.columns)[1:] if 'B' in e]
    normalised_B_loci = pd.DataFrame(softmax(input_df[B_alleles],axis=1))
    normalised_B_loci.columns=B_alleles
    # 3. normalise the C alleles
    C_alleles = [e for e in list(input_df.columns)[1:] if 'C' in e]
    normalised_C_loci = pd.DataFrame(softmax(input_df[C_alleles],axis=1))
    normalised_C_loci.columns=C_alleles
    # 4. normalise the DR alleles
    DR_alleles = [e for e in list(input_df.columns)[1:] if 'DR' in e]
    normalised_DR_loci = pd.DataFrame(softmax(input_df[DR_alleles],axis=1))
    normalised_DR_loci.columns=DR_alleles
    # 5. normalise the DR alleles
    DQ_alleles = [e for e in list(input_df.columns)[1:] if 'DQ' in e]
    normalised_DQ_loci = pd.DataFrame(softmax(input_df[DQ_alleles],axis=1))
    normalised_DQ_loci.columns=DQ_alleles
    # 6. normalise the DR alleles
    DP_alleles = [e for e in list(input_df.columns)[1:] if 'DP' in e]
    normalised_DP_loci = pd.DataFrame(softmax(input_df[DP_alleles],axis=1))
    normalised_DP_loci.columns = DP_alleles
    # get the normalised table 
    normalized_probs_per_loci = pd.concat([
        input_df[['sample_name']], 
        normalised_A_loci,
        normalised_B_loci,
        normalised_C_loci,
        normalised_DR_loci,
        normalised_DQ_loci,
        normalised_DP_loci
    ],axis=1)
    return normalized_probs_per_loci

def translate_probabilities_to_tsv(input_dict:Dict[str,Dict[str,np.ndarray]])->pd.DataFrame:
    sample_names = sorted(list(input_dict.keys()))
    column_names = sorted(list(input_dict[sample_names[0]].keys()))

    # add a sanity check to make sure the same column names
    for sname, res in input_dict.items(): 
        if sorted(list(res.keys())) != column_names:
            raise ValueError(f"Something went wron with decoding the model's output, please contact the developer at: h.elabd@ikmb.uni-kiel.de")
    
    columns = [list() for _ in range(len(column_names)+1)]
    for sname in sample_names:
        columns[0].append(sname)
        for idx, col_name in enumerate(column_names):
            columns[idx+1].append(input_dict[sname][col_name][True])
        #print(columns)
    database_as_table = pd.DataFrame({cn:c for (cn, c) in zip(['sample_name']+column_names, columns) })
    return database_as_table