"""
@author: Hesham ElAbd
@brief: performing HLA imputation on a given TCR repertoire or a collection of repertoires
"""
## Load the modules 
#-------------------
import os
import time
import pickle
import pandas as pd
import numpy as np
from typing import Dict
from tqdm import tqdm

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
                with open(f'{path2models}/TRA_weights.pickle','rb') as reader:
                    self.weighted_HLA_associated_clonotypes = pickle.load(reader)
            
            case 'TRB':
                with open(f'{path2models}/TRB_models.pickle','rb') as reader:
                    self.hla_models = pickle.load(reader)
                with open(f'{path2models}/TRB_weights.pickle','rb') as reader:
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



def translate_to_tsv(input_directory:Dict[str,Dict[str,np.ndarray]]):pass