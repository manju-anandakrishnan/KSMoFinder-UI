import service
from service import KinaseSubstrateService, KinaseInferenceService
import streamlit as st
import pandas as pd
import os

PRED_DIR = 'data/ksmo_predictions'
def get_predictions(substrate,motif,prob):
    file_path = os.path.join(PRED_DIR,f'{substrate}.feather')
    pred_df = pd.read_feather(file_path)
    pred_df = pred_df[(pred_df['substrate_motif']==substrate+'_'+motif)
                      & (pred_df['ksf2_pred']>=prob)].copy()
    pred_df = pred_df[['kinase','ksf2_pred']].copy()
    pred_df['Kinase Gene'] =  pred_df['kinase'].apply(lambda x:get_protein_gene_dict()[x])

    pred_df.rename({'kinase':'Kinase UniProt ID',
                    'ksf2_pred':'Prediction Probability'}, inplace=True, axis=1)
    pred_df = pred_df.reset_index(drop=True)
    pred_df.drop_duplicates(inplace=True)
    pred_df.sort_values(by='Prediction Probability',ascending=False,inplace=True)
    pred_df = pred_df[['Kinase Gene','Kinase UniProt ID','Prediction Probability']]
    return pred_df

@st.cache_data
def get_substrate_genes():
    substrate_genes = ['']
    substrate_genes.extend(service.get_substrate_genes())
    substrate_genes.sort()
    return substrate_genes

@st.cache_data
def get_protein_gene_dict():
    df = pd.read_csv('data/gene_protein.csv',sep='|')
    protein_gene_dict = dict(zip(df['protein'].to_list(),df['gene'].to_list()))
    return protein_gene_dict

class Controller:

    def __init__(self):
        self.ks_service = KinaseSubstrateService()
        self.ki_service = KinaseInferenceService()

    def get_proteins(self,gene):
        protein_phosphosites = service.get_protein_phosphosites()
        protein_phosphosites = protein_phosphosites[gene]
        proteins = protein_phosphosites.keys()
        return proteins

    def get_phosphosite(self,gene, protein):
        phosphosites = service.get_protein_phosphosites()
        protein_phosphosites = phosphosites[gene][protein]
        site_motif_pairs = {site:m for site,m in protein_phosphosites}
        return site_motif_pairs

    def get_predicted_kinases(self,substrate_protein,motif,prediction_prob_cutoff=0):
        return get_predictions(substrate_protein,motif,prediction_prob_cutoff)
    
    def get_kinase_substrate_links(self,kinase,substrate,intermediate_max):
        return self.ks_service.get_kinase_substrate_links(kinase,substrate,intermediate_max)
    
    def get_dysregulated_kinases(self, df, log2FC=None, pval=None):
        input_fg_df_shape = df.shape[0]
        enrichment_results, fg_df_shape = self.ki_service.get_dysregulated_kinases(df,log2FC,pval)
        result_df = pd.DataFrame.from_dict(enrichment_results, orient='index')
        result_df = result_df.reset_index()
        result_df = result_df.rename(columns={'index': 'Kinase UniProt ID'})
        result_df['Kinase Gene'] =  result_df['Kinase UniProt ID'].apply(lambda x:get_protein_gene_dict()[x])
        result_df.sort_values(by=['adj_pval'],inplace=True)
        print(result_df.columns)
        result_df = result_df[['Kinase Gene','Kinase UniProt ID','predicted_sites','pval','adj_pval']]
        return input_fg_df_shape, result_df, fg_df_shape
    
