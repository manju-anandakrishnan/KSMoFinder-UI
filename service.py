import dao
from core import KinaseInferenceCore
from dao import DatabaseAccess as db_access, KinaseSubstrateDAO as ks_dao
import streamlit as st
import pandas as pd

# @st.cache_data
# def get_phosphoproteins():
#     kg_substrates = dao.get_kg_substrates()
#     substrate_proteins = {x[0] for x in kg_substrates}
#     return substrate_proteins

@st.cache_data
def get_substrate_genes():
    kg_substrates = dao.get_kg_substrates()
    substrate_genes = {x[0] for x in kg_substrates}
    return substrate_genes


@st.cache_data
def get_protein_phosphosites():
    kg_substrates = dao.get_kg_substrates()
    protein_sites = {}
    for gene,p,site,m in kg_substrates:
        protein_sites.setdefault(gene,{})
        protein_sites[gene].setdefault(p,[])
        protein_sites[gene][p].append((site,m))
    return protein_sites

class KinaseSubstrateService:

    def __init__(self):
        self.dao = dao.KinaseSubstrateDAO(db_access.get_driver())

    def get_kinase_substrate_links(self,kinase,substrate,intermediate_max):
        return self.dao.get_ks_links(kinase,substrate,intermediate_max)
    
class KinaseInferenceService:

    def __init__(self):
        pass

    def get_dysregulated_kinases(self, df, log2FC=None, pval=None):
        ki_core = KinaseInferenceCore(df,log2FC,pval)
        return ki_core.get_enrichment_results()



