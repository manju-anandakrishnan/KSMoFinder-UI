import streamlit as st
import pandas as pd
from scipy.stats import fisher_exact, false_discovery_control

class KinaseInferenceCore:

    '''
    This method loads the background kinase-substrate dataset in the cache
    '''
    @st.cache_data
    def get_background():
        bg_kss_df = pd.read_feather('data/bg_kss.feather')
        bg_sites = bg_kss_df.copy()
        bg_sites.drop(columns=['kinase'],inplace=True,axis=1)
        bg_sites.drop_duplicates(inplace=True)
        bg_sites_count = bg_sites.shape[0]
        return bg_kss_df, bg_sites_count

    def __init__(self,fg_df,log2FC,pval):
        if log2FC:
            log2FC = float(log2FC)
            fg_df = fg_df[fg_df['log2FC'].abs() > log2FC]
        if pval:
            pval = float(pval)
            fg_df = fg_df[fg_df['pval'] < pval]
        total_fg_sites = fg_df.shape[0]
        fg_df.drop_duplicates(inplace=True)
        fg_sites_count = fg_df.shape[0]
        print(f'Count of duplicates in the input::{total_fg_sites-fg_sites_count}')
        print(f'Remaining total sites in the input::{fg_sites_count}')
        self.fg_sites_count = fg_sites_count
        self.fg_df = fg_df

    def _merge_fg_bg(self,bg_kss_df):

        fg_merged = self.fg_df.merge(bg_kss_df,how='outer',left_on=['protein','motif'],right_on=['substrate','motif'],indicator=True)
        fg_df = fg_merged[fg_merged._merge.isin(['both'])].copy()
        bg_df = fg_merged[fg_merged._merge=='right_only'].copy()

        self.fg_df_kinase_cnt_dict = fg_df['kinase'].value_counts().to_dict()
        self.bg_df_kinase_cnt_dict = bg_df['kinase'].value_counts().to_dict()

    def get_contingency_tbl_cnts(self,kinase):
        a_cnt = self.fg_df_kinase_cnt_dict.get(kinase,0)
        c_cnt = self.fg_sites_count - a_cnt
        b_cnt = self.bg_df_kinase_cnt_dict.get(kinase,0)
        d_cnt = self.bg_sites_count - a_cnt - c_cnt - b_cnt
        return a_cnt, b_cnt, c_cnt, d_cnt


    def get_enrichment_results(self):

        bg_kss_df, self.bg_sites_count = KinaseInferenceCore.get_background()
        self._merge_fg_bg(bg_kss_df)

        enrich_results = {}

        kinases = []
        kinases_pval = []

        for bg_kinase in bg_kss_df['kinase'].unique():

            a_cnt, b_cnt, c_cnt, d_cnt = self.get_contingency_tbl_cnts(bg_kinase)

            if any(x < 0 for x in (a_cnt, b_cnt, c_cnt, d_cnt)) or any(x ==0 for x in (b_cnt+d_cnt, a_cnt+c_cnt)):continue
            res = fisher_exact([[a_cnt,b_cnt],[c_cnt,d_cnt]], alternative='greater')
            p_value = res.pvalue
            enrich_results[bg_kinase] = {'pval':p_value,'adj_pval':None,
                                        'predicted_sites':a_cnt}
            kinases.append(bg_kinase)
            kinases_pval.append(p_value) 

        adjusted_p_values = false_discovery_control(kinases_pval, method='bh')
        for idx,adj_pval in enumerate(adjusted_p_values):
            kinase = kinases[idx]
            enrich_results[kinase]['adj_pval'] = adj_pval

        return enrich_results, self.fg_df.shape[0]


