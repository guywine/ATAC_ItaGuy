from mRNA_gonads import Table_mRNA
from Ahringer import Ahringer
import utilities as ut 
import pandas as pd 
import matplotlib.pyplot as plt 
from gene_sets import Gene_sets
from atac_signal import ATAC_signal

def create_rna_and_atac_df():

    at = ATAC_signal('exp1')
    atac_sig_gfp_df = at.df_list_to_calc(at.cond1)
    atac_sig_oma1_df = at.df_list_to_calc(at.cond2)

    atac_gfp_norm = ut.normalize_df_cols(atac_sig_gfp_df, (0,100))
    atac_oma1_norm = ut.normalize_df_cols(atac_sig_oma1_df, (0,100))

    atac_GFP_mean = atac_gfp_norm.mean(axis=1)
    atac_oma1_mean = atac_oma1_norm.mean(axis=1)

    rna_all = ut.load_gene_expression_df()
    gfp_df = pd.DataFrame({'atac_gfp_mean':atac_GFP_mean})
    oma1_df = pd.DataFrame({'atac_oma1_mean':atac_oma1_mean})
    atac_means = pd.concat([gfp_df, oma1_df], axis=1)

    atac_means = atac_means[~atac_means.index.duplicated(keep='first')]

    rna_atac_table = pd.concat([rna_all, atac_means], axis=1)
    
    # still need to add ahringer atac
    ar_atac = Ahringer().atac.loc[:,'Germline']

if __name__=='__main__':
    # ut.print_gene_expression('oma-1')
    # ut.print_gene_expression('oma-2')
    ut.print_gene_expression('WBGene00000001')

    # gs = Gene_sets()
    # hrde1_kennedy = gs.big_table.loc[:,'hrde-1-Kennedy']

    m = Table_mRNA()
    fc_over4 = 


    








