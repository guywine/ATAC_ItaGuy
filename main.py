import pandas as pd 

from gene_id import Gene_IDs
import plotting as my_plots
import utilities as ut
from atac_signal import ATAC_signal
from gene_sets import Gene_sets
import calc_signals as cas

def where_are_the_mutation():
    david_df = pd.read_csv('tables/MutAccum_Number_Mutations_per_Gene.csv')
    david_list = list(david_df['Gene'])
    dic_david = {"David's mutated genes":david_list}
    my_plots.plot_groups_signals(exp1, dic_david, mean_flag=True, bootstrap=True, boot_size=len(david_list), boot_iters=2000)

def get_hrde_regulated():
    hrde1_kennedy = gs.get_list('hrde-1-Kennedy')
    hrde_FC_sig = gs.get_list('mRNA_isSig')
    hrde_up = gs.get_list('mRNA_log2_FC', thresh=0)
    hrde_up_sig = ut.intersect_lists(hrde_FC_sig, hrde_up)
    hrde_regulated = ut.intersect_lists(hrde_up_sig, hrde1_kennedy)
    return hrde_regulated


if __name__=='__main__':
    oma_wbid = 'WBGene00003864'
    gid = Gene_IDs()

    if "gs" not in locals():
        gs = Gene_sets()

    if "exp1" not in locals():
        exp1 = ATAC_signal("exp1")

    # if "exp_mss" not in locals():
    #     exp_mss = ATAC_signal("exp_metsetset")
    
    if "exp_hrde1" not in locals():
        exp_hrde1 = ATAC_signal("exp_hrde_guy")
    
    hrde_regulated = get_hrde_regulated()

    ### 15 / 151 missing
    hrde_reg_intersected = ut.intersect_lists(hrde_regulated, exp_hrde1.scores1.index)

    my_plots.scatter_genes_both_conds(exp_hrde1, marked_list=hrde_reg_intersected, shown_value='score')

    

    

