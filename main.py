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

    if "gs" not in locals():
        gs = Gene_sets()

    if "exp1" not in locals():
        exp1 = ATAC_signal("exp1")

    if "exp_mss" not in locals():
        exp_mss = ATAC_signal("exp_metsetset")
    
    if "exp_hrde1" not in locals():
        exp_hrde1 = ATAC_signal("exp_hrde_guy")
    
    hrde_regulated = get_hrde_regulated()

    ################ hrde-1 nearbys and bootstrap
    hrde1_kennedy = gs.get_list('hrde-1-Kennedy')
    hrde_FC_sig = gs.get_list('mRNA_isSig')
    hrde_up = gs.get_list('mRNA_log2_FC', thresh=0)
    hrde_up_sig = ut.intersect_lists(hrde_FC_sig, hrde_up)
    hrde_down = gs.get_list('mRNA_log2_FC', thresh=0, bottom=True)
    hrde_down_sig = ut.intersect_lists(hrde_FC_sig, hrde_down)
    hrde_regulated = ut.intersect_lists(hrde_up_sig, hrde1_kennedy)

    hrde1_nearby_up, hrde1_nearby_down = ut.get_nearby_genes_list(hrde_regulated, 2000) # len 75, len 28
    hrde1_nearby_up_1200, hrde1_nearby_down_1200 = ut.get_nearby_genes_list(hrde_regulated, 1200) # len 49, len 7

    hrde_dic = {'hrde1_kennedy':hrde1_kennedy, 'hrde_reg':hrde_regulated, 'hrde-1 upstream':hrde1_nearby_up, 'hrde down sig':hrde_down_sig}
    
    my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde-1 upstream':hrde1_nearby_up}, mean_flag=True)

    


    my_plots.plot_fc_groups_dots(exp_hrde1, hrde_dic, mean_flag=True)
    my_plots.plot_fc_groups_dots(exp_hrde1, hrde_dic, mean_flag=False)

    # my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde-1 upstream':hrde1_nearby_up}, mean_flag=False, var_type='none')

    my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde-1 upstream':hrde1_nearby_up}, mean_flag=True, bootstrap=True, boot_size=75)
    # my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde-1 upstream':hrde1_nearby_up}, mean_flag=False, var_type='none', bootstrap=True, boot_size=75)

    ###### testing z-score
    # my_plots.plot_groups_signals(exp1, mean_flag=False, var_type='none')
    # my_plots.plot_groups_signals(exp1, mean_flag=False, var_type='none', zscore_signal=True)
    # wbid_oma1 = exp1.gid.to_wbid('oma-1')
    # dic_oma_gfp = {'oma-1':[wbid_oma1],'GFP':['GFP']}
    # my_plots.plot_groups_signals(exp1, groups_dic = dic_oma_gfp, mean_flag=True, zscore_signal=True)
    # my_plots.plot_groups_signals(exp1, groups_dic = dic_oma_gfp, mean_flag=True, zscore_signal=True, plot_range=(-800,600))

    ######## hrde-1 zscore upstream 2000
    # my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde-1 upstream':hrde1_nearby_up}, mean_flag=True, zscore_signal=True)
    # my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde-1 upstream':hrde1_nearby_up}, mean_flag=False, var_type='none', zscore_signal=True)


    ######## hrde-1 zscore upstream 1200
    # my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde-1 upstream 1200':hrde1_nearby_up_1200}, mean_flag=True, zscore_signal=True)
    # my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde-1 upstream 1200':hrde1_nearby_up_1200}, mean_flag=False, var_type='none', zscore_signal=True)
    
    

    

