import pandas as pd 
from gene_id import Gene_IDs
import plotting as my_plots
import utilities as ut
from atac_signal import ATAC_signal
from gene_sets import Gene_sets

if __name__=='__main__':
    # david_df = pd.read_csv('tables/MutAccum_Number_Mutations_per_Gene.csv')
    # david_list = list(david_df['Gene'])
    # dic_david = {"David's mutated genes":david_list}

    if "gs" not in locals():
        gs = Gene_sets()

    if "exp1" not in locals():
        exp1 = ATAC_signal("exp1")

    if "exp_mss" not in locals():
        exp_mss = ATAC_signal("exp_metsetset")
    
    if "exp_hrde1" not in locals():
        exp_hrde1 = ATAC_signal("exp_hrde_guy")
    
    # my_plots.plot_groups_signals(exp1, dic_david, mean_flag=True, bootstrap=True, boot_size=len(david_list), boot_iters=2000)
    # my_plots.plot_gene_atac_signal_distribution(exp1, 'oma-1', mean_flag=False)
    
    # my_plots.plot_signal_gene(exp1, 'oma-1', drop_rep=0)
    # my_plots.plot_groups_signals(exp1, mean_flag=True)
    # my_plots.plot_gene_atac_signal_distribution(exp1, 'oma-1')

    ################ hrde-1 nearbys and bootstrap
    # hrde1_kennedy = gs.get_list('hrde-1-Kennedy')
    # hrde_FC_sig = gs.get_list('mRNA_isSig')
    # hrde_up = gs.get_list('mRNA_log2_FC', thresh=0)
    # hrde_up_sig = ut.intersect_lists(hrde_FC_sig, hrde_up)
    # hrde_regulated = ut.intersect_lists(hrde_up_sig, hrde1_kennedy)

    # hrde1_nearby_up, hrde1_nearby_down = ut.get_nearby_genes_list(hrde_regulated, 2000) # len 82, len 33

    # dic_hrde = {'hrde-1 regulated':hrde_regulated, 'hrde-1 upstream':hrde1_nearby_up, 'hrde-1 downstream':hrde1_nearby_down}

    # my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde-1 upstream':hrde1_nearby_up}, mean_flag=True)
    # my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde-1 upstream':hrde1_nearby_up}, mean_flag=False, var_type='none')

    # my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde-1 upstream':hrde1_nearby_up}, mean_flag=False, bootstrap=True, boot_size=75, var_type='none')

    ###### testing z-score
    my_plots.plot_groups_signals(exp1, mean_flag=False, var_type='none')
    my_plots.plot_groups_signals(exp1, mean_flag=False, var_type='none', zscore_signal=True)


    

