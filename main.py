import pandas as pd 
import random
import matplotlib.pyplot as plt

from gene_id import Gene_IDs
import plotting as my_plots
import utilities as ut
from atac_signal import ATAC_signal
from gene_sets import Gene_sets
import calc_signals as cas
from mRNA_gonads import Table_mRNA


def analyze_single_gene(ATAC_exp, gene_name: str, plot_range=(-1000,1000)):
    '''
    Plots both signal of gene and FC distribution.
    Both seperately and meaned reps.
    '''
    print(f'Analyze gene {gene_name} for experiment {ATAC_exp.exp_name}')
    my_plots.plot_signal_gene(ATAC_exp, gene_name, mean_flag=False, var_type='none', plot_range=plot_range)
    my_plots.plot_signal_gene(ATAC_exp, gene_name, mean_flag=True, var_type='sem', plot_range=plot_range)
    my_plots.plot_gene_atac_signal_distribution(ATAC_exp, gene_name, mean_flag=False, plot_type='violin') 
    my_plots.plot_gene_atac_signal_distribution(ATAC_exp, gene_name, mean_flag=True, plot_type='violin')

def get_grde1_kennedy_fc_tables(hrde1_kennedy):
    '''
    Creates table, each row is a gene:
    columns: mRNA-SX-mean, mRNA-hrde1-mean, ATAC-FC, ATAC-score-SX, ATAC-score-hrde1
    '''
    hrde1_k_fc_exp1 = exp1.fc.loc[hrde1_kennedy_intersected,:]


if __name__=='__main__':
    oma_wbid = 'WBGene00003864'
    gid = Gene_IDs()
    m = Table_mRNA()

    if "gs" not in locals():
        gs = Gene_sets()

    if "exp1" not in locals():
        exp1 = ATAC_signal("exp1")

    if "exp_mss" not in locals():
        exp_mss = ATAC_signal("exp_metsetset")
    
    if "exp_hrde1" not in locals():
        exp_hrde1 = ATAC_signal("exp_hrde_guy")
    
    
    #### get hrde-1 lists
    hrde1_regulated = ut.get_hrde_regulated(gs)
    hrde1_reg_intersected = ut.intersect_lists(hrde1_regulated, exp1.scores1.index) # later: 15 / 151 missing

    hrde1_kennedy = gs.get_list('hrde-1-Kennedy')
    hrde1_kennedy_intersected = ut.intersect_lists(hrde1_kennedy, exp1.fc.index) # later: 68 / 1527 missing
    hrde1_kennedy_mRNA_intersected = ut.intersect_lists(hrde1_kennedy, m.table.index) # later: 47 / 1527 missing

    hrde1_dic = {'hrde1 kennedy':hrde1_kennedy, 'hrde1 regulated':hrde1_regulated}

    highly, lowly = ut.get_highly_lowly()

    hrde1_kennedy_highly = ut.intersect_lists(hrde1_kennedy_intersected, highly)
    ### check the missing genes
    # hrde1_kennedy_missing = ut.intersect_lists(hrde1_kennedy, exp1.scores1.index, 'only first')
    # x = pd.DataFrame({'missing hrde1-Kennedy genes:' : hrde1_kennedy_missing})
    # x.to_csv('missing hrde1 genes', index=False)


    if False:#### verify normalization ######
        my_plots.plot_groups_signals(exp1)
        my_plots.plot_groups_signals(exp_mss)
        my_plots.plot_groups_signals(exp_hrde1)


    if False:
        #### 1.A - Signal along GFP gene in exp1: [# later - rep1 is problematic]
        print('1.A - GFP')
        analyze_single_gene(exp1, 'GFP', plot_range=(-1000,700)) # later - ranks changed

        #### 1.B - Signal along oma-1 gene in exp1: [# later - rep0 is problematic]
        print('1.B - oma-1')
        analyze_single_gene(exp1, 'oma-1', plot_range=(-1000,700)) # later - ranks changed

        #### 1.c [Sup.] - Variability across gene location
        print('1.C - construction')


        #### 2.A - Signal along GFP gene in exp_mss:
        print('2.A')
        analyze_single_gene(exp_mss, 'GFP', plot_range=(-1000,700))

        #### 2.B - Signal along oma-1 gene in exp_mss:
        print('2.B')
        analyze_single_gene(exp_mss, 'oma-1', plot_range=(-1000,700))

        #### 2.C - 2.C. (sup ?? ) ATAC signal of H3K9 target groups WT versus met set set (H3K9//mRNA changing)
        print('2.C - construction')


        #### 3.A - 'define HRDE-1 regulated group'  - scatter plot
        print('3.A - not done')
        my_plots.scatter_genes_both_conds(exp_hrde1, marked_list=hrde1_reg_intersected, shown_value='score', log_flag=True)
        my_plots.scatter_mRNA_both_conds(log_flag=True, marked_list=hrde1_reg_intersected)

        ### scatter atac-score
        my_plots.scatter_genes_both_conds(exp_hrde1, marked_list=hrde1_kennedy_intersected, shown_value='score', log_flag=True)
        my_plots.scatter_mRNA_both_conds(log_flag=True, marked_list=hrde1_kennedy_intersected)

        
        ### exp1:
        my_plots.scatter_genes_both_conds(exp1, marked_list=hrde1_reg_intersected, shown_value='score', log_flag=True)
        my_plots.scatter_genes_both_conds(exp1, marked_list=hrde1_kennedy_intersected, shown_value='score', log_flag=True)



        #### 3.B - HRDE-1 regulated targets get open:
        print('3.B - hrde-1')
        
        my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde1 kennedy':hrde1_kennedy})
        my_plots.plot_groups_signals(exp_hrde1, groups_dic=hrde1_dic)
        my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde1 kennedy':hrde1_kennedy}, mean_flag=True, var_type='sem')


        ## bootstrap - Mean FC score of the hrde-1 kennedy group (FC = SX / hrde-1)
        ##### bootstrap atac-FC
        print('exp hrde1 - hrde kennedy, bootstrap atac-FC')
        print('r1:')
        cas.bootstrap_group_score_fc_histogram(exp_hrde1.fc.iloc[:,0], hrde1_kennedy) # later - WTF??? Some biological thing
        print('r2:')
        cas.bootstrap_group_score_fc_histogram(exp_hrde1.fc.iloc[:,1], hrde1_kennedy)# later - WTF???
        print('r3:')
        cas.bootstrap_group_score_fc_histogram(exp_hrde1.fc.iloc[:,2], hrde1_kennedy) # later - correct!
        print('exp hrde1 - hrde kennedy mean')
        cas.bootstrap_group_score_fc_histogram(exp_hrde1.fc.mean(axis=1), hrde1_kennedy) # later - opposite!


        ##### bootstrap mRNA FC

        ####### bootstrap for mRNA levels
        print('exp hrde1 - hrde kennedy, bootstrap mRNA-FC')
        print('r1:')
        cas.bootstrap_group_score_fc_histogram(m.fc.iloc[:,0], hrde1_kennedy)
        print('r2:')
        cas.bootstrap_group_score_fc_histogram(m.fc.iloc[:,1], hrde1_kennedy) 
        print('r3:')
        cas.bootstrap_group_score_fc_histogram(m.fc.iloc[:,2], hrde1_kennedy) 

        ## bootstrap - Mean FC score of the hrde-1 regulated group (FC = SX / hrde-1)
        print('exp hrde1 - hrde regulated')
        cas.bootstrap_group_score_fc_histogram(exp_hrde1.fc.iloc[:,0], hrde1_regulated) # ok
        cas.bootstrap_group_score_fc_histogram(exp_hrde1.fc.iloc[:,1], hrde1_regulated) # ok
        cas.bootstrap_group_score_fc_histogram(exp_hrde1.fc.iloc[:,2], hrde1_regulated) # ok
        print('exp hrde1 - hrde regulated mean')
        cas.bootstrap_group_score_fc_histogram(exp_hrde1.fc.mean(axis=1), hrde1_regulated) # done


        ## bootstrap - negative control (wrong experiment) ### later - all wrong
        print('exp1 - hrde kennedy')
        cas.bootstrap_group_score_fc_histogram(exp1.fc.iloc[:,0], hrde1_kennedy)
        cas.bootstrap_group_score_fc_histogram(exp1.fc.iloc[:,1], hrde1_kennedy)
        cas.bootstrap_group_score_fc_histogram(exp1.fc.iloc[:,2], hrde1_kennedy)
        cas.bootstrap_group_score_fc_histogram(exp1.fc.iloc[:,3], hrde1_kennedy)
        print('exp1 - hrde kennedy mean')
        cas.bootstrap_group_score_fc_histogram(exp1.fc.mean(axis=1), hrde1_kennedy)

        wbid_list = list(exp1.fc.index)
        rand_genes_1500 = random.sample(wbid_list, 1500)
        cas.bootstrap_group_score_fc_histogram(exp1.fc.mean(axis=1), rand_genes_1500)


        print('exp1 - hrde regulated')
        cas.bootstrap_group_score_fc_histogram(exp1.fc.iloc[:,0], hrde1_regulated)
        cas.bootstrap_group_score_fc_histogram(exp1.fc.iloc[:,1], hrde1_regulated)
        cas.bootstrap_group_score_fc_histogram(exp1.fc.iloc[:,2], hrde1_regulated)
        cas.bootstrap_group_score_fc_histogram(exp1.fc.iloc[:,3], hrde1_regulated)
        print('exp1 - hrde regulated mean')
        cas.bootstrap_group_score_fc_histogram(exp1.fc.mean(axis=1), hrde1_regulated)


        ####################################################
        ############   4. by standers   ####################
        ####################################################
        
        #### 4.A - Proximity
        print('\n\n\............n\n\n')
        print('4.A - Proximity')

        ### C09G9.5
        analyze_single_gene(exp1, 'C09G9.5')
        analyze_single_gene(exp_mss, 'C09G9.5')

        ### spr-2
        analyze_single_gene(exp1, 'spr-2')
        analyze_single_gene(exp_mss, 'spr-2')

        ### F14E5.8
        analyze_single_gene(exp1, 'F14E5.8')
        analyze_single_gene(exp_mss, 'F14E5.8')

        ### F14E5.7
        analyze_single_gene(exp1, 'F14E5.7')
        analyze_single_gene(exp_mss, 'F14E5.7')

        ### F14E5.1
        analyze_single_gene(exp1, 'F14E5.1')
        analyze_single_gene(exp_mss, 'F14E5.1')

        ### efl-3
        analyze_single_gene(exp1, 'efl-3')
        analyze_single_gene(exp_mss, 'efl-3')

        ### unc-119
        analyze_single_gene(exp1, 'unc-119')
        analyze_single_gene(exp_mss, 'unc-119')

        ### rad-26
        analyze_single_gene(exp1, 'rad-26')
        analyze_single_gene(exp_mss, 'rad-26')

        ### C27B7.2
        analyze_single_gene(exp1, 'C27B7.2')
        analyze_single_gene(exp_mss, 'C27B7.2')

        ### npax-4
        analyze_single_gene(exp1, 'npax-4')
        analyze_single_gene(exp_mss, 'npax-4')
        
        ### C09G9.8
        analyze_single_gene(exp1, 'C09G9.8')
        analyze_single_gene(exp_mss, 'C09G9.8')


        #### 4.B - Sequence similarity
        print('\n\n\............n\n\n')
        print('4.B')

        ### oma-2
        analyze_single_gene(exp1, 'oma-2') # later: rep-0 is a problem
        analyze_single_gene(exp_mss, 'oma-2') 


        #### 4.D - hrde-1 by-standers
        ds = [1000, 1500, 2000, 3000, 5000, 10000, 15_000, 25_000]

        ### genes up:
        for distance in ds:
            genes_up, genes_down = ut.get_nearby_genes_list(hrde1_regulated, distance)
            print(f'for distance: {distance}, num  of genes:{len(genes_up)}')
            cas.bootstrap_group_score_fc_histogram(exp_hrde1.fc.mean(axis=1), genes_up)
            plt.show()

        ### genes up:
        for distance in ds:
            genes_up, genes_down = ut.get_nearby_genes_list(hrde1_regulated, distance)
            print(f'for distance: {distance}, num  of genes:{len(genes_down)}')
            cas.bootstrap_group_score_fc_histogram(exp_hrde1.fc.mean(axis=1), genes_down)
            plt.show()
    
    print('end')






    







    

