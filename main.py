import pandas as pd
import random
import numpy as np
import matplotlib.pyplot as plt

from gene_id import Gene_IDs
import plotting as my_plots
import utilities as ut
from atac_signal import ATAC_signal
from gene_sets import Gene_sets
import calc_signals as cas
from mRNA_gonads import Table_mRNA

import seaborn as sns

'''
same graph of "plot_group_signal" but with grey-bootstrap:
- for regulated
- for kennedy

two conditions, 3 reps.
'''

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

def get_grde1_kennedy_fc_tables():
    '''
    Creates table, each row is a gene:
    columns: mRNA-SX-mean, mRNA-hrde1-mean, ATAC-FC, ATAC-score-SX, ATAC-score-hrde1
    '''
    hrde1_k_fc_exp1 = exp1.fc.loc[hrde1_kennedy_intersected,:]
    hrde1_k_fc_exp_hrde1 = exp_hrde1.fc.loc[hrde1_kennedy_intersected,:]

    hrde1_k_mrna = m.mRNA.loc[hrde1_kennedy_intersected, 'hrde-1 mean':'sx mean']

    cat_plot_df(hrde1_k_fc_exp1, 'hrde-1 kennedy genes: ATAC FC, exp1')
    cat_plot_df(hrde1_k_fc_exp_hrde1, 'hrde-1 kennedy genes: ATAC FC, exp_hrde1')

    cat_plot_df(hrde1_k_mrna, 'hrde-1 kennedy genes: hrde1 mRNA')


    cat_plot_df(exp1.fc, 'exp1 FC', hue_hrde1=True)




def cat_plot_df(df, title:str='XXX', log_flag=False, hue_hrde1:bool=False):
    '''
    '''
    long_df = df.melt(ignore_index=False)
    if hue_hrde1:
        add_col_if_hrde1_kennedy(long_df)
        g = sns.catplot(x="variable", y="value", hue='hrde-1_k_flag', data=long_df, s=3)
    else:
        g = sns.catplot(x="variable", y="value", data=long_df, s=3)
    
    # if log_flag:
    #     g.set(yscale='log')

    g.fig.suptitle(title)

    plt.show()



def add_col_if_hrde1_kennedy(df):
    '''
    '''
    global gs
    hrde1_kennedy = gs.get_list('hrde-1-Kennedy')
    df['hrde-1_k_flag'] = [1 if ind in hrde1_kennedy else 0 for ind in df.index]


def plot_hist_values(df, xlim=5, bins=15, logy=False):
    df_new = df[df<xlim]
    df_new.stack().plot.hist(bins=bins, logx=False, logy=logy, density=1)


def list_fold_change_of_scores(ATAC_exp, thresh1=2, thresh2=2):
    list_fc_df = []
    df_cond1 = ATAC_exp.scores1 + 1e-7 # later
    df_cond2 = ATAC_exp.scores2 + 1e-7 # later

    for rep_i in range(ATAC_exp.num_of_reps):
        cond1 = ATAC_exp.scores1[f'rep {rep_i}']
        cond2 = ATAC_exp.scores2[f'rep {rep_i}']
        inds_above_thresh_1 = cond1[cond1>thresh1].index
        inds_above_thresh_2 = cond2[cond2>thresh2].index

        inds_intersect = ut.intersect_lists(inds_above_thresh_1, inds_above_thresh_2)

        cond1_masked = cond1.loc[inds_intersect]
        cond2_masked = cond2.loc[inds_intersect]

        rep_fc = cond2_masked / cond1_masked

        list_fc_df.append(rep_fc)
    
    return list_fc_df


#####
def analyze_single_gene_hist(ATAC_exp, gene_name: str, plot_range=(-1000,1000)):
    '''
    Plots both signal of gene and FC distribution.
    Both seperately and meaned reps.
    '''
    print(f'Analyze gene {gene_name} for experiment {ATAC_exp.exp_name}')
    my_plots.plot_signal_gene(ATAC_exp, gene_name, mean_flag=True, var_type='sem', plot_range=plot_range)
    my_plots.plot_gene_atac_signal_distribution(ATAC_exp, gene_name, mean_flag=False, plot_type='hist') 
    my_plots.plot_gene_atac_signal_distribution(ATAC_exp, gene_name, mean_flag=True, plot_type='hist')
#####


if __name__=='__main__':
    oma_wbid = 'WBGene00003864'
    gid = Gene_IDs()
    m = Table_mRNA()

    if "gs" not in locals():
        gs = Gene_sets()

    if "exp1" not in locals():
        exp1 = ATAC_signal("exp1")
        # exp1_pc = ATAC_signal('exp1', pc_only_flag=True)

    if "exp_mss" not in locals():
        exp_mss = ATAC_signal("exp_metsetset")
        # exp_mss_pc = ATAC_signal("exp_metsetset", pc_only_flag=True)
    
    if "exp_hrde1" not in locals():
        exp_hrde1 = ATAC_signal("exp_hrde_guy")
        # exp_hrde1_pc = ATAC_signal("exp_hrde_guy", pc_only_flag=True)
    
    
    #### get hrde-1 lists
    hrde1_regulated = ut.get_hrde_regulated(gs)
    hrde1_reg_intersected = ut.intersect_lists(hrde1_regulated, exp1.scores1.index) # later: 15 / 151 missing

    hrde1_kennedy = gs.get_list('hrde-1-Kennedy')
    hrde1_kennedy_intersected = ut.intersect_lists(hrde1_kennedy, exp1.fc.index) # later: 68 / 1527 missing
    hrde1_kennedy_mRNA_intersected = ut.intersect_lists(hrde1_kennedy, m.table.index) # later: 47 / 1527 missing

    hrde1_dic = {'hrde1 kennedy':hrde1_kennedy, 'hrde1 regulated':hrde1_regulated}

    highly, lowly = ut.get_highly_lowly() # 1018 each

    highly_no_hrde = ut.intersect_lists(highly, hrde1_kennedy, 'only first') # 930
    lowly_no_hrde = ut.intersect_lists(lowly, hrde1_kennedy, 'only first') # 1011

    hrde1_kennedy_highly = ut.intersect_lists(hrde1_kennedy_intersected, highly)
    ### check the missing genes
    # hrde1_kennedy_missing = ut.intersect_lists(hrde1_kennedy, exp1.scores1.index, 'only first')
    # x = pd.DataFrame({'missing hrde1-Kennedy genes:' : hrde1_kennedy_missing})
    # x.to_csv('missing hrde1 genes', index=False)


    if False:#### verify normalization ######
        print('normalization - exp1')
        my_plots.plot_groups_signals(exp1)
        print('normalization - exp1 pc_only')
        my_plots.plot_groups_signals(exp1_pc)


        my_plots.plot_groups_signals(exp_mss)
        my_plots.plot_groups_signals(exp_mss_pc)


        my_plots.plot_groups_signals(exp_hrde1,groups_dic={'highly w.o hrde-genes':highly_no_hrde, 'lowly w.o hrde-genes':lowly_no_hrde}, add_highly_lowly=True)
        my_plots.plot_groups_signals(exp_hrde1_pc)


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

    if False:
        #### 3.A - 'define HRDE-1 regulated group'  - scatter plot
        print('3.A - not done')
        my_plots.scatter_genes_both_conds(exp_hrde1, marked_list=hrde1_reg_intersected, shown_value='score', log_flag=True)
        my_plots.scatter_mRNA_both_conds(log_flag=True, marked_list=hrde1_reg_intersected)

        ### scatter atac-score
        my_plots.scatter_genes_both_conds(exp_hrde1, marked_list=hrde1_kennedy_intersected, shown_value='score', log_flag=True)
        my_plots.scatter_mRNA_both_conds(log_flag=True, marked_list=hrde1_kennedy_intersected)

        
        ### exp1:
        my_plots.scatter_genes_both_conds(exp1, marked_list=hrde1_reg_intersected, shown_value='score', log_flag=True)

##### Itamar
        highly_intersected = ut.intersect_lists(highly, exp1.fc.index)
        my_plots.scatter_genes_both_conds(exp_hrde1, marked_list=highly_intersected, shown_value='score', log_flag=True)

        score_higher_than_16 = list(exp_hrde1_pc.scores2[exp_hrde1_pc.scores2['rep 0']>16].index)

##### Itamaer End


        #### 3.B - HRDE-1 regulated targets get open:
        print('3.B - hrde-1')
        
        my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde1 kennedy':hrde1_kennedy})
        my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde1 regulated':hrde1_regulated})
        my_plots.plot_groups_signals(exp_hrde1, groups_dic=hrde1_dic)
        my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde1 kennedy':hrde1_kennedy}, mean_flag=True, var_type='sem')
        my_plots.plot_groups_signals(exp_hrde1, groups_dic={'hrde1 regulated':hrde1_regulated}, mean_flag=True, var_type='sem')
        my_plots.plot_groups_signals(exp_hrde1, groups_dic=hrde1_dic, mean_flag=True, var_type='sem')


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



def plot_groups_signals(
    ATAC_exp,
    groups_dic: dict = {},
    mean_flag: bool = False,
    var_type: str = "none",
    add_highly_lowly: bool = True,
    bootstrap: bool = False,
    boot_size: int = 2315,
    boot_iters: int = 1000,
    drop_rep: int = 10,
    zscore_signal: bool = False,
    plot_range: tuple = (0, 0),
    drop_gene_list: list = []
):
    """
    Takes an experiment, a dictionary with groups, plots panel with two axes.
    - one ax for each condition: each group is a line with std.

    If mean_flag: variance is between reps. if seperate: variance is between genes.

    Parameters
    --------
    - ATAC_exp: ATAC_signal object.
    - groups_dic: keys - group_name, values - list of wbids.
    """
    if add_highly_lowly:
        my_plots.add_highly_lowly_to_dic(groups_dic)

    if bootstrap:
        print(f"bootstrapping {boot_iters} iterations...")

    if not mean_flag:
        for rep_i in range(ATAC_exp.num_of_reps):
            fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
            fig.suptitle(
                f"{ATAC_exp.exp_name}, Replicate {rep_i+1}, (variance between genes)",
                fontsize=14,
            )
            for cond_i in [
                0,
                1,
            ]:  # for each condition (sample) (get for this sample a df of means, and df of vars):
                sample_df = ATAC_exp.exp_df.iloc[rep_i, cond_i]
                if zscore_signal:  # later
                    sample_df = ut.normalize_zscore_df(sample_df)
                means_df, vars_df = my_plots.groups_df_mean_and_var_dfs_for_sample(
                    sample_df, groups_dic, var_type
                )

                if bootstrap:
                    (
                        means_df[f"bootstrap ({boot_size} genes)"],
                        bootstrap_var,
                    ) = cas.bootstrap_atac_signal(
                        sample_df, group_size=boot_size, num_of_iters=boot_iters
                    )  # later
                    if not isinstance(vars_df, int):
                        vars_df[f"bootstrap ({boot_size} genes)"] = bootstrap_var

                if plot_range.count(0) != 2:  # if range was given:
                    means_df = my_plots.narrow_to_range(means_df, plot_range[0], plot_range[1])
                    vars_df = my_plots.narrow_to_range(vars_df, plot_range[0], plot_range[1])

                axes[cond_i].set_title(f"{ATAC_exp.condition_names[cond_i]}")
                legend_flag = cond_i  # 0 / 1 [only legend on right ax]
                my_plots.plot_ax(axes[cond_i], means_df, vars_df, legend_flag)

    else:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
        fig.suptitle(
            f"{ATAC_exp.exp_name}, Mean of All Replicates (variance between replicates)",
            fontsize=14,
        )

        for cond_i in [0, 1]:
            means_df_list = []
            for rep_i in range(ATAC_exp.num_of_reps):
                if rep_i == drop_rep:
                    print(f"dropped rep {rep_i}")
                    continue
                sample_df = ATAC_exp.exp_df.iloc[rep_i, cond_i]

                ##### added:
                if len(drop_gene_list)>0:
                    sample_df.drop(drop_gene_list, inplace=True)

                if zscore_signal:  # later
                    sample_df = ut.normalize_zscore_df(sample_df)
                means_df, _ = my_plots.groups_df_mean_and_var_dfs_for_sample(
                    sample_df, groups_dic, var_type="none"
                )

                if bootstrap:
                    (
                        means_df[f"bootstrap ({boot_size} genes)"],
                        _,
                    ) = cas.bootstrap_atac_signal(
                        sample_df, group_size=boot_size, num_of_iters=boot_iters
                    )  # later
                means_df_list.append(means_df)

            df_means_all_reps, df_vars = my_plots.get_mean_variance_of_df_list(
                means_df_list, var_type
            )

            if plot_range.count(0) != 2:  # if range was given:
                df_means_all_reps = my_plots.narrow_to_range(
                    df_means_all_reps, plot_range[0], plot_range[1]
                )
                df_vars = my_plots.narrow_to_range(df_vars, plot_range[0], plot_range[1])

            axes[cond_i].set_title(f"{ATAC_exp.condition_names[cond_i]}")
            legend_flag = cond_i  # 0 / 1 [only legend on right ax]

            my_plots.plot_ax(axes[cond_i], df_means_all_reps, df_vars, legend_flag)






    







    

