import read_tables as rdt
import matplotlib.pyplot as plt
import pandas as pd
from gene_sets import Gene_sets
import seaborn as sns
import calc_signals as cas 

def plot_experiment_dfs(
    results_df_dict: dict,
    conditions: tuple = ("c1", "c2"),
    mean_all_reps: bool = False,
    compare_conditions: bool = False,
    variance_type: str = "none",
):
    """
    Gets a dict with results, each sample is a df.
    Plots all results with required variance.
    - If the results_df_dict has only one df for each repeat (e.g. FC between conditions was calculated), Than "conditions" is the calculation type.
    - Variance type is plotted only if mean_all_reps is True.


    Parameters
    ------------
    - results_df_dict: dict. {#rep_num : #list_of_dfs_for_this_rep}. Each df row is a gene_name/group_name. Each rep is a list of length 2 or one (Depends if conditions were merged in calculation).
    - conditions: tuple of str. Either two condition names or one name for the calculation (e.g. "ATAC FC - GFP/OMA-1")
    - mean_all_reps: bool=False. If true, means all replicates an plots on the same panel.
    - compare_conditions: bool=False. If True, plots both condition on the same axe.
    - variance_type: str='none'. Can be ['none','std','sem']. Only relevant if mean_all_reps = True.
    """

    #### verify input
    # num_of_conds = # get 1 or 2
    # if num_of_conds==1:
        # assert compare_conditions=False
    # if mean_all_reps==False:
        # assert variance_type=='none'
    # else:
        # assert variance_type.lower() in ['none', 'std', 'sem']
    # if compare_conditions == True:
        # assert that the rows are the same in all dfs


    # if mean_all:
        ### calculate mean and variance
        # initiate mean_dfs_list, var_df_list
        # for each condition:
            # df_mean, df_var = get_mean_variance(df_list, variance_type)
            # add them to mean_dfs_list, var_df_list

        # plot_panel(mean_dfs_list, var_df_list, conditions, compare_conditions)
    
    # If not mean all:
        # for each rep_num in results_df_dict:
            # print: replicate {rep_num}
            # plot_panel (dfs_list, conditions, compare_conditions)


def plot_panel(dfs_list, conditions, compare_conditions, var_df_list=0):
    '''
    Gets 1/2 dfs and plots them on 1/2 axes in same panel.
    If plotted on the same panel, there will be two different styles.

    If var_df_list is [0,0,0,0] it will not 
    '''
    # create color pallete?

    # if dfs_list has only 1:
        # create a panel of 1 ax, title "conditions"
        # legend!!!
        # loop on cols of df (location)
            # plot_vector(df_col, variance_col, ax, color, style)
    
    # if dfs_list has 2:
        # if compare==False:
            # create a panel of 2 axes.
            # title 1 and 2 from conditions
            # loop on cols:
                # plot_vector(df1_col, variance_col, ax1, color, style)
                # plot_vector(df2_col, variance_col, ax2, color, style)
        # if compare==True:
            # create a panel of 1 ax.
            # legend!!!
            # loop on cols:
                # plot_vector(df1_col, variance_col, ax, color, style1)
                # plot_vector(df2_col, variance_col, ax, color, style2)

def plot_vector(vec, variance_vec, ax, color, style):
    '''
    Plots a line on the graph, with variance as a band.
    if variance_vec==0, doesn't plot variance.
    '''
    # plot line on ax
    # if variance_vec!=0:
        # plot variance band


if __name__ == "__main__":
    ## generate dictionary of all samples of experiment
    exp_dic = rdt.read_experiment_to_dic(exp_name="exp1")

    ## generate lists of genes to be plotted
    gs = Gene_sets()
    dic_list = {
        "hrde-1": ["hrde-1-Kennedy"],
        "pol-2": ["isPol2"],
        "highly 10%": ["expression_mean", 10]
    }
    dic_groups = gs.get_multiple_lists(dic_list)

    # create means of these groups, one for each condition:
    df_means_list_a = cas.get_df_means_list(exp_dic, dic_groups, cond_num=0)
    df_means_list_b = cas.get_df_means_list(exp_dic, dic_groups, cond_num=0)




