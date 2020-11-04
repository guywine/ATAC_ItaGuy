import read_tables as rdt
import matplotlib.pyplot as plt
import pandas as pd
from gene_sets import Gene_sets
import seaborn as sns
import calc_signals as cas 

def plot_experiment_dfs(
    dict_conditions: dict,
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
    - dict_conditions: dict. "condition_name":list_of_sample_dfs_from_all_reps. Can have 1/2 conditions
    - mean_all_reps: bool=False. If true, means all replicates an plots on the same panel.
    - compare_conditions: bool=False. If True, plots both condition on the same axe.
    - variance_type: str='none'. Can be ['none','std','sem']. Only relevant if mean_all_reps = True.
    """

    #### verify input
    conditions = list(dict_conditions.keys())
    num_of_conds = len(conditions)

    assert num_of_conds in {1,2}
    if num_of_conds==1:
        assert compare_conditions == False, 'Cannot comapre conditions, only one condition given'
    
    if mean_all_reps == False:
        assert variance_type=='none', 'Cannot plot variance for independent repeats'
    else:
        assert variance_type.lower() in {'none', 'std', 'sem'}, 'Variance type should be "none"/"std"/"sem"'
    
    if compare_conditions == True:
        assert assert_conditions_are_comparable(dict_conditions[conditions[0]], dict_conditions[conditions[1]])

    #### plot
    if mean_all:
        ### calculate mean and variance
        mean_dfs_list = [] 
        var_df_list = []
        for condition in conditions:
            df_mean, df_var = cas.get_mean_variance(dict_conditions[condition], variance_type)
            mean_dfs_list.append(df_mean)
            var_df_list.append(df_var)

        plot_panel(mean_dfs_list, conditions, compare_conditions, var_df_list)
    
    else: # if not mean all:
        for rep_num in len(df_list_a):
            dfs_list = [df_list[rep_num] for df_list in dict_conditions.values()]
            print(f'Replicate #{rep_num}:')
            plot_panel(dfs_list, conditions, compare_conditions)


def plot_panel(dfs_list:list, conditions: list, compare_conditions: bool, var_df_list=[0]):
    '''
    Gets 1/2 dfs and plots them on 1/2 axes in same panel.
    If plotted on the same panel, there will be two different styles.

    If var_df_list is list of zeros, it will not plot.

    '''
    return 0
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


def assert_conditions_are_comparable(df_list_a: pd.DataFrame, df_list_b: pd.DataFrame):
    '''
    Takes two conditions, each with its df_list, and verifies that the structures of the data are the same.
    Verifies the shape, indices, and columns.

    * Checks only first df in list, assuming inside the list the dfs are the same.

    Parameters
    ----------
    - df_list_a: pd.DataFrame
    - df_list_b: pd.DataFrame

    Return
    ---------
    - bool. True if the dfs are equal.
    '''
    a0 = df_list_a[0]
    b0 = df_list_b[0]

    if a0.shape!=b0.shape:
        return False
    
    inds_same = a0.index==b0.index
    if not all(inds_same):
        return False
    
    cols_same = a0.columns==b0.columns
    if not all(cols_same):
        return False
    
    return True

    


    


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
    df_means_list_a = cas.get_group_means_df_list(exp_dic, dic_groups, cond_num=0)
    df_means_list_b = cas.get_group_means_df_list(exp_dic, dic_groups, cond_num=1)

    dict_conditions = {'group a':df_means_list_a, 'group b':df_means_list_b}

    plot_experiment_dfs(dict_conditions, mean_all_reps=True, compare_conditions=False, variance_type='std')
    




