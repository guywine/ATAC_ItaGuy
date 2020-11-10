import read_tables as rdt
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
from gene_sets import Gene_sets
import seaborn as sns
import calc_signals as cas


def plot_experiment_dfs(
    dict_conditions: dict,
    mean_all_reps: bool = True,
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
    - dict_conditions: dict. "condition_name":list_of_dfs_from_all_reps. Can have 1/2 conditions
    - mean_all_reps: bool=False. If true, means all replicates an plots on the same panel.
    - compare_conditions: bool=False. If True, plots both condition on the same axe.
    - variance_type: str='none'. Can be ['none','std','sem']. Only relevant if mean_all_reps = True.
    """

    #### verify input
    conditions = list(dict_conditions.keys())
    num_of_conds = len(conditions)

    assert num_of_conds in {1, 2}
    if num_of_conds == 1:
        assert (
            compare_conditions == False
        ), "Cannot comapre conditions, only one condition given"
        df_list_a = dict_conditions[conditions[0]]
    else: 
        print('here1') ### later
        df_list_a = dict_conditions[conditions[0]]
        df_list_b = dict_conditions[conditions[1]]

    if mean_all_reps == False:
        assert variance_type == "none", "Cannot plot variance for independent repeats"
    else:
        assert variance_type.lower() in {
            "none",
            "std",
            "sem",
        }, 'Variance type should be "none"/"std"/"sem"'

    if compare_conditions == True:
        assert assert_conditions_are_comparable(df_list_a, df_list_b)

    #### plot
    if mean_all_reps:
        ### calculate mean and variance
        mean_dfs_list = []
        var_df_list = []
        print(f'var type: {variance_type}') ### later
        for condition in conditions:
            df_mean, df_var = cas.get_mean_variance(
                dict_conditions[condition], variance_type
            )
            mean_dfs_list.append(df_mean)
            var_df_list.append(df_var)

        # print(f'df_var 1: {df_var}') ### later
        plot_panel(mean_dfs_list, conditions, compare_conditions, var_df_list)

    else:  # if not mean all:
        for rep_num in range(len(df_list_a)):
            dfs_list = [df_list[rep_num] for df_list in dict_conditions.values()]
            print(f"Replicate #{rep_num+1}:")
            plot_panel(dfs_list, conditions, compare_conditions)
        
   
def plot_panel(
    dfs_list: list, conditions: list, compare_conditions: bool, var_df_list=[0]
):
    """
    Gets 1/2 dfs and plots them on 1/2 axes in same panel.
    If plotted on the same panel, there will be two different styles.

    If var_df_list is list of zeros, it will not plot.

    """
    # generate list of zeros if no variance is given
    if isinstance(var_df_list[0],int):
        var_df_list = [0] * len(dfs_list)

    if len(dfs_list) == 1:
        df = dfs_list[0]
        var_df = var_df_list[0]
        # print(f'var_df 2: {var_df}') ### later
        fig, axes = plt.subplots()
        plt.title(f"{conditions[0]}", fontsize=14)
        plt.xlabel('Location relative to TSS')
        plt.ylabel('ATAC-seq signal (norm.)')

        plot_ax(df, axes, var_df)


    else: # if dfs_list has 2:
        df_a = dfs_list[0]
        var_df_a = var_df_list[0]
        df_b = dfs_list[1]
        var_df_b = var_df_list[1]

        if compare_conditions==False:
            fig, axes = plt.subplots(1,2, figsize=(12,5), sharey=True)
            axes[0].set_title(conditions[0], fontsize=16)
            axes[1].set_title(conditions[1], fontsize=16)

            plot_ax(df_a, axes[0], var_df_a, legend_flag=False)
            plot_ax(df_b, axes[1], var_df_b)
        
        else: # if compare==True:
            fig, axes = plt.subplots()
            plt.title(f"{conditions[0]} and {conditions[1]}", fontsize=16)

            plot_ax(df_a, axes, var_df_a, style='solid')
            plot_ax(df_b, axes, var_df_b, style='dashed', legend_flag=False)

            ### add legend of conditions:
            black_line = mlines.Line2D([], [], color='black', label=conditions[0])
            dashed_line = mlines.Line2D([], [], color='black', label=conditions[1], linestyle='dashed')
            plt.legend(handles=[black_line, dashed_line], loc='right')

    plt.show()

def plot_ax(vec_df, ax, var_df=0, style='solid', legend_flag:bool = True):
    '''
    Gets:
    - vec_df: df. 
    - ax: ax to plot on.
    - colors: color pallete
    - style: str. linestyle.
    - var_df: either 0 or given pd.DataFrame with same columns as vec_df.
    - legend_flag: bool. If true, plot legend on this axis.
    '''
    lines = []
    colors = plt.get_cmap('Accent')

    # print(f'var_df: {var_df}') ### later

    if not isinstance(var_df, int): # if a variance vector was given
        for col_i in range(vec_df.shape[1]):
            line = plot_vector(vec=vec_df.iloc[:,col_i], ax=ax, color=colors(col_i), style=style, var_vec=var_df.iloc[:,col_i])
            lines.append(line)
    else:
        for col_i in range(vec_df.shape[1]):
            line = plot_vector(vec=vec_df.iloc[:,col_i], ax=ax, color=colors(col_i), style=style)
            lines.append(line)
    
    if legend_flag:
        first_leg = plt.legend(lines, vec_df.columns)
        ax = plt.gca().add_artist(first_leg)
        


def plot_vector(vec, ax, color, style="solid", var_vec=0):
    """
    Plots a line on the graph, with variance as a band.
    if variance_vec==0, doesn't plot variance.
    """
    ## print(f'variance_vec: {variance_vec}') ### later

    line = ax.plot(vec, c=color, linestyle=style)
    if not isinstance(var_vec, int):
        ax.fill_between(vec.index, vec-var_vec, vec+var_vec, fc=color, alpha=0.2)
    return line[0]


def assert_conditions_are_comparable(df_list_a: pd.DataFrame, df_list_b: pd.DataFrame):
    """
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
    """
    a0 = df_list_a[0]
    b0 = df_list_b[0]

    if a0.shape != b0.shape:
        return False

    inds_same = a0.index == b0.index
    if not all(inds_same):
        return False

    cols_same = a0.columns == b0.columns
    if not all(cols_same):
        return False

    return True


if __name__ == "__main__":
    if "dic_groups" not in locals():
        ## generate dictionary of all samples of experiment
        exp_dic = rdt.read_experiment_to_dic(exp_name="exp1")

        ## generate lists of genes to be plotted
        gs = Gene_sets()
        dic_names = {
            "hrde-1": ["hrde-1-Kennedy"],
            "pol-2": ["isPol2"],
            "highly 10%": ["expression_mean", 10],
        }
        dic_groups = gs.get_multiple_lists(dic_names)

        # dic_groups = {'oma-1':['WBGene00003864']}

    # create means of these groups, one for each condition:
    df_means_list_a = cas.get_group_means_df_list(exp_dic, dic_groups, cond_num=0)
    df_means_list_b = cas.get_group_means_df_list(exp_dic, dic_groups, cond_num=1)

    dict_conditions = {"GFP": df_means_list_a, "OMA-1": df_means_list_b}
    dict_condition = {"Calculated": df_means_list_b}


    print('mean, 1 condition')
    plot_experiment_dfs(
        dict_condition,
        compare_conditions=False,
        variance_type="std",
    )
    plot_experiment_dfs(
        dict_condition,
        compare_conditions=False,
        variance_type="sem",
    )
    plot_experiment_dfs(
        dict_condition,
        compare_conditions=False,
        variance_type="none",
    )

    # print('\n\ndo not mean, one condition')
    # plot_experiment_dfs(dict_condition, mean_all_reps=False)

    # print('\n\ndo not mean, 2 conditions')
    # plot_experiment_dfs(dict_conditions, mean_all_reps=False)

    print('\n\nmean, 2 conditions')
    plot_experiment_dfs(dict_conditions, mean_all_reps=True)
    plot_experiment_dfs(dict_conditions, mean_all_reps=True, variance_type='std')
    plot_experiment_dfs(dict_conditions, mean_all_reps=True, variance_type='sem')

    print('\n\nmean, 2 conditions, compare')
    plot_experiment_dfs(dict_conditions, mean_all_reps=True, compare_conditions=True)
    plot_experiment_dfs(dict_conditions, mean_all_reps=True, compare_conditions=True, variance_type='sem')

    # print('\n\nno mean, 2 conditions, compare')
    # plot_experiment_dfs(dict_conditions, mean_all_reps=False, compare_conditions=True)


    df_mean, df_std = cas.get_mean_variance(df_means_list_a, 'std')
    df_col = df_mean.iloc[:,0]
    std_col = df_std.iloc[:,0]

