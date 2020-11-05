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

    # generate color pallete
    colors = plt.get_cmap('Accent')
    lines = []

    if len(dfs_list) == 1:
        df = dfs_list[0]
        var_df = var_df_list[0]
        fig, axes = plt.subplots()
        plt.title(f"{conditions[0]}", fontsize=14)
        plt.xlabel('Location relative to TSS')
        plt.ylabel('ATAC-seq signal (norm.)')

        if not isinstance(var_df, int):
            for col_i in range(df.shape[1]):
                line = plot_vector(vec=df.iloc[:,col_i], ax=axes, color=colors(col_i), variance_vec=var_df.iloc[:,col_i])
                lines.append(line)
        else:
            for col_i in range(df.shape[1]):
                line = plot_vector(vec=df.iloc[:,col_i], ax=axes, color=colors(col_i))
                lines.append(line)

        plt.legend(lines, df.columns)

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

    plt.show()


def plot_vector(vec, ax, color, style="solid", variance_vec=0):
    """
    Plots a line on the graph, with variance as a band.
    if variance_vec==0, doesn't plot variance.
    """
    ## print(f'variance_vec: {variance_vec}') ### later

    line = ax.plot(vec, c=color, linestyle=style)
    if not isinstance(variance_vec, int):
        ax.fill_between(vec.index, vec-variance_vec, vec+variance_vec, fc=color, alpha=0.2)
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
        dic_list = {
            "hrde-1": ["hrde-1-Kennedy"],
            "pol-2": ["isPol2"],
            "highly 10%": ["expression_mean", 10],
        }
        dic_groups = gs.get_multiple_lists(dic_list)

    # create means of these groups, one for each condition:
    df_means_list_a = cas.get_group_means_df_list(exp_dic, dic_groups, cond_num=0)
    df_means_list_b = cas.get_group_means_df_list(exp_dic, dic_groups, cond_num=1)

    dict_conditions = {"group a": df_means_list_a, "group b": df_means_list_b}
    dict_condition = {"group a": df_means_list_a}


    print('do not mean, one condition')
    plot_experiment_dfs(
        dict_condition,
        mean_all_reps=True,
        compare_conditions=False,
        variance_type="std",
    )

    print('do not mean, one condition')
    plot_experiment_dfs(
        dict_condition,
        mean_all_reps=False,
        compare_conditions=False,
        variance_type="none",
    )

    df_mean, df_std = cas.get_mean_variance(df_means_list_a, 'std')
    df_col = df_mean.iloc[:,0]
    std_col = df_std.iloc[:,0]

    # fig, axes = plt.subplots()
    # line = axes.plot(df_col, color='red')
    # y1 = df_col+std_col
    # y2 = df_col-std_col
    # axes.fill_between(df_col.index, y1, y2, fc='red', alpha=0.2)
    # plt.legend(line,[df_col.name])

