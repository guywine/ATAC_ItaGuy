"""
Plot ATAC-seq signal. Choose:
- Experiment: exp1 / exp_gonads / exp_hrde_gonads / exp_hrde_guy / exp_metsetset
- set of genes: list of sets ([hrde-1 set, csr-1 set])
- Normalize: highly_lowly / none
- Mean all replicates or display them seperately?
"""
import read_tables as rdt
import matplotlib.pyplot as plt
import pandas as pd
from gene_sets import Gene_sets
import seaborn as sns


def mean_gene_groups_of_sample(sample_df: pd.DataFrame, dic_groups: dict):
    """
    Gets df of sample and groups and returns a dict with mean of these groups.

    Parameters
    -----------
    - sample_df: pd.DataFrame of ATAC-seq signal
    - dic_groups: dict, group_name : list_of_ids

    Return
    -----------
    - df_groups_means: pd.DataFrame, group_name:mean_of_all_genes_in_group
    """
    df_groups_means = pd.DataFrame([])
    for group in dic_groups:
        group_ids = dic_groups[group]
        # take only genes that appear in the sample df:
        intersected_list = list(set(sample_df.index) & set(group_ids))
        df_groups_means[group] = sample_df.loc[intersected_list, :].mean()
    return df_groups_means

def plot_new(exp_dic: dict, dic_groups: dict, conditions: tuple = ("c1", "c2"), mean_all:bool=False, compare:bool=False):
    '''
    '''
    ## get a list of dfs for conditions, each with means of all groups
    df_means_list_a = get_df_means_list(exp_dic, dic_groups, cond_num=0)
    df_means_list_b = get_df_means_list(exp_dic, dic_groups, cond_num=1)

    if mean_all:
        print('Mean of all replicates:')
        # concat
        df_means_concat_a = pd.concat(df_means_list_a)
        df_means_concat_b = pd.concat(df_means_list_b)
        
        # mean all:
        df_means_a = df_means_concat_a.groupby(level=0).mean()
        df_means_b = df_means_concat_b.groupby(level=0).mean()

        # plot
        plot_panel(df_means_a, df_means_b, conditions, compare)

    else:
        for rep_i in range(len(df_means_list_a)): ## later fix
            print(f"Replicate {rep_i}:")
            # plot
            plot_panel(df_means_list_a[rep_i], df_means_list_b[rep_i], conditions, compare)


# def plot_replicates(
#     reps_dic: dict, dic_groups: dict, conditions: tuple = ("condition 1", "condition 2")
# ):
#     """
#     Gets a dict with all samples (dfs) ordered by replicates,
#     and the desired groups of genes, and plots them seperately.
#     """
#     for rep_num in reps_dic:
#         print(f"Replicate {rep_num} :")
    
#         ## get mean of groups for first condition
#         cond_a_df = reps_dic[rep_num][0]
#         df_means_a = mean_gene_groups_of_sample(
#             sample_df=cond_a_df, dic_groups=dic_groups
#         )

#         ## get mean of groups for second condition
#         cond_b_df = reps_dic[rep_num][1]
#         df_means_b = mean_gene_groups_of_sample(
#             sample_df=cond_b_df, dic_groups=dic_groups
#         )

#         ## plot
#         plot_panel(df_means_a, df_means_b, conditions)


def plot_panel(df_means_a, df_means_b, conditions: tuple=('cond1', 'cond2'), compare:bool=False):
    '''
    Can plot panel with:
    - 2 axes, one for each condition
    - 1 ax, with both conditions, when the condition name is added to the group name
    '''
    if compare:
        fig, axes = plt.subplots(nrows=1, ncols=1, sharey=True, figsize=(4, 4))

        add_cond_to_col_names(df=df_means_a, condition=conditions[0])
        add_cond_to_col_names(df=df_means_b, condition=conditions[1])
        df_both_conds = pd.concat([df_means_a, df_means_b])
        df_both_conds.plot(ax=axes)

    else:
        fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(9, 4))

        # plot both conditions
        df_means_a.plot(ax=axes[0])
        df_means_b.plot(ax=axes[1])
        
        # set titles
        axes[0].set_title(conditions[0], fontsize=18)
        axes[1].set_title(conditions[1], fontsize=18)

    plt.show()

    return fig, axes


def get_mean_of_groups_all_replicates(
    reps_dic: dict, dic_groups: dict, cond_num:int=0
):
    """
    """
    df_means_list = get_df_means_list(reps_dic, dic_groups, cond_num)

    ### create single df average across replicates:
    df_mean_all_reps = pd.concat(df_means_list).groupby(level=0).mean()
    df_std = pd.concat(df_means_list).groupby(level=0).std()

    return df_mean_all_reps, df_std


def get_df_means_list(reps_dic: dict, dic_groups: dict, cond_num:int):
    """
    Returns a list of all samples from this condition, means of groups.
    --------
    """
    df_means_list = []
    for rep_num in reps_dic:
        sample_df = reps_dic[rep_num][cond_num]
        df_means = mean_gene_groups_of_sample(
            sample_df=sample_df, dic_groups=dic_groups
        )
        df_means_list.append(df_means)
    return df_means_list

# def plot_all_replicates_mean(reps_dic: dict, dic_groups: dict, conditions: tuple = ("condition 1", "condition 2")):
#     '''
#     ### compute mean of samples and std (as shadow) between samples.
#     ### https://seaborn.pydata.org/generated/seaborn.lineplot.html
#     '''
#     print('Mean of all repeats:')
#     df_mean_all_reps_a, df_std_a = get_mean_of_groups_all_replicates(reps_dic, dic_groups, cond_num=0)
#     df_mean_all_reps_b, df_std_b = get_mean_of_groups_all_replicates(reps_dic, dic_groups, cond_num=1)
#     plot_panel(df_mean_all_reps_a, df_mean_all_reps_b, conditions)


def add_cond_to_col_names(df:pd.DataFrame, condition:str='treated'):
    '''
    Changes inplace.
    '''
    new_cols = [f'{name} {condition}' for name in df.columns]
    cols_dict = dict(zip(df.columns, new_cols))
    df.rename(columns=cols_dict, inplace=True)



if __name__ == "__main__":
    ## generate dictionary of all samples of experiment
    exp1_dic = rdt.read_experiment_to_dic(exp_name="exp1")

    ## generate lists of genes to be plotted
    gs = Gene_sets()

    dic_list = {
        "hrde-1": ["isHrde1"],
        "pol-2": ["isPol2"],
        "highly": ["R1-SX_S14", 10],
        "all genes": ["ALL"],
    }

    dic_groups = gs.get_multiple_lists(dic_list)

    # oma_1_gene = ['WBGene00003864']
    # dic_groups = {'oma-1 gene': oma_1_gene}

    # plot_replicates(
    #     exp1_dic, dic_groups, conditions=("anti-gfp RNAi", "anti-oma-1 RNAi")
    # )

    plot_new(exp1_dic, dic_groups, conditions=("anti-gfp RNAi", "anti-oma-1 RNAi"), mean_all=False)
    plot_new(exp1_dic, dic_groups, conditions=("anti-gfp RNAi", "anti-oma-1 RNAi"), mean_all=True)
    plot_new(exp1_dic, dic_groups, conditions=("anti-gfp RNAi", "anti-oma-1 RNAi"), mean_all=False, compare=True)



    # plot_all_replicates_mean(exp1_dic, dic_groups, conditions=('gfp RNAi','oma-1 RNAi'))

    ### work with four concat tables to create STD?
    # df_means_list = get_df_means_list(exp1_dic, dic_groups, cond_num=0)
    # means_concat = pd.concat(df_means_list)
    # means_computed_a = means_concat.groupby(level=0).mean()

    # df_means_list_b = get_df_means_list(exp1_dic, dic_groups, cond_num=1)
    # means_concat_b = pd.concat(df_means_list_b)
    # means_computed_b = means_concat_b.groupby(level=0).mean()

    # means_computed_a['cond']='a'
    # means_computed_b['cond']='b'



    
