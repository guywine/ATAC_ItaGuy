"""
Gets df of samples and can calculate:
- mean of each gene-group
- ATAC_FC for each gene between two replicates
- Normalize to highly and lowly?
"""

import pandas as pd
import random
import matplotlib.pyplot as plt
import utilities as ut

def bootstrap_group_score_fc_histogram(gene_table, wbid_list, num_of_iters: int=1_000):
    '''
    - gene_table: one_column df - each gene has only a single value.
    '''
    i_bootstrap_means, group_mean = ut.bootstrap_group_score(gene_table, wbid_list, num_of_iters)

    fig, ax = plt.subplots(1, 1)
    ax.hist(i_bootstrap_means.iloc[:,0], bins=20, zorder=0)
    mark_x = group_mean  ### later
    mark_y = 5  ### later
    hand = ax.scatter(mark_x, mark_y, c="red", marker=7, zorder=5)

    plt.show()




def bootstrap_atac_signal(df_sample: pd.DataFrame, group_size: int, num_of_iters: int=1_000):
    '''
    '''
    gene_pool = pd.read_csv('tables/protein_coding_wbids.csv')

    df_boot_means = pd.DataFrame([])

    for i in range(num_of_iters):
        boot_group = random.sample(list(gene_pool['genes']), group_size)
        intersected_list = list(set(df_sample.index) & set(boot_group))

        df_boot_means[i] = df_sample.loc[intersected_list, :].mean(axis=0)
    
    boot_mean_series = df_boot_means.mean(axis=1)
    boot_std_series = df_boot_means.std(axis=1)
    
    return boot_mean_series, boot_std_series





def get_group_means_df_list(reps_dic: dict, dic_groups: dict, cond_num: int):
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


def mean_gene_groups_of_sample(sample_df: pd.DataFrame, dic_groups: dict):
    """
    Gets df of sample and dic_groups and returns a df with mean for each of these groups.

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



def get_mean_variance(df_list: list, variance_type: str = "none"):
    """
    Gets a list of dfs with identicle structure, and calculates for each cell the mean and variance across all dfs in list.
    Returns a df containing means, and a seperate df containing variance.

    Parameters
    ----------
    - df_list: list of dfs to mean (each df a sample of rep)
    - variance_type: str. ['std' / 'sem' / 'none']
    """
    ### create single df average across replicates:
    df_mean_all_reps = pd.concat(df_list).groupby(level=0).mean()
    if variance_type.lower() == "std":
        df_variance = pd.concat(df_list).groupby(level=0).std()
    elif variance_type.lower() == "sem":
        df_variance = pd.concat(df_list).groupby(level=0).sem()
    elif variance_type.lower() == "none":
        df_variance = 0

    return df_mean_all_reps, df_variance