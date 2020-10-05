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


def generate_mock_groups():
    gene_list_1 = ["WBGene00007063", "WBGene00007064", "WBGene00007067"]
    gene_list_2 = [
        "WBGene00017071",
        "WBGene00019895",
        "WBGene00009583",
        "WBGene00018682",
    ]

    dic_groups = {"group 1": gene_list_1, "group 2": gene_list_2}
    return dic_groups


def mean_gene_groups_of_sample(sample_df: pd.DataFrame, dic_groups: dict):
    """
    Gets df of sample and groups and returns a dict with mean of these groups.

    Parameters
    -----------
    - sample_df: pd.DataFrame of ATAC-seq signal
    - dic_groups: dict, group_name : list_of_ids

    Return
    -----------
    - dic_means: dict of pd.Series, group_name : mean_of_all_genes_in_group
    """
    dic_means = {}
    for group in dic_groups:
        group_ids = dic_groups[group]
        # take only genes that appear in the sample df:
        intersected_list = list(set(sample_df.index) & set(group_ids))
        dic_means[group] = sample_df.loc[intersected_list, :].mean()
    return dic_means


def plot_groups(dic_means, ax_place):
    """
    Gets a dict of means and plots.
    """
    df_groups = pd.DataFrame(dic_means)
    df_groups.plot(ax=ax_place)


def plot_replicates(reps_dic: dict, dic_groups: dict, conditions: tuple=('gfp','oma-1')):
    """
    Gets a dict with all samples (dfs) ordered by replicates, 
    and the desired groups of genes, and plots them seperately.
    """
    for rep in reps_dic:
        group_a = reps_dic[rep][0]
        group_b = reps_dic[rep][1]

        fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(9,4))

        print(f"Replicate {rep} :")
        dic_means_a = mean_gene_groups_of_sample(sample_df=group_a, dic_groups=dic_groups)
        ax_a = axes[0]
        plot_groups(dic_means_a, ax_a)
        axes[0].set_title(conditions[0], fontsize=18)

        dic_means_b = mean_gene_groups_of_sample(sample_df=group_b, dic_groups=dic_groups)
        ax_b = axes[1]
        plot_groups(dic_means_b, ax_b)
        axes[1].set_title(conditions[1], fontsize=18)
        plt.show()

def plot_mean_of_replicates(reps_dic: dict, dic_groups: dict):
    '''
    '''
    pass
    ### compute mean of samples and std (as shadow) between samples.
    ### https://seaborn.pydata.org/generated/seaborn.lineplot.html


if __name__ == "__main__":
    ## generate dictionary of all samples of experiment
    exp1_dfs_list = rdt.read_experiment(exp_name="exp1")
    gfp_0, gfp_1, gfp_2, gfp_4, oma1_0, oma1_1, oma1_2, oma1_4 = exp1_dfs_list
    reps_dic = {
        0: [gfp_0, oma1_0],
        1: [gfp_1, oma1_1],
        2: [gfp_2, oma1_2],
        4: [gfp_4, oma1_4],
    }

    ## generate lists of genes to be plotted
    gs = Gene_sets()

    dic_list = {'hrde-1':['isHrde1', 10], 
        'pol-2':['isPol2'], 
        'highly':['R1-SX_S14', 10],
        'all genes':['ALL']}

    dic_groups = gs.get_multiple_lists(dic_list)

    # oma_1_gene = ['WBGene00003864']
    # dic_groups = {'oma-1 gene': oma_1_gene}

    plot_replicates(reps_dic, dic_groups, conditions=('anti-gfp RNAi', 'anti-oma-1 RNAi'))
