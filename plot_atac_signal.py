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



def plot_replicates(
    reps_dic: dict, dic_groups: dict, conditions: tuple = ("condition 1", "condition 2")
):
    """
    Gets a dict with all samples (dfs) ordered by replicates,
    and the desired groups of genes, and plots them seperately.
    """
    for rep_num in reps_dic:
        print(f"Replicate {rep_num} :")
    
        ## get mean of groups for first condition
        cond_a_df = reps_dic[rep_num][0]
        df_means_a = mean_gene_groups_of_sample(
            sample_df=cond_a_df, dic_groups=dic_groups
        )

        ## get mean of groups for second condition
        cond_b_df = reps_dic[rep_num][1]
        df_means_b = mean_gene_groups_of_sample(
            sample_df=cond_b_df, dic_groups=dic_groups
        )

        ## plot
        plot_panel(df_means_a, df_means_b, conditions)

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

def plot_panel(df_means_a, df_means_b, conditions: tuple=('cond1', 'cond2')):
    fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(9, 4))
    
    # plot first condition
    df_means_a.plot(ax=axes[0])
    axes[0].set_title(conditions[0], fontsize=18)

    # plot second condition
    df_means_b.plot(ax=axes[1])
    axes[1].set_title(conditions[1], fontsize=18)

    plt.show()

    return fig, axes


def get_mean_of_replicates(
    reps_dic: dict, dic_groups: dict, conditions: tuple = ("condition 1", "condition 2")
):
    """

    """
    for rep_num in reps_dic:

        cond_a_df = reps_dic[rep_num][0]
        dic_means_a = mean_gene_groups_of_sample(
            sample_df=cond_a_df, dic_groups=dic_groups
        )

    ### compute mean of samples and std (as shadow) between samples.
    ### https://seaborn.pydata.org/generated/seaborn.lineplot.html


if __name__ == "__main__":
    ## generate dictionary of all samples of experiment / 
    ## namedtuple of [conditions=('a', 'b'), reps_dic = {0:}]
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

    dic_list = {
        "hrde-1": ["isHrde1"],
        "pol-2": ["isPol2"],
        "highly": ["R1-SX_S14", 10],
        "all genes": ["ALL"],
    }

    dic_groups = gs.get_multiple_lists(dic_list)

    # oma_1_gene = ['WBGene00003864']
    # dic_groups = {'oma-1 gene': oma_1_gene}

    plot_replicates(
        reps_dic, dic_groups, conditions=("anti-gfp RNAi", "anti-oma-1 RNAi")
    )
