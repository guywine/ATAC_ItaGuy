'''
Plot ATAC-seq signal. Choose:
- Experiment: exp1 / exp_gonads / exp_hrde_gonads / exp_hrde_guy / exp_metsetset
- set of genes: list of sets ([hrde-1 set, csr-1 set])
- Normalize: highly_lowly / none
- Mean all replicates or display them seperately?
'''
import read_tables as rdt 
import matplotlib.pyplot as plt
import pandas as pd


def generate_mock_groups():
    gene_list_1 = ['WBGene00007063', 'WBGene00007064', 'WBGene00007067']
    gene_list_2 = ['WBGene00017071', 'WBGene00019895', 'WBGene00009583', 'WBGene00018682']

    dic_groups = {'group 1':gene_list_1, 'group 2':gene_list_2}
    return dic_groups

def mean_gene_groups_of_sample(sample_df: pd.DataFrame, dic_groups: dict):
    '''
    Gets df of sample and groups and returns a dict with mean of these groups.

    Parameters
    -----------
    - sample_df: pd.DataFrame of ATAC-seq signal
    - dic_groups: dict, group_name : list_of_ids

    Return
    -----------
    - dic_means: dict of pd.Series, group_name : mean_of_all_genes_in_group
    '''
    dic_means = {}
    for group in dic_groups:
        group_ids = dic_groups[group]
        dic_means[group] = gfp_0.loc[group_ids].mean()
    return dic_means
    
def plot_groups(dic_means):
    '''
    Gets a dict of means and plots.
    '''
    df_groups = pd.DataFrame(dic_means)
    df_groups.plot()


if __name__=='__main__':
    exp1_dfs_list = rdt.read_experiment(exp_name='exp1')
    gfp_0, gfp_1, gfp_2, gfp_4, oma1_0, oma1_1, oma1_2, oma1_4 = exp1_dfs_list

    dic_groups = generate_mock_groups()

    dic_means = mean_gene_groups_of_sample(sample_df=gfp_0, dic_groups=dic_groups)
    plot_groups(dic_means)








