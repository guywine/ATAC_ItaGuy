import matplotlib.pyplot as plt
import pandas as pd 
from gene_id import Gene_IDs


#### to change to a function that plots a gene on a histogram:

# func1: from series, name, and gene list - return ax with histogram and markings.
# func2: from df with multiple series, use func 1 to create multiple panels, and a legend.


def plot_reps_hist_mark_gene(df_reps: pd.DataFrame, genes_to_mark):
    """
    Plots histogram for every column, marking the desired genes.

    Parameters
    ----------
    - df_reps: pd.DataFrame. Row:Gene. Column: Replicate. Values can be any calculated value for gene.
    - genes_to_mark: [string / list of strings]. Either gene names / Wbid. (e.g "oma-1" / ["WBGene00003864"])
    """
    genes_list = [genes_to_mark]  ### later: make_strings_a_list 'oma-1' -> ['oma-1']
    gid = Gene_IDs()  ### later
    list_of_wbids = [gid.to_wbid(gene) for gene in genes_list]  ### later
    num_of_reps = df_reps.shape[1]
    fig, axes = plt.subplots(1, num_of_reps, figsize=(num_of_reps * 5, 5))
    if num_of_reps == 1:
        axes.set_title(f"{df_reps.columns[0]}")
        axes.hist(df_reps.iloc[:, 0], bins=20, zorder=0)
        genes_x = df_reps.loc[list_of_wbids[0]][0]  ### later
        genes_y = 10  ### later
        hand = axes.scatter(genes_x, genes_y, c="red", marker=7, zorder=5)
    else:
        for rep_i in range(num_of_reps):
            # list_of_points = plot_values_for_genes(ax = axes[rep_i], value_series = df_reps.iloc[:,rep_i], list_of_indices =list_of_wbids)
            axes[rep_i].set_title(f"{df_reps.columns[rep_i]}")
            axes[rep_i].hist(df_reps.iloc[:, rep_i], bins=20, zorder=0)
            genes_x = df_reps.loc[list_of_wbids[0]][rep_i]  ### later
            genes_y = 10  ### later
            hand = axes[rep_i].scatter(genes_x, genes_y, c="red", marker=7, zorder=5)

    ## legend
    plt.show()


# def plot_values_for_genes(ax: plt.axis, value_series: pd.Series, list_of_indices: list):
#     '''
#     Gets a list of indices and plots them on the axis. also adds legend.
#     '''
#     value_list = [value_series[ind] for ind in list_of_indices]
#     y_values = [0]*len(value_list)

#     dic_points = {'oma-1':value_list[0], 'oma-2':value_list[1]}
#     df_dict = pd.DataFrame([dic_points])