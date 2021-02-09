import matplotlib.pyplot as plt
import pandas as pd 
from gene_id import Gene_IDs
import utilities as ut
import seaborn as sns


def plot_signal_gene(ATAC_exp, gene_name: str, mean_flag: bool = True, var_type: str = 'none'):
    '''
    Plots for a single gene one of two options:
    - Single panel with mean and variance of this gene for all replicates (line per condition)
    - Panel for each replicate, with the gene's signal (line per condition)

    Parameters
    --------
    - ATAC_exp: ATAC_signal object. Contains ATAC-seq data of chosen experiment.
    - gene_name: str. 
    - mean_flag: bool. If True, plots all replicates meaned on the same panel.
    - var_type: str. ['std'/'sem'/'none']. Which variance to plot. Only relevant if mean_flag True.
    '''

    if mean_flag:
        gene_means, gene_vars = ATAC_exp.get_gene_mean_and_var_both_conditions(gene_name, var_type)
        fig, ax0 = plt.subplots()
        plt.title(f"Signal for gene: {gene_name},\t{ATAC_exp.exp_name}", fontsize=14)
        plt.xlabel('Location relative to TSS')
        plt.ylabel('ATAC-seq signal (norm.)')

        plot_ax(ax0, gene_means, gene_vars)

    else:
        wbid = ATAC_exp.gid.to_wbid(gene_name)
        num_reps = ATAC_exp.num_of_reps

        list_of_rep_dfs = []

        for rep_i in range(num_reps):
            gene_rep_df = pd.DataFrame([])
            gene_rep_df[ATAC_exp.condition_names[0]] = ATAC_exp.cond1[rep_i].loc[wbid,:]
            gene_rep_df[ATAC_exp.condition_names[1]] = ATAC_exp.cond2[rep_i].loc[wbid,:]

            list_of_rep_dfs.append(gene_rep_df)

        fig, axes = plt.subplots(1,num_reps, figsize=(num_reps*6,5), sharey=True)
        for rep_i in range(num_reps):
            axes[rep_i].set_title(f'replicate {rep_i+1}')
            ## later - add y title and x title
            legend_flag = False
            if rep_i==num_reps-1:
                legend_flag=True
            plot_ax(axes[rep_i], list_of_rep_dfs[rep_i], legend_flag)


def plot_ax(ax, vec_df, var_df=0, legend_flag:bool = True):
    '''
    Plots all lines in vec_df columns (with legend)
    Adds fill_between of var_df if not 0.
    Adds legend if specified.

    Parameters
    ----------
    - ax: ax to plot on.
    - vec_df: pd.DataFrame. Comtains vectors. Col names - vectors names.
    - var_df: pd.DataFrame. Equal in structure to "vec_df". Contains variance data.
    - legend_flag: bool. If true, plot legend on this axis.

    * No return - plots on ax.
    '''
    lines = []
    colors = plt.get_cmap('Accent') # later (define outside?)

    if isinstance(var_df, int): # if no varince df given
        for col_i in range(vec_df.shape[1]):
            line = plot_vector(vec=vec_df.iloc[:,col_i], ax=ax, color=colors(col_i))
            lines.append(line)
    else:   # if svariance given
        for col_i in range(vec_df.shape[1]):
            line = plot_vector(vec=vec_df.iloc[:,col_i], ax=ax, color=colors(col_i), var_vec=var_df.iloc[:,col_i])
            lines.append(line)
    
    if legend_flag:
        first_leg = plt.legend(lines, vec_df.columns)
        ax = plt.gca().add_artist(first_leg)
    

def plot_vector(vec, ax, color, var_vec=0):
    """
    Plots a line on the graph, with variance as a band.
    if variance_vec==0, doesn't plot variance.
    """
    ## print(f'variance_vec: {variance_vec}') ### later

    line = ax.plot(vec, c=color)
    if not isinstance(var_vec, int):
        ax.fill_between(vec.index, vec-var_vec, vec+var_vec, fc=color, alpha=0.3)
    return line[0]





####################################################################
#### to change to a function that plots a gene on a histogram:

# func1: from series, name, and gene list - return ax with histogram and markings.
# func2: from df with multiple series, use func 1 to create multiple panels, and a legend.

def plot_gene_atac_signal_histogram(ATAC_exp, gene_to_mark: str, mean_flag: bool = True):
    if mean_flag:
        mean_df = pd.DataFrame({"mean_FC": ATAC_exp.fc.mean(axis=1)})
        plot_reps_hist_mark_gene(df_reps=mean_df, genes_to_mark=gene_to_mark)
        ut.get_gene_rank(mean_df.iloc[:, 0], gene_to_mark)
    else:
        print(f"gene {gene_to_mark}:")
        plot_reps_hist_mark_gene(df_reps=ATAC_exp.fc, genes_to_mark=gene_to_mark)


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


######################################