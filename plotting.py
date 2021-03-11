import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
from gene_id import Gene_IDs
import utilities as ut
import seaborn as sns
import calc_signals as cas


def plot_fc_gene(
    ATAC_exp1,
    ATAC_exp2,
    gene_name: str,
    std_flag: bool = True,
    log2_flag: bool = True,
    plot_range: tuple = (0, 0),
):
    """
    Plots for a single gene panels, each line is the FC across the gene in one experiment.
    - Single panel with FC (then mean) of all reps.
    - Panel for each replicate.

    Parameters
    --------
    - ATAC_exp1: ATAC_signal object. Contains ATAC-seq data of chosen experiment.
    - ATAC_exp2: "     "       "
    - gene_name: str.
    - var_type: str. ['std'/'sem'/'none']. Which variance to plot. Only relevant if mean_flag True.
    """
    gene_fc_1 = ATAC_exp1.get_gean_fc_df(gene_name, log2=log2_flag)
    gene_fc_2 = ATAC_exp2.get_gean_fc_df(gene_name, log2=log2_flag)

    gene_fc_means = pd.DataFrame([])
    gene_fc_means[ATAC_exp1.exp_name] = gene_fc_1.mean(axis=1)
    gene_fc_means[ATAC_exp2.exp_name] = gene_fc_2.mean(axis=1)

    if plot_range.count(0) != 2:  # if range was given:
        gene_fc_means = narrow_to_range(gene_fc_means, plot_range[0], plot_range[1])

    if std_flag:  # later - not working
        gene_fc_std = pd.DataFrame([])
        gene_fc_std[ATAC_exp1.exp_name] = gene_fc_1.std(axis=1)
        gene_fc_std[ATAC_exp2.exp_name] = gene_fc_2.std(axis=1)
        if plot_range.count(0) != 2:  # if range was given:
            gene_fc_std = narrow_to_range(gene_fc_std, plot_range[0], plot_range[1])
    else:
        gene_fc_std = 0

    fig, ax0 = plt.subplots()
    plt.title(f"Fold-Change of signal for gene: {gene_name}", fontsize=14)
    plt.xlabel("Location relative to TSS")
    plt.ylabel("FC of ATAC-seq signal (norm.)")

    plot_ax(ax0, gene_fc_means, gene_fc_std)


def plot_signal_gene(
    ATAC_exp,
    gene_name: str,
    mean_flag: bool = True,
    var_type: str = "none",
    plot_range: tuple = (0, 0),
):
    """
    Plots for a single gene one of two options:
    - Single panel with mean and variance of this gene for all replicates (line per condition)
    - Panel for each replicate, with the gene's signal (line per condition)

    Parameters
    --------
    - ATAC_exp: ATAC_signal object. Contains ATAC-seq data of chosen experiment.
    - gene_name: str.
    - mean_flag: bool. If True, plots all replicates meaned on the same panel.
    - var_type: str. ['std'/'sem'/'none']. Which variance to plot. Only relevant if mean_flag True.
    """

    if mean_flag:
        gene_means, gene_vars = ATAC_exp.get_gene_mean_and_var_both_conditions(
            gene_name, var_type
        )
        fig, ax0 = plt.subplots()
        plt.title(f"Signal for gene: {gene_name}, {ATAC_exp.exp_name}", fontsize=14)
        plt.xlabel("Location relative to TSS")
        plt.ylabel("ATAC-seq signal (norm.)")

        if plot_range.count(0) != 2:  # if range was given:
            gene_means = narrow_to_range(gene_means, plot_range[0], plot_range[1])
            gene_vars = narrow_to_range(gene_vars, plot_range[0], plot_range[1])

        plot_ax(ax0, gene_means, gene_vars)

    else:
        wbid = ATAC_exp.gid.to_wbid(gene_name)
        num_reps = ATAC_exp.num_of_reps

        list_of_rep_dfs = []

        for rep_i in range(num_reps):
            gene_rep_df = pd.DataFrame([])
            gene_rep_df[ATAC_exp.condition_names[0]] = ATAC_exp.cond1[rep_i].loc[
                wbid, :
            ]
            gene_rep_df[ATAC_exp.condition_names[1]] = ATAC_exp.cond2[rep_i].loc[
                wbid, :
            ]

            if plot_range.count(0) != 2:  # if range was given:
                gene_rep_df = narrow_to_range(gene_rep_df, plot_range[0], plot_range[1])

            list_of_rep_dfs.append(gene_rep_df)

        fig, axes = plt.subplots(1, num_reps, figsize=(num_reps * 6, 5), sharey=True)
        for rep_i in range(num_reps):
            axes[rep_i].set_title(f"replicate {rep_i+1}")
            ## later - add y title and x title
            legend_flag = False
            if rep_i == num_reps - 1:
                legend_flag = True
            plot_ax(axes[rep_i], list_of_rep_dfs[rep_i], legend_flag)


def plot_groups_signals(
    ATAC_exp,
    groups_dic: dict,
    mean_flag: bool = False,
    var_type: str = "std",
    add_highly_lowly: bool = True,
    bootstrap: bool = False,
    boot_size: int = 2315,
    boot_iters: int = 1000
):
    """
    Takes an experiment, a dictionary with groups, plots panel with two axes.
    - one ax for each condition: each group is a line with std.

    If mean_flag: variance is between reps. if seperate: variance is between genes.

    Parameters
    --------
    - ATAC_exp: ATAC_signal object.
    - groups_dic: keys - group_name, values - list of wbids.
    """
    if add_highly_lowly:
        add_highly_lowly_to_dic(groups_dic)

    if bootstrap:
        print(f'bootstrapping {boot_iters} iterations...')

    if not mean_flag:
        for rep_i in range(ATAC_exp.num_of_reps):
            fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
            fig.suptitle(f"{ATAC_exp.exp_name}, Replicate {rep_i+1}, (variance between genes)", fontsize=14)
            for cond_i in [0,1]: # for each condition (sample) (get for this sample a df of means, and df of vars):
                sample_df = ATAC_exp.exp_df.iloc[rep_i,cond_i]
                means_df, vars_df = groups_df_mean_and_var_dfs_for_sample(sample_df, groups_dic, var_type)

                if bootstrap:
                    means_df[f'bootstrap ({boot_size} genes)'], vars_df[f'bootstrap ({boot_size} genes)'] = cas.bootstrap_atac_signal(sample_df, group_size=boot_size, num_of_iters=boot_iters) # later

                axes[cond_i].set_title(f"{ATAC_exp.condition_names[cond_i]}")
                legend_flag = cond_i # 0 / 1 [only legend on right ax]
                plot_ax(axes[cond_i], means_df, vars_df, legend_flag)

    else:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
        fig.suptitle(f"{ATAC_exp.exp_name}, Mean of All Replicates (variance between replicates)", fontsize=14)

        for cond_i in [0,1]:
            means_df_list = []
            for rep_i in range(ATAC_exp.num_of_reps):
                sample_df = ATAC_exp.exp_df.iloc[rep_i, cond_i]
                means_df, _ = groups_df_mean_and_var_dfs_for_sample(sample_df, groups_dic, var_type='none')
                if bootstrap:
                    means_df[f'bootstrap ({boot_size} genes)'], _ = cas.bootstrap_atac_signal(sample_df, group_size=boot_size, num_of_iters=boot_iters) # later
                means_df_list.append(means_df)
            
            df_means_all_reps, df_vars = get_mean_variance_of_df_list(means_df_list, var_type)

            axes[cond_i].set_title(f"{ATAC_exp.condition_names[cond_i]}")
            legend_flag = cond_i # 0 / 1 [only legend on right ax]
            plot_ax(axes[cond_i], df_means_all_reps, df_vars, legend_flag)


def groups_df_mean_and_var_dfs_for_sample(df_sample, group_dic: dict, var_type="std"):
    '''
    Gets a dictionary of gene groups. Gets for this sample:
    - means_df: col - group_name, row - location
    - vars_df: col - group_name, row - location. If var_type='none', returns 0.
    '''
    means_df = pd.DataFrame([])
    vars_df = pd.DataFrame([])
    
    for group_name in group_dic:
        means_df[group_name], vars_df[group_name] = group_mean_and_var_for_sample(df_sample, group_dic[group_name], var_type)

    if var_type=='none':
        vars_df=0
    
    return means_df, vars_df


def group_mean_and_var_for_sample(df_sample, wbid_list, var_type="std"):
    """
    Gets a list of genes, and a df sample, returns the means and vars.

    return
    --------
    - group_mean: pd.Series, mean along gene
    - group_var: pd.Series, var along gene. If var_type='none', returns 0.
    """
    # take only genes that appear in the sample df:
    intersected_list = list(set(df_sample.index) & set(wbid_list))

    df_of_genes = df_sample.loc[intersected_list,:]
    group_mean = df_of_genes.mean()
    if var_type.lower() == "std":
        group_var = df_of_genes.std()
    elif var_type.lower() == "sem":
        group_var = df_of_genes.sem()
    elif var_type.lower() == "none":
        group_var = 0
    
    return group_mean, group_var


def add_highly_lowly_to_dic(dic_groups: dict):
    '''
    Adds to the dic two groups: "highly expressed", "lowly expressed"
    changes inplace?
    '''
    highly, lowly = ut.get_highly_lowly(prcnt=5)
    dic_groups['highly expressed (top 5%}'] = highly
    dic_groups['lowly expressed (bottom 5%}'] = lowly

def narrow_to_range(df, first_row, last_row):
    """
    Narrows df of signal by locations (relative to tss).
    Inclusive.
    """
    df.index = df.index.astype("int64")

    return df.loc[first_row:last_row, :]


def plot_ax(ax, vec_df, var_df=0, legend_flag: bool = True):
    """
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
    """
    lines = []
    colors = plt.get_cmap("Set1")  # later (define outside?)

    if isinstance(var_df, int):  # if no varince df given
        for col_i in range(vec_df.shape[1]):
            line = plot_vector(vec=vec_df.iloc[:, col_i], ax=ax, color=colors(col_i))
            lines.append(line)
    else:  # if svariance given
        for col_i in range(vec_df.shape[1]):
            line = plot_vector(
                vec=vec_df.iloc[:, col_i],
                ax=ax,
                color=colors(col_i),
                var_vec=var_df.iloc[:, col_i],
            )
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
        ax.fill_between(
            vec.index.astype("int64"), vec - var_vec, vec + var_vec, fc=color, alpha=0.3
        )
    return line[0]


def get_mean_variance_of_df_list(df_list: list, var_type: str = "none"):
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
    if var_type.lower() == "std":
        df_variance = pd.concat(df_list).groupby(level=0).std()
    elif var_type.lower() == "sem":
        df_variance = pd.concat(df_list).groupby(level=0).sem()
    elif var_type.lower() == "none":
        df_variance = 0

    return df_mean_all_reps, df_variance


####################################################################
#### to change to a function that plots a gene on a histogram:

# func1: from series, name, and gene list - return ax with histogram and markings.
# func2: from df with multiple series, use func 1 to create multiple panels, and a legend.



def plot_gene_atac_signal_distribution(
    ATAC_exp, gene_to_mark: str, mean_flag: bool = True, plot_type: str = 'violin', zscore: bool=False
):
    '''
    - plot_type: str ['violin' / 'hist']
    '''
    gid = Gene_IDs()
    wbid = gid.to_wbid(gene_to_mark)

    if zscore:
        df_ready = normalize_zscore_df(ATAC_exp.fc)
    else:
        df_ready = ATAC_exp.fc

    # df_ready = ATAC_exp.fc

    if mean_flag:
        mean_df = pd.DataFrame({"mean_FC": df_ready.mean(axis=1)})
        
        gene_mean = mean_df.loc[wbid]
        gene_std = ATAC_exp.fc.loc[wbid].std()

        plot_df_cols_mark_gene(df_reps=mean_df, gene_to_mark=gene_to_mark, plot_type=plot_type, std_whiskers=gene_std)
        ut.get_gene_rank(mean_df.iloc[:, 0], gene_to_mark)
    
    else:
        print(f"gene {gene_to_mark}:")
        plot_df_cols_mark_gene(df_reps=df_ready, gene_to_mark=gene_to_mark, plot_type=plot_type)


def plot_df_cols_mark_gene(df_reps: pd.DataFrame, gene_to_mark, plot_type: str = 'violin', std_whiskers:float=0):
    """
    Plots histogram / violinplot for every column, marking the desired genes.

    Parameters
    ----------
    - df_reps: pd.DataFrame. Row:Gene. Column: Replicate. Values can be any calculated value for gene.
    - genes_to_mark: [string / list of strings]. Either gene names / Wbid. (e.g "oma-1" / ["WBGene00003864"])
    - plot: str ['hist' / 'violin']
    """
    gid = Gene_IDs()  ### later
    wbid = gid.to_wbid(gene_to_mark)  ### later
    num_of_reps = df_reps.shape[1]
    fig, axes = plt.subplots(1, num_of_reps, figsize=(num_of_reps * 5, 5))
    ax_now = axes
    for rep_i in range(num_of_reps):
        if num_of_reps>1:
            ax_now = axes[rep_i]
        ax_now.set_title(f"{df_reps.columns[rep_i]}")

        if plot_type=='hist':
            ax_now.hist(df_reps.iloc[:, rep_i], bins=20, zorder=0)
            genes_x = df_reps.loc[wbid][rep_i]  ### later
            genes_y = 10  ### later
            hand = ax_now.scatter(genes_x, genes_y, c="red", marker=7, zorder=5)

        elif plot_type=='violin':
            ax_now.violinplot(df_reps.iloc[:, rep_i], showextrema=False, quantiles=[0.05,0.5,0.95])
            genes_y = df_reps.loc[wbid][rep_i]  ### later
            genes_x = 1  ### later
            hand = ax_now.scatter(genes_x, genes_y, c="red", marker='.', zorder=5)

            if std_whiskers:
                ax_now.add_patch(Rectangle((0.95, genes_y-std_whiskers),0.1, std_whiskers*2, alpha=0.3, facecolor='purple'))

    plt.show()
        
