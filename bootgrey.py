def bootstrap_signal_df(
    df_sample: pd.DataFrame, group_size: int, num_of_iters: int = 1_000
):
    """"""
    gene_pool = pd.read_csv("tables/protein_coding_wbids.csv")

    df_boot_means = pd.DataFrame([])

    for i in range(num_of_iters):
        boot_group = random.sample(list(gene_pool["genes"]), group_size)
        intersected_list = list(set(df_sample.index) & set(boot_group))

        df_boot_means[i] = df_sample.loc[intersected_list, :].mean(axis=0)

    return df_boot_means


def plot_groups_signals_bootgrey(
    ATAC_exp,
    groups_dic: dict = {},
    mean_flag: bool = False,
    var_type: str = "none",
    add_highly_lowly: bool = True,
    bootgrey: bool = False,
    boot_size: int = 2315,
    boot_iters: int = 1000,
    drop_rep: int = 10,
    zscore_signal: bool = False,
    plot_range: tuple = (0, 0),
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

    if bootgrey:
        print(f"bootstrapping {boot_iters} iterations for group size {boot_size}...")

    if not mean_flag:
        for rep_i in range(ATAC_exp.num_of_reps):
            fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
            fig.suptitle(
                f"{ATAC_exp.exp_name}, Replicate {rep_i+1}, (variance between genes)",
                fontsize=14,
            )
            for cond_i in [
                0,
                1,
            ]:  # for each condition (sample) (get for this sample a df of means, and df of vars):
                sample_df = ATAC_exp.exp_df.iloc[rep_i, cond_i]
                if zscore_signal:  # later
                    sample_df = ut.normalize_zscore_df(sample_df)
                means_df, vars_df = groups_df_mean_and_var_dfs_for_sample(
                    sample_df, groups_dic, var_type
                )

                if plot_range.count(0) != 2:  # if range was given:
                    means_df = narrow_to_range(means_df, plot_range[0], plot_range[1])
                    vars_df = narrow_to_range(vars_df, plot_range[0], plot_range[1])

                axes[cond_i].set_title(f"{ATAC_exp.condition_names[cond_i]}")
                legend_flag = cond_i  # 0 / 1 [only legend on right ax]

                if bootgrey:  # add greys to ax and show
                    boot_dfs = bootstrap_signal_df(
                        sample_df, group_size=boot_size, num_of_iters=boot_iters
                    )
                    for boot_i in range(boot_iters):
                        axes[cond_i].plot(boot_dfs.iloc[:, boot_i], c="grey")

                return_plotted_ax(axes[cond_i], means_df, vars_df, legend_flag)

    else:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
        fig.suptitle(
            f"{ATAC_exp.exp_name}, Mean of All Replicates (variance between replicates)",
            fontsize=14,
        )

        for cond_i in [0, 1]:
            means_df_list = []
            for rep_i in range(ATAC_exp.num_of_reps):
                if rep_i == drop_rep:
                    print(f"dropped rep {rep_i}")
                    continue
                sample_df = ATAC_exp.exp_df.iloc[rep_i, cond_i]
                
                if zscore_signal:  # later
                    sample_df = ut.normalize_zscore_df(sample_df)
                
                means_df, _ = groups_df_mean_and_var_dfs_for_sample(
                    sample_df, groups_dic, var_type="none"
                )
                means_df_list.append(means_df)

            df_means_all_reps, df_vars = get_mean_variance_of_df_list(
                means_df_list, var_type
            )

            if plot_range.count(0) != 2:  # if range was given:
                df_means_all_reps = narrow_to_range(
                    df_means_all_reps, plot_range[0], plot_range[1]
                )
                df_vars = narrow_to_range(df_vars, plot_range[0], plot_range[1])

            axes[cond_i].set_title(f"{ATAC_exp.condition_names[cond_i]}")
            legend_flag = cond_i  # 0 / 1 [only legend on right ax]

            if bootgrey:  # add greys to ax and show
                if cond_i==0:
                    sample_mean_repeats, _ = get_mean_variance_of_df_list(ATAC_exp.cond1)
                elif cond_i==1:
                    sample_mean_repeats, _ = get_mean_variance_of_df_list(ATAC_exp.cond2)

                boot_dfs = bootstrap_signal_df(
                    sample_mean_repeats, group_size=boot_size, num_of_iters=boot_iters
                )
                for boot_i in range(boot_iters):
                    axes[cond_i].plot(boot_dfs.iloc[:, boot_i], c="grey")

            return_plotted_ax(axes[cond_i], df_means_all_reps, df_vars, legend_flag)
    
    plt.show()


def groups_df_mean_and_var_dfs_for_sample(df_sample, group_dic: dict, var_type="none"):
    """
    Gets a dictionary of gene groups. Gets for this sample:
    - means_df: col - group_name, row - location
    - vars_df: col - group_name, row - location.

    If var_type 'none', returns 0.
    """
    means_df = pd.DataFrame([])
    vars_df = pd.DataFrame([])

    for group_name in group_dic:
        means_df[group_name], vars_df[group_name] = group_mean_and_var_for_sample(
            df_sample, group_dic[group_name], var_type
        )

    if var_type == "none":
        vars_df = 0

    return means_df, vars_df


def group_mean_and_var_for_sample(df_sample, wbid_list, var_type="none"):
    """
    Gets a list of genes, and a df sample, returns the means and vars.

    return
    --------
    - group_mean: pd.Series, mean along gene
    - group_var: pd.Series, var along gene. If var_type='none', returns 0.
    """
    # take only genes that appear in the sample df:
    df_of_genes = ut.get_wbid_group_df(df_sample, wbid_list)
    group_mean = df_of_genes.mean()
    if var_type.lower() == "std":
        group_var = df_of_genes.std()
    elif var_type.lower() == "sem":
        group_var = df_of_genes.sem()
    elif var_type.lower() == "none":
        group_var = 0

    return group_mean, group_var


def add_highly_lowly_to_dic(dic_groups: dict):
    """
    Adds to the dic two groups: "highly expressed", "lowly expressed"
    changes inplace?
    """
    highly, lowly = ut.get_highly_lowly(prcnt=5)
    dic_groups["highly expressed (top 5%}"] = highly
    dic_groups["lowly expressed (bottom 5%}"] = lowly


def narrow_to_range(df, first_row, last_row):
    """
    Narrows df of signal by locations (relative to tss).
    Inclusive.
    """
    df.index = df.index.astype("int64")

    return df.loc[first_row:last_row, :]


def return_plotted_ax(ax, vec_df, var_df=0, legend_flag: bool = True):
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

    * returns ax - plots on ax.
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

    return ax


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



if __name__ == "__main__":
    ## load experiments first

    hrde1_regulated = ut.get_hrde_regulated(gs)
    hrde1_reg_intersected = ut.intersect_lists(
        hrde1_regulated, exp1.scores1.index
    )  # later: 15 / 151 missing

    hrde1_kennedy = gs.get_list("hrde-1-Kennedy")
    hrde1_kennedy_intersected = ut.intersect_lists(
        hrde1_kennedy, exp1.fc.index
    )  # later: 68 / 1527 missing

    plot_groups_signals_bootgrey(
        exp1,
        groups_dic={"regulated": hrde1_reg_intersected},
        bootgrey=True,
        boot_size=150,
    )

    plot_groups_signals_bootgrey(
        exp_hrde1,
        groups_dic={"regulated": hrde1_reg_intersected},
        bootgrey=True,
        boot_size=150,
    )

    plot_groups_signals_bootgrey(
        exp_hrde1,
        groups_dic={"regulated": hrde1_reg_intersected},
        mean_flag=True,
        bootgrey=True,
        boot_size=len(hrde1_reg_intersected),
    )

    plot_groups_signals_bootgrey(
        exp_hrde1,
        groups_dic={"kennedy": hrde1_kennedy_intersected},
        mean_flag=True,
        bootgrey=True,
        boot_size=len(hrde1_kennedy_intersected),
    )




