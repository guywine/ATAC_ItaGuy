def scatter_genes_both_conds(
    ATAC_exp, marked_list: list = [], shown_value: str = "sum"
):
    """
    - shown_value: str ['sum' / 'score']
    Plots a panel with scatter for each rep:
    - Each dot is a gene.
    - X: sum in cond1
    - Y: sum in cond2

    If marked_list given, all the genes in list are colored in red.
    """
    if not set(marked_list).issubset(ATAC_exp.scores1.index):
        raise KeyError(
            'some of the genes in the "marked" list were not found in the tables'
        )

    if shown_value == "sum":
        df_reps_cond1, df_reps_cond2 = ATAC_exp.gene_sum_tables
        suptitle = f"Gene sums, exp {ATAC_exp.exp_name}"

    if shown_value == "score":
        df_reps_cond1 = ATAC_exp.scores1
        df_reps_cond2 = ATAC_exp.scores2
        suptitle = f"Gene scores, exp {ATAC_exp.exp_name}"

    _scatter_reps(
        ATAC_exp, suptitle, df_reps_cond1, df_reps_cond2, marked_list
    )


def _scatter_reps(
    ATAC_exp, suptitle: str, df_reps_cond1, df_reps_cond2, marked_list: list = []
):
    """
    Gets two dfs with rep value for each gene.

    Plots a panel with scatter for each rep:
    - Each dot is a gene.
    - X: some value of cond1
    - Y: some value of cond1

    If marked_list given, all the genes in list are colored in red.
    """
    num_of_reps = df_reps_cond1.shape[1]

    fig, axes = plt.subplots(1, num_of_reps, figsize=(num_of_reps * 5, 5))

    # titles and labels:
    fig.suptitle(suptitle, fontsize=16)
    fig.text(0.5, 0.04, f"{ATAC_exp.condition_names[0]}", ha="center")
    axes[0].set_ylabel(f"{ATAC_exp.condition_names[1]}")

    for rep_i in range(ATAC_exp.num_of_reps):
        ax_now = axes[rep_i]
        ax_now.set_title(f"{df_reps_cond1.columns[rep_i]}")

        genes_x = df_reps_cond1.iloc[:, rep_i]
        genes_y = df_reps_cond2.iloc[:, rep_i]
        ax_now.scatter(genes_x, genes_y, s=2, c="cornflowerblue")

        if marked_list:  # if list not empty

            marked_xs = df_reps_cond1.loc[marked_list, f"rep {rep_i}"]
            marked_ys = df_reps_cond2.loc[marked_list, f"rep {rep_i}"]
            ax_now.scatter(marked_xs, marked_ys, s=5, c="r")

    plt.show()
