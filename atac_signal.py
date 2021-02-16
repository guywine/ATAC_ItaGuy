import utilities as ut
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt  # later: remove
import read_tables as rt
from gene_id import Gene_IDs
import plotting as my_plots

# from sklearn import preprocessing

class ATAC_signal:
    def __init__(self, exp_name: str = "exp1", var_type: str = "std"):
        self.hotspot = (-500, -100)  # user to define
        self.add_to_avoid_zero_division = 1  # user to define

        self.cond1, self.cond2 = rt.read_experiment(exp_name)
        self.exp_df = rt.create_exp_df(self.cond1, self.cond2, exp_name)

        self.mean1, self.var1 = ut.calc_mean_variance_of_dfs(self.cond1, var_type)
        self.mean2, self.var2 = ut.calc_mean_variance_of_dfs(self.cond2, var_type)

        self.scores1 = self.calc_hotspot_scores(
            self.cond1
        )  # median signal of all genes, in all reps
        self.scores2 = self.calc_hotspot_scores(
            self.cond2
        )  # median signal of all genes, in all reps
        
        self.num_of_reps = self.exp_df.shape[0]
        self.exp_name = exp_name
        self.condition_names = tuple(self.exp_df.columns)

        self.fc = self.generate_FC_median_df()

        self.gid = Gene_IDs()

    def calc_hotspot_scores(self, df_list, calc_type: str = "median"):
        """
        Produces a df: row - gene, column - rep_num.
        The value of Mean / Median is calculated only for the "hotspot" defined in the Class object.

        Parameters
        ----------
        - df_list: list of df. Each df is a sample - its ATAC-signal results (-1000:1000). All samples must have the same indices.
        - calc_type: str ['median' / 'mean']

        return
        ----------
        - df_calc: pd.DataFrame. Row - gene. Column - rep_num. Value - as calculated.
        """
        num_of_reps = len(df_list)
        df_calc = pd.DataFrame([])
        for df_i in range(num_of_reps):
            rep_series = self.calc_node_hotspot_of_sample(df_list[df_i], calc_type)
            df_calc[f"rep {df_i}"] = rep_series
        return df_calc

    def calc_node_hotspot_of_sample(self, signal_df: pd.DataFrame, calc_type: str):
        """
        For a single sample df, produces a Series, with the calculated value for all genes (in this rep).
        The value of Mean / Median is calculated only for the "hotspot" defined in the Class object.

        Parameters
        ----------
        - sample_df: pd.DataFrame. A sample - Row:gene. Columns: ATAC-signal results (-1000:1000).
        - calc_type: str ['median' / 'mean']

        return
        ----------
        - rep_series: pd.Series. The calculated value for each gene.
        """
        hot_inds = (self.hotspot[0] + 1000, self.hotspot[1] + 1000)
        if calc_type == "median":
            rep_series = signal_df.iloc[:, hot_inds[0] : hot_inds[1]].median(axis=1)
        elif calc_type == "mean":
            rep_series = signal_df.iloc[:, hot_inds[0] : hot_inds[1]].mean(axis=1)
        return rep_series

    def generate_FC_median_df(self, div_2_by_1: bool = True, log2: bool = True):
        """
        Generats a df. Row: gene, Column: replicate. 
        Value is the fc_median_parameter, calculated so:
        * Fold change of signals
        * median of all values in hotspot
        Values calculated for "hotspot" defined in the Class object.

        Parameters
        ----------
        - div_2_by_1: bool. If True, divide condition 2 by condition 1. If False, other way around. (default True)
        - log2: bool. If True, return FC results in log2 scale. (default True)

        return
        ----------
        - df_FC_median: pd.DataFrame. Row: gene, Column: replicate. Value is the fc_median_parameter.
        """
        df_FC_median = pd.DataFrame([])
        for rep_i in range(self.num_of_reps):
            median_FC_series = self.calc_median_fc_hotspot_of_sample(
                rep_i, div_2_by_1, log2
            )
            df_FC_median[f"rep {rep_i}"] = median_FC_series
        return df_FC_median

    def calc_median_fc_hotspot_of_sample(
        self, rep_num: int, div_2_by_1: bool = True, log2: bool = True
    ):
        """
        Given a replicate number, calculates the FC between conditions and then median (log2) for every gene.
        Calculated for the "hotspot" defined in the ATAC_signal object.

        Parameters
        ----------
        - rep_num: int. number of replicate to do the operation for.
        - div_2_by_1: bool. If True, divide condition 2 by condition 1. If False, other way around.  (default True)
        - log2: bool. If True, return results in log2 scale. (default True)

        return
        ----------
        - median_FC_log2_series: pd.Series. For each gene calculated log2(FC-median).
        """
        FC_df = self._fc_signal_of_rep(rep_num, div_2_by_1, log2)

        median_FC_series = self.calc_node_hotspot_of_sample(
            signal_df=FC_df, calc_type="median"
        )

        return median_FC_series

    def _fc_signal_of_rep(
        self, rep_num: int, div_2_by_1: bool = True, log2: bool = True
    ):
        """
        For a given rep_num, calculates FC for every spot in every gene.
        (Divides B / A)
        Also outputs as log2.

        Parameters
        ---------
        - rep_num: int. number of replicate to do the operation for.
        - div_2_by_1: bool. If True, divide condition 2 by condition 1. If False, other way around.  (default True)

        return
        ---------
        - df_FC_along_gene: Row - gene. Col - position. Value: FC bewteen conditions.
        """
        ## choose direction:
        if div_2_by_1:
            mone = 1
        else:
            mone = 0

        ## add to avoid zero:
        df_a = self.exp_df.iloc[rep_num, 1 - mone] + self.add_to_avoid_zero_division
        df_b = self.exp_df.iloc[rep_num, mone] + self.add_to_avoid_zero_division

        ## divide:
        df_FC_along_gene = df_b.div(df_a)

        ## log2:
        if log2:
            df_FC_along_gene = np.log2(df_FC_along_gene)

        return df_FC_along_gene

    def _update_mean_and_var_of_exp(self, var_type: str = "std"):
        """
        Updates for this object the mean and var of all repeats.
        """
        self.mean1, self.var1 = ut.calc_mean_variance_of_dfs(self.cond1, var_type)
        self.mean2, self.var2 = ut.calc_mean_variance_of_dfs(self.cond2, var_type)

    ###############################################################
    #### from here it is getting the data needed to plot lines ####
    ###############################################################

    def get_gene_mean_and_var_both_conditions(
        self, gene_name: str, var_type: str = "std"
    ):
        """
        """
        wbid = self.gid.to_wbid(gene_name)
        mean_1, var_1 = self.get_gene_mean_and_var_for_cond(
            cond_num=1, wbid=wbid, var_type=var_type
        )
        mean_2, var_2 = self.get_gene_mean_and_var_for_cond(
            cond_num=2, wbid=wbid, var_type=var_type
        )

        cond_names = self.exp_df.columns

        gene_means = pd.DataFrame([])
        gene_means[cond_names[0]] = mean_1
        gene_means[cond_names[1]] = mean_2

        gene_vars = pd.DataFrame([])
        gene_vars[cond_names[0]] = var_1
        gene_vars[cond_names[1]] = var_2

        return gene_means, gene_vars

    def get_gene_mean_and_var_for_cond(
        self, cond_num: int, wbid: str, var_type: str = "std"
    ):
        """
        Get for a gene the mean and variance for this condition.

        Parameters
        --------
        - cond_num: int. [1 / 2]
        - gene_wbid: str.
        - var_type: str ['std' / 'sem' / 'none']

        Return
        --------
        - mean_series: pd.Series.
        - var_series: pd.Series.
        """
        if cond_num == 1:
            gene_reps = get_gene_replicates(
                self.cond1, wbid
            )  # later where is this function
        elif cond_num == 2:
            gene_reps = get_gene_replicates(
                self.cond2, wbid
            )  # later where is this function

        mean_series = gene_reps.mean()
        if var_type.lower() == "std":
            var_series = gene_reps.std()

        elif var_type.lower() == "sem":
            var_series = gene_reps.sem()

        elif var_type.lower() == "none":
            var_series = 0

        return mean_series, var_series
    
    def get_gean_fc_df(self, gene_name: str, div_2_by_1: bool = True, log2: bool = True):
        '''
        '''
        wbid = self.gid.to_wbid(gene_name)

        gene_fc_df = pd.DataFrame([])
        for rep_i in range(self.num_of_reps):
            cond1_sig = self.cond1[rep_i].loc[wbid,:]+self.add_to_avoid_zero_division
            cond2_sig = self.cond2[rep_i].loc[wbid,:]+self.add_to_avoid_zero_division
            if div_2_by_1:
                gene_fc_df[rep_i] = cond2_sig / cond1_sig
            else:
                gene_fc_df[rep_i] = cond1_sig / cond2_sig
        
        if log2:
            gene_fc_df = np.log2(gene_fc_df)
        
        return gene_fc_df



def get_gene_replicates(list_of_dfs: list, gene_wbid: str):
    """
    Gets a df for: Rows: rep_num. Cols: location.
    """
    gene_reps_df = pd.DataFrame([])
    for rep_num in range(len(list_of_dfs)):
        gene_reps_df = gene_reps_df.append(list_of_dfs[rep_num].loc[gene_wbid, :])

    gene_reps_df.reset_index(inplace=True, drop=True)

    return gene_reps_df


def mean_and_var_gene_list_for_signal_df(
    signal_df: pd.DataFrame, wbid_list: list, var_type: str = "std"
):
    """
    Caclculates the mean and variance of ATAC-signal of a wbid_list for the desired sample,
    along the gene.

    Parameters
    ----------
    - signal_df: pd.DataFrame. Rows: wbids, Cols (2001): location. Value - ATAC-signal.
    - wbid_list: list of wbids to calculate for.
    - var_type: std.

    Return
    ---------
    - mean_series: pd.Series (len 2001)
    - var_series: pd.Series (len 2001)
    """
    # only genes that appear in the sample:
    intersected_list = list(set(signal_df.index) & set(wbid_list))

    our_genes_df = signal_df.loc[intersected_list, :]

    mean_series = our_genes_df.mean()

    if var_type.lower() == "std":
        var_series = our_genes_df.std()
    elif var_type.lower() == "sem":
        var_series = our_genes_df.sem()
    elif var_type.lower() == "none":
        var_series = 0

    return mean_series, var_series


if __name__ == "__main__":
    import utilities as ut

    if "exp1" not in locals():
        exp1 = ATAC_signal("exp1")

    if "exp_mss" not in locals():
        exp_mss = ATAC_signal("exp_metsetset")

    # my_plots.plot_signal_gene(exp1, 'oma-1', var_type='std')

    # gene_fc_df_log2 = exp1.get_gean_fc_df('GFP')
    # gene_fc_df = exp1.get_gean_fc_df('GFP', log2=False)

    # my_plots.plot_fc_gene(exp1, exp_mss, 'oma-1')

    # genes_to_test = [
    #     "oma-1",
    #     "oma-2",
    #     "GFP",
    #     "C09G9.5",
    #     "spr-2",
    #     "F14E5.8",
    #     "F14E5.1",
    #     "efl-3",
    #     "unc-119",
    #     "rad-26",
    #     "C27B7.2",
    #     "npax-4",
    #     "C09G9.8",
    # ]

    # my_plots.plot_gene_atac_signal_histogram(exp1, "oma-1", False)

    # for gene in genes_to_test:
    #     print("Regular:")
    #     my_plots.plot_gene_atac_signal_histogram(exp1, gene)
        # print('met;set;set mutant:')
        # plot_gene_atac_signal(exp_mss, gene)
