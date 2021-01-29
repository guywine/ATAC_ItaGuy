import utilities as ut
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt # later: remove
import read_tables as rt
from gene_id import Gene_IDs
import plotting as my_plots
# from sklearn import preprocessing


class ATAC_signal:
    def __init__(self, exp_name: str = "exp1"):
        self.hotspot = (-500, -100)  # user to define
        self.add_to_avoid_zero_division = 1  # user to define

        self.cond1, self.cond2 = rt.read_experiment(exp_name)
        self.exp_df = rt.create_exp_df(self.cond1, self.cond2, exp_name)

        self.mean1, _ = ut.calc_mean_variance_of_dfs(self.cond1)
        self.mean2, _ = ut.calc_mean_variance_of_dfs(self.cond2)

        self.scores1 = self.calc_hotspot_scores(
            self.cond1
        )  # median signal of all genes, in all reps
        self.scores2 = self.calc_hotspot_scores(
            self.cond2
        )  # median signal of all genes, in all reps

        self.fc = self.generate_FC_median_df()

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
        Generats a df. Row: gene, Column: replicate. Value is the fc_median_parameter.
        Values calculated for "hotspot" defined in the Class object.

        Parameters
        ----------
        - div_2_by_1: bool. If True, divide condition 2 by condition 1. If False, other way around. (default True)
        - log2: bool. If True, return FC results in log2 scale. (default True)

        return
        ----------
        - df_FC_median: pd.DataFrame. Row: gene, Column: replicate. Value is the fc_median_parameter.
        """
        num_of_reps = self.exp_df.shape[0]

        df_FC_median = pd.DataFrame([])
        for rep_i in range(num_of_reps):
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
        
        # if log2:
        #     median_FC_log2_series = np.log2(median_FC_series)
        #     return median_FC_log2_series
        
        return median_FC_series

    def _fc_signal_of_rep(self, rep_num: int, div_2_by_1: bool = True, log2: bool = True):
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
        df_a = (
            self.exp_df.iloc[rep_num, 1 - mone] + self.add_to_avoid_zero_division
        ) 
        df_b = (
            self.exp_df.iloc[rep_num, mone] + self.add_to_avoid_zero_division
        ) 

        ## divide:
        df_FC_along_gene = df_b.div(df_a)

        ## log2:
        if log2:
            df_FC_along_gene = np.log2(df_FC_along_gene)

        return df_FC_along_gene




if __name__ == "__main__":
    import utilities as ut

    if "exp1" not in locals():
        exp1 = ATAC_signal("exp1")

    a1 = exp1.exp_df.loc[1, "anti gfp"]
    b1 = exp1.exp_df.loc[1, "anti OMA-1"]

    fc_oma1_by_gfp = exp1.generate_FC_median_df()
    fc_oma1_by_gfp_no_log = exp1.generate_FC_median_df(log2=False)
    fc_gfp_by_oma1 = exp1.generate_FC_median_df(div_2_by_1=False)

    gid = Gene_IDs()
    oma1 = gid.to_wbid("oma-1")
    oma2 = gid.to_wbid("oma-2")

    my_plots.plot_reps_hist_mark_gene(df_reps=fc_gfp_by_oma1, genes_to_mark="oma-1")

    gfp_by_oma1_mean = pd.DataFrame({"mean_FC": fc_gfp_by_oma1.mean(axis=1)})
    gfp_by_oma1_mean_123 = pd.DataFrame(
        {"mean_FC": fc_gfp_by_oma1.iloc[:, 1:4].mean(axis=1)}
    )

    print("oma-1, FC normalized from 4 replicates")
    my_plots.plot_reps_hist_mark_gene(df_reps=gfp_by_oma1_mean, genes_to_mark="oma-1")
    print("oma-2, FC normalized from 4 replicates")
    my_plots.plot_reps_hist_mark_gene(df_reps=gfp_by_oma1_mean_123, genes_to_mark="oma-2")

    ut.get_gene_rank(gfp_by_oma1_mean_123.iloc[:, 0], "oma-2")

    #### more genes:
    my_plots.plot_reps_hist_mark_gene(df_reps=gfp_by_oma1_mean_123, genes_to_mark="unc-119")
    ut.get_gene_rank(gfp_by_oma1_mean_123.iloc[:, 0], "unc-119")
