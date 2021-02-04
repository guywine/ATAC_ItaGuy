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


def plot_gene_atac_signal(exp: ATAC_signal, gene_to_mark: str, mean_flag: bool=True):
    if mean_flag:
        mean_df = pd.DataFrame({"mean_FC": exp.fc.mean(axis=1)})
        my_plots.plot_reps_hist_mark_gene(df_reps=mean_df, genes_to_mark=gene_to_mark)
        ut.get_gene_rank(mean_df.iloc[:, 0], gene_to_mark)
    else:
        print(f'gene {gene_to_mark}:')
        my_plots.plot_reps_hist_mark_gene(df_reps=exp.fc, genes_to_mark=gene_to_mark)


if __name__ == "__main__":
    import utilities as ut

    if "exp1" not in locals():
        exp1 = ATAC_signal("exp1")
        exp_mss = ATAC_signal('exp_metsetset')
    


    # plot_gene_atac_signal(exp1, 'oma-1', mean_flag=False)

    genes_to_test = ['oma-1','oma-2', 'GFP','C09G9.5','spr-2','F14E5.8','F14E5.1','efl-3','unc-119','rad-26','C27B7.2','npax-4','C09G9.8']

    
    for gene in genes_to_test:
        print('Regular:')
        plot_gene_atac_signal(exp1, gene)
        # print('met;set;set mutant:')
        # plot_gene_atac_signal(exp_mss, gene)

    

