import pandas as pd 
import numpy as np
import read_tables as rt
import calc_signals as cas

class ATAC_signal():

    def __init__(self, exp_name: str ='exp1'):
        self.hotspot = (-500,-100)
        self.cond1, self.cond2 = rt.read_experiment(exp_name)
        self.exp_df = rt.create_exp_df(self.cond1, self.cond2, exp_name)
        self.mean1, _ = cas.get_mean_variance(self.cond1)
        self.mean2, _ = cas.get_mean_variance(self.cond2)

    
    def df_list_to_calc(self, df_list, calc_type: str='median'):
        '''
        Produces a df, each column is the calculated values for this replicate. 
        The value of Mean / Median is calculated only for the "hotspot" defined in the ATAC_signal object.

        Parameters
        ----------
        - df_list: list of df. Each df is a sample - its ATAC-signal results (-1000:1000). All samples must have the same indices.
        - calc_type: str ['median' / 'mean']

        return
        ----------
        - df_calc: pd.DataFrame. Row: gene. Column: calculated value of this replicate.
        '''
        num_of_reps = len(df_list)
        df_calc = pd.DataFrame([])
        for df_i in range(num_of_reps):
            rep_series = self.calc_hotspot_of_sample(df_list[df_i], calc_type)
            df_calc[f'rep {df_i}']=rep_series
        return df_calc

    def calc_hotspot_of_sample(self, signal_df: pd.DataFrame, calc_type: str):
        '''
        Produces a Series, with the calculated values for all genes.
        The value of Mean / Median is calculated only for the "hotspot" defined in the ATAC_signal object.

        Parameters
        ----------
        - sample_df: pd.DataFrame. A sample - Row:gene. Columns: ATAC-signal results (-1000:1000).
        - calc_type: str ['median' / 'mean']

        return
        ----------
        - rep_series: pd.Series. The calculated value for each gene.
        '''
        hot_inds = (self.hotspot[0]+1000, self.hotspot[1]+1000)
        if calc_type=='median':
            rep_series = signal_df.iloc[:,hot_inds[0]:hot_inds[1]].median(axis=1)
        elif calc_type=='mean':
            rep_series = signal_df.iloc[:,hot_inds[0]:hot_inds[1]].mean(axis=1)
        return rep_series
    
    def calc_fc_hotspot(self, rep_num: int, div_2_by_1: bool=True, log2: bool=True):
        '''
        Given a replicate number, calculates the FC-median (log2) parameter for every gene.
        Calculated for the "hotspot" defined in the ATAC_signal object.

        Parameters
        ----------
        - rep_num: int. number of replicate to do the operation for.
        - div_2_by_1: bool. If True, divide condition 2 by condition 1. If False, other way around.  (default True)
        - log2: bool. If True, return results in log2 scale. (default True)
        
        return
        ----------
        - median_FC_log2_series: pd.Series. For each gene calculated log2(FC-median). 
        '''
        df_FC_along_gene = self._fc_signal_of_rep(rep_num, div_2_by_1)
        median_FC_series = self.calc_hotspot_of_sample(signal_df=df_FC_along_gene, calc_type='median')
        if log2:
            median_FC_log2_series = np.log2(median_FC_series)
            return median_FC_log2_series
        return median_FC_series

    def _fc_signal_of_rep(self, rep_num: int, div_2_by_1: bool=True):
        '''
        Gets both conditions of replicates. Calculates log2_FC for every spot in every gene. 
        (Divides B / A)

        Parameters
        ---------
        - rep_num: int. number of replicate to do the operation for.
        - div_2_by_1: bool. If True, divide condition 2 by condition 1. If False, other way around.  (default True)

        return
        ---------
        - df_FC
        '''
        if div_2_by_1:
            mone = 1
        else:
            mone = 0
        
        df_a = self.exp_df.iloc[rep_num, 1-mone]
        df_b = self.exp_df.iloc[rep_num, mone]

        df_FC_along_gene = df_b.div(df_a)
        return df_FC_along_gene





if __name__=='__main__':
    exp1_atac = ATAC_signal('exp1')

    a1 = exp1_atac.exp_df.loc[1,'anti gfp']
    b1 = exp1_atac.exp_df.loc[1,'anti OMA-1']

    gfp_med = exp1_atac.df_list_to_calc(exp1_atac.cond1)
    gfp_mean = exp1_atac.df_list_to_calc(exp1_atac.cond1, 'mean')
    oma1_med = exp1_atac.df_list_to_calc(exp1_atac.cond2)
    oma1_mean = exp1_atac.df_list_to_calc(exp1_atac.cond2, 'mean')

    log2_fc_rep1 = exp1_atac._fc_signal_of_rep(1)


