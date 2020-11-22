import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import read_tables as rt
import calc_signals as cas
from gene_id import Gene_IDs
from sklearn import preprocessing


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
            rep_series = self.calc_node_hotspot_of_sample(df_list[df_i], calc_type)
            df_calc[f'rep {df_i}']=rep_series
        return df_calc
    
    def calc_node_hotspot_of_sample(self, signal_df: pd.DataFrame, calc_type: str):
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
    
    def generate_FC_median_df(self, div_2_by_1: bool=True, log2: bool=True):
        '''
        Generats a df. Row: gene, Column: replicate. Value is the fc_median_parameter.
        Values calculated for "hotspot" defined in the object.

        Parameters
        ----------
        - div_2_by_1: bool. If True, divide condition 2 by condition 1. If False, other way around.  (default True)
        - log2: bool. If True, return results in log2 scale. (default True)

        return
        ----------
        - df_FC_median: pd.DataFrame. Row: gene, Column: replicate. Value is the fc_median_parameter.
        '''
        num_of_reps = self.exp_df.shape[0]
        df_FC_median = pd.DataFrame([])
        for rep_i in range(num_of_reps):
            median_FC_series = self.calc_median_fc_hotspot_of_sample(rep_i, div_2_by_1, log2)
            df_FC_median[f'rep {rep_i}']=median_FC_series
        return df_FC_median

    
    def calc_median_fc_hotspot_of_sample(self, rep_num: int, div_2_by_1: bool=True, log2: bool=True):
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
        median_FC_series = self.calc_node_hotspot_of_sample(signal_df=df_FC_along_gene, calc_type='median')
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
        add_to_avoid_zero = 1 # define by user

        if div_2_by_1:
            mone = 1
        else:
            mone = 0
        
        df_a = self.exp_df.iloc[rep_num, 1-mone] + add_to_avoid_zero
        df_b = self.exp_df.iloc[rep_num, mone] + add_to_avoid_zero

        df_FC_along_gene = df_b.div(df_a)
        return df_FC_along_gene

def plot_reps_hist_mark_gene(df_reps: pd.DataFrame, genes_to_mark):
    '''
    Plots histogram for every columns, marking the asked genes.

    Parameters
    ----------
    - df_reps: pd.DataFrame. Row:Gene. Column: Replicate. Values can be any calculated value for gene.
    - genes_to_mark: [string / list of strings]. Either gene names / Wbid. (e.g "oma-1" / ["WBGene00003864"])
    '''
    genes_list = [genes_to_mark] ### later: make_strings_a_list 'oma-1' -> ['oma-1']
    gid = Gene_IDs() ### later
    list_of_wbids = [gid.to_wbid(gene) for gene in genes_list] ### later 
    num_of_reps = df_reps.shape[1]
    fig, axes = plt.subplots(1, num_of_reps, figsize = (num_of_reps*5,5))
    if num_of_reps==1:
        axes.set_title(f'{df_reps.columns[0]}')
        axes.hist(df_reps.iloc[:,0], bins=20, zorder=0)
        genes_x = df_reps.loc[list_of_wbids[0]][0] ### later
        genes_y = 10 ### later
        hand = axes.scatter(genes_x, genes_y, c='red', marker=7, zorder=5)
    else:
        for rep_i in range(num_of_reps):
            # list_of_points = plot_values_for_genes(ax = axes[rep_i], value_series = df_reps.iloc[:,rep_i], list_of_indices =list_of_wbids)
            axes[rep_i].set_title(f'{df_reps.columns[rep_i]}')
            axes[rep_i].hist(df_reps.iloc[:,rep_i], bins=20, zorder=0)
            genes_x = df_reps.loc[list_of_wbids[0]][rep_i] ### later
            genes_y = 10 ### later
            hand = axes[rep_i].scatter(genes_x, genes_y, c='red', marker=7, zorder=5)

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




if __name__=='__main__':
    import utilities as ut 

    if 'exp1' not in locals():
        exp1 = ATAC_signal('exp1')

    a1 = exp1.exp_df.loc[1,'anti gfp']
    b1 = exp1.exp_df.loc[1,'anti OMA-1']

    fc_oma1_by_gfp = exp1.generate_FC_median_df()
    fc_oma1_by_gfp_no_log = exp1.generate_FC_median_df(log2=False)
    fc_gfp_by_oma1 = exp1.generate_FC_median_df(div_2_by_1=False)

    gid = Gene_IDs()
    oma1 = gid.to_wbid('oma-1')
    oma2 = gid.to_wbid('oma-2')

    plot_reps_hist_mark_gene(df_reps=fc_gfp_by_oma1, genes_to_mark='oma-1')

    gfp_by_oma1_mean = pd.DataFrame({'mean_FC':fc_gfp_by_oma1.mean(axis=1)})
    gfp_by_oma1_mean_123 = pd.DataFrame({'mean_FC':fc_gfp_by_oma1.iloc[:,1:4].mean(axis=1)})

    print('oma-1, FC normalized from 4 replicates')
    plot_reps_hist_mark_gene(df_reps=gfp_by_oma1_mean, genes_to_mark='oma-1')
    print('oma-2, FC normalized from 4 replicates')
    plot_reps_hist_mark_gene(df_reps=gfp_by_oma1_mean, genes_to_mark='oma-2')



    ut.get_gene_rank(gfp_by_oma1_mean.iloc[:,0],'oma-2')
    ut.get_gene_rank(gfp_by_oma1_mean_123.iloc[:,0],'GFP')


    #### atac scores:
    if False:
        gfp_medians = exp1.df_list_to_calc(exp1.cond1)
        oma1_medians = exp1.df_list_to_calc(exp1.cond2)

        gfp_score = gfp_medians.mean(axis=1)
        oma1_score = oma1_medians.mean(axis=1)

        my_atac_scores = pd.DataFrame({'anti-OMA-1':oma1_score, 'anti-GFP':gfp_score})

        ut.print_gene_ranks_in_df(my_atac_scores, 'oma-1', print_res=True)











