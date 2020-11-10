import pandas as pd 
import read_tables as rt
import calc_signals as cas

class ATAC_signal():

    def __init__(self, exp_name: str ='exp1'):
        self.hotspot = (-500,-100)
        self.cond1, self.cond2 = rt.read_experiment(exp_name)
        self.df = rt.create_exp_df(self.cond1, self.cond2, exp_name)
        self.mean1, _ = cas.get_mean_variance(self.cond1)
        self.mean2, _ = cas.get_mean_variance(self.cond2)


    def df_list_to_calc(self, df_list, calc_type: str='median'):
        '''
        '''
        num_of_reps = len(df_list)
        df_meds = pd.DataFrame([])
        for df_i in range(num_of_reps):
            rep_series = self.calc_genes_of_sample(df_list[df_i], calc_type)
            df_calc[f'rep {df_i}']=rep_series
        return df_calc

    def calc_genes_of_sample(self, sample_df: pd.DataFrame, calc_type: str):
        hot_inds = (self.hotspot[0]+1000, self.hotspot[1]+1000)
        if calc_type=='median':
            rep_series = sample_df.iloc[:,hot_inds[0]:hot_inds[1]].median(axis=1)
        elif calc_type=='mean':
            rep_series = sample_df.iloc[:,hot_inds[0]:hot_inds[1]].mean(axis=1)
        return rep_series





if __name__=='__main__':
    exp1_atac = ATAC_signal('exp1')

