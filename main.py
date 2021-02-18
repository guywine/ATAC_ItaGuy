import pandas as pd 
from gene_id import Gene_IDs
import plotting as my_plots
import utilities as ut
from atac_signal import ATAC_signal

if __name__=='__main__':
    david_df = pd.read_csv('tables/MutAccum_Number_Mutations_per_Gene.csv')
    david_list = list(david_df['Gene'])
    dic_david = {"David's mutated genes":david_list}

    if "exp1" not in locals():
        exp1 = ATAC_signal("exp1")

    if "exp_mss" not in locals():
        exp_mss = ATAC_signal("exp_metsetset")

    my_plots.plot_groups_signals(exp1, dic_david, mean_flag=True)
    my_plots.plot_groups_signals(exp_mss, dic_david)

