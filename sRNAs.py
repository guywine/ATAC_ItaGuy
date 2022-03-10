import pandas as pd 
import matplotlib.pyplot as plt

# class Small_RNA():

#     def __init__():


def read_sRNAs_and_calc_rpm(direction: str):
    '''
    direction in ['reverse', 'stranded']
    '''
    special_rows = ['__no_feature','__ambiguous','__too_low_aQual','__not_aligned','__alignment_not_unique']

    all_reps = []
    
    for rep_i in range(1,4):
        rpm_df = read_rep_and_calc_rpm(rep_i, direction, special_rows)
        all_reps.append(rpm_df)
    
    df_all_reps = pd.concat(all_reps, axis=1)
    return df_all_reps

def read_rep_and_calc_rpm(rep_num: int, direction: str, special_rows: list):
    '''
    '''
    fname = f'tables/secondary_sRNAs_Itai/wt_{direction}_R{rep_num}.txt'
    rep_df = pd.read_csv(fname, sep='\t', names=['wbid',f'R{rep_num}'], index_col=['wbid'])
    total_reads = rep_df.iloc[0:-3,0].sum() # skip last three rows

    rep_df.drop(special_rows, inplace=True)

    rpm_factor = 1e6/total_reads
    print(f'rep num {rep_num}, direction {direction}, factor = {rpm_factor}')

    rep_df[f'R{rep_num} rpm'] = rep_df[f'R{rep_num}'] * rpm_factor

    return rep_df

def get_df_all_reps_both_dirs():
    df_reverse = read_sRNAs_and_calc_rpm('reverse')
    df_forward = read_sRNAs_and_calc_rpm('stranded')
    both_dirs = df_reverse + df_forward
    return both_dirs


def get_rpm_df():
    df_all_reps = get_df_all_reps_both_dirs()
    rpm_df = df_all_reps.drop(['R1','R2','R3'], axis=1)
    rpm_means = rpm_df.mean(axis=1)
    return rpm_df

def thresh_by_rep(rpm_df: pd.DataFrame, thresh: int):
    '''
    '''
    rows_to_drop_i = []
    for row_i in range(rpm_df.shape[0]):
        if (rpm_df.iloc[row_i,:]<thresh).all():
            rows_to_drop_i.append(row_i)
    
    cutoff_df = rpm_df.drop(rpm_df.index[rows_to_drop_i])
    return cutoff_df



if __name__=='__main__':
    if True:
        rpm_df = get_rpm_df()
        over_5 = thresh_by_rep(rpm_df, 5)
        over_5['mean'] = over_5.mean(axis=1)

        mean_over_5_genes = over_5[over_5['mean']>=5]

        genes_some_over_5 = over_5.index.tolist()

        over_5.to_csv('sRNA_over_5_rpm.csv')
    
    x = pd.read_csv('tables/sRNA_over_5_rpm.csv', index_col='wbid')
    





