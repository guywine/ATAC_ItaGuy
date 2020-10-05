'''
This will provide the user to choose which gene set to get.
- groups from table: by ATAC_bigtable (hrde-1 score, H3K9, mRNA level...)
- a set of the most / least 10% of ATAC_FC for a desired condition (choose percentage)
'''
import pandas as pd
import numpy as np

def read_big_table():
    '''
    Reads big_table of ATAC project, removes first header row, sets wbid to index.
    Replaces 0 with NaN in the desired columns.

    Return
    ---------
    - big_table: pd.DataFrame. index is wbid.
    '''
    big_table = pd.read_excel('ATAC_bigtable_1.xlsx', header=1)
    big_table.set_index('gene', inplace=True)

    # clean zeros from Ketting and Kennedy
    cols_to_drop_zero = ['isHrde1', 'ketting_score_all', 'ketting_score_20_23']
    big_table[cols_to_drop_zero] = big_table[cols_to_drop_zero].replace(0,np.nan)
    return big_table


def get_list(big_table, col_name:str, prcnt:float=0, upper:bool = True):
    '''
    Parameters
    ----------
    - col_name: str. name of columns

    Return
    ----------
    - wbid_list: list.
    '''
    col_orig = big_table[col_name].dropna()
    
    if prcnt:
        if upper:
            # get upper percnage
            quantile = col_orig.quantile(1-(prcnt/100))
            new_col = col_orig[col_orig>=quantile]    
        else:
            # get lower percentage
            quantile = col_orig.quantile((prcnt/100))
            new_col = col_orig[col_orig<=quantile]
    else:
        new_col = col_orig
    
    wbid_list = list(new_col.index)
    
    return wbid_list


if __name__=='__main__':
    big_table = read_big_table()

    ## add option of all genes
    col_name = 'isHrde1'
    ids_list = get_list(big_table, col_name)
    ids_list_2 = get_list(big_table, col_name, prcnt=10)
    ids_list_3 = get_list(big_table, col_name, prcnt=20, upper=False)

    print('done')