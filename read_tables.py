import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Union
import pathlib
import seaborn as sns


def read_atac_table(f_name: Union[pathlib.Path, str]):
    """
    Gets path to table and reads it to pandas.
    Adds correct columns

    Return
    ---------
    - atac_table: pd.DataFrame
    
    """
    f_path = pathlib.Path(f_name)
    if not f_path.exists():
        raise ValueError("File does not exist")
    col_names = ["wbid", "feature_type", "symbol"]
    col_names.extend(list(range(-1000, 1001)))
    atac_table = pd.read_csv(f_path, header=None, index_col=False, names=col_names)
    return atac_table

def drop_unneeded_columns(atac_df:pd.DataFrame):
    '''
    Drops 'feature_type' and 'symbol' columns, inplace
    '''
    atac_df.drop(columns=["feature_type", "symbol"], inplace=True)

def index_wbid(atac_df:pd.DataFrame):
    '''
    Sets the wbid column as the index, inplace
    '''
    atac_df.set_index('wbid', inplace=True)

def read_and_format_atac_table(f_name: Union[pathlib.Path, str]):
    '''
    Reads table from file (if exists)
    Removes un needed columns
    Sets wbid as index

    Return
    ---------
    - atac_table: pd.DataFrame
    '''
    atac_table = read_atac_table(f_name)
    drop_unneeded_columns(atac_table)
    index_wbid(atac_table)
    return atac_table

def read_experiment(exp_name='exp1': str):
    '''
    Gets an experiment name, reads all of the relevant tables to a list.

    exp_names can be: 'exp1' / 'exp_gonads' / 'exp_hrde_gonads' / 'exp_hrde_guy' / 'exp_metsetset'
    
    Return
    -------
    - list_of_tables
    '''
    exp1_files = ["DATA/ATAC_R0-iGFP.csv", "DATA/ATAC_R1-iGFP.csv", "DATA/ATAC_R2-iGFP.csv", "DATA/ATAC_R4-iGFP.csv", 'DATA/ATAC_R0-iOMA-1.csv', 'DATA/ATAC_R1-iOMA-1.csv', 'DATA/ATAC_R2-iOMA-1.csv', 'DATA/ATAC_R4-iOMA-1.csv']
    exp_gonads_files = [] # to fill later
    exp_hrde_gonads_files = [] # to fill later
    exp_hrde_guy_files = [] # to fill later
    exp_metsetset_files = [] # to fill later

    exps_dict = {'exp1':exp1_files, 'exp_gonads':exp_gonads_files, 'exp_hrde_gonads':exp_hrde_gonads_files, 'exp_hrde_guy':exp_hrde_guy_files, 'exp_metsetset':exp_metsetset_files}

if __name__=='__main__':
    gfp_0 = read_and_format_atac_table("DATA/ATAC_R0-iGFP.csv")
    gfp_1 = read_and_format_atac_table("DATA/ATAC_R1-iGFP.csv")
    gfp_2 = read_and_format_atac_table("DATA/ATAC_R2-iGFP.csv")
    gfp_4 = read_and_format_atac_table("DATA/ATAC_R4-iGFP.csv")
    oma1_0 = read_and_format_atac_table('DATA/ATAC_R0-iOMA-1.csv')
    oma1_1 = read_and_format_atac_table('DATA/ATAC_R1-iOMA-1.csv')
    oma1_2 = read_and_format_atac_table('DATA/ATAC_R2-iOMA-1.csv')
    oma1_4 = read_and_format_atac_table('DATA/ATAC_R4-iOMA-1.csv')
    
    plt.plot(gfp_0.iloc[0])
# shape: 20083 * 2004
# type:
# 19906 - protein_coding
# 13 - pseudogenes
# rest - NaN
