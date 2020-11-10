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


def drop_unneeded_columns(atac_df: pd.DataFrame):
    """
    Drops 'feature_type' and 'symbol' columns, inplace
    """
    atac_df.drop(columns=["feature_type", "symbol"], inplace=True)


def index_wbid(atac_df: pd.DataFrame):
    """
    Sets the wbid column as the index, inplace
    """
    atac_df.set_index("wbid", inplace=True)


def read_and_format_atac_table(f_name: Union[pathlib.Path, str]):
    """
    Reads table from file (if exists)
    Removes un needed columns
    Sets wbid as index

    Return
    ---------
    - atac_table: pd.DataFrame
    """
    atac_table = read_atac_table(f_name)
    drop_unneeded_columns(atac_table)
    index_wbid(atac_table)
    return atac_table


def read_experiment(exp_name: str = "exp1"):
    """
    Gets an experiment name, reads all of the relevant tables to two lists, one for each condition.

    exp_names can be: 'exp1' / 'exp_gonads' / 'exp_hrde_gonads' / 'exp_hrde_guy' / 'exp_metsetset'

    Return
    -------
    - exp_dfs_list: list of dfs, each is a sample of the desired experiment
    """
    #### file_names
    exp1_files = [
        ["DATA/ATAC_R0-iGFP.csv",
        "DATA/ATAC_R1-iGFP.csv",
        "DATA/ATAC_R2-iGFP.csv",
        "DATA/ATAC_R4-iGFP.csv"],
        ["DATA/ATAC_R0-iOMA-1.csv",
        "DATA/ATAC_R1-iOMA-1.csv",
        "DATA/ATAC_R2-iOMA-1.csv",
        "DATA/ATAC_R4-iOMA-1.csv",]
    ]
    exp_gonads_files = []  # to fill later
    exp_hrde_gonads_files = []  # to fill later
    exp_hrde_guy_files = []  # to fill later
    exp_metsetset_files = []  # to fill later

    exps_dict = {
        "exp1": exp1_files,
        "exp_gonads": exp_gonads_files,
        "exp_hrde_gonads": exp_hrde_gonads_files,
        "exp_hrde_guy": exp_hrde_guy_files,
        "exp_metsetset": exp_metsetset_files,
    }

    exp_dfs_cond1 = []
    exp_dfs_cond2 = []
    num_of_reps = len(exps_dict[exp_name][0])
    for rep_i in range(num_of_reps):
        exp_dfs_cond1.append(read_and_format_atac_table(exps_dict[exp_name][0][rep_i]))
        exp_dfs_cond2.append(read_and_format_atac_table(exps_dict[exp_name][1][rep_i]))

    return exp_dfs_cond1, exp_dfs_cond2


def create_exp_dic(exp_dfs_cond1: list, exp_dfs_cond2: list, exp_name: str):
    """"""
    exp_dic = {}
    if exp_name == "exp1": # first GFP, than OMA-1
        exp_dic[0] = [exp_dfs_cond1[0], exp_dfs_cond2[0]]
        exp_dic[1] = [exp_dfs_cond1[1], exp_dfs_cond2[1]]
        exp_dic[2] = [exp_dfs_cond1[2], exp_dfs_cond2[2]]
        exp_dic[4] = [exp_dfs_cond1[3], exp_dfs_cond2[3]]
    else:  # if other experiments than the order is the same:
        for rep_num in range(3):
            exp_dic[rep_num] = [exp_dfs_cond1[rep_num], exp_dfs_cond2[rep_num]]
    return exp_dic

def create_exp_df(exp_dfs_cond1: list, exp_dfs_cond2: list, exp_name: str):
    '''
    '''
    columns_dic = {'cond1':exp_dfs_cond1, 'cond2':exp_dfs_cond2}
    exp_df = pd.DataFrame(columns_dic)
    
    if exp_name in ['exp1', 'exp_gonads', 'exp_metsetset']:
        conds = {'cond1':'anti gfp','cond2':'anti OMA-1'}
    elif exp_name in ['exp_hrde_gonads', 'exp_hrde_guy']:
        conds = {'cond1':'SX','cond2':'hrde-1;SX'}
    else:
        NameError, 'exp name not recognized'
    
    exp_df.rename(columns=conds, inplace=True)
    exp_df.index.name = 'rep'

    return exp_df


def read_experiment_to_dic(exp_name: str = "exp1"):
    """"""
    exp_dfs_cond1, exp_dfs_cond2 = read_experiment(exp_name)
    exp_dic = create_exp_dic(exp_dfs_cond1, exp_dfs_cond2, exp_name)
    return exp_dic

def read_experiment_to_df(exp_name: str = "exp1"):
    exp_dfs_cond1, exp_dfs_cond2 = read_experiment(exp_name)
    exp_df = create_exp_df(exp_dfs_cond1, exp_dfs_cond2, exp_name)
    return exp_df


if __name__ == "__main__":
    # exp_dfs_cond1, exp_dfs_cond2 = read_experiment(exp_name="exp1")

    # exp1_dic = read_experiment_to_dic(exp_name="exp1")
    exp1_df = read_experiment_to_df()

    # plt.plot(gfp_0.iloc[0])
