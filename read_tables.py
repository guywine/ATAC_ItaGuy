import pandas as pd
from typing import Union
import pathlib
import seaborn as sns


def read_experiment_to_df(exp_name: str = "exp1", pc_only_flag:bool=False):
    '''
    Gets an experiment name ['exp1' / 'exp_gonads' / 'exp_hrde_gonads' / 'exp_hrde_guy' / 'exp_metsetset'].\n
    Reads all its samples to df, and puts all of the dfs into a big df:
    - rows: rep_num (starts from 0)
    - cols: condition name [anti gfp / anti OMA-1] or [hrde-1;SX / SX]
    '''
    exp_dfs_cond1, exp_dfs_cond2 = read_experiment(exp_name, pc_only_flag)
    exp_df = create_exp_df(exp_dfs_cond1, exp_dfs_cond2, exp_name)

    return exp_df


def read_experiment(exp_name: str = "exp1", pc_only_flag:bool=False):
    """
    Gets an experiment name, reads all of the relevant tables to two lists, one for each condition.

    exp_names can be: 'exp1' / 'exp_gonads' / 'exp_hrde_gonads' / 'exp_hrde_guy' / 'exp_metsetset'

    Return
    -------
    - exp_dfs_cond1: list of dfs, each is a sample of the first condition (gfp / hrde-1;SX)
    - exp_dfs_cond1: list of dfs, each is a sample of the second condition (oma-1 / SX)
    """
    exps_dict = get_exps_f_names_dict()

    exp_dfs_cond1 = []
    exp_dfs_cond2 = []
    num_of_reps = len(exps_dict[exp_name][0])
    for rep_i in range(num_of_reps):
        # now:
        sample_df_0 = read_and_format_atac_table(exps_dict[exp_name][0][rep_i])
        sample_0_norm = normalize_sample(sample_df_0, exp_name, rep_i, 0, pc_only_flag)
        exp_dfs_cond1.append(sample_0_norm)
        ## before: exp_dfs_cond1.append(read_and_format_atac_table(exps_dict[exp_name][0][rep_i]))

        sample_df_1 = read_and_format_atac_table(exps_dict[exp_name][1][rep_i])
        sample_1_norm = normalize_sample(sample_df_1, exp_name, rep_i, 1, pc_only_flag)
        exp_dfs_cond2.append(sample_1_norm)
        ## exp_dfs_cond2.append(read_and_format_atac_table(exps_dict[exp_name][1][rep_i]))

    return exp_dfs_cond1, exp_dfs_cond2


def get_exps_f_names_dict():
    '''
    '''
    exp1_files = [
        [
            "DATA/exp1/ATAC_R0-iGFP.csv.gz",
            "DATA/exp1/ATAC_R1-iGFP.csv.gz",
            "DATA/exp1/ATAC_R2-iGFP.csv.gz",
            "DATA/exp1/ATAC_R4-iGFP.csv.gz",
        ],
        [
            "DATA/exp1/ATAC_R0-iOMA-1.csv.gz",
            "DATA/exp1/ATAC_R1-iOMA-1.csv.gz",
            "DATA/exp1/ATAC_R2-iOMA-1.csv.gz",
            "DATA/exp1/ATAC_R4-iOMA-1.csv.gz",
        ],
    ]

    exp_hrde_guy_files = [
        [
            'DATA/exp_hrde_guy/hrde1-sx-R1.zip',
            'DATA/exp_hrde_guy/hrde1-sx-R2.zip',
            'DATA/exp_hrde_guy/hrde1-sx-R3.zip'
        ],
        [
            'DATA/exp_hrde_guy/sx-R1.zip',
            'DATA/exp_hrde_guy/sx-R2.zip',
            'DATA/exp_hrde_guy/sx-R3.zip'
        ]
    ]
    exp_metsetset_files = [
        [
            'DATA/exp_metsetset/mss-gfp-R1.zip',
            'DATA/exp_metsetset/mss-gfp-R2.zip',
            'DATA/exp_metsetset/mss-gfp-R3.zip'
        ],
        [
            'DATA/exp_metsetset/mss-oma1-R1.zip',
            'DATA/exp_metsetset/mss-oma1-R2.zip',
            'DATA/exp_metsetset/mss-oma1-R3.zip'
        ]
    ]

    # exp_gonads_files = []  # to fill later
    # exp_hrde_gonads_files = []  # to fill later

    exps_dict = {
        "exp1": exp1_files,
        "exp_hrde_guy": exp_hrde_guy_files,
        "exp_metsetset": exp_metsetset_files,
    #    "exp_gonads": exp_gonads_files,
    #    "exp_hrde_gonads": exp_hrde_gonads_files,
    }

    return exps_dict



def create_exp_df(exp_dfs_cond1: list, exp_dfs_cond2: list, exp_name: str):
    """
    Gets two lists of dfs, each list is samples of a condition.
    Puts all df_samples in a single df of dfs:
        - Row: rep_number
        - Col: condition name
            * ['anti gfp' and 'anti OMA-1'] or [hrde-1;SX' and 'SX']
    """
    columns_dic = {"cond1": exp_dfs_cond1, "cond2": exp_dfs_cond2}
    exp_df = pd.DataFrame(columns_dic)

    if exp_name in ["exp1", "exp_gonads", "exp_metsetset"]:
        conds = {"cond1": "anti gfp", "cond2": "anti OMA-1"}
    elif exp_name in ["exp_hrde_gonads", "exp_hrde_guy"]:
        conds = {"cond1": "hrde-1;SX", "cond2": "SX"}
    else:
        NameError, "exp name not recognized"

    exp_df.rename(columns=conds, inplace=True)
    exp_df.index.name = "rep"

    ### add to avoid zero:
    exp_add_to_avoid_zero(exp_df)

    return exp_df


def read_and_format_atac_table(f_name: Union[pathlib.Path, str]):
    """
    Reads table of a single sample from file (if exists)
    Removes un needed columns.
    Sets wbid as index. [Removes rows with NaN in wbid]

    Return
    ---------
    - atac_table: pd.DataFrame
    """
    atac_table = read_atac_table(f_name)
    add_GFP_wbid(atac_table)
    drop_unneeded_columns(atac_table)
    atac_table.dropna(subset=["wbid"], inplace=True)
    index_wbid(atac_table)

    ## add_to_avoid_zero
    ## atac_table += 0.5

    return atac_table


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


def add_GFP_wbid(atac_df: pd.DataFrame):
    """
    Adds to the GFP transgene a wbid "GFP".

    * inplace.
    """
    GFP_ind = atac_df[atac_df["symbol"] == "GFP"].index
    atac_df.loc[GFP_ind, "wbid"] = "GFP"


def drop_unneeded_columns(atac_df: pd.DataFrame):
    """
    Drops 'feature_type' and 'symbol' columns.

    * inplace
    """
    atac_df.drop(columns=["feature_type", "symbol"], inplace=True)


def index_wbid(atac_df: pd.DataFrame):
    """
    Sets the wbid column as the index.

    * inplace
    """
    atac_df.set_index("wbid", inplace=True)


def read_sam_size_dic(pc_only_flag):
    '''
    Reads "all_aligned" table
    '''
    sam_size_dic = {}
    if pc_only_flag:
        sam_size_dic['exp1'] = pd.read_excel('DATA/exp1_pc_only.xlsx', index_col=0)
        sam_size_dic['exp_hrde_guy'] = pd.read_excel('DATA/exp_hrde_guy_pc_only.xlsx', index_col=0)
        sam_size_dic['exp_metsetset'] = pd.read_excel('DATA/exp_metsetset_pc_only.xlsx', index_col=0)
    else:
        sam_size_dic['exp1'] = pd.read_excel('DATA/exp1.xlsx', index_col=0)
        sam_size_dic['exp_hrde_guy'] = pd.read_excel('DATA/exp_hrde_guy.xlsx', index_col=0)
        sam_size_dic['exp_metsetset'] = pd.read_excel('DATA/exp_metsetset.xlsx', index_col=0)
    return sam_size_dic


def normalize_sample(sample_df, exp_name, rep_i, cond_i, pc_only_flag):
    '''
    '''
    sam_size_dic = read_sam_size_dic(pc_only_flag)
    all_aligned_num = sam_size_dic[exp_name].iloc[rep_i, cond_i]
    norm_sample = (sample_df/all_aligned_num)*1e6
    return norm_sample


def exp_add_to_avoid_zero(exp_df):
    '''
    Checks for each rep what is the minimum value of both conditions, and adds it to both dfs of rep.
    '''
    for rep_i in range(exp_df.shape[0]):
        min0 = find_min_after_zero(exp_df.iloc[rep_i,0])
        min1 = find_min_after_zero(exp_df.iloc[rep_i,1])
        min_of_both = min(min0, min1)

        a = exp_df.iloc[rep_i,0]
        b = exp_df.iloc[rep_i,1]
        a+=min_of_both
        b+=min_of_both


def find_min_after_zero(df):
    '''
    finds min of whole df after zero (assuming not all zero. and positive)
    '''
    return (df[df>0]).min().min()


def plot_hist_values(df, xlim=5, bins=15, logy=False, sns_flag=False):
    df_new = df[df<xlim]

    if sns_flag:
        sns.histplot(df_new,stat="probability")
    else:
        df_new.stack().plot.hist(bins=bins, logx=False, logy=logy, density=1)


def find_smallest_sample(df, n=1000):
    df_new = df[df>0]
    stacked_df = df_new.stack().reset_index(drop=True).drop_duplicates()
    small = stacked_df.nsmallest(n)
    print(small)
    return small






if __name__ == "__main__":
    if False:
        f_name = "DATA/exp1/ATAC_R0-iGFP.csv.gz"
        atac_table_initial = read_and_format_atac_table(f_name)

        exp1_df = read_experiment_to_df()

        exp1_df_pc_only = read_experiment_to_df(pc_only_flag=True)
    
    exps_dict = get_exps_f_names_dict()
    e1_g0 = read_and_format_atac_table(exps_dict['exp1'][0][0])
    e1_g1 = read_and_format_atac_table(exps_dict['exp1'][0][1])
    e1_g2 = read_and_format_atac_table(exps_dict['exp1'][0][2])
    e1_g3 = read_and_format_atac_table(exps_dict['exp1'][0][3])

    e1_o0 = read_and_format_atac_table(exps_dict['exp1'][1][0])
    e1_o1 = read_and_format_atac_table(exps_dict['exp1'][1][1])
    e1_o2 = read_and_format_atac_table(exps_dict['exp1'][1][2])
    e1_o3 = read_and_format_atac_table(exps_dict['exp1'][1][3])

    h1_h1 = read_and_format_atac_table(exps_dict['exp_hrde_guy'][0][0])
    h1_h2 = read_and_format_atac_table(exps_dict['exp_hrde_guy'][0][1])
    h1_h3 = read_and_format_atac_table(exps_dict['exp_hrde_guy'][0][2])

    h1_wt1 = read_and_format_atac_table(exps_dict['exp_hrde_guy'][1][0])
    h1_wt2 = read_and_format_atac_table(exps_dict['exp_hrde_guy'][1][1])
    h1_wt3 = read_and_format_atac_table(exps_dict['exp_hrde_guy'][1][2])

    mss1_g1 = read_and_format_atac_table(exps_dict['exp_metsetset'][0][0])
    mss1_g2 = read_and_format_atac_table(exps_dict['exp_metsetset'][0][1])
    mss1_g3 = read_and_format_atac_table(exps_dict['exp_metsetset'][0][2])

    mss1_o1 = read_and_format_atac_table(exps_dict['exp_metsetset'][1][0])
    mss1_o2 = read_and_format_atac_table(exps_dict['exp_metsetset'][1][1])
    mss1_o3 = read_and_format_atac_table(exps_dict['exp_metsetset'][1][2])

    s = find_smallest_sample(h1_h1)

    s.to_csv('hrde1_hrde_R1_smallest.csv', index=False)



    
    
    # exp_hrde_df = read_experiment_to_df('exp_hrde_guy')
    # exp_mss_df = read_experiment_to_df('exp_metsetset')

    # gfp0 = exp1_df.iloc[0,0]
    # sx_2 = exp_hrde_df.iloc[2,1]

    # sam_size_dic = read_sam_aize_dic()

    # gfp0_norm = normalize_sample(gfp0, 'exp1',0,0)
    # sx_2_norm = normalize_sample(sx_2, 'exp_hrde_guy', 2,1)





    

