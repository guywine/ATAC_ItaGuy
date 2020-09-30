import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Union
import pathlib


def read_atac_table(f_name: Union[pathlib.Path, str]):
    """
    Gets path to table and reads it to pandas.
    Adds correct columns
    """
    f_path = pathlib.Path(f_name)
    if not f_path.exists():
        raise ValueError("File does not exist")
    col_names = ["wbid", "feature_type", "symbol"]
    col_names.extend(list(range(-1000, 1001)))
    atac_table = pd.read_csv(f_path, header=None, index_col=False, names=col_names)
    return atac_table


gfp_0 = read_atac_table("DATA/ATAC_R0-iGFP.csv")
gfp_1 = read_atac_table("DATA/ATAC_R1-iGFP.csv")
gfp_2 = read_atac_table("DATA/ATAC_R2-iGFP.csv")
gfp_4 = read_atac_table("DATA/ATAC_R4-iGFP.csv")
oma1_0 = read_atac_table('DATA/ATAC_R0-iOMA-1.csv')
oma1_1 = read_atac_table('DATA/ATAC_R1-iOMA-1.csv')
oma1_2 = read_atac_table('DATA/ATAC_R2-iOMA-1.csv')
oma1_4 = read_atac_table('DATA/ATAC_R4-iOMA-1.csv')


# shape: 20083 * 2004
# type:
# 19906 - protein_coding
# 13 - pseudogenes
# rest - NaN
