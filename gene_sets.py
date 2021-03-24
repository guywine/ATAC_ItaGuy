"""
This will provide the user to choose which gene set to get.
- groups from table: by ATAC_bigtable (hrde-1 score, H3K9, mRNA level...)
- a set of the most / least 10% of ATAC_FC for a desired condition (choose percentage)
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

from gene_id import Gene_IDs

class Gene_sets:
    def __init__(self):
        self.big_table = self._read_big_table()
        self.all_wbids = list(self.big_table.index)

    def _read_big_table(self):
        """
        Reads big_table of ATAC project, removes first header row, sets wbid to index.
        Replaces 0 with NaN in the desired columns.

        Return
        ---------
        - big_table: pd.DataFrame. index is wbid.
        """
        big_table = pd.read_excel("tables/ATAC_bigtable_1.xlsx", header=1)
        big_table.set_index("gene", inplace=True)

        big_table_with_exp = self._add_expression_levels(big_table)

        # clean zeros from Ketting and Kennedy
        cols_to_drop_zero = [
            "hrde-1-Kennedy",
            "hrde-1-Ketting-all",
            "hrde-1-Ketting-20-23",
        ]
        big_table_with_exp[cols_to_drop_zero] = big_table_with_exp[
            cols_to_drop_zero
        ].replace(0, np.nan)
        return big_table_with_exp

    def _add_expression_levels(self, big_table):
        exp_df = self.get_expression_df()
        big_table_with_exp = big_table.join(exp_df, how="outer")
        return big_table_with_exp
    
    @staticmethod
    def get_expression_df():
        f_name = "tables/geneExprTable_YA.csv"
        exp_df = pd.read_csv(f_name)
        exp_df.set_index("geneID", inplace=True)
        exp_df.rename(
            columns={
                "median_adult Ce": "expression_median",
                "mean_adult Ce": "expression_mean",
            },
            inplace=True,
        )
        return exp_df


    def get_list(self, col_name: str, prcnt: float = 0, bottom: bool = False, thresh=None):
        """
        Parameters
        ----------
        - col_name: str. name of columns

        Return
        ----------
        - wbid_list: list.
        """
        if col_name.lower() == "all":
            return self.all_wbids

        col_orig = self.big_table[col_name].dropna()
        if prcnt:
            if not bottom:
                # get upper percnage
                quantile = col_orig.quantile(1 - (prcnt / 100))
                print(f"quantile upper of {col_name} = {quantile}")
                new_col = col_orig[col_orig > quantile]
            else:
                # get lower percentage
                quantile = col_orig.quantile((prcnt / 100))
                print(f"quantile lower of {col_name} = {quantile}")
                new_col = col_orig[col_orig < quantile]
        else:
            new_col = col_orig
        
        if thresh is not None:
            if not bottom:
                new_col = new_col[new_col > thresh]
            else:
                new_col = new_col[new_col < thresh]

        wbid_list = list(new_col.index)

        return wbid_list

    def get_multiple_lists(self, dic_list):
        """
        dic_list = {'hrde-1':['isHrde1', 10],
                'pol-2':['isPol2'],
                'lowly':['R1-SX_S14', 10, False],
                'all genes':['ALL']}

        return
        ---------
        - dic_groups: dic.
        """
        dic_groups = {}
        for lst_name in dic_list:
            dic_groups[lst_name] = self.get_list(*dic_list[lst_name])

        return dic_groups


def plot_venn_from_dic(dic_groups: dict, list_of_names: list):
    """
    Gets dic_groups, and list of lists names (len 2 or 3).
    Plots venn diagram of the matching lists
    """
    for name in list_of_names:
        if name not in dic_groups.keys():
            KeyError, f"group {name} is not found in the dic_groups"

    list_of_sets = [set(dic_groups[lst_name]) for lst_name in list_of_names]
    plot_ven(list_of_sets, list_of_names)


def plot_ven(list_of_sets: list, list_of_names: list):
    """
    Plots venn diagram for 2/3 sets.

    list_of_sets: list of lists or of sets
    """
    assert len(list_of_names) in [2, 3], "Venn diagram only works for 2/3 sets"
    assert len(list_of_names) == len(
        list_of_sets
    ), "Num of names does not match num of groups"

    if not all(isinstance(elem, set) for elem in list_of_sets):  # if some are not sets
        list_of_sets = [set(group) for group in list_of_sets]

    plt.figure()
    plt.title("Gene sets", fontsize=16)
    if len(list_of_names) == 2:
        venn2(subsets=(list_of_sets), set_labels=(list_of_names))
    if len(list_of_names) == 3:
        venn3(subsets=(list_of_sets), set_labels=(list_of_names))


def add_intersect(
    dic_groups: dict,
    group1: str,
    group2: str,
    inter_type: str = "intersection",
    replace: bool = False,
):
    """
    Takes a dictionary containing lists of genes, and operates on the two groups.
    Adds the resulted group to the dictionary with an automatic name.
    Can also replace the two original groups if chosen.

    Parameters
    ----------
    - dic_groups: dict. 'group_name':list_of_wbids
    - group1: str. name of group.
    - group2: str. name of group.
    - inter_type: str. 'intersection' / 'only first' / 'only second' / 'union'
    - replace: bool. Default: False. If True, removes the 2 original groups from the dict.

    Return
    ----------
    **changes dic_groups in place.**
    """
    ### Verify all group names appear in dictionary
    if group1 not in dic_groups.keys():
        KeyError, f"group {group1} is not found in the dic_groups"
    if group2 not in dic_groups.keys():
        KeyError, f"group {group2} is not found in the dic_groups"

    ### generate list
    new_list = intersect_lists(dic_groups[group1], dic_groups[group2], inter_type)

    ### generate new name
    if inter_type.lower() == "union":
        new_name = f"{group1} union {group2}"
    elif inter_type.lower() == "intersection":
        new_name = f"{group1} and {group2}"
    elif inter_type.lower() == "only first":
        new_name = f"{group1} (without {group2})"
    elif inter_type.lower() == "only second":
        new_name = f"{group2} (without {group1})"

    ### edit dictionary
    dic_groups[new_name] = new_list
    if replace:
        del dic_groups[group1], dic_groups[group2]


def intersect_lists(lst1, lst2, inter_type: str = "intersection"):
    """
    Takes two lists, returns a list of required relation.

    Parameters
    ----------
    - lst1: list
    - lst2: list
    - inter_type: str. 'intersection' / 'only first' / 'only second' / 'union'

    Return
    ------
    - list
    """
    assert inter_type.lower() in [
        "intersection",
        "only first",
        "only second",
        "union",
    ], f'intersection type "{inter_type}" not recognized'

    st1 = set(lst1)
    st2 = set(lst2)

    assert st1 != st2, "Lists sent to intersect contain the same set of elements"

    if inter_type.lower() == "union":
        return list(st1.union(st2))
    else:
        assert (
            len(st1.intersection(st2)) != 0
        ), "Lists sent to intersect have no intersection"
        if inter_type.lower() == "intersection":
            return list(st1.intersection(st2))
        elif inter_type.lower() == "only first":
            return list(st1.difference(st2))
        elif inter_type.lower() == "only second":
            return list(st2.difference(st1))

def wbid_list_to_names(wbid_list: list):
    '''
    '''
    gid = Gene_IDs()
    gene_names = [gid.to_name(wbid) for wbid in wbid_list if str(gid.to_name(wbid)) != 'nan']
    return gene_names

def gene_is_in_list(wbid_list:list , gene_name:str, print_flag: bool=False):
    '''
    '''
    gene_names = wbid_list_to_names(wbid_list)
    if gene_name in gene_names:
        if print_flag:    
            gid = Gene_IDs()
            wbid = gid.to_wbid(gene_name)
            print(f'gene {wbid} found')
        return True
    else:
        if print_flag:
            print('gene {gene_name} not found.')
        return False

if __name__ == "__main__":
    if "gs" not in locals():
        gs = Gene_sets()

    # dic_list = {
    #     "hrde-1": ["hrde-1-Kennedy"],
    #     "pol-2": ["isPol2"],
    #     "highly 10%": ["expression_mean", 10],
    #     "lowly 10%": ["expression_mean", 10, True],
    #     "all genes": ["ALL"],
    # }

    # dic_groups = gs.get_multiple_lists(dic_list)


    # name1 = 'WBGene00268208'
    # plot_venn_from_dic(dic_groups, list_of_names=["highly 10%", "hrde-1"])
    # add_intersect(dic_groups, 'highly 10%', 'hrde-1', inter_type='only second')


    ### create hrde-1 regulated list:
    hrde1_kennedy = gs.get_list('hrde-1-Kennedy')
    hrde_FC_sig = gs.get_list('mRNA_isSig')
    hrde_up = gs.get_list('mRNA_log2_FC', thresh=0)
    hrde_up_sig = intersect_lists(hrde_FC_sig, hrde_up)
    hrde_regulated = intersect_lists(hrde_up_sig, hrde1_kennedy)
    # df_mRNA_and_hrde1 = gs.big_table.loc[hrde_regulated,['mRNA_log2_FC','hrde-1-Kennedy']]

    # ### get strongest genes on both parameters
    # hrde1_strongest = df_mRNA_and_hrde1[df_mRNA_and_hrde1['hrde-1-Kennedy']>1200]
    # mRNA_FC_strongest = df_mRNA_and_hrde1[df_mRNA_and_hrde1['mRNA_log2_FC']>3]
    # hrde1_top = ut.intersect_lists(list(hrde1_strongest.index), list(mRNA_FC_strongest.index))
    # hrde1_top_names = ut.list_to_name(hrde1_top)

    # ### 

    # hrde1_top_germ_under5 = ut.intersect_lists(germline_under5, hrde1_top)



    

    # import utilities as ut 
    # ut.print_gene_ranks_in_df(df_mRNA_and_hrde1, 'oma-1', print_res=True)



    


