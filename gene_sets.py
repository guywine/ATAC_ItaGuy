"""
This will provide the user to choose which gene set to get.
- groups from table: by ATAC_bigtable (hrde-1 score, H3K9, mRNA level...)
- a set of the most / least 10% of ATAC_FC for a desired condition (choose percentage)
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3


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
        big_table = pd.read_excel("ATAC_bigtable_1.xlsx", header=1)
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
        f_name = "geneExprTable_YA.csv"
        exp_table = pd.read_csv(f_name)
        exp_table.set_index("geneID", inplace=True)
        exp_table.rename(
            columns={
                "median_adult Ce": "expression_median",
                "mean_adult Ce": "expression_mean",
            },
            inplace=True,
        )

        big_table_with_exp = big_table.join(exp_table, how="outer")
        return big_table_with_exp

    def get_list(self, col_name: str, prcnt: float = 0, bottom: bool = False):
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


if __name__ == "__main__":
    gs = Gene_sets()

    ## add option of all genes
    col_name = "isHrde1"
    # ids_list = at.get_list(col_name, 10, False)

    dic_list = {
        "hrde-1 all": ["hrde-1-Kennedy"],
        "pol-2": ["isPol2"],
        "highly 10%": ["expression_mean", 10],
        "lowly 10%": ["expression_mean", 10, True],
        "all genes": ["ALL"],
    }

    dic_groups = gs.get_multiple_lists(dic_list)

    print("done")

    # name1 = 'WBGene00268208'
    plot_venn_from_dic(dic_groups, list_of_names=["highly 10%", "hrde-1 all"])
