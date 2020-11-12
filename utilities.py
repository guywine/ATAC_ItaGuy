import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3


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