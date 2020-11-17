import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from sklearn import preprocessing

from gene_id import Gene_IDs
from gene_sets import Gene_sets
from mRNA_gonads import Table_mRNA
from Ahringer import Ahringer



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

def print_gene_ranks_in_df(gene_df: pd.DataFrame, gene:str, print_res: bool=False):
    '''
    '''
    gene_name = Gene_IDs().to_name(gene)
    print(f'gene - {gene_name}')
    for column in gene_df.columns:
        if print_res:
            print(f'{column}')
        perc = get_gene_rank(gene_df[column], gene, print_res, print_gene=False)
        if not print_res:
            print(f'{column}:\t{perc:.2f}%')
    print('\n')


def get_gene_rank(gene_series: pd.Series, gene:str, print_res: bool=True, print_gene: bool=True):
    '''
    Series where index is gene_name.

    Return
    --------
    - percentile: float. 
    '''
    # get gene name and wbid:
    if gene.lower()=='gfp':
        wbid='GFP'
        name = wbid
    else:
        gid = Gene_IDs()
        wbid = gid.to_wbid(gene)
        name = gid.to_name(gene)
        if str(name)=='nan':
            name = wbid
    
    gene_rank_df = pd.DataFrame({'value':gene_series, 'rank':gene_series.rank()})

    if ('GFP' in gene_series.index[0]) or ('WBG' in gene_series.index[0]):
        gene_ind = wbid
    else:
        gene_ind = name

    value = gene_rank_df.loc[gene_ind,'value']
    if value == 0:
        print(f'Value for gene {name} is 0')
        return 0
    else:
        rank = gene_rank_df.loc[gene_ind,'rank']
        percentile = (rank/gene_rank_df.shape[0])*100
        if print_gene:
            print(f'Gene - {name}\n')
        
        if print_res:
            print(f'Value:\t{value:.2f}\nRank: {rank} ({percentile:.2f}%)\n')

        return percentile

def print_gene_expression(gene:str, print_res: bool=False):
    rna_all = load_gene_expression_df()
    print_gene_ranks_in_df(rna_all, gene, print_res)


def load_gene_expression_df():
    '''
    Reads the gene expression values of all protein coding genes from 3 sources:
    1) 'ours': our mRNA in gonads
    2) 'Ahringer'
    3) 'expression_median': From Hila's table of YA
    4) 'expression_mean': From Hila's table of YA
    '''
    our_col = Table_mRNA().mean_exp
    our_df = pd.DataFrame({'ours (Gonads)':our_col})

    ahri_df = Ahringer().rna.loc[:,'Germline':'Muscle']
    ahri_df.columns = ['Ahringer (Germline)','Ahringer (Muscle)','Ahringer (Neurons)']

    exp_df = Gene_sets.get_expression_df()

    rna_all = pd.concat([our_df, ahri_df, exp_df], axis=1)
    rna_pc = protein_coding_only(rna_all)
    return rna_pc


def normalize_df_cols(df: pd.DataFrame, min_max: tuple=(-1,1)):
    '''
    '''
    x = df.values #returns a numpy array
    min_max_scaler = preprocessing.MinMaxScaler(feature_range=min_max)
    x_scaled = min_max_scaler.fit_transform(x)
    new_df=df.copy()
    new_df.loc[:,:] = x_scaled
    return new_df


def protein_coding_only(df: pd.DataFrame):
    '''
    Returns df with only protein coding elements. (Assumes wbid indices)
    '''
    pc_wbids = pd.read_csv('protein_coding_wbids.csv')

    ### add gfp
    pc_list = list(pc_wbids['genes'])
    pc_list.append('GFP')

    pc_series = pd.Series(pc_list)

    new_df = df[df.index.isin(pc_series)]
    return new_df

def get_top_somatic_rna_genes(prcnt: int):
    '''
    Uses Ahringer's data to find genes that are higly expressed in:
    - Muscle
    - Neurons

    Intersects lists and returns list with wbids.
    '''
    ar = Ahringer()
    top_neuro = get_list_of_column(ar.rna['Neurons'],prcnt)
    top_musc = get_list_of_column(ar.rna['Muscle'],prcnt)
    top_soma = intersect_lists(top_neuro, top_musc)
    return top_soma


def get_list_of_column(
        gene_series: pd.Series,
        prcnt: float = 0,
        bottom: bool = False,
        thresh=None,
        under_thresh: bool = False,
    ):
    """
    Parameters
    ----------
    - gene_series: pd.Series. Inds are wbis/gene_ids

    Return
    ----------
    - gene_id_list: list.
    """
    col_orig = gene_series.dropna()
    if prcnt:
        if not bottom: # get upper percnage
            quantile = col_orig.quantile(1 - (prcnt / 100))
            new_col = col_orig[col_orig > quantile]
        
        else: # get lower percentage
            quantile = col_orig.quantile((prcnt / 100))
            new_col = col_orig[col_orig <= quantile]
    else:
        new_col = col_orig

    if thresh is not None:
        if under_thresh:
            new_col = new_col[new_col < thresh]
        else:
            new_col = new_col[new_col > thresh]

    index_list = list(new_col.index)

    if ('GFP' in index_list[0]) or ('WBG' in index_list[0]):
        wbid_list = index_list
    else:
        wbid_list = list_to_wbids(index_list)

    return wbid_list


def list_to_wbids(gene_list: list):
    '''
    Translate list of genes to wbid list
    '''
    gid = Gene_IDs()
    wbid_list = [gid.to_wbid(gene) for gene in gene_list if gid.to_wbid(gene)]
    return wbid_list
        



