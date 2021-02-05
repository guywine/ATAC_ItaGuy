import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from sklearn import preprocessing

from gene_id import Gene_IDs
from gene_sets import Gene_sets
from mRNA_gonads import Table_mRNA
from Ahringer import Ahringer
# from atac_signal import ATAC_signal


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
    """"""
    gid = Gene_IDs()
    gene_names = [
        gid.to_name(wbid) for wbid in wbid_list if str(gid.to_name(wbid)) != "nan"
    ]
    return gene_names


def gene_is_in_list(wbid_list: list, gene_name: str, print_flag: bool = False):
    """"""
    gene_names = wbid_list_to_names(wbid_list)
    if gene_name in gene_names:
        if print_flag:
            gid = Gene_IDs()
            wbid = gid.to_wbid(gene_name)
            print(f"gene {wbid} found")
        return True
    else:
        if print_flag:
            print("gene {gene_name} not found.")
        return False


def print_gene_ranks_in_df(gene_df: pd.DataFrame, gene: str, print_res: bool = False):
    """"""
    gene_name = Gene_IDs().to_name(gene)
    print(f"gene - {gene_name}")
    for column in gene_df.columns:
        if print_res:
            print(f"{column}")
        perc = get_gene_rank(gene_df[column], gene, print_res, print_gene=False)
        if not print_res:
            print(f"{column}:\t{perc:.2f}%")
    print("\n")


def calc_mean_variance_of_dfs(df_list: list, variance_type: str = "none"):
    """
    Gets a list of dfs with identicle structure, and calculates for each cell the mean and variance across all dfs in list.
    Returns a df containing means, and a seperate df containing variance.

    Parameters
    ----------
    - df_list: list of dfs to mean (each df a sample of rep)
    - variance_type: str. ['std' / 'sem' / 'none']
    """
    ### create single df average across replicates:
    df_mean_all_reps = pd.concat(df_list).groupby(level=0).mean()
    if variance_type.lower() == "std":
        df_variance = pd.concat(df_list).groupby(level=0).std()
    elif variance_type.lower() == "sem":
        df_variance = pd.concat(df_list).groupby(level=0).sem()
    elif variance_type.lower() == "none":
        df_variance = 0

    return df_mean_all_reps, df_variance


def print_gene_atac(atac_sig, gene: str):
    """
    Gets an ATAC_signal object and a gene.
    Prints its ranks in ours (GFP + oma-1) and ahringer.
    """
    gfp_score = atac_sig.scores1.mean(axis=1)
    oma1_score = atac_sig.scores2.mean(axis=1)
    ahri_germline = Ahringer().atac["Germline"]

    gene_name = Gene_IDs().to_name(gene)
    print(f"ATAC for gene - {gene_name}")

    print("Ours (anti-GFP)")
    get_gene_rank(gfp_score, gene_name, print_gene=False)

    print("Ours (anti-OMA-1)")
    get_gene_rank(oma1_score, gene_name, print_gene=False)

    print("Ahringer Germline")
    get_gene_rank(ahri_germline, gene_name, print_gene=False)


def get_gene_rank(
    gene_series: pd.Series, gene: str, print_res: bool = True, print_gene: bool = True
):
    """
    Series where index is gene_name.

    Return
    --------
    - percentile: float.
    """
    ### later remove:
    if gene.lower() == "gfp":
        wbid = "GFP"
        name = wbid
    else:
        gid = Gene_IDs()
        wbid = gid.to_wbid(gene)
        name = gid.to_name(gene)
        if str(name) == "nan":
            name = wbid
    #####

    gene_rank_df = pd.DataFrame({"value": gene_series, "rank": gene_series.rank()})

    if ("GFP" in gene_series.index[0]) or ("WBG" in gene_series.index[0]):
        gene_ind = wbid
    else:
        gene_ind = name

    value = gene_rank_df.loc[gene_ind, "value"]
    # if value == 0:
    #     print(f"Value for gene {name} is 0")
    #     return 0
    # else:
    rank = gene_rank_df.loc[gene_ind, "rank"]
    percentile = (rank / gene_rank_df.shape[0]) * 100
    if print_gene:
        print(f"Gene - {name}\n")

    if print_res:
        print(f"Value:\t{value:.2f}\nRank: {rank} ({percentile:.2f}%)\n")

    return percentile


def find_rank_by_value(gene_series: pd.Series, val: float, print_flag: bool = True):
    """
    Series where index is gene_name.

    Return
    --------
    - percentile: float.
    """
    smaller = (gene_series <= val).sum()
    percentile = 100 * (smaller / gene_series.size)
    if print_flag:
        print(f"value: {val}\tpercentile : {percentile:.2f}%")
    return percentile


def print_gene_expression(gene: str, print_res: bool = False):
    rna_all = load_gene_expression_df()
    print_gene_ranks_in_df(rna_all, gene, print_res)


def load_gene_expression_df():
    """
    Reads the gene expression values of all protein coding genes from 3 sources:
    1) 'ours': our mRNA in gonads
    2) 'Ahringer'
    3) 'expression_median': From Hila's table of YA
    4) 'expression_mean': From Hila's table of YA
    """
    our_col = Table_mRNA().mean_exp
    our_df = pd.DataFrame({"ours (Gonads)": our_col})

    ahri_df = Ahringer().rna.loc[:, "Germline":"Muscle"]
    ahri_df.columns = ["Ahringer (Germline)", "Ahringer (Muscle)", "Ahringer (Neurons)"]

    exp_df = Gene_sets.get_expression_df()

    rna_all = pd.concat([our_df, ahri_df, exp_df], axis=1)
    rna_pc = protein_coding_only(rna_all)
    return rna_pc


def normalize_df_cols(df: pd.DataFrame, min_max: tuple = (-1, 1)):
    """"""
    x = df.values  # returns a numpy array
    min_max_scaler = preprocessing.MinMaxScaler(feature_range=min_max)
    x_scaled = min_max_scaler.fit_transform(x)
    new_df = df.copy()
    new_df.loc[:, :] = x_scaled
    return new_df


def protein_coding_only(df: pd.DataFrame):
    """
    Returns df with only protein coding elements. (Assumes wbid indices)
    """
    pc_wbids = pd.read_csv("tables/protein_coding_wbids.csv")

    ### add gfp
    pc_list = list(pc_wbids["genes"])
    pc_list.append("GFP")

    pc_series = pd.Series(pc_list)

    new_df = df[df.index.isin(pc_series)]
    return new_df


def get_top_somatic_rna_genes(prcnt: int):
    """
    Uses Ahringer's data to find genes that are higly expressed in:
    - Muscle
    - Neurons

    Intersects lists and returns list with wbids.
    """
    ar = Ahringer()
    top_neuro = get_list_of_column(ar.rna["Neurons"], prcnt)
    top_musc = get_list_of_column(ar.rna["Muscle"], prcnt)
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
        if not bottom:  # get upper percnage
            quantile = col_orig.quantile(1 - (prcnt / 100))
            new_col = col_orig[col_orig > quantile]

        else:  # get lower percentage
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

    if ("GFP" in index_list[0]) or ("WBG" in index_list[0]):
        wbid_list = index_list
    else:
        wbid_list = list_to_wbids(index_list)

    return wbid_list


def list_to_wbids(gene_list: list):
    """
    Translate list of genes to wbid list
    """
    gid = Gene_IDs()
    wbid_list = [gid.to_wbid(gene) for gene in gene_list if gid.to_wbid(gene)]
    return wbid_list


def list_to_name(gene_list: list):
    """"""
    gid = Gene_IDs()
    name_list = [gid.to_name(gene) for gene in gene_list if gid.to_name(gene)]
    return name_list


def screen_chromosome(gene_list: list, chr_str: str, to_names: bool = False):

    chrom_df = pd.read_csv("tables/gene_locs.csv", index_col="Wbid")

    wbid_list = list_to_wbids(gene_list)

    chr_list = [gene for gene in wbid_list if chrom_df.loc[gene, "chr"] == chr_str]

    if to_names:
        chr_list = list_to_name(chr_list)

    return chr_list


def screen_gene_location(gene_list: list, to_names: bool = False):
    '''
    Screens only genes that are in the "middle" of the chromosome.
    Borders defined by paper publishing mutation rate change.
    '''
    chrom_df = pd.read_csv("tables/gene_locs.csv", index_col="Wbid")

    wbid_list = list_to_wbids(gene_list)

    arms_dict = {
        "I": (3858000, 11040000),
        "II": (4879000, 12020000),
        "III": (3722000, 10340000),
        "IV": (3896000, 12970000),
        "V": (5897000, 16550000),
        "X": (6137000, 12480000),
    }

    new_list = []
    for wbid in wbid_list:
        chrom = chrom_df.loc[wbid, "chr"]
        range_tuple = arms_dict[chrom]
        if gene_within_range(chrom_df, wbid, range_tuple):
            new_list.append(wbid)

    if to_names:
        new_list = list_to_name(new_list)

    return new_list

def print_gene_loc(gene: str):
    chrom_df = pd.read_csv("tables/gene_locs.csv", index_col="Wbid")
    
    gid = Gene_IDs()
    wbid = gid.to_wbid(gene)
    name = gid.to_name(gene)

    chrom = chrom_df.loc[wbid,'chr']
    start = chrom_df.loc[wbid,'start']
    stop = chrom_df.loc[wbid,'stop']
    strand = chrom_df.loc[wbid,'strand']

    print(f'gene {name}\nchrom:\t{chrom} , {strand} , {start}-{stop}')


def get_gene_size(gene: str):
    chrom_df = pd.read_csv("tables/gene_locs.csv", index_col="Wbid")
    
    gid = Gene_IDs()
    wbid = gid.to_wbid(gene)
    name = gid.to_name(gene)

    start = chrom_df.loc[wbid,'start']
    stop = chrom_df.loc[wbid,'stop']
    
    return abs(stop-start)

def screen_by_size(gene_list: str, thresh, max: bool=True):
    new_list = []
    for gene in gene_list:
        size = get_gene_size(gene)
        if max:
            if size<=thresh:
                new_list.append(gene)
        else:
            if size>=thresh:
                new_list.append(gene)
    return new_list



def gene_within_range(chrom_df: pd.DataFrame, wbid: str, range_tuple: tuple):
    start = chrom_df.loc[wbid, "start"]
    stop = chrom_df.loc[wbid, "stop"]
    if start >= range_tuple[0] and stop <= range_tuple[1]:
        return True
    else:
        return False

def screen_ahringer_1_prom(gene_name_list: list):
    ar = Ahringer()
    # gene_counts = ar.atac.index.value_counts()
    # genes_1_prom = gene_counts[gene_counts==1].index.tolist()
    # genes_1_prom_sep = list_seperate_commas(genes_1_prom)

    gene_inds = ar.atac.index.tolist()
    gene_inds_sep = list_seperate_commas(gene_inds)
    
    all_genes_with_prom = list(set(gene_inds_sep))
    genes_1_prom = screen_only_singles(gene_inds_sep)
    genes_more_than_1_prom = ut.intersect_lists(all_genes_with_prom, genes_1_prom, 'only first')


def list_seperate_commas(l:list):
    new_list = []
    for item in l:
        seperated = item.split(',')
        new_list.extend(seperated)
    
    return new_list


def screen_only_singles(l:list):
    '''
    Returns list of only items appearing 1 time in list.
    
    '''
    singles_list = []
    item_counter = collections.Counter(l)
    for item in item_counter:
        if item_counter[item]==1:
            singles_list.append(item)
    
    return singles_list

def print_gene_sRNAs(gene: str):
    gid = Gene_IDs()
    wbid = gid.to_wbid(gene)
    sRNAs_df_full = pd.read_csv('tables/sRNAs_rpm_df.csv', index_col='wbid')
    try:
        print(f'gene {gid.to_name(gene)} sRNAs:')
        print(sRNAs_df_full.loc[wbid,:])
    except KeyError:
        print(f'Gene {gene} is under the threshold')

def get_sRNAs_list(mean_above_thresh: bool=False):
    sRNAs_df = pd.read_csv('tables/sRNA_over_5_rpm.csv', index_col='wbid')
    if mean_above_thresh:
        sRNAs_df = sRNAs_df[sRNAs_df['mean']>=5]
    return sRNAs_df.index.tolist()


def screen_genes_with_name(gene_list: list):
    gid = Gene_IDs()
    names_only = [gid.to_name(gene) for gene in gene_list if '-' in gid.to_name(gene)]
    return names_only


def intersect_all(*lists_args):
    last_set = set(lists_args[0])
    for lst_i in range(1, len(lists_args)):
         st_new = set(lists_args[lst_i])
         last_set = last_set.intersection(st_new)
         print(last_set)
    
    return list(last_set)

# def thresh_df_by_col(df: pd.DataFrame, thresh: int=5):
#     '''
#     '''
#     rows_to_drop_i = []
#     for row_i in range(df.shape[0]):
#         if (df.iloc[row_i,:]<thresh).all():
#             rows_to_drop_i.append(row_i)
    
#     cutoff_df = df.drop(df.index[rows_to_drop_i])
#     return cutoff_df

