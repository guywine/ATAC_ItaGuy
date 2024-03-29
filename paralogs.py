import pandas as pd
from gene_id import Gene_IDs
import utilities as ut 
from gene_id import Gene_IDs
import matplotlib.pyplot as plt

class Paralogs:
    def __init__(self):
        self.table = pd.read_csv(
            "tables/caenorhabditis_elegans.PRJNA13758.WBPS14.paralogs.tsv", delimiter="\t"
        )
        self.table_ids = self.table.set_index('gene_id')
        self.gid = Gene_IDs()
        self.cutoff_default = 98

    def get_cutoff(self, cutoff: int = 0):
        """
        cutoff preformed on 'target_identity'
        """
        if cutoff == 0:
            cutoff = self.cutoff_default

        table_cutoff = self.table[self.table["target_identity"] > cutoff]
        return table_cutoff

    def get_cutoff_wbids(self, cutoff: int = 0):
        if cutoff == 0:
            cutoff = self.cutoff_default

        table_cutoff = self.get_cutoff(cutoff)
        wbids_list = list(table_cutoff["gene_id"])
        wbids_list.extend(list(table_cutoff["paralog_gene_id"]))
        wbids_list = list(set(wbids_list))
        return wbids_list

    def get_paralogs_of_gene(self, gene: str, cutoff: int = 0):
        """
        Returns False if no paralogs for this cutoff
        """
        if cutoff == 0:
            cutoff = self.cutoff_default

        table_cutoff = self.get_cutoff(cutoff)
        wbid = self.gid.to_wbid(gene)
        gene_rows_mask = table_cutoff["gene_id"] == wbid
        if gene_rows_mask.sum() == 0:
            return False
        paralogs_list = list(table_cutoff["paralog_gene_id"][gene_rows_mask])
        return paralogs_list


    def are_paralogs(self, gene_lst: list, cutoff: int = 0):
        if cutoff == 0:
            cutoff = self.cutoff_default

        paralogs = []
        for gene in gene_lst:
            if self.is_paralog(gene, cutoff):
                paralogs.append(gene)

        return paralogs


    def is_paralog(self, gene: str, cutoff: int = 0):
        if cutoff == 0:
            cutoff = self.cutoff_default

        table_cutoff = self.get_cutoff(cutoff)
        wbid = self.gid.to_wbid(gene)
        wbids_list = self.get_cutoff_wbids(cutoff)
        if wbid in wbids_list:
            return True
        else:
            return False

    def different_trait(self, wbids1: list, wbids2: list, to_names: bool = False, cutoff: int = 0):
        """
        Takes 2 lists and checks if theres a couple of paralogs in which
        one is in first list and the other on the second list.
        For instance: one appears in "list_highly" and other in "list_lowly"
        """
        if cutoff == 0:
            cutoff = self.cutoff_default

        table_cutoff = self.get_cutoff(cutoff)
        couples = []
        for row_i in range(table_cutoff.shape[0]):
            gene1 = table_cutoff.iloc[row_i, :]["gene_id"]
            gene2 = table_cutoff.iloc[row_i, :]["paralog_gene_id"]
            if gene1 in wbids1:
                if gene2 in wbids2:
                    couples.append((gene1, gene2))
            elif gene1 in wbids2:
                if gene2 in wbids1:
                    couples.append((gene1, gene2))

        if not couples:
            print("No paralogs matching these lists found")
            return False
        
        if to_names:
            couples = self.couples_to_name(couples)

        return couples
    
    def couples_to_name(self, couples: list):
        couples_names = []
        for couple in couples:
            name1 = self.gid.to_name(couple[0])
            name2 = self.gid.to_name(couple[1])
            couples_names.append((name1, name2))
        return couples_names
    
    def get_best_para(self, gene: str, print_flag: bool = False):
        wbid = self.gid.to_wbid(gene)
        name = self.gid.to_name(gene)
        if wbid not in self.table_ids.index:
            return False, False
        paras_df = self.table_ids.loc[wbid,:]
        row_i = paras_df['target_identity'].argmax()
        best_para_score = paras_df.iat[row_i, 1]
        best_para_wbid = paras_df.iat[row_i, 0]

        best_para_name = self.gid.to_name(best_para_wbid)
        if print_flag:
            print(f'Best paralog for {name}:\t {best_para_name} ({best_para_score})\n')

        return best_para_score, best_para_wbid




if __name__ == "__main__":
    gid = Gene_IDs()

    para = Paralogs()

    para.table_ids['target_identity'].hist()

    # ut.find_rank_by_value(para.table_ids['target_identity'], 40)


    panel_genes = ['WBGene00003864','WBGene00001595',
                    'WBGene00000216','WBGene00003093','WBGene00000751',
                    'WBGene00000468','WBGene00000714','WBGene00000905',
                    'WBGene00009602', 'WBGene00008010', 'WBGene00016639']

    for gene in panel_genes:
        _,_ = para.get_best_para(gene, print_flag=True)

    '''

for gene in panel_genes:...
for gene WBGene00003864:
score: 59.0331, gene: oma-2
for gene WBGene00001595:
score: 33.4646, gene: K07H8.9
for gene WBGene00000216:
score: 38.7387, gene: asp-4
for gene WBGene00003093:
score: 84.5794, gene: lys-6
for gene WBGene00000751:
score: 93.2862, gene: col-179
for gene WBGene00000468:
score: 38.342, gene: snai-1
for gene WBGene00000714:
score: 71.0769, gene: col-142
for gene WBGene00000905:
score: 26.3158, gene: cyp-33E1
for gene WBGene00009602:
score: 46.6321, gene: Y61A9LA.12
for gene WBGene00008010:
score: 60.6925, gene: F15D4.5
    '''

    # germline_high_10_ours = ut.get_list_of_column(rna_all['ours (Gonads)'], prcnt=10)
    # germline_high_10_ahri = ut.get_list_of_column(rna_all['Ahringer (Germline)'], prcnt=10)
    # germline_high_10_both = ut.intersect_lists(germline_high_10_ours, germline_high_10_ahri) # 1085

    # para_high_10 = para.are_paralogs(germline_high_10_both)

    # germline_below_10_ours = ut.get_list_of_column(rna_all['ours (Gonads)'], thresh=10, under_thresh=True)
    # germline_below_10_ahri = ut.get_list_of_column(rna_all['Ahringer (Germline)'], thresh=10, under_thresh=True)
    # germline_below_10_both = ut.intersect_lists(germline_below_10_ours, germline_below_10_ahri) # 11129

    # couples_high_low = para.different_trait(germline_high_10_both, germline_below_10_both, to_names=True) # 67
    # ## all are his-x genes!

    # germline_below_80_ours = ut.get_list_of_column(rna_all['ours (Gonads)'], thresh=80, under_thresh=True) 
    # germline_below_80_ahri = ut.get_list_of_column(rna_all['Ahringer (Germline)'], thresh=80, under_thresh=True)
    # germline_below_80_both = ut.intersect_lists(germline_below_80_ours, germline_below_80_ahri) # ???

    # couples_high_low80 = para.different_trait(germline_high_10_both, germline_below_80_both, to_names=True) # 77
    # new_couples = ut.intersect_lists(couples_high_low, couples_high_low80, 'only second')

    # couples_high_low80_cutoff_90 = para.different_trait(germline_high_10_both, germline_below_80_both, to_names=True, cutoff=90) # 77

    # germline_bottom_70_ours = ut.get_list_of_column(rna_all['ours (Gonads)'], prcnt=70, bottom=True) # 14199
    # germline_bottom_70_ahri = ut.get_list_of_column(rna_all['Ahringer (Germline)'], prcnt=70, bottom=True) # 14102
    # germline_bottom_70_both = ut.intersect_lists(germline_bottom_70_ours, germline_bottom_70_ahri) # 12899

    # couples_bottom_70_cutoff_90 = para.different_trait(germline_high_10_both, germline_bottom_70_both, to_names=True, cutoff=90) # 70 couples




    

