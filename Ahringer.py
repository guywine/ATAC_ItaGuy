import pandas as pd
import gene_sets as gsets
import regex as re
from mRNA_gonads import Table_mRNA
from gene_id import Gene_IDs
from atac_signal import ATAC_signal

# import utilities as ut
import matplotlib.pyplot as plt


class Ahringer:
    '''
    This class reads Ahringer two types of data:
    1) RNA
    2) ATAC - Some genes appear multiple times.
    Both tables will be indexed by wbid.
    '''
    def __init__(self):
        self.rna = self._read_RNA()
        self.atac_all = self._read_ATAC()
        self.atac = self._filter_atac_all()

    @staticmethod  # makes an instance method NOT set first arg as 'self'
    def _read_RNA():
        rna = pd.read_csv(
            "DATA/Ahringer/tissue-specific.RNA-seq.dataset.txt", delimiter="\t"
        )
        rna.rename(columns={"WormBaseID": "Wbid"}, inplace=True)
        rna.index = rna["Wbid"]
        return rna

    @staticmethod
    def _read_ATAC():
        atac = pd.read_csv(
            "DATA/Ahringer/tissue-specific.ATAC-seq.dataset.txt", delimiter="\t"
        )
        cols = list(atac.columns)
        new_cols = []
        for col in cols:
            m = re.search("_TPM", col)
            if m:
                new_cols.append(col[: m.start()])
            else:
                new_cols.append(col)
        atac.columns = new_cols
        atac.index = atac["geneID"]
        return atac

    def _filter_atac_all(self):
        """
        The Ahringer table consists of many peaks which are not assigned to any gene.
        This function filters two things:
        a) Only peaks assigned to a gene
        b) only peaks that are marked as "coding promoter"
        """
        atac_new = self.atac_all[
            self.atac_all["geneID"] != "."
        ]  # only peaks related to genes
        atac_new = atac_new[atac_new["regulatory_class"] == "coding_promoter"]
        return atac_new

    def screen_genes_for_all_peaks(
        self, gene_list: list, thresh: float, below_thresh: bool = False
    ):
        """
        Some genes in the table have multiple peaks assigned to them as promoters.
        This function returns a list of gene names of which all of their atac peaks in the Germline match the desired threshold.
        """
        gid = Gene_IDs()
        name_list = [gid.to_name(gene) for gene in gene_list if gid.to_name(gene)]

        new_list = []
        for gene in name_list:
            try:
                atac_values = self.atac.loc[gene, "Germline"]
                if below_thresh:
                    if (atac_values <= thresh).all():
                        new_list.append(gene)
                else:
                    if (atac_values >= thresh).all():
                        new_list.append(gene)
            except KeyError:
                pass

        return new_list


if __name__ == "__main__":

    import utilities as ut

    ar = Ahringer()

    # m = Table_mRNA()

    # germline_ahri_under5 = ut.get_list_of_column(
    #     ar.rna["Germline"], thresh=5, under_thresh=True
    # )  # 9762
    # germline_ours_under5 = ut.get_list_of_column(
    #     m.mean_exp, thresh=5, under_thresh=True
    # )  # 13954
