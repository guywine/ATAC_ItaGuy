import pandas as pd
import gene_sets as gsets
import regex as re
from mRNA_gonads import Table_mRNA
from gene_id import Gene_IDs
from atac_signal import ATAC_signal

# import utilities as ut
import matplotlib.pyplot as plt


class Ahringer:
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
        atac_new = self.atac_all[
            self.atac_all["geneID"] != "."
        ]  # only peaks related to genes
        atac_new = atac_new[atac_new["regulatory_class"] == "coding_promoter"]
        return atac_new

    def get_list(
        self,
        data_type: str,
        col_name: str,
        prcnt: float = 0,
        bottom: bool = False,
        thresh=None,
        under_thresh: bool = False,
    ):
        """
        Parameters
        ----------
        - col_name: str. name of columns
        - data_type: str. ['rna' / 'atac']

        Return
        ----------
        - gene_id_list: list.
        """
        if data_type.lower() == "rna":
            data = self.rna
        if data_type.lower() == "atac":
            data = self.atac

        col_orig = data[col_name].dropna()
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
                new_col = col_orig[col_orig <= quantile]
        else:
            new_col = col_orig

        if thresh is not None:
            if under_thresh:
                new_col = new_col[new_col < thresh]
            else:
                new_col = new_col[new_col > thresh]

        index_list = list(new_col.index)

        gene_id_list = list(data.loc[index_list, "geneID"])
        gene_id_list = list(set(gene_id_list))

        return gene_id_list


if __name__ == "__main__":
    import utilities as ut

    ar = Ahringer()
    m = Table_mRNA()
    gid = Gene_IDs()
    ats = ATAC_signal()


    ### new
    # germline_ahri_top10 = ut.get_list_of_column(ar.rna['Germline'],prcnt=10) # 2023
    # germline_ours_top10 = ut.get_list_of_column(m.mean_exp, prcnt=10) # 2029
    # germline_top_10 = ut.intersect_lists(germline_ahri_top10, germline_ours_top10) # 1087

    # atac_ahri_low = ut.get_list_of_column(ar.atac['Germline'],thresh=10, under_thresh=True) # 4478

    # igfp_mean_scores = ats.scores1.mean(axis=1)
    # ioma1_mean_scores = ats.scores2.mean(axis=1)
    # atac_ours_gfp_low = ut.get_list_of_column(igfp_mean_scores, prcnt=70, bottom=True) # 14184
    # atac_ours_oma1_low = ut.get_list_of_column(ioma1_mean_scores, prcnt=70, bottom=True) # 14221

    # atac_ours_low = ut.intersect_lists(atac_ours_gfp_low, atac_ours_oma1_low) # 13883
    # atac_all_low = ut.intersect_lists(atac_ours_low, atac_ahri_low) # 2892

    # high_rna_low_atac = ut.intersect_lists(germline_top_10, atac_all_low)


