import pandas as pd
import gene_sets as gsets
import regex as re
from mRNA_gonads import Table_mRNA
from gene_id import Gene_IDs

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

    rna_all = ut.load_gene_expression_df()

    germline_ahri_top5 = ut.get_list_of_column(ar.rna['Germline'],prcnt=5) # len 1012
    germline_ours_top5 = ut.get_list_of_column(m.mean_exp, prcnt=5) # len 1015

    germline_top_5 = ut.intersect_lists(germline_ahri_top5, germline_ours_top5) # len 450

    germ_top_5_df = ar.rna.loc[germline_top_5,:]
    germ_top5_1st = germ_top_5_df[germ_top_5_df['1st.max.tissue']=='Germline'] # 275

    germ_top5_over400_1st = germ_top5_1st[germ_top5_1st['Germline']>400] # 98

    soma_under_150_mask = (germ_top5_over400_1st.loc[:,'Neurons':'Intest.']<150).all(axis=1)

    germ_top5_over400_1st_soma_under_150 = germ_top5_over400_1st[soma_under_150_mask]

    strong_germline_rna = list(germ_top5_over400_1st_soma_under_150.index)

    atac_top_10 = ut.get_list_of_column(ar.atac['Germline'],prcnt=10)

    germline_atac = ut.intersect_lists(strong_germline_rna, atac_top_10)

