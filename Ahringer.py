import pandas as pd 
import gene_sets as gsets 
import regex as re
import utilities as ut 
import matplotlib.pyplot as plt

class Ahringer():  
    def __init__(self):
        self.rna = self._read_RNA()
        self.atac_all = self._read_ATAC()
        self.atac = self.atac_all[self.atac_all['geneID']!='.'] # only peaks related to genes
    
    def _read_RNA(self):
        rna = pd.read_csv('DATA/Ahringer/tissue-specific.RNA-seq.dataset.txt', delimiter='\t')
        rna.rename(columns={'WormBaseID':'Wbid'}, inplace=True)
        return rna
    
    def _read_ATAC(self):
        atac = pd.read_csv('DATA/Ahringer/tissue-specific.ATAC-seq.dataset.txt', delimiter='\t')
        cols = list(atac.columns)
        new_cols = []
        for col in cols:
            m = re.search('_TPM',col)
            if m:
                new_cols.append(col[:m.start()])
            else:
                new_cols.append(col)
        atac.columns = new_cols
        return atac
        
        


    def get_list(self, data_type: str, col_name: str, prcnt: float = 0, bottom: bool = False, thresh=None, under_thresh: bool=False):
        """
        Parameters
        ----------
        - col_name: str. name of columns
        - data_type: str. ['rna' / 'atac']

        Return
        ----------
        - gene_id_list: list.
        """
        if data_type.lower()=='rna':
            data = self.rna
        if data_type.lower()=='atac':
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
        
        gene_id_list = list(data.loc[index_list,'geneID'])
        gene_id_list = list(set(gene_id_list))

        return gene_id_list
    

if __name__=='__main__':
    ar = Ahringer()

    ### high rna and low atac: ###
    rna_above_1000 = ar.get_list('rna','Germline',thresh=1000) # len 116
    rna_top_5p = ar.get_list('rna','Germline',prcnt=5) # len 1012

    atac_below_5 = ar.get_list('atac','Germline',thresh=5, under_thresh=True) # len 733
    atac_bottom_5p = ar.get_list('atac','Germline',prcnt=5, bottom=True) # len 1568

    # no intersection --- high_ex_low_atac = ut.intersect_lists(rna_above_1000, atac_below_5)
    high_ex_low_atac = ut.intersect_lists(rna_top_5p, atac_below_5)

    ### High RNA in soma (Muscles, Neurons), low RNA in germline ###
    # only zeros: rna_germline_bottom_5p = ar.get_list('rna','Germline',prcnt=5,bottom=True)
    # only zeros: rna_germline_bottom_10p = ar.get_list('rna','Germline',prcnt=10,bottom=True)
    rna_germline_below_10 = ar.get_list('rna', 'Germline', thresh=10, under_thresh=True) # len 11,376


    rna_muscles_top_5p = ar.get_list('rna','Muscle',prcnt=5) # len 1012
    rna_neurons_top_5p = ar.get_list('rna','Neurons',prcnt=5) # len 1012
    top_muscle_neurons = ut.intersect_lists(rna_muscles_top_5p, rna_neurons_top_5p) # len 789

    only_somatic = ut.plot_ven([rna_germline_below_10, top_muscle_neurons],['germline low','somatic high'])
    only_somatic_names = ut.intersect_lists(rna_germline_below_10, top_muscle_neurons)




