import pandas as pd 
import pathlib

class Table_mRNA():

    def __init__(self):
        self.table = self._read_table()
        self.mRNA = self._create_mean_ranks()
    
    def _read_table(self):
        mRNA_table = '/Users/guyweintraub/Desktop/Google Drive/Rechavi Lab/Project ATAC/hrde-1/hrde1_mRNA/mRNA_rpm.csv'
        mRNA_path = pathlib.Path(mRNA_table)
        mRNA = pd.read_csv(mRNA_path, index_col='Wbid')
        mRNA.drop('Unnamed: 0', axis=1, inplace=True)
        return mRNA

    def _create_mean_ranks(self):
        means = self._add_mean_cols()
        means_ranked = self._add_ranks(means)
        return means_ranked 

    def _add_mean_cols(self):
        means = pd.DataFrame([])
        hrde1 = self.table.iloc[:,:3]
        means['hrde-1 mean'] = hrde1.mean(axis=1)

        sx = self.table.iloc[:,3:]
        means['sx mean'] = sx.mean(axis=1)
        return means
    
    def _add_ranks(self, means):
        means['sx rank'] = means['sx mean'].rank()
        means['hrde-1 rank'] = means['hrde-1 mean'].rank()
        return means

    def get_gene_rank(self, gene_name:str):
        exp_level = self.mRNA.loc[gene_name,'sx mean']
        if exp_level == 0:
            print(f'Expression for gene {gene_name} is 0')
        else:
            sx_rank = self.mRNA.loc[gene_name,'sx rank']
            percentile = (sx_rank/self.mRNA.shape[0])*100
            print(f'Gene - {gene_name}\nExp:\t{exp_level:.2f}\nRank: {sx_rank} ({percentile:.2f}%)')
        
            mutant_rank = self.mRNA.loc[gene_name,'hrde-1 rank']
            if abs(sx_rank-mutant_rank)>1000:
                perc_mutant = (mutant_rank/self.mRNA.shape[0])*100
                print(f'in the mutant, the rank is {mutant_rank} ({perc_mutant:.2f}%)')
        
        print('\n')




if __name__=='__main__':
    m = Table_mRNA()

    t = m.mRNA
    m.get_gene_rank('GFP')
    m.get_gene_rank('WBGene00000002')


