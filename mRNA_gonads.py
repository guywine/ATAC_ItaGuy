import pandas as pd 
import pathlib
from gene_id import Gene_IDs

class Table_mRNA():

    def __init__(self):
        self.table = self._read_table()
        self.mRNA = self._create_mean_ranks()
        self.mean_exp = self.mRNA.loc[:,'sx mean']
        self.fc = self._create_FC()
    
    def _read_table(self):
        mRNA_table = '/Users/guyweintraub/Desktop/Google Drive/Rechavi Lab/Project ATAC/hrde-1/hrde1_mRNA/mRNA_rpm.csv'
        mRNA_path = pathlib.Path(mRNA_table)
        mRNA = pd.read_csv(mRNA_path, index_col='Wbid')
        mRNA.drop('Unnamed: 0', axis=1, inplace=True)

        mRNA = self._protein_coding_only(mRNA)

        return mRNA
    
    def _protein_coding_only(self, table):
        '''
        '''
        pc_wbids = pd.read_csv('protein_coding_wbids.csv')

        ### add gfp
        pc_list = list(pc_wbids['genes'])
        pc_list.append('GFP')

        pc_series = pd.Series(pc_list)

        table = table[table.index.isin(pc_series)]
        return table

    def _create_FC(self):
        add_to_avoid_zero = 1

        FC_table = pd.DataFrame([])
        for rep_i in range(1,4):
            hrde1_col = self.table[f'hrde1_r{rep_i}_rpm']+add_to_avoid_zero
            sx_col = self.table[f'sx_r{rep_i}_rpm']+add_to_avoid_zero
            FC_table[f'FC rep {rep_i}'] = hrde1_col / sx_col
        
        FC_table['mean fc']=FC_table.mean(axis=1)
        
        return FC_table

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

    def get_gene_rank(self, gene:str):
        if gene.lower()=='gfp':
            wbid='GFP'
            name = wbid
        else:
            gid = Gene_IDs()
            wbid = gid.to_wbid(gene)
            name = gid.to_name(gene)
            if str(name)=='nan':
                name = wbid
        exp_level = self.mRNA.loc[wbid,'sx mean']
        if exp_level == 0:
            print(f'Expression for gene {name} is 0')
        else:
            sx_rank = self.mRNA.loc[wbid,'sx rank']
            percentile = (sx_rank/self.mRNA.shape[0])*100
            print(f'Gene - {name}\nExp:\t{exp_level:.2f}\nRank: {sx_rank} ({percentile:.2f}%)')
        
            mutant_rank = self.mRNA.loc[wbid,'hrde-1 rank']
            if abs(sx_rank-mutant_rank)>1000:
                perc_mutant = (mutant_rank/self.mRNA.shape[0])*100
                print(f'in the mutant, the rank is {mutant_rank} ({perc_mutant:.2f}%)')
        
        print('\n')




if __name__=='__main__':
    m = Table_mRNA()

    t = m.mRNA
    # m.get_gene_rank('GFP')
    # m.get_gene_rank('WBGene00000002')
    # m.get_gene_rank('oma-1')
    # m.get_gene_rank('oma-2')
    # m.get_gene_rank('ppw-1')
    # m.get_gene_rank('sago-2')
    # m.get_gene_rank('pptr-1')
    # m.get_gene_rank('pptr-2')
    # m.get_gene_rank('lin-15a')


