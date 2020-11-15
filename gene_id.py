import pandas as pd 
import math

class Gene_IDs():

    def __init__(self):
        self.table = self.read_table()
        # self.name_table = self.reverse_table()

    def read_table(self):
        big_table = pd.read_csv('extendedCelegansIDs_bigTable_39col_2019format_edit.csv', index_col='gene')

        id_table = big_table[['gene ID','Sequence ID','Other IDs']].copy()
        id_table.rename(columns={'gene ID':'name'}, inplace=True)
        id_table['name'] = id_table['name'].str.lower()
        id_table.index.names = ['wbid']
        id_table['wbid'] = id_table.index
        return id_table
    
    # def reverse_table(self):
    #     rev = self.table.dropna(subset=['name'])
    #     rev['wbid'] = rev.index
    #     rev.set_index('name',inplace=True)
    #     return rev

    # def to_name_old(self, gene: str):
    #     '''
    #     '''
    #     if 'WBGene' not in gene:
    #         if gene in set(self.table['name']):
    #             return gene
    #         else: 
    #             return f'gene "{gene}" not found'
        
    #     ### later add assertion this gene exists
    #     try:
    #         name = self.table.loc[gene,'name']
    #     except KeyError:
    #         name = f'gene "{gene}" not found' 
    #     return name
    
    def to_name(self, gene):
        '''
        Returns name, if no name, returns other ID, if no ID, returns 'nan'.
        Excepts name, wbid, and sequence ID.

        Returns false if gene not found.

        Parameters
        ----------
        - gene: str. wbid / name / seqID
        '''
        gene_wbid = self.to_wbid(gene)
        if not gene_wbid:
            return False
        
        if str(self.table.loc[gene_wbid,'name'])!='nan':
            return self.table.loc[gene_wbid,'name']
        elif str(self.table.loc[gene_wbid,'Sequence ID'])!='nan':
            return self.table.loc[gene_wbid,'Sequence ID']
        elif str(self.table.loc[gene_wbid,'Other IDs'])!='nan':
            return self.table.loc[gene_wbid,'Other IDs']
        else:
            return 'nan'

    
    def to_wbid(self, gene: str):
        '''
        Returns wbid of gene.
        Excepts name, wbid, and sequence ID.

        Returns false if gene not found.

        Parameters
        ----------
        - gene: str. wbid / name / seqID
        '''
        row = self.table[self.table.isin([gene]).any(axis=1)]
        if row.shape[0]==0:
            print(f'gene {gene} not found')
            return False
        else:
            wbid = row.iloc[0]['wbid']
            return wbid

    
        
if __name__=='__main__':
    gid = Gene_IDs()
    # gid.to_wbid('oma-1')

