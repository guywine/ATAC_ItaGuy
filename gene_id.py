import pandas as pd 
import math

class Gene_IDs():

    def __init__(self):
        self.table = self.read_table()
        self.name_table = self.reverse_table()

    def read_table(self):
        big_table = pd.read_csv('extendedCelegansIDs_bigTable_39col_2019format_edit.csv', index_col='gene')

        id_table = big_table[['gene ID','Sequence ID','Other IDs']].copy()
        id_table.rename(columns={'gene ID':'name'}, inplace=True)
        id_table['name'] = id_table['name'].str.lower()
        id_table.index.names = ['wbid']
        return id_table
    
    def reverse_table(self):
        rev = self.table.dropna(subset=['name'])
        rev['wbid'] = rev.index
        rev.set_index('name',inplace=True)
        return rev

    def to_name(self, gene: str):
        '''
        '''
        if 'WBGene' not in gene:
            if gene in self.table['name']:
                return gene
            else: 
                return f'gene "{gene}" not found'
        
        ### later add assertion this gene exists
        try:
            name = self.table.loc[gene,'name']
        except KeyError:
            name = f'gene "{gene}" not found' 
        return name
    
    def to_wbid(self, gene: str):
        '''
        '''
        if 'WBGene' in gene:
            if gene in self.table['wbid']:
                return gene
            else: 
                return f'gene "{gene}" not found'            
        
        gene = gene.strip().lower()

        try:
            wbid = self.name_table.loc[gene, 'wbid']
        except KeyError:
            wbid = f'gene "{gene}" not found' 
        return wbid
    
        
if __name__=='__main__':
    gid = Gene_IDs()
    gid.to_wbid('oma-1')

