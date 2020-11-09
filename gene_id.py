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
        id_table.index.names = ['wbid']
        return id_table
    
    def reverse_table(self):
        rev = self.table.dropna(subset=['name'])
        rev['wbid'] = rev.index
        rev.set_index('name',inplace=True)
        return rev

    def to_name(self, wbid: str):
        '''
        '''
        name = self.table.loc[wbid,'name']
        return name
    
    def to_wbid(self, name: str):
        '''
        '''
        gene_name = name.strip().lower()
        wbid = self.name_table.loc[gene_name, 'wbid']
        return wbid
        
if __name__=='__main__':
    gid = Gene_IDs()

