import pandas as pd 
import math

class Gene_IDs():

    def __init__(self):
        self.table = self.read_table()

    def read_table(self):
        big_table = pd.read_csv('extendedCelegansIDs_bigTable_39col_2019format_edit.csv', index_col='gene')
        id_table = big_table[['gene ID','Sequence ID','Other IDs']].copy()
        return id_table

    def to_name(self, wbid: str):
        '''
        '''
        name = self.table.loc[wbid,'gene ID']
        return name

        
if __name__=='__main__':
    gid = Gene_IDs()

