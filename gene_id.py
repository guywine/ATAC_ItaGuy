import pandas as pd 
import math

class Gene_IDs():

    def __init__(self):
        self.table = self.read_gene_IDs_table()

    @staticmethod
    def read_gene_IDs_table():
        id_table = pd.read_csv('gene_IDs_table.txt',sep='\t')
        id_table.drop('Your Input', axis=1, inplace=True)
        id_table.columns = ['wbid','gene name', 'sequence ID','other name']
        id_table.set_index('wbid', drop=False, inplace=True)
        return id_table

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
        
        if str(self.table.loc[gene_wbid,'gene name'])!='nan':
            return self.table.loc[gene_wbid,'gene name']
        else:
            return 'nan'

    
        
if __name__=='__main__':
    gid = Gene_IDs()
    # gid.to_wbid('oma-1')

