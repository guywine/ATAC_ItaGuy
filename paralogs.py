import pandas as pd 
from gene_id import Gene_IDs

class Paralogs:
    def __init__(self):
        self.table = pd.read_csv('caenorhabditis_elegans.PRJNA13758.WBPS14.paralogs.tsv', delimiter='\t')
        self.gid = Gene_IDs()
        self.cutoff_default = 98


    def get_cutoff(self, cutoff: int = 0):
        '''
        cutoff preformed on 'target_identity'
        '''
        if cutoff==0:
            cutoff = self.cutoff_default
        
        table_cutoff = self.table[self.table['target_identity']>cutoff]
        return table_cutoff


    def get_cutoff_wbids(self, cutoff: int = 0):
        if cutoff==0:
            cutoff = self.cutoff_default

        table_cutoff = self.get_cutoff(cutoff)
        wbids_list = list(table_cutoff['gene_id'])
        wbids_list.extend(list(table_cutoff['paralog_gene_id']))
        wbids_list = list(set(wbids_list))
        return wbids_list
    

    def get_paralogs_of_gene(self, gene: str, cutoff: int = 0):
        '''
        Returns False if no paralogs for this cutoff
        '''
        if cutoff==0:
            cutoff = self.cutoff_default

        table_cutoff = self.get_cutoff(cutoff)
        wbid = self.gid.to_wbid(gene)
        gene_rows_mask = table_cutoff['gene_id']==wbid
        if gene_rows_mask.sum()==0:
            return False
        paralogs_list = list(table_cutoff['paralog_gene_id'][gene_rows_mask])
        return paralogs_list

    def is_paralog(self, gene:str, cutoff: int = 0):
        if cutoff==0:
            cutoff = self.cutoff_default
        
        table_cutoff = self.get_cutoff(cutoff)
        wbid = self.gid.to_wbid(gene)
        wbids_list = self.get_cutoff_wbids(cutoff)
        if wbid in wbids_list:
            return True
        else:
            return False


    def different_trait(self, wbids1:list, wbids2:list, cutoff: int = 0):
        '''
        Takes 2 lists and checks if theres a couple of paralogs in which 
        one is in first list and the other on the second list.
        For instance: one appears in "list_highly" and other in "list_lowly"
        '''
        if cutoff==0:
            cutoff = self.cutoff_default

        table_cutoff = self._get_cutoff(cutoff)
        couples = []
        for row_i in range(table_cutoff.shape[0]):
            gene1 = table_cutoff.iloc[row_i,:]['gene_id']
            gene2 = table_cutoff.iloc[row_i,:]['genparalog_gene_ide_id']
            if gene1 in wbids1:
                if gene2 in wbids2:
                    couples.append((gene1, gene2))
            elif gene1 in wbids2:
                if gene2 in wbids1:
                    couples.append((gene1, gene2))
        
        if not couples:
            print('No paralogs matching these lists found')
            return False
        
        return couples
        

if __name__=='__main__':
    para = Paralogs()
    para_ppw1 = para.get_paralogs_of_gene('ppw-1')

