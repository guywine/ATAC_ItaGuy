'''
This will provide the user to choose which gene set to get.
- groups from table: by ATAC_bigtable (hrde-1 score, H3K9, mRNA level...)
- a set of the most / least 10% of ATAC_FC for a desired condition (choose percentage)
'''
import pandas as pd
import numpy as np

class Big_table():
    def __init__(self):
        self.big_table = self._read_big_table()
        self.all_wbids = list(self.big_table.index)

    def _read_big_table(self):
        '''
        Reads big_table of ATAC project, removes first header row, sets wbid to index.
        Replaces 0 with NaN in the desired columns.

        Return
        ---------
        - big_table: pd.DataFrame. index is wbid.
        '''
        big_table = pd.read_excel('ATAC_bigtable_1.xlsx', header=1)
        big_table.set_index('gene', inplace=True)

        # clean zeros from Ketting and Kennedy
        cols_to_drop_zero = ['isHrde1', 'ketting_score_all', 'ketting_score_20_23']
        big_table[cols_to_drop_zero] = big_table[cols_to_drop_zero].replace(0,np.nan)
        return big_table


    def get_list(self, col_name:str, prcnt:float=0, upper:bool = True):
        '''
        Parameters
        ----------
        - col_name: str. name of columns

        Return
        ----------
        - wbid_list: list.
        '''
        if col_name.lower() == 'all':
            return self.all_wbids

        col_orig = self.big_table[col_name].dropna()
        if prcnt:
            if upper:
                # get upper percnage
                quantile = col_orig.quantile(1-(prcnt/100))
                print(f'quantile upper of {col_name} = {quantile}')
                new_col = col_orig[col_orig>quantile]    
            else:
                # get lower percentage
                quantile = col_orig.quantile((prcnt/100))
                print(f'quantile lower of {col_name} = {quantile}')
                new_col = col_orig[col_orig<quantile]
        else:
            new_col = col_orig
        
        wbid_list = list(new_col.index)
        
        return wbid_list

    def get_multiple_lists(self, dic_list):
        '''
        dic_list = {'hrde-1':['isHrde1', 10], 
                'pol-2':['isPol2'], 
                'lowly':['R1-SX_S14', 10, False],
                'all genes':['ALL']}

        return
        ---------
        - dic_groups: dic.
        '''
        dic_groups = {}
        for lst_name in dic_list:
            dic_groups[lst_name] = self.get_list(*dic_list[lst_name])
   
        return dic_groups


if __name__=='__main__':
    at = Big_table()

    ## add option of all genes
    col_name = 'isHrde1'
    # ids_list = at.get_list(col_name, 10, False)

    dic_list = {'hrde-1':['isHrde1', 10], 
        'pol-2':['isPol2'], 
        'highly':['R1-SX_S14', 10],
        'all genes':['ALL']}

    dic_groups = at.get_multiple_lists(dic_list)

    print('done')