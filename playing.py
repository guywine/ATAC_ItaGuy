import pandas as pd 
import matplotlib.pyplot as plt 

import utilities as ut 
from mRNA_gonads import Table_mRNA
from Ahringer import Ahringer
from gene_sets import Gene_sets
from atac_signal import ATAC_signal
from gene_id import Gene_IDs


if __name__=='__main__':
    ### create instances of all classes
    ar = Ahringer()
    gid = Gene_IDs()
    
    if "m" not in locals():
        m = Table_mRNA()

    if "ats" not in locals():
        ats = ATAC_signal()

    if "gs" not in locals():
        gs = Gene_sets()


    ## df of hrde-1 up regulated genes
    if True:
        hrde1_kennedy = gs.get_list('hrde-1-Kennedy')
        hrde_FC_sig = gs.get_list('mRNA_isSig')
        hrde_up = gs.get_list('mRNA_log2_FC', thresh=0)
        hrde_up_sig = ut.intersect_lists(hrde_FC_sig, hrde_up)
        hrde_regulated = ut.intersect_lists(hrde_up_sig, hrde1_kennedy)
        df_mRNA_and_hrde1 = gs.big_table.loc[hrde_regulated,['mRNA_log2_FC','hrde-1-Kennedy']]

    ## low rna:
    if True: 
        germline_ahri_low = ut.get_list_of_column(
            ar.rna["Germline"], thresh=10, under_thresh=True
        )  # 11,380
        germline_ours_low = ut.get_list_of_column(
            m.mean_exp, thresh=5, under_thresh=True
        )  # 13,954
        rna_low = ut.intersect_lists(
            germline_ahri_low, germline_ours_low
        )  # 10,697

    ## low atac:
    if True:
        igfp_mean_scores = ats.scores1.mean(axis=1)
        ioma1_mean_scores = ats.scores2.mean(axis=1)

        atac_ours_gfp_low = ut.get_list_of_column(
            igfp_mean_scores, prcnt=40, bottom=True
        )  #
        atac_ours_oma1_low = ut.get_list_of_column(
            ioma1_mean_scores, prcnt=40, bottom=True
        )  #
        atac_ours_low = ut.intersect_lists(atac_ours_gfp_low, atac_ours_oma1_low)  # 7577 

        atac_ahri_low = ut.get_list_of_column(
            ar.atac["Germline"], prcnt=40, bottom=True
        )  # 5585

        atac_low = ut.intersect_lists(atac_ahri_low, atac_ours_low) # 1537

    #### strong hrde-1, low rna, low atac
    if True:
        low_rna_low_atac = ut.intersect_lists(rna_low, atac_ours_low) # 1212
        low_closed_hrde1 = ut.intersect_lists(low_rna_low_atac, hrde_regulated)

    if True:
        low_rna_low_atacscreen_chromosome

    








