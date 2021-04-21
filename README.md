# ATAC_ItaGuy
## Hackers of the ATAC project

**To Do**
- Add GFP on the read_tables.py
- Add to gene_sets the ability to get genes by atac signal.
- Add function "ATAC_FC" for a single gene
- Add function "ATAC_FC_Group"

---
**Done**
- Fix column names in big_table file
- Add option to intersect two gene sets (replace them or add them to groups_dic)
- ADD "std/sem/none" parameter in plotting functions


# More
DATA:
- ATAC-seq signal
- ALL kinds of features ((hrde-1 score, H3K9, mRNA level...)


Function to calculate ATAC_contin_param for all genes:
- choose range: -500:-100
- choose experiment: exp1 / exp_gonads / exp_hrde_gonads / exp_hrde_guy / exp_metsetset
- Normalize: to highly_lowly / none
- choose "treated" condition: gfp / oma-1 / hrde-1 / sx / metsetset / 
- Mean 3 all replicates? : yes or seperately?
- log2 FC
- median
returns a table/dic, for each gene the value of paramaters

** To test the range: it can be interesting to plot a correlation between the range of 500:100 to a range of 150-450 for instance. 
** To test the normalization: plot a correlation between results with and without normalization


Function to get list of desired genes:
- groups from table: by ATAC_bigtable (hrde-1 score, H3K9, mRNA level...)
- a set of the most / least 10% of ATAC_contin_param for a desired condition (choose percentage)


Function plot ATAC-seq signal. Choose:
- Experiment: exp1 / exp_gonads / exp_hrde_gonads / exp_hrde_guy / exp_metsetset
- set of genes: list
- Normalize: to highly_lowly / none
- Mean 3 all replicates? : yes or seperately?
- STD size?
- Plot


Test models:
- Can we predict "strange genes" from one of the features?





## to do
- ATAC-signal + violin FC for:
* WBGene00015351
* WBGene00016177
* WBGene00000224
* WBGene00013100

- median instead of means
- for every group for which we want to show an ATAC-fold change as a group:
    - add bootstrap dotplot:
        * each iteration takes n genes and mean/median of the fc_score of all of them
        * we plot the vector which is the size of i_iteration, each dot is a mean/median
        * test - how many of these bootstrap values in this vector are more extreme than our median.


- ATAC_Signal z-score? no

Talk to hila:
    * change "fold_change" to difference in integral

- Nearbys: What is counted as "close"?
    * What is the score of difference?
    * Maybe we need to score the closest 2000, then after 2000, then after... and see what is the radius of "nearby" 

- Hrde-1 Nearbys: Are hrde-1 targets isolated? Because genes will escape from the "bystander zone".

** Kolmogorov-Smirnoff: aqqumulative distribution


- change highly and lowly in hrde-1 plotting to hrde-1 highly and lowly


## done
- hrde-1 nearbys:
    - get all "fc_scores" of :
        * hrde-1 regulated 
        * hrde-1 targeted
        * hrde-1 upstream
        * hrde-1 downstream
    - plot violin plot / dotplot

* zscore all columns? Or does the negative screws it up?
    * plot groups to see if normalized (highly, lowly, ...)

- fix bootstrap bug and look at "hrde-1 nearbys upstream"

- Make it possible to remove a replicate (for all plots)
    * plot gene signal
    * plot group signal (mean_flag)
    * plot gene distribution (mean_flag)

- scatter plot of fold-changes (with all genes, all replicates, our gene dotted in color with std whiskers) + dotted line of 95% [can also use to compare two experiments]
    - plot not the FC_param itself but its Z-score / rank instead.
- add function "bootstrap_gene_list(df_sample, list_size, gene_pool: default[protein_coding], num_of_iters)"
- test functions of fold-change and signal calculating (spr-2 fold change)
- Add possibility to plot only between range
- plot lines of groups of genes (default = highly + lowly)