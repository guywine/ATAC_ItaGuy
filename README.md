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
- problem with rep 2/4 of exp1
- as groups, add function "bootstrap_gene_list(df_sample, list_size, gene_pool: default[protein_coding], num_of_iters)"
'''
Each iteration, select random group, mean signal.
Than returns the mean of all means, and variance of all means.
'''
- bystanders function
- scatter plot of fold-changes (with all genes, all replicates, our gene dotted in color with std whiskers) + dotted line of 95% [can also use to compare two experiments]
    - plot not the FC_param itself but its Z-score / rank instead.


** Normalize Z-score. (x_value - mean_all)/std_all
** Kolmogorov-Smirnoff: aqqumulative distribution

## done
- test functions of fold-change and signal calculating (spr-2 fold change)
- Add possibility to plot only between range
- plot lines of groups of genes (default = highly + lowly)