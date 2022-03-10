# ATAC_ItaGuy
## Hackers of the ATAC project




** To test the range: it can be interesting to plot a correlation between the range of 500:100 to a range of 150-450 for instance. 
** To test the normalization: plot a correlation between results with and without normalization

- Hrde-1 Nearbys: Are hrde-1 targets isolated? Because genes will escape from the "bystander zone".

** Kolmogorov-Smirnoff: aqqumulative distribution

- change highly and lowly in hrde-1 plotting to hrde-1 highly and lowly




## to do

"hrde-1-Kennedy" strange normalization and behaviour
- Look at the overlap between the kennedy list (1459) to "highly expressed" (1018): only 88 genes

- Look for outliers: Cat-plot the FC values.
- Look for outliers: Get a table of atac and mRNA values.

"hrde-1 regulated" bystander effect of downstream genes
- Plot group mean atac-signal in the background of a bootstrap atac signal.
- Plot bootstrap of group mean mRNA-FC (similar to what I did in the ATAC-FC)

by-stander of oma-1/GFP:
- Look for pattern of effect in regards to several factors.


- normalization: 
    - add_to_avoid_zero lowest value in both conditions after normalization
    - if not working: filter low values

- violinplot of each condition: lowly, highly, kennedy, hrde-reg

- test "bootstrap" function




- by-stander hrde-1: do it only where the hrde-1 regulated get really condensed.

- Graph that shows that the least variable is the hotspot (most conserved)

- nearbys: all values are fold-change (groups that add + groups segregated [0-2500, 2500-5000])
    - violin plot of every group
    +
    - bootstrap of mean for every group





## done

- run bootstrap of hrde-1 groups on mRNA data. 
    - plot scatter of both conditions and mark hrde-groups.
    - bootstrap mRNA-FC to see if it's the same as ATAC-FC

- If normaization doesn't work: Plot, scatter: each dot is a gene
    - x axis: value of hotspot in WT (*mean and median)
    - y axis: value of hotspot in hrde-1 mutant (*mean and median)
- Color hrde-regulated genes in red (the rest in grey).

- Normalization:
    - send hila of all_aligned to see if I messed up
    - Check for mistakes

- missing genes in hrde-1 lists: send Hila the missing genes

- for every group for which we want to show an ATAC-fold change as a group:
    - add bootstrap dotplot:
        * each iteration takes n genes and mean/median of the fc_score of all of them
        * we plot the vector which is the size of i_iteration, each dot is a mean/median
        * test - how many of these bootstrap values in this vector are more extreme than our median.


- Nearbys: What is counted as "close"?
    * What is the score of difference?
    * Maybe we need to score the closest 2000, then after 2000, then after... and see what is the radius of "nearby" 


- ATAC-signal + violin FC for:
* WBGene00015351
* WBGene00016177
* WBGene00000224
* WBGene00013100


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



** to think about **
