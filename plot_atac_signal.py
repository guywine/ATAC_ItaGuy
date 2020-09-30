'''
Plot ATAC-seq signal. Choose:
- Experiment: exp1 / exp_gonads / exp_hrde_gonads / exp_hrde_guy / exp_metsetset
- set of genes: list
- Normalize: highly_lowly / none
- Mean all replicates or display them seperately?
'''
import read_tables as rdt 
import matplotlib.pyplot as plt
import pandas as pd



exp1_dfs_list = rdt.read_experiment(exp_name='exp1')
gfp_0, gfp_1, gfp_2, gfp_4, oma1_0, oma1_1, oma1_2, oma1_4 = exp1_dfs_list

## example list of genes wbid:
genes = list(gfp_0.index[:3])

only_genes = gfp_0.loc[genes]
mean_of_genes = only_genes.mean()

plt.plot(mean_of_genes)
x = list(range(-1000, 1001))
a = gfp_0.iloc[0]
b = gfp_0.iloc[1]
c = gfp_0.iloc[2]

plt.plot(x, a, x, b, x, c)
plt.show()

df = pd.DataFrame({'a':a,'b':b,'c':c})
df.plot()




