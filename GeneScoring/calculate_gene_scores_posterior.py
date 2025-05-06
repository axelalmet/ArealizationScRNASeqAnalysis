import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
# plt.rcParams['font.family'] = 'Arial'

gene_lists = pd.read_csv('../scrnaseqdata/TemporalOccipitalGenes.csv')

genes_to_score = gene_lists['gene'].unique().tolist()

adata = sc.read_h5ad('../scrnaseqdata/cortical_organoid_merged.h5ad')

sc.tl.score_genes(adata, genes_to_score, score_name='Temporal/Occipital')

samples = adata.obs['samples'].unique().tolist()

plt.rcParams['figure.figsize'] = (4, 3)
days_and_celltypes = {'D35': ['aRG', 'IP', 'PP/CR'],
                      'D56': ['aRG', 'DLN', 'IP', 'newborn neurons', 'CR'],
                      'D97': ['aRG/oRG', 'CR', 'DLN', 'IP', 'ULN', 'newborn neurons']}

for day, celltypes in days_and_celltypes.items():

    for celltype in celltypes:
        adata_sample = adata[(adata.obs['day'] == day)&(adata.obs['celltype'] == celltype)]

        sc.pl.violin(adata_sample,
                    keys='Temporal/Occipital',
                    groupby='condition',
                    order=['FGF', 'ctr', 'BMP'],
                    rotation=90, 
                    save='cortical_organoid_' + day + '_' + celltype.replace(' ', '-').replace('/', '-') + '_TemporalOccipital.pdf')

adata.write_h5ad('../scrnaseqdata/cortical_organoid_merged.h5ad', compression='gzip')