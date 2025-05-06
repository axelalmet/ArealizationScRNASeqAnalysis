import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
# plt.rcParams['font.family'] = 'Arial'

gene_lists = pd.read_csv('../scrnaseqdata/ArealGenes.csv')

areal_identities = gene_lists['area'].dropna().unique().tolist()

areal_gene_lists = {}

for area in areal_identities:
    areal_gene_lists[area] = gene_lists.loc[gene_lists['area'] == area, 'gene'].unique().tolist()

adata = sc.read_h5ad('../scrnaseqdata/cortical_organoid_merged.h5ad')

for area in areal_gene_lists:
    sc.tl.score_genes(adata, areal_gene_lists[area], score_name=area)

samples = adata.obs['samples'].unique().tolist()

plt.rcParams['figure.figsize'] = (4, 3)
days_and_celltypes = {'D35': ['aRG', 'PP/CR', 'IP'],
                      'D56': ['aRG', 'DLN', 'IP', 'newborn neurons', 'CR'],
                      'D97': ['aRG/oRG', 'CR', 'IP', 'DLN', 'ULN', 'newborn neurons']}

for day, celltypes in days_and_celltypes.items():

    for celltype in celltypes:
        adata_sample = adata[(adata.obs['day'] == day)&(adata.obs['celltype'] == celltype)]

        for area in areal_gene_lists:
            sc.pl.violin(adata_sample,
                         keys=area,
                         groupby='condition',
                         order=['FGF', 'ctr', 'BMP'],
                        rotation=90, 
                        save='cortical_organoid_' + day + '_' + celltype.replace(' ', '-').replace('/', '-') + '_' + area.replace(' ', '-').replace('/', '+') + '.pdf')

adata.write_h5ad('../scrnaseqdata/cortical_organoid_merged.h5ad', compression='gzip')