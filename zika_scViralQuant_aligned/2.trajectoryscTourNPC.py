# activate newer version of python 3.8 source /share/lwhitmo/pythonvenvs/py3.8venv/bin/activate
import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

countmatrix='/share/lwhitmo/projects/Stokes_HFB_iNPC_Analysis/zika_scViralQuant_aligned/inpc_results/Convert4SCtour/counts_npc.txt'
metadata='/share/lwhitmo/projects/Stokes_HFB_iNPC_Analysis/zika_scViralQuant_aligned/inpc_results/Convert4SCtour/metadata_npc.csv'
results_folder='2.TrajectoryAnalysisSCtouriNPC_NPC/'

adata= sc.read(countmatrix).T
info = pd.read_csv(metadata, sep=';', index_col=0)
cells = adata.obs_names.intersection(info.index)
adata.obs['nFeature_RNA'] = info.loc[cells, 'nFeature_RNA'].copy()
adata.obs['orig.ident'] = info.loc[cells, 'orig.ident'].copy()
adata.obs['garnett_cluster_extend_lw'] = info.loc[cells, 'garnett_cluster_extend_lw'].copy()
adata.obs['ZIKA'] = info.loc[cells, 'ZIKA'].copy()
adata.obs['condition'] = info.loc[cells, 'condition'].copy()

#adata = adata.raw
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
sc.pp.filter_genes(adata, min_cells=20)
sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000, subset=True)

#By default the percentage of cells used to train the model is set to 0.9 
#when the total number of cells is less than 10,000 and 0.2 when greater than 10,000. 
#Users can adjust the percentage by using the parameter percent (for example percent=0.6).
tnode = sct.train.Trainer(adata, loss_mode='nb')
tnode.train()

#Infer the developmental pseudotime based on the trained model.
adata.obs['ptime'] = tnode.get_time()

#Larger alpha_z skews the latent space towards the intrinsic transcriptomic 
#structure while larger alpha_predz is more representative of the extrinsic
#pseudotime ordering. Users can adjust the two parameters according to their purposes.
#Currently using default
mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.2, alpha_predz=0.8)
adata.obsm['X_TNODE'] = mix_zs

#Infer the transcriptomic vector field.
adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])

# Generate a UMAP embedding based on the inferred latent space.
# Optionally, you can order the cells according to their pseudotime 
# before this step, which is demonstrated to yield a better trajectory for some datasets.
adata = adata[np.argsort(adata.obs['ptime'].values), :]
sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=15)
sc.tl.umap(adata, min_dist=0.1)

fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(18, 5))
sc.pl.umap(adata, color='ptime', size=20, ax=axs[0], show=False)
sc.pl.umap(adata, color='garnett_cluster_extend_lw', size=20, ax=axs[1], legend_loc='on data', show=False)
sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', 
                        use_rep_neigh='X_TNODE', color='garnett_cluster_extend_lw', ax=axs[2], 
                        legend_loc='on data', frameon=False, size=100, alpha=0.2)
plt.savefig(results_folder+"trajectory.png")

fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(18, 5))
sc.pl.umap(adata, color='ptime', size=20, ax=axs[0], show=False)
sc.pl.umap(adata, color='orig.ident', size=20, ax=axs[1], legend_loc='on data', show=False)
sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', 
                        use_rep_neigh='X_TNODE', color='orig.ident', ax=axs[2], 
                        legend_loc='on data', frameon=False, size=100, alpha=0.2)
plt.savefig(results_folder+"trajectory_orig.ident.png")

fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(18, 5))
sc.pl.umap(adata, color='ptime', size=20, ax=axs[0], show=False)
sc.pl.umap(adata, color='ZIKA', size=20, ax=axs[1], legend_loc='on data', show=False)
sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', 
                        use_rep_neigh='X_TNODE', color='ZIKA', ax=axs[2], 
                        legend_loc='on data', frameon=False, size=100, alpha=0.2)
plt.savefig(results_folder+"trajectory_ZIKA.png")

adata.obs['ptime_reverse'] = sct.train.reverse_time(adata.obs['ptime'].values)
sc.pl.umap(adata, color=['orig.ident', 'ptime_reverse'], legend_loc='on data')
plt.savefig(results_folder+"reversetrajectory.png")

adata.var.to_csv(results_folder+"variablegeneinfo.csv")
df = pd.DataFrame(adata.obsm['X_umap'], columns = ['UMAP_1','UMAP_2'])
df.index=adata.obs.index

df['nFeature_RNA'] = adata.obs['nFeature_RNA'].copy()
df['orig.ident'] = adata.obs['orig.ident'].copy()
df['garnett_cluster_extend_lw'] =  adata.obs['garnett_cluster_extend_lw'].copy()
df['ZIKA'] =  adata.obs['ZIKA'].copy()
df['ptime'] =  adata.obs[ 'ptime'].copy()
df['condition'] =  adata.obs[ 'condition'].copy()

df.to_csv(results_folder+"umapcoordinates.csv")
