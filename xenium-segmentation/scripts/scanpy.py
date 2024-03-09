#!/usr/bin/env python
# coding: utf-8

# Scanpy

import pandas as pd
import tifffile
import os
import json
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import anndata as ad
import warnings

warnings.filterwarnings('ignore')

# Import the transcripts CSV

# Add path to the Baysor ouput

output_path = 'data'
transcripts = pd.read_csv(os.path.join(output_path, 'transcripts.csv'))

# Create a cell x gene table

cross_tab = pd.crosstab(index=transcripts["cell"].values,
                        columns=transcripts['gene'].values)

# Get the spatial position of the cells. Here we just take the mean of x and y.

spatial = transcripts[~pd.isna(transcripts.cell)]
spatial = spatial.groupby("cell")[['x', 'y']].mean()
spatial = spatial.reindex(cross_tab.index)

# Put it together in an anndata object. This is also saved.

adata = ad.AnnData(
    X=cross_tab,
    obs=pd.DataFrame(
        index=cross_tab.index.values,
        data={
        'cell':cross_tab.index.values
    }),
    var=pd.DataFrame(
        index=cross_tab.columns,
        data={
            'gene':cross_tab.columns.values
        }
    )
)

adata.layers['raw'] = adata.X
adata.obsm['spatial'] = spatial.to_numpy()

adata.write("anndata.h5ad")

# Preprocessing

sc.pl.highest_expr_genes(adata, n_top=20)

# Basic filtering

sc.pp.filter_cells(adata, min_genes=20)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

# Principal component analysis

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata)
sc.pl.pca_variance_ratio(adata, log=True)

# Computing the neighborhood graph

sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.pl.umap(adata, color=['leiden'])
sc.pl.embedding(adata, basis='spatial', color=['leiden'])

# Subset to a smaller FOV

max_width = int(os.getenv("WIDTH"))
max_height = int(os.getenv("HEIGHT"))

x_offset = int(os.getenv("X_OFFSET"))
y_offset = int(os.getenv("Y_OFFSET"))

def import_image(path: str):
    file = os.path.join(path, "morphology_mip.ome.tif")
    img = tifffile.imread(file)
    return img

img = import_image("data/xenium")

def get_pixel_size(path: str) -> float:
    file = open(os.path.join(path, "experiment.xenium"))
    experiment = json.load(file)
    pixel_size = experiment['pixel_size']
    return pixel_size

pixel_size = get_pixel_size("data/xenium")

# check boundaries

if max_width > img.shape[1]:
    max_width = img.shape[1]
    x_offset = 0
if max_height > img.shape[0]:
    max_height = img.shape[0]
    y_offset = 0

if (x_offset < 0) and (img.shape[1] > max_width):
    x_offset = round(img.shape[1] /2 - max_width /2)
if (max_width + x_offset) > img.shape[1]:
    x_offset = img.shape[1] - max_width

if (y_offset < 0) and (img.shape[0] > max_height):
    y_offset = round(img.shape[0] /2 - max_height /2)
if (max_height + y_offset) > img.shape[0]:
    y_offset = img.shape[0] - max_height

ax = sc.pl.embedding(adata, basis='spatial', color=['leiden'], show=False, size=75)
ax.set_xlim((x_offset * pixel_size,  (x_offset + max_width) * pixel_size  ))
ax.set_ylim(((y_offset + max_height) * pixel_size, y_offset * pixel_size  ))

# Differential expression

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
marker_genes = np.unique(pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5).to_numpy().flatten())
sc.pl.dotplot(adata, marker_genes, groupby='leiden')
