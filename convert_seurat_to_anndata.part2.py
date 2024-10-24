import scanpy as sc
import numpy as np
import os
import anndata2ri
import pathlib
from scipy import io
import anndata#
import pandas as pd
from tqdm import tqdm
import argparse
import sys

# Activate the anndata2ri conversion between SingleCellExperiment and AnnData
# anndata2ri.activate()

#Loading the rpy2 extension enables cell magic to be used
#This runs R code in jupyter notebook cells
# %load_ext rpy2.ipython

sc.settings.verbosity = 3
# sc.logging.print_versions()

import warnings
warnings.filterwarnings("ignore")

#####---------------------------------------------------------------------------------------------------------#####
##### CONFIGURATIONS
outdir = "/home/hieunguyen/CRC1382/outdir"
PROJECT = "SAlBounny_full.filter_contaminated_cells.clusterRes_0.5"
output.version = "20241021"

path_to_main_output = os.path.join(outdir, "20231018_SAlBounny", "data_analysis")
path_to_01_output = os.path.join(path_to_main_output, "01_output")
path_to_02_output = os.path.join(path_to_main_output, "02_output")
path_to_03_output = os.path.join(path_to_main_output, "03_output")
path_to_08_output = os.path.join(path_to_main_output, "08_output")

path_to_seurat2anndata = os.path.join(path_to_08_output, "seurat2anndata")
#####---------------------------------------------------------------------------------------------------------#####
# load sparse matrix:
print(os.path.join(path_to_seurat2anndata, "counts_{}.mtx".format(PROJECT)))
X = io.mmread(os.path.join(path_to_seurat2anndata, "counts_{}.mtx".format(PROJECT)))

# create anndata object
adata = anndata.AnnData(X=X.transpose().tocsr())

# load cell metadata:
cell_meta = pd.read_csv(os.path.join(path_to_seurat2anndata, "metadata_{}.csv".format(PROJECT)))

# load gene names:
with open(os.path.join(path_to_seurat2anndata, "gene_names_{}.csv".format(PROJECT)), 'r') as f:
  gene_names = f.read().splitlines()
# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv(os.path.join(path_to_seurat2anndata, "pca_{}.csv".format(PROJECT)))
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# save dataset as anndata format
adata.write(os.path.join(path_to_seurat2anndata, '{}.h5ad'.format(PROJECT)))

# reload dataset
adata = sc.read_h5ad(os.path.join(path_to_seurat2anndata, '{}.h5ad'.format(PROJECT)))
