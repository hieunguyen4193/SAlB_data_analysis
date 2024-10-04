#####-----------------------------------------------------------------#####
# PATHS AND CONFIGURATIONS
#####-----------------------------------------------------------------#####
gc()
rm(list = ls())

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/CRC1382_SAlBounny_projects/CRC1382_SAlBounny_fulldata/RNA_velocity_end2end"

source(file.path(path.to.main.src, "R00_helper_functions.R"))
source(file.path(path.to.main.src, "R00_import_libraries.R"))

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "20231018_SAlBounny"
cluster.resolution <- 0.5
umap_assay <- "INTE_UMAP"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path_to_seurat2anndata <- file.path(path.to.main.output, "seurat2anndata")
dir.create(path_to_seurat2anndata, showWarnings = FALSE, recursive = TRUE)

path.to.input.seurat <- file.path(path.to.main.output, 
                                  "03_output", 
                                  sprintf("20231018_SAlBounny.filter_contaminated_cells.clusterRes_%s.subClusterFoxp3_cluster6.rds", cluster.resolution))
object.name <- sprintf("%s_03_output", PROJECT)

s.obj <- readRDS(path.to.input.seurat)

s.obj$barcode <- colnames(s.obj)

s.obj$UMAP_1 <- s.obj@reductions[[umap_assay]]@cell.embeddings[,1]
s.obj$UMAP_2 <- s.obj@reductions[[umap_assay]]@cell.embeddings[,2]

write.csv(s.obj@meta.data, file=file.path(path_to_seurat2anndata, sprintf('metadata_%s.csv', object.name)), quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(s.obj, assay='RNA', slot='data')
writeMM(counts_matrix, file=file.path(path_to_seurat2anndata, sprintf('counts_%s.mtx', object.name)))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(s.obj@reductions[[str_replace(umap_assay, "UMAP", "PCA")]]@cell.embeddings, file=file.path(path_to_seurat2anndata, sprintf('pca_%s.csv', object.name)), quote=F, row.names=F)

# write gene names
write.table( data.frame('gene'=rownames(counts_matrix)),file=file.path(path_to_seurat2anndata, sprintf('gene_names_%s.csv', object.name)),
             quote=F,row.names=F,col.names=F)

################################################
DimPlot(object = s.obj, reduction = umap_assay, label = TRUE, label.box = TRUE, repel = TRUE, pt.size = 1)