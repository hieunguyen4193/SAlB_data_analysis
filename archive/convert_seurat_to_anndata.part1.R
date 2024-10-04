gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

cluster.resolution <- 0.5

#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
path.to.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "20231018_SAlBounny"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")

path.to.08.output <- file.path(path.to.main.output, "08_output")
path.to.seurat2anndata <- file.path(path.to.08.output, "seurat2anndata")
dir.create(path.to.seurat2anndata, showWarnings = FALSE, recursive = TRUE)

path.to.seurat.obj <- file.path(path.to.03.output, sprintf("%s.filter_contaminated_cells.clusterRes_%s.rds", PROJECT, cluster.resolution))
s.obj <- readRDS(path.to.seurat.obj)

object.name <- str_replace(basename(path.to.seurat.obj), ".rds", "")

s.obj$barcode <- colnames(s.obj)

s.obj$UMAP_1 <- s.obj@reductions$INTE_UMAP@cell.embeddings[,1]
s.obj$UMAP_2 <- s.obj@reductions$INTE_UMAP@cell.embeddings[,2]
write.csv(s.obj@reductions$INTE_UMAP@cell.embeddings, file=file.path(path.to.seurat2anndata, sprintf('pca_%s.csv', object.name)), quote=F, row.names=F)

write.csv(s.obj@meta.data, file=file.path(path.to.seurat2anndata, sprintf('metadata_%s.csv', object.name)), quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(s.obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=file.path(path.to.seurat2anndata, sprintf('counts_%s.mtx', object.name)))

# write gene names
write.table( data.frame('gene'=rownames(counts_matrix)),file=file.path(path.to.seurat2anndata, sprintf('gene_names_%s.csv', object.name)),
             quote=F,row.names=F,col.names=F)