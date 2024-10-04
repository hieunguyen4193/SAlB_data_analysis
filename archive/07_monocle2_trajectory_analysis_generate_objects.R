#####----------------------------------------------------------------------#####
#
# trnguyen@ukaachen.de
#
#####----------------------------------------------------------------------#####

##### clean up #####
# gc()
# rm(list = ls())

my_random_seed <- 42
set.seed(my_random_seed)

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

library(devtools)

if ("monocle3" %in% installed.packages()){
  remove.packages("monocle3")
  devtools::install_github("cysouw/qlcMatrix")
  install.packages("DDRTree")
  install.packages("densityClust")
  BiocManager::install("monocle.objSingleCell", update = FALSE)
  install.packages("fastICA")
  BiocManager::install("biocViews", update = FALSE)
  # remove.packages("BiocGenerics")
  BiocManager::install("HSMMSingleCell", update = FALSE)
  # install.packages("/media/hieunguyen/HD0/storage/offline_pkgs/monocle_2.30.0.tar.gz", type = "source", repos = NULL)
  install.packages("/media/hieunguyen/HD0/storage/offline_pkgs/monocle_2.30.0/monocle", type = "source", repos = NULL)
  BiocManager::install("tradeSeq", update = FALSE)
}

library(monocle)

outdir <- "/media/hieunguyen/CRC1382H/CRC1382/outdir"
PROJECT <- "20231018_SAlBounny"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.07.output <- file.path(path.to.main.output, "07_output")
dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)

cluster.resolution <- 0.5
path.to.input.s.obj <- file.path(path.to.03.output, sprintf("20231018_SAlBounny.filter_contaminated_cells.clusterRes_%s.rds", cluster.resolution))
s.obj <- readRDS(path.to.input.s.obj)

if (file.exists(file.path(path.to.07.output, sprintf("monocle_obj_cluster_res_%s", cluster.resolution))) == FALSE){
  
  ##### CREATE MONOCLE OBJECT
  data <- GetAssayData(s.obj, slot = "count", assay = "RNA")
  
  pd <- new('AnnotatedDataFrame', data = s.obj@meta.data)
  
  fd <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fd)
  
  monocle.obj <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  monocle.obj <- estimateSizeFactors(monocle.obj)
  monocle.obj <- estimateDispersions(monocle.obj)
  
  monocle.obj <- detectGenes(monocle.obj, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(monocle.obj),
                                      num_cells_expressed >= 100))
  
  fData(monocle.obj)$use_for_ordering <-
    fData(monocle.obj)$num_cells_expressed > 0.05 * ncol(monocle.obj)
  
  ordering.genes <- subset(fData(monocle.obj), fData(monocle.obj)$use_for_ordering == TRUE)$gene_short_name
  ordering.genes <- intersect(ordering.genes, expressed_genes)
  monocle.obj <- monocle.obj[ordering.genes,]
  
  set.seed(my_random_seed)
  monocle.obj <- reduceDimension(monocle.obj,
                                 max_components = 2,
                                 norm_method = 'log',
                                 num_dim = 3,
                                 reduction_method = 'tSNE',
                                 verbose = T, 
                                 random_seed = my_random_seed)
  set.seed(my_random_seed)
  monocle.obj <- clusterCells(monocle.obj, verbose = F, random_seed = my_random_seed)
  
  clustering_DEG_genes <-
    differentialGeneTest(monocle.obj[ordering.genes,],
                         fullModelFormulaStr = '~Cluster',
                         cores = 20)
  
  HSMM_ordering_genes <-
    row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
  
  monocle.obj <-
    setOrderingFilter(monocle.obj,
                      ordering_genes = HSMM_ordering_genes)
  set.seed(my_random_seed)
  monocle.obj <-
    reduceDimension(monocle.obj, method = 'DDRTree', random_seed = my_random_seed)
  set.seed(my_random_seed)
  monocle.obj <- orderCells(monocle.obj)
  
  p <- plot_cell_trajectory(monocle.obj)
  
  saveRDS(monocle.obj, file.path(path.to.07.output, sprintf("monocle_obj_cluster_res_%s", cluster.resolution)) )
  saveRDS(data.frame(status = c(sprintf("finished case_cluster_res %s", cluster.resolution))), file.path(path.to.07.output, sprintf("finished case_cluster_res_%s.csv", cluster.resolution)))
}
