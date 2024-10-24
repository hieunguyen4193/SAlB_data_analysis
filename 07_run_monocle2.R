##### clean up #####
gc()
rm(list = ls())

my_random_seed <- 42
set.seed(my_random_seed)
# to solve unable to access index for repository https://mran.microsoft.com/snapshot/2020-07-16/src/contrib
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org"
options(repos=r)})

if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

library(devtools)
library(monocle)

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "SAlBounny_full"
output.version <- "20241021"

path.to.main.output <- file.path(outdir, PROJECT, output.version, "data_analysis")
path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.07.output <- file.path(path.to.main.output, "07_output")
dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)

path.to.monocle2.input <- file.path(path.to.07.output, "monocle2_input")
dir.create(path.to.monocle2.input, showWarnings = FALSE, recursive = TRUE)

cluster.resolution <- 0.5
path.to.input.s.obj <- file.path(path.to.03.output, sprintf("SAlBounny_full.filter_contaminated_cells.clusterRes_%s.rds", cluster.resolution))
s.obj <- readRDS(path.to.input.s.obj)

#####----------------------------------------------------------------------#####
##### RUN MONOCLE2 FROM S.OBJ
#####----------------------------------------------------------------------#####
run_monocle2_from_presave_obj <- function(monocle.obj, path.to.save.monocle.obj){
  library(monocle)
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
  
  saveRDS(monocle.obj, file.path(path.to.save.monocle.obj, sprintf("monocle_obj.rds")))
  return(monocle.obj)
}

#####----------------------------------------------------------------------#####
##### MAIN RUN
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.07.output, "monocle_obj.rds")) == FALSE){
  print("reading in monocle object prepared from seurat object data...")
  monocle.obj <- readRDS(file.path(path.to.monocle2.input, 
                                   sprintf("SAlBounny_full.filter_contaminated_cells.clusterRes_%s.monocle2.rds", cluster.resolution)))
  monocle.obj <- run_monocle2_from_presave_obj(monocle.obj, path.to.07.output)  
} else {
  monocle.obj <- readRDS(file.path(path.to.07.output, "monocle_obj.rds"))
}


##### plot cell trajectory, color by seurat clusters
p <- plot_cell_trajectory(monocle.obj, color_by = "seurat_clusters")
ggsave(plot = p, filename = sprintf("cell_trajectory.seurat_clsuters.svg"), path = path.to.07.output, device = "svg", dpi = 300, width = 14, height = 10)

##### plot cell trajectory, color by monocle2 states
p <- plot_cell_trajectory(monocle.obj, color_by = "State")
ggsave(plot = p, filename = sprintf("cell_trajectory.State.svg"), path = path.to.07.output, device = "svg", dpi = 300, width = 14, height = 10)

##### plot cell trajectory, color by pseudotime
p <- plot_cell_trajectory(monocle.obj, color_by = "Pseudotime")
ggsave(plot = p, filename = sprintf("cell_trajectory.pseudotime.svg"), path = path.to.07.output, device = "svg", dpi = 300, width = 14, height = 10)

##### plot cell trajectory, color by pseudotime
monocle.obj.reverse <- orderCells(monocle.obj, reverse = TRUE)
p <- plot_cell_trajectory(monocle.obj.reverse, color_by = "Pseudotime")
ggsave(plot = p, filename = sprintf("cell_trajectory.rev_Pseudotime.svg"), path = path.to.07.output, device = "svg", dpi = 300, width = 14, height = 10)

##### save monocle data to csv file
monocledf <- data.frame(
  barcode = colnames(monocle.obj),
  state = monocle.obj$State,
  pseudotime = monocle.obj$Pseudotime
)
monocle.reversedf <- data.frame(
  barcode = colnames(monocle.obj.reverse),
  state = monocle.obj.reverse$State,
  pseudotime = monocle.obj.reverse$Pseudotime
)
write.csv(monocledf, file.path(path.to.07.output, "monocledf.csv"))
write.csv(monocle.reversedf, file.path(path.to.07.output, "monocledf.rev.csv"))