##### clean up #####
gc()
rm(list = ls())

my_random_seed <- 42
set.seed(my_random_seed)

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

library(devtools)
if ("monocle" %in% installed.packages() == FALSE){
  BiocManager::install("monocle", update = FALSE)
}
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

data <- GetAssayData(s.obj, slot = "count", assay = "RNA")

pd <- new('AnnotatedDataFrame', data = s.obj@meta.data)

fd <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fd)

library(monocle)
monocle.obj <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
saveRDS(monocle.obj, file.path(path.to.monocle2.input, sprintf("SAlBounny_full.filter_contaminated_cells.clusterRes_%s.monocle2.rds", cluster.resolution)))