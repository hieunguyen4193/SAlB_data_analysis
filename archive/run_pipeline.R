gc()
rm(list = ls())
my_random_seed <- 42

set.seed(my_random_seed)
# __________VDJ DATA ANYLYSIS PIPELINE__________
PROJECT <- "20231018_SAlBounny"

# __________GEX DATA ANALYSIS PIPELINE__________
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline"
path2src <- file.path(path.to.pipeline.src, "processes_src")


source(file.path(path2src, "import_libraries.R"))
source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))

path2input <- file.path("/home/hieunguyen/CRC1382/storage", "SAlBounny_full")

# _____stage lst for single sample_____
stage_lst <- list()

stage_lst = c(
  A_48h = "A_48h",
  N_48h =   "N_48h",
  A_CD4 = "A_CD4",
  N_CD4 =   "N_CD4",
  Ad7_2 =   "Ad7_2",
  Nd7_2 = "Nd7_2",
  I     =   "I"
)

MINCELLS  <- 5
MINGENES  <- 50

save.RDS <- list(s1 = TRUE,
                 s2 = TRUE,
                 s3 = TRUE,
                 s4 = TRUE,
                 s5 = TRUE,
                 s6 = TRUE,
                 s7 = FALSE,
                 s8 = TRUE,
                 s8a = TRUE,
                 s9 = FALSE)

sw <- list(s1 = "on",
           s2 = "on",
           s3 = "on",
           s4 = "on",
           s5 = "on",
           s6 = "on",
           s7 = "off",
           s8 = "on",
           s8a = "off",
           s9 = "off")

rerun <- list(s1 = FALSE, 
              s2 = FALSE,
              s3 = FALSE,
              s4 = FALSE,
              s5 = FALSE,
              s6 = FALSE,
              s7 = FALSE,
              s8 = FALSE,
              s8a = FALSE,
              s9 = FALSE)


filter.thresholds <- list(nFeatureRNAfloor = 300,
                          nFeatureRNAceiling = NULL,
                          nCountRNAfloor = 500, 
                          nCountRNAceiling = 40000,
                          pct_mitofloor = NULL, 
                          pct_mitoceiling = 25,
                          pct_ribofloor = NULL, 
                          pct_riboceiling = NULL,
                          ambientRNA_thres = 0.5,
                          log10GenesPerUMI = 0.8)

remove_doublet <- FALSE
path.to.10X.doublet.estimation <- "/home/hieunguyen/CRC1382/storage/DoubletEstimation10X.csv"


num.PCA <- 30
num.PC.used.in.UMAP <- 30
num.PC.used.in.Clustering <- 30
num.dim.integration <- 30
num.dim.cluster <- 30
cluster.resolution <- 0.5

path.to.output <- file.path("/home/hieunguyen/CRC1382/outdir/SAlBounny", PROJECT)
dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)

filtered.barcodes <- NULL

s.obj <- run_pipeline_GEX(path2src = path2src,
                          path2input = path2input,
                          path.to.logfile.dir = file.path(path.to.output, "logs"),
                          stage_lst = stage_lst,
                          path.to.10X.doublet.estimation = path.to.10X.doublet.estimation,
                          MINCELLS = MINCELLS,
                          MINGENES = MINGENES,
                          PROJECT = PROJECT,
                          remove_doublet = remove_doublet,
                          save.RDS = save.RDS,
                          path.to.output = path.to.output,
                          rerun = rerun, 
                          DE.test = "wilcox",
                          num.PCA = num.PCA,
                          num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                          num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                          use.sctransform = FALSE,
                          filtered.barcodes = filtered.barcodes,
                          filter.thresholds = filter.thresholds,
                          input.method = "normal",
                          cluster.resolution = cluster.resolution,
                          num.dim.integration = num.dim.integration,
                          inte_pca_reduction_name = "INTE_PCA",
                          inte_umap_reduction_name = "INTE_UMAP",
                          with.VDJ = FALSE)

#### ALWAYS REMEMBER TO SAVE SESSIONINFO !!!!!!
writeLines(capture.output(sessionInfo()), file.path(path.to.output, sprintf("%s_sessionInfo.txt", PROJECT)))

