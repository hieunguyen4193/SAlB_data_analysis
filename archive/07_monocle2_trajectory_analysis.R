
my_random_seed <- 42
set.seed(my_random_seed)

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

#####-----------------------------------------------------------------------#####
##### install monocle
#####-----------------------------------------------------------------------#####
# if ("monocle3" %in% installed.packages()){
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
# }
#####-----------------------------------------------------------------------#####
library(devtools)
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
monocle.obj <- readRDS(file.path(path.to.07.output, sprintf("monocle_obj_cluster_res_%s", cluster.resolution)))


# Plot

DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE)


## UMAP

plot_cell_trajectory(monocle.obj, color_by = "seurat_clusters")


## Cell trajectory from `monocle2`

plot_cell_trajectory(monocle.obj)




plot_cell_trajectory(monocle.obj, color_by = "Pseudotime")




##### Add pseudotime information to the main seurat object
seurat.metadata <- s.obj@meta.data %>% rownames_to_column("barcode")
monocle.metadata <- pData(monocle.obj) %>% rownames_to_column("barcode") %>%
  subset(select = c(barcode, Pseudotime, State))

meta.data <- merge(seurat.metadata, monocle.metadata, by.x = "barcode", by.y = "barcode") %>%
  column_to_rownames("barcode")
meta.data <- meta.data[row.names(s.obj@meta.data), ]
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$Pseudotime, col.name = "Pseudotime")
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$State, col.name = "State")


## UMAP: grouped by trajectory STATE

DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "State")




FeaturePlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, features = c("Pseudotime"))


# Differential gene expression tests (trajectories)

## State vs State

if (file.exists(file.path(path.to.07.output, sprintf("pairwise_state_markers"))) == FALSE){
  num.states <- length(unique(s.obj$State))
  i.s <- seq(1, num.states - 1)
  j.s <- seq(2, num.states)
  
  pairwise.state.markers <- hash()
  for (item in seq(1, length(i.s))){
    state1 <- i.s[[item]]
    state2 <- j.s[[item]]
    tmp.state.markers <- FindMarkers(object = s.obj, assay = "RNA", test.use = "wilcox", group.by = "State", ident.1 = state1, ident.2 = state2)
    pairwise.state.markers[[sprintf("%s_%s", state1, state2)]] <- tmp.state.markers %>% subset(p_val_adj <= 0.05)
  } 
  saveRDS(pairwise.state.markers, file.path(path.to.07.output, sprintf("pairwise_state_markers")))
} else {
  pairwise.state.markers <- readRDS(file.path(path.to.07.output, sprintf("pairwise_state_markers")))
}


## Differential gene expression test along trajectories (VGAM)

pseudotime.markers <- differentialGeneTest(monocle.obj, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 16, verbose = TRUE)




sig_gene_names <- row.names(subset(pseudotime.markers, qval < 0.1))
plot_pseudotime_heatmap(monocle.obj[head(sig_gene_names, 10),],
                num_clusters = 3,
                cores = 1,
                show_rownames = T)


## Differential gene expression test (TradeSeq)

devtools::install_github("statOmics/tradeSeq")
library(tradeSeq)
extract_monocle_info <- function(cds) {
  if (cds@dim_reduce_type != "DDRTree") {
    stop(paste0("For now tradeSeq only support Monocle with DDRTree",
                "reduction. If you want to use another type",
                "please use another format for tradeSeq inputs."))
  }
  # Get the reduced dimension of DDRT
  rd <- t(monocle::reducedDimS(cds)) %>% as.data.frame()
  
  # Get the various lineages info for weights and pseudotime
  y_to_cells <- cds@auxOrderingData[["DDRTree"]]
  y_to_cells <- y_to_cells$pr_graph_cell_proj_closest_vertex %>%
    as.data.frame()
  y_to_cells$cells <- rownames(y_to_cells)
  y_to_cells$Y <- y_to_cells$V1
  root <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  root <- y_to_cells$Y[y_to_cells$cells == root]
  mst <- monocle::minSpanningTree(cds)
  endpoints <- names(which(igraph::degree(mst) == 1))
  endpoints <- endpoints[endpoints != paste0("Y_", root)]
  cellWeights <- lapply(endpoints, function(endpoint) {
    path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
    path <- as.character(path)
    df <- y_to_cells[y_to_cells$Y %in% path, ]
    df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
    colnames(df) <- endpoint
    return(df)
  }) %>% do.call(what = 'cbind', args = .)
  pseudotime <- sapply(cellWeights, function(w) cds$Pseudotime)
  rownames(cellWeights) <- rownames(pseudotime) <- colnames(cds)
  # Get the lineages representation
  edges_rd <- t(monocle::reducedDimK(cds)) %>% as.data.frame()
  rd_lineages <- lapply(endpoints, function(endpoint){
    path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
    path <- as.character(path)
    path <- paste("Y", path, sep = "_")
    return(edges_rd[path, ])
  })
  return(list("pseudotime" = pseudotime,
              "cellWeights" = as.matrix(cellWeights)))
}
info <- extract_monocle_info(monocle.obj)
sce <- fitGAM(counts = Biobase::exprs(monocle.obj),
              cellWeights = info$cellWeights,
              pseudotime = info$pseudotime)
saveRDS(sce, file.path(path.to.07.output, "sce.rds"))


# Branch analysis (monocle-BEAM)

