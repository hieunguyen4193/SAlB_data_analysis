library(Seurat)
library(dplyr)
s.obj.old <- readRDS("/media/hieunguyen/HD01/outdir/CRC1382/20231018_SAlBounny/s8_output/20231018_SAlBounny.output.s8.rds")
s.obj.new <- readRDS("/home/hieunguyen/CRC1382/outdir/SAlBounny_full/test1/s8_output/SAlBounny_full.output.s8.rds")

s.obj <- RunUMAP(s.obj.new, 
                 reduction = "INTE_PCA", dims = 1:30, 
                 reduction.name="INTE_UMAP_2", 
                 seed.use = 0, 
                 umap.method = "uwot")
DimPlot(object = s.obj, reduction = "INTE_UMAP_2")

# s.obj.old <- readRDS("/media/hieunguyen/HD01/outdir/CRC1382/20231018_SAlBounny/s6_output/20231018_SAlBounny.output.s6.rds")
# s.obj.new <- readRDS("/home/hieunguyen/CRC1382/outdir/SAlBounny_20241019/20231018_SAlBounny/s6_output/20231018_SAlBounny.output.s6.rds")

count.old.integrated <- GetAssayData(object = s.obj.old, slot = "data", assay = "integrated")
count.new.integrated <- GetAssayData(object = s.obj.new, slot = "data", assay = "integrated")

count.old <- GetAssayData(object = s.obj.old, slot = "data", assay = "RNA")
count.new <- GetAssayData(object = s.obj.new, slot = "data", assay = "RNA")

# count.old.scale <- GetAssayData(object = s.obj.old, slot = "scale", assay = "RNA")
# count.new.scale <- GetAssayData(object = s.obj.new, slot = "scale", assay = "RNA")

count.old.scale <- GetAssayData(object = s.obj.old, slot = "scale", assay = "integrated")
count.new.scale <- GetAssayData(object = s.obj.new, slot = "scale", assay = "integrated")

path.to.save.output <- "/home/hieunguyen/CRC1382/outdir/check"
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
# count.old <- GetAssayData(object = s.obj.old, slot = "data", assay = "RNA")
# count.new <- GetAssayData(object = s.obj, slot = "data", assay = "RNA")

# sum(count.old[, "A_48h_A_48h_AAACCTGAGAACTCGG-1"] - count.new[, "A_48h_A_48h_AAACCTGAGAACTCGG-1"] )

old.commands <- s.obj.old@commands
new.commands <- s.obj.new@commands

for (f in names(old.commands)){
  # fileConn<-file(file.path(path.to.save.output, sprintf("old.s.obj.%s.txt", f)))
  # writeLines(old.commands[[f]], fileConn)
  # close(fileConn)
  
  # fileConn<-file(file.path(path.to.save.output, sprintf("new.s.obj.%s.txt", f)))
  # writeLines(new.commands[[f]], fileConn)
  # close(fileConn)  
  
  sink(file = file.path(path.to.save.output, sprintf("old.s.obj.%s.txt", f)))
  old.commands[[f]] %>% print()
  sink()
  sink(file = file.path(path.to.save.output, sprintf("new.s.obj.%s.txt", f)))
  new.commands[[f]] %>% print()
  sink()
}



#########
num.dim.integration <- 30
s.obj <- s.obj.new
DefaultAssay(s.obj) <- "RNA"
data.list <- SplitObject(s.obj, split.by = "name")

data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})

k.filter <- 200

anchors <- FindIntegrationAnchors(object.list = data.list, dims = 1:num.dim.integration, scale=F,
                                  k.filter = k.filter)## THIS IS CCA DIMENSIONS

s.obj_inte <- IntegrateData(anchorset = anchors, dims = 1:num.dim.integration, k.weight = k.filter) ## THIS IS PCA DIMENSION

## keep the order of integration obj
s.obj_inte <- s.obj_inte[, colnames(s.obj)]

s.obj[['integrated']] <- s.obj_inte[['integrated']]

s.obj@commands <- c(s.obj@commands, s.obj_inte@commands)

s.obj@tools <- c(s.obj@tools, s.obj_inte@tools)

DefaultAssay(s.obj) <- "integrated"

s.obj <- ScaleData(s.obj, verbose = FALSE, features = row.names(s.obj))

s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = FALSE, reduction.name="INTE_PCA")

count.old.scale <- GetAssayData(object = s.obj.old, slot = "scale", assay = "RNA")
count.new.scale <- GetAssayData(object = s.obj.new, slot = "scale", assay = "RNA")
sum(count.old.scale - count.new.scale)
