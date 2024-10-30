library("Seurat")

s.obj.old <- readRDS("/home/hieunguyen/CRC1382/outdir/SAlBounny_full/20241021/data_analysis/03_output/SAlBounny_full.filter_contaminated_cells.clusterRes_0.5.subClusterFoxp3_cluster6.rds")
s.obj.new <- readRDS("/home/hieunguyen/CRC1382/outdir/SAlBounny_full/20241021/data_analysis/06_output/SAlBounny_full.filter_contaminated_cells.clusterRes_0.5_subcluster_gene_Izumo1r.rds")

counts_matrix.old <- GetAssayData(s.obj.old, assay='RNA', slot='counts')
counts_matrix.new <- GetAssayData(s.obj.new, assay='RNA', slot='counts')
