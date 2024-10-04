gc()
rm(list = ls())
library(Seurat)

run_QC_on_s.obj <- function(s.obj){
  all.QC <- list()
  
  # Number of cells obtained per experiment/sample
  all.QC$cell.counts.plot <- ggplot(s.obj@meta.data, 
                                    aes(x=name , fill=name)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 14)) +
    ggtitle("Number of cells in each dataset")
  
  # distribution of number of UMI in each sample
  all.QC$nCountRNA.distribution <- ggplot(s.obj@meta.data,
                                          aes(color=name, x=nCount_RNA, fill = name)) + 
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 14)) +
    ylab("Cell density") +
    geom_vline(xintercept = 500, color = "red") +
    ggtitle("Distribution of read depths in each sample")
  
  # distribution of number of features (genes detected) in each sample
  all.QC$nFeature_RNA.distribution <- ggplot(s.obj@meta.data,
                                             aes(color=name, x=nFeature_RNA, fill = name, y = ..scaled..)) + 
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 14)) +
    ylab("Cell density") +
    geom_vline(xintercept = 1000, color = "red") +
    xlim(1000, 10000) +
    ggtitle("Distribution of number of detected genes in each sample")
  
  
  # scatter plot showing the relation between cell read-depth and number of genes detected.
  ## with Mitochondrial percentage
  all.QC$nCount.vs.nFeature.MT <- ggplot(s.obj@meta.data, 
                                         aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm, formula = y ~ x) + # apply a linear regression to show the relation, if existed.
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    facet_wrap(~name) +
    ggtitle("Scatter plot: nCount_RNA vs. nFeature_RNA, cmap % Mitochondrial genes")
  
  ## with Ribosome percentage
  all.QC$nCount.vs.nFeature.Ribo <- ggplot(s.obj@meta.data, 
                                           aes(x=nCount_RNA, y=nFeature_RNA, color=percent.ribo)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm, formula = y ~ x) + # apply a linear regression to show the relation, if existed.
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    facet_wrap(~name) +
    ggtitle("Scatter plot: nCount_RNA vs. nFeature_RNA, cmap % Ribosome genes")
  
  # Complexity: 
  
  ## We can see the samples where we sequenced each cell less have a higher overall complexity, 
  # that is because we have not started saturating the sequencing for any given gene for these samples. 
  # Outlier cells in these samples might be cells that have a less complex RNA species than other cells. 
  # Sometimes we can detect contamination with low complexity cell types like red blood cells via this metric. 
  # Generally, we expect the novelty score to be above 0.80. 
  
  ## More expanations: they are looking for cells that have a low number of genes with a high number of UMI counts. 
  # This likely means that you only captured transcripts from a low number of genes, and simply sequenced transcripts 
  # from those lower number of genes over and over again. This could be because of the cell type 
  # (such as a red blood cell having little to no RNA as they mentioned), or some other strange artifact.
  
  # Compute the complexity and add it to the s.obj@meta.data
  s.obj@meta.data <- s.obj@meta.data %>% 
    mutate(log10GenesPerUMI = log10(nFeature_RNA) / log10(nCount_RNA))
  
  all.QC$complexity <- ggplot(s.obj@meta.data,
                              aes(x=log10GenesPerUMI, color = name, fill=name)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 14)) +
    geom_vline(xintercept = 0.8) +
    ggtitle("Complexity: Log10(nCount_RNA) / log10(nFeature_RNA)")
  
  # add new slot for all.QC into the existed SEURAT OBJECT. 
  s.obj@misc$new.all.QC <- all.QC  
  return(s.obj)
}

################################################################################
path.to.s.obj <- "/home/hieunguyen/CRC1382/outdir/SAlBounny/SAlBounny_first_analysis/OUTPUT/R01_output/preprocessed_seurat_objects.rds"

s.obj.prev <- readRDS(path.to.s.obj)

s.obj.prev <- run_QC_on_s.obj(s.obj.prev)
################################################################################

outdir <- "/home/hieunguyen/CRC1382/outdir/SAlBounny"
PROJECT <- "SAlBounny_16102023_rerunQC"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

path.to.s.obj <- file.path(outdir, PROJECT, "s8_output", sprintf("%s.output.s8.rds", PROJECT))

s.obj.new <- readRDS(path.to.s.obj)
s.obj.new.noCD4 <- subset(s.obj.new, name %in% c("A_CD4", "N_CD4") == FALSE)
s.obj.new.noCD4 <- subset(s.obj.new.noCD4, nFeature_RNA >= 300 & log10GenesPerUMI >= 0.8 & nCount_RNA <= 40000)

s.obj.new.noCD4 <- run_QC_on_s.obj(s.obj.new.noCD4)

s.obj.Nd7 <- subset(s.obj.new, name == "Nd7_2")

s.obj.Nd7 <- run_QC_on_s.obj(s.obj.Nd7)
