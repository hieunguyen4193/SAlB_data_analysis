gc()
rm(list = ls())

# remove.packages("DOSE")
# remove.packages("GOSemSim")
# remove.packages("yulab.utils")
# remove.packages("clusterProfiler")
# devtools::install_github("YuLab-SMU/yulab.utils", upgrade = "never")
# remotes::install_github("GuangchuangYu/GOSemSim", upgrade = "never")
# remotes::install_github("GuangchuangYu/clusterProfiler", upgrade = "never")
# install.packages("/home/hieunguyen/CRC1382/storage/offline_pkgs/org.Mm.eg.db_3.18.0.tar.gz", type = "sources", repos = NULL)
# install.packages("/home/hieunguyen/CRC1382/storage/offline_pkgs/org.Mm.eg.db_3.18.0.tar.gz", type = "sources", repos = NULL)
# BiocManager::install("org.Mm.eg.db", update = FALSE )
# install.packages("heatmaply")

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

library(scales)
library(scatterpie)

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "SAlBounny_full"
output.version <- "20241021"
s.obj.name <- "SAlBounny_full.filter_contaminated_cells.clusterRes_0.5.subClusterFoxp3_cluster6.rds"

path.to.main.output <- file.path(outdir, PROJECT, output.version, "data_analysis")
path.to.html.output <- file.path(path.to.main.output, "html_output", "05_output")
dir.create(path.to.html.output, showWarnings = FALSE, recursive = TRUE)
path.to.rmd <- "/home/hieunguyen/CRC1382/src_2023/SAlB_data_analysis/05_DEG.Rmd"

pair.DEA <- list(
  `pair1` = list(
    sample1 = "Ad7_2",
    sample2 = "A_CD4"
  ),
  `pair2` = list(
    sample1 = "N_48h",
    sample2 = "A_48h"
  ),
  `pair3` = list(
    sample1 = "Nd7_2",
    sample2 = "Ad7_2"
  ),
  `pair4` = list(
    sample1 = "Nd7_2",
    sample2 = "N_CD4"
  )
)

for (i in names(pair.DEA)){
  sample1 <- pair.DEA[[i]]$sample1
  sample2 <- pair.DEA[[i]]$sample2
  html.filename <- sprintf("DEA_%s_vs_%s.html", sample1, sample2)
  rmarkdown::render(input = path.to.rmd,
                    params = list(
                      sample1 = sample1,
                      sample2 = sample2
                    ),
                    output_dir = path.to.html.output,
                    output_file = html.filename)
}
