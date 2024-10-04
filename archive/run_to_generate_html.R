gc()
rm(list = ls())

set.seed(42)

path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline"

source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))
source(file.path(path.to.pipeline.src, "processes_src", "s9_findMarkers.R"))
outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "20231018_SAlBounny"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output")
path.to.03.output <- file.path(path.to.main.output, "03_output")

path.to.save.html <- file.path(path.to.main.output, "HTML_reports")
dir.create(path.to.save.html, recursive = TRUE, showWarnings = FALSE)

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/CRC1382_SAlBounny_projects/CRC1382_SAlBounny_fulldata"
  
for (cluster.resolution in c(0.3, 0.5, 0.8, 1, 1.25, 1.5)){
  path.to.Rmd.file <- file.path(path.to.project.src, "03_remove_contaminated_cells.Rmd")
  html_name <- sprintf(str_replace(basename(path.to.Rmd.file), "Rmd", sprintf("clsuterRes_%s.html", cluster.resolution)))
  rmarkdown::render(input = path.to.Rmd.file,
                    output_file = html_name,
                    output_dir = path.to.save.html,
                    params = list(cluster.resolution = cluster.resolution))   
}
