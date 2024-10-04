gc()
rm(list = ls())

# remove.packages("DOSE")
# remove.packages("GOSemSim")
# remove.packages("yulab.utils")
# remove.packages("clusterProfiler")
# 
# remotes::install_github("GuangchuangYu/GOSemSim", upgrade = "never")
# remotes::install_github("GuangchuangYu/clusterProfiler", upgrade = "never")

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/CRC1382_SAlBounny_projects/CRC1382_SAlBounny_fulldata"

sample.pairs <- list(`1` = c("N_48h", "A_48h"),
                     `2` = c("Nd7_2", "Ad7_2"),
                     `3` = c("Ad7_2", "A_CD4"),
                     `3` = c("Nd7_2", "N_CD4")
)
path.to.main.output <- "/media/hieunguyen/HD0/outdir/CRC1382/20231018_SAlBounny/data_analysis"

for (pair in sample.pairs){
  input.sample1 <- pair[[1]]
  input.sample2 <- pair[[2]]
  path.to.rmd.file <- file.path(path.to.main.src, "05_DEG.Rmd")
  save.html.name <- sprintf("DEA_%s_vs_%s.html", input.sample1, input.sample2)
  path.to.save.html <- file.path(path.to.main.output, "HTML_reports")
  if (file.exists(file.path(path.to.save.html, save.html.name)) == FALSE){
    rmarkdown::render(input = path.to.rmd.file,
                      params = list(sample1 = input.sample1, 
                                    sample2 = input.sample2),
                      output_file = save.html.name,
                      output_dir = path.to.save.html)
  }
}
