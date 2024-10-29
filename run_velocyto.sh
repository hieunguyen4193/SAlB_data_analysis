path="/hpcwork/uk104163/storage/tmp";
path_to_masked_regions_gtf="/hpcwork/uk104163/storage/mm10_rmsk.gtf";
path_to_modified_gtf="/hpcwork/uk104163/storage/build-mm10/mm10-2020-A_build/gencode.vM23.primary_assembly.annotation.gtf.filtered";
for sample in SAlBounny_200709_A_CD4 SAlBounny_200709_N_CD4 SAlBounny_merged_I;do \
samtools_threads=20;
path_to_cellranger_output=/hpcwork/uk104163/storage/tmp/${sample};
velocyto run10x -m ${path_to_masked_regions_gtf} ${path_to_cellranger_output} ${path_to_modified_gtf} -@ ${samtools_threads};done
