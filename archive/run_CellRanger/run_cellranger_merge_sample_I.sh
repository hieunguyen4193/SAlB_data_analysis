path_to_cellranger="/home/hieunguyen/CRC1382/src/src_pipeline/CellRanger/cellranger-7.1.0"
export PATH=${path_to_cellranger}:$PATH

path_to_fastq="/media/hieunguyen/HNSD_MBPro/raw_data/200403_200709_merged_sample_I";
path_to_outputdir=${path_to_fastq}/cellranger_outputs;

for sample in I;do \
cellranger count --id=SAlBounny_merged_I_${sample} \
                   --transcriptome=/media/hieunguyen/HD0/ext_HDD/raw_data/build-ref-genome-10x/build-mm10/mm10 \
                   --fastqs=${path_to_fastq}/${sample} \
                   --sample=${sample} \
                   --localcores=20;

mkdir -p ${path_to_outputdir}/${sample};

rsync -avh --progress ./SAlBounny_merged_I_${sample} ${path_to_outputdir}/${sample};

touch SAlBounny_merged_I_${sample}.finished.txt;
rm -rf ./SAlBounny_merged_I_${sample};done;
