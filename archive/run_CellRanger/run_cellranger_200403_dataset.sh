path_to_cellranger="/home/hieunguyen/CRC1382/src/src_pipeline/CellRanger/cellranger-7.1.0"
export PATH=${path_to_cellranger}:$PATH

path_to_fastq="/media/hieunguyen/HNSD_MBPro/raw_data/200403_NB501289_0317_AHKKFHBGXC/genomics.rwth-aachen.de/data/230830_Shahed_old_data/200403_NB501289_0317_AHKKFHBGXC/compressed_tars/NT200402/NT200402_";
path_to_outputdir="/media/hieunguyen/HNSD_MBPro/raw_data/200403_NB501289_0317_AHKKFHBGXC/genomics.rwth-aachen.de/data/230830_Shahed_old_data/200403_NB501289_0317_AHKKFHBGXC/cellranger_outputs";

for sample in Ad7 I Nd7;do \
cellranger count --id=SAlBounny_200403_${sample} \
                   --transcriptome=/media/hieunguyen/HD0/ext_HDD/raw_data/build-ref-genome-10x/build-mm10/mm10 \
                   --fastqs=${path_to_fastq}/${sample} \
                   --sample=${sample} \
                   --localcores=20;

mkdir -p ${path_to_outputdir}/${sample};

rsync -avh --progress ./SAlBounny_200403_${sample} ${path_to_outputdir}/${sample};

touch SAlBounny_200403_${sample}.finished.txt;
rm -rf ./SAlBounny_200403_${sample};done;
