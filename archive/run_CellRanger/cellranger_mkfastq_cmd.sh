inputdir="/media/hieunguyen/HNSD_MBPro/raw_data/200709_NB501289_0341_AHLKC5BGXF/genomics.rwth-aachen.de/data/230830_Shahed_old_data/200709_NB501289_0341_AHLKC5BGXF";

path_to_cellranger="/home/hieunguyen/CRC1382/src/src_pipeline/CellRanger/cellranger-7.1.0"
export PATH=${path_to_cellranger}:$PATH
cellranger mkfastq --run ${inputdir} --sample-sheet ${inputdir}/SampleSheet.csv;
