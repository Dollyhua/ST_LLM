# ST
step1. Download ST data, both raw and metadata. You can download it from GEO_ID_download.sh or GEO_list_download.sh

step2. The unified processing format is h5 file, if it is "barcodes.tsv", "genes.tsv" and "matrix.mtx", convert it with MAESTRO mtx-to-h5, if it is "counts.tsv", "mat.txt", "Counts. CSV", use MAESTRO count-to-h5 conversion, see https://github.com/liulab-dfci/MAESTRO, or on the command line, enter MAESTRO count-to-h5 --help to ask for help.

step3. ST data preprocessing: After being unified into h5 format, R language "SEURAT"package is used for processing, the specific code is in ST_pre_processing.R

step4.
