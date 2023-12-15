# Workflow to preprocess ST data

## 1. Collect and download ST data 
We usually download ST data from GEO database. Related codes can be found at `GEO_ID_download.sh` and `GEO_list_download.sh`.

## 2. Standardize file format for batch processing
To run this pipeline aotumatically, we need standardize file format.
### 2.1 Standardize ST gene expression file
ST expression data may have multiple file formats, such as HDF5, plain text, and mtx format. This pipeline mainly uses HDF5 format. Expression files in other formats can be converted into HDF5 files through MAESTRO.

The following commands can be used for format conversion.

    MAESTRO mtx-to-h5	      #Convert mtx format to HDF5 format
    MAESTRO count-to-h5       #Convert plain text format to HDF5 format

See https://github.com/liulab-dfci/MAESTRO for more information, or enter 'MAESTRO count-to-h5 --help' on the command line to ask for help.
### 2.2 Standardize spatial image file


## 3. Preprocess ST data
R language "Seurat" package is used to preprocess ST data, the specific code is in `ST_pre_processing.R`.

