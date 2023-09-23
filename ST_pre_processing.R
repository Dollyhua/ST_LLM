####ST data pre_processing
rm(list = ls())
library(dplyr)
library(hdf5r)
library(Seurat)
library(rhdf5)

####writeh5 function:Write the seurat object to the HDF5 file
writeh5 <- function(h5path, count_data){
  h5createFile(h5path)
  
  h5createGroup(h5path,"matrix")
  h5write(count_data@Dimnames[[2]] , h5path, "matrix/barcodes")
  h5write(count_data@x, h5path, "matrix/data")
  h5createGroup(h5path,"matrix/features")
  h5write("genome", h5path, "matrix/features/_all_tag_keys")
  Genes <- rep("Gene Expression", length(count_data@Dimnames[[1]]))
  h5write(Genes,h5path, "matrix/features/feature_type")
  Genome <- rep("GRCh38", length(count_data@Dimnames[[1]]))
  h5write(Genome, h5path,"matrix/features/genome")
  
  h5write(count_data@Dimnames[[1]],h5path, "matrix/features/id")
  h5write(count_data@Dimnames[[1]],h5path, "matrix/features/name")
  h5write(count_data@i, h5path, "matrix/indices") 
  h5write(count_data@p, h5path, "matrix/indptr")
  h5write(count_data@Dim, h5path, "matrix/shape")
  
  h5closeAll()
  
}

####ST_preprocess function
ST_preprocess <- function(ID, suffix){
    path <- list.files(paste0("/fs/home/dongzhonghua/STARDUST/Data/", study_ID, "/suppl/"))[grep(file_suffix, list.files(paste0("/fs/home/dongzhonghua/STARDUST/Data/", study_ID, "/suppl/")))]
    samples <- unlist(lapply(strsplit(path,file_suffix),function(x) x[1]))
    print(samples)

    for (sample in samples){
        print(sample)
        ####read h5 file
        h5_file <- Read10X_h5(filename = paste0("/fs/home/dongzhonghua/STARDUST/Data/", study_ID, "/suppl/", sample, file_suffix))
        h5_object <- CreateSeuratObject(counts=h5_file, project=paste0(sample))
        #####QC，Discard low quality Spots
        h5_object[["percent.mt"]] <- PercentageFeatureSet(object = h5_object, pattern = "^MT-")
        h5_object_sub <- subset(h5_object, subset = nFeature_RNA > 200 & percent.mt < 10)

        h5_object <- SCTransform(h5_object, assay = "RNA", verbose = FALSE)

        h5_object <- FindVariableFeatures(h5_object, selection.method = "vst", nfeatures = 2000)

        all.genes <- rownames(h5_object)

        h5_object <- RunPCA(h5_object,features = VariableFeatures(object=h5_object))

        h5_object <- FindNeighbors(h5_object, dims = 1:30) 
        h5_object <- FindClusters(h5_object, resolution=0.5)

        h5_object <- RunUMAP(h5_object, dims = 1: 30)

        #save.rds file 
        saveRDS(h5_object, file = paste0("/fs/home/dongzhonghua/STARDUST/Data/", study_ID, "/suppl/", sample, ".rds"))

        #######Save the ST gene count file after QC filtering for STRIDE deconvolution
        writeh5(h5path = paste0("/fs/home/dongzhonghua/STARDUST/Data/", study_ID, "/suppl/", sample, '_gene_count_QC.h5'), count_data = h5_object@assays$RNA@data)

    }
####platforms:10x Genomics
####GSE158803 for example
#####study_ID <- "GSE158803"
####file_suffix <- "_gene_count.h5"  ######The filename suffix is "_gene_count.h5";" _filtered_feature_bc_matrix.h5";" .h5" etc

####star working
#study_ID <- "GSE202740"
#file_suffix <- "_gene_count.h5"
#ST_preprocess(ID = study_ID, suffix = file_suffix)

