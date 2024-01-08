# Load packages
library(dplyr)
library(Seurat)
library(Matrix)
library(png)
library(tiff)
library(jpeg)
library(jsonlite)

expr_preprocess <- function(resDir, sample, samplePath, genome){
    # Preprocess spatial gene expression
    #
    # Args:
    #   resDir: Path to store results.
    #   sample: Name of sample.
    #   samplePath: Path to the HDF5 format gene expression file of sample.
    #   genome: Genome version. Choose from: 'GRCh38', 'GRCh37', 'GRCm38', 'GRCm37'.
    #
    # Return:
    #   QC-filtered .rds file.
    #   QC-filtered .h5 file.

    print("Preprocessing ST gene expression...")
    # Read in h5 data
    h5_file <- Read10X_h5(filename = samplePath)
    h5_object <- CreateSeuratObject(counts = h5_file, project = sample)
    
    # Filter out low quality spots
    h5_object[["percent.mt"]] <- PercentageFeatureSet(object = h5_object, pattern = "^MT-")
    h5_object <- subset(h5_object, subset = nFeature_RNA > 200 & percent.mt < 20)
    h5_object <- SCTransform(h5_object, assay = "RNA", verbose = FALSE)
    h5_object <- FindVariableFeatures(h5_object, selection.method = "vst", nfeatures = 2000)
    h5_object <- RunPCA(h5_object, features = VariableFeatures(object = h5_object))
    h5_object <- FindNeighbors(h5_object, dims = 1:30)
    h5_object <- FindClusters(h5_object, resolution = 0.5)
    h5_object <- RunUMAP(h5_object, dims = 1:30)

    # Save the .rds file for future analysis.
    saveRDS(h5_object, file = file.path(resDir, paste0(sample, "_filtered.rds")))

    # Save the QC-filtered .h5 file for STRIDE deconvolution.
    writeh5(h5path = file.path(resDir, paste0(sample, "_gene_count_QC.h5")), 
            count_data = h5_object@assays$RNA@data, 
            genome_assembly = genome)

    return(h5_object)
    
}


# To import spatial image information, create a folder named "sample_spatial."
# This folder should contain both images and spot location information.
# Acceptable image formats include tissue_hires_image.png, tissue_lowres_image.png, sample_HE.tif, sample_HE.png and sample_HE.jpeg(sample_HE.jpg).
# Acceptable location information formats include tissue_positions_list.csv and sample_spot_location.txt.
# Additionally, you can also read in scalefactors_json.json for additional image information if available.

Read_image <- function(sample, spatial_dir){
    # Read in spatial image and location
    # 
    # Args:
    #   sample: Name of the sample.
    #   saptial_dir: Path to directory containing images and spot location information.
    #
    # Return:
    #   A 'VisiumV1' class object with image, scale factors, coordinates, and spot radius.

    print("Reading image...")
    spatial_files <- list.files(spatial_dir)

    # Read in image
    if ("tissue_hires_image.png" %in% spatial_files) {
        image.path <- file.path(spatial_dir, "tissue_hires_image.png")
        img <- readPNG(image.path)
    } else if ("tissue_lowres_image.png" %in% spatial_files) {
        image.path <- file.path(spatial_dir, "tissue_lowres_image.png")
        img <- readPNG(image.path)
    } else if (paste0(sample,"_HE.tif") %in% spatial_files) {
        image.path <- file.path(saptial_dir, paste0(sample,"_HE.tif"))
        img <- readTIFF(image.path)
    } else if(paste0(sample,"_HE.png") %in% spatial_files) {
        image.path <- file.path(saptial_dir, paste0(sample,"_HE.png"))
        img <- readPNG(image.path)
    } else if(paste0(sample,"_HE.jpeg") %in% spatial_files) {
        image.path <- file.path(saptial_dir, paste0(sample,"_HE.jpeg"))
        img <- readJPEG(image.path)
    } else if(paste0(sample,"_HE.jpg") %in% spatial_files) {
        image.path <- file.path(saptial_dir, paste0(sample,"_HE.jpg"))
        img <- readJPEG(image.path)
    } else {
        image.path <- NA
        img <- NA
        warning("No acceptable image found. Please check the directory for image files or rename them to acceptable formats.")
    }

    # Read in scalefactors_json.json
    if ("scalefactors_json.json" %in% spatial_files){
        scale.factors <- fromJSON(txt = file.path(spatial_dir, 'scalefactors_json.json')) #library(jsonlite)  
        if (basename(image.path) == "tissue_hires_image.png"){
            unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_hires_scalef
            spot.radius <-  unnormalized.radius / max(dim(x = img))
        } else if (basename(image.path) == "tissue_lowres_image.png"){
            unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
            spot.radius <-  unnormalized.radius / max(dim(x = img))
        } else{
            unnormalized.radius <- NA
            spot.radius <- NA
        }
        scale.factors = scalefactors(
        spot = scale.factors$spot_diameter_fullres,
        fiducial = scale.factors$fiducial_diameter_fullres,
        hires = scale.factors$tissue_hires_scalef,
        lowres = scale.factors$tissue_lowres_scalef
        )
    } else{
        scale.factors <- scalefactors(
        spot = NA,
        fiducial = NA,
        hires = NA, 
        lowres = NA
        )
        unnormalized.radius <- NA
        spot.radius <- NA
    }

    # Read in tissue positions
    if ("tissue_positions_list.csv" %in% spatial_files){
        tissue_positions_df <- read.csv(file.path(spatial_dir, "tissue_positions_list.csv"), 
                                    col.names = c('barcodes', 'tissue', 'row', 'col', 'imagerow', 'imagecol'),
                                    quote = "", header = FALSE, row.names = 1)
    } else if(paste0(sample, "_spot_location.txt") %in% spatial_files){
        tissue_positions_df <- read.table(file.path(spatial_dir, paste0(sample, "_spot_location.txt")), 
                                        sep = "\t", quote = "", header = TRUE, 
                                        col.names = c("barcodes","row","col"), row.names = 1)
    } else{
        tissue_positions_df <- NA
        warning("No acceptable tissue positions file found. Please check the directory for positions files or rename them to acceptable formats.")
    }
        
    # Return a VisiumV1 object with image, scale factors, coordinates, and spot radius
    return(new(
    Class = 'VisiumV1',
    image = img,
    scale.factors = scale.factors,
    coordinates = tissue_positions_df,
    spot.radius = spot.radius
    ))
}

  

ST_preprocess <- function(outDir, study_ID, file_suffix, genome) {
  # Preprocess ST data
  #
  # Args:
  #   outDir: Diretory for data saving. Choose from: '/fs/home/dongzhonghua/STARDUST/Data/' or '/fs/home/dongzhonghua/ST_LLM/Data/' or '/fs/home/shixiaoying/Project/ST_LLM/Data/'.
  #   study_ID: GSE number. 
  #   file_suffix: Suffix of files. Example: '_gene_count.h5', '_filtered_feature_bc_matrix.h5', '.h5' and so on.
  #   genome: Genome version. Choose from: 'GRCh38', 'GRCh37', 'GRCm38', 'GRCm37'.
  # 
  # Return:
  #   QC-filtered .rds file.
  #   QC-filtered .h5 file.
  
  resDir <- file.path(outDir, study_ID, "suppl")
  allh5 <- list.files(resDir, pattern = file_suffix, recursive = T, full.names = T)
  samples <- gsub(file_suffix, '', basename(allh5))
  print(samples)
  
  for (sample in samples) {
    print(sample)
    samplePath = allh5[grep(sample, allh5)]
    
    # Perform quality control and basic analysis of spatial gene expression.
    ST_object <- expr_preprocess(resDir = resDir, sample = sample, samplePath = samplePath, genome = genome)

    # Add image information to ST_object.
    ST_object <- readRDS(file.path(resDir, paste0(sample, "_filtered.rds"))) 
    spatial_dir <- file.path(resDir, paste0(sample,"_spatial"))
    image <- Read_image(sample = sample, spatial_dir = spatial_dir)
    DefaultAssay(object = image) <- 'RNA'
    image <- image[colnames(x = ST_object)]
    ST_object[[sample]] <- image

    # Save the .rds file for future analysis.
    saveRDS(ST_object, file = file.path(resDir, paste0(sample, "_filtered.rds")))
  }
}

writeh5 <- function(h5path, count_data, genome_assembly){
  library(rhdf5)
  library(hdf5r)
  
  h5createFile(h5path)
  h5createGroup(h5path,"matrix")
  h5write(count_data@Dimnames[[2]] , h5path, "matrix/barcodes")
  h5write(as.integer(count_data@x), h5path, "matrix/data")
  h5createGroup(h5path,"matrix/features")
  h5write("genome", h5path, "matrix/features/_all_tag_keys")
  Genes <- rep("Gene Expression", length(count_data@Dimnames[[1]]))
  h5write(Genes,h5path, "matrix/features/feature_type")
  Genome <- rep(genome_assembly, length(count_data@Dimnames[[1]])) 
  h5write(Genome, h5path,"matrix/features/genome")
  h5write(count_data@Dimnames[[1]],h5path, "matrix/features/id")
  h5write(count_data@Dimnames[[1]],h5path, "matrix/features/name")
  h5write(count_data@i, h5path, "matrix/indices") 
  h5write(count_data@p, h5path, "matrix/indptr")
  h5write(count_data@Dim, h5path, "matrix/shape")
  h5closeAll()
}

# Usage
outDir <- "/fs/home/dongzhonghua/STARDUST/Data"
study_ID <- "GSE210041"
file_suffix <- "_filtered_feature_bc_matrix.h5"
genome <- "GRCh38"
ST_preprocess(outDir, study_ID, file_suffix, genome)
