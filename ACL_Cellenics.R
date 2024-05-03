
# Libraries ---------------------------------------------------------------
BiocManager::install("DropletUtils")
alibrary(DropletUtils)
library(Seurat)

# Load Data ---------------------------------------------------------------

setwd("~/Desktop/UCSF/ACL/Data/")
acl <- readRDS("acl_firstfive_recluster.rds")

Seurat::DefaultAssay(acl) <- "RNA"
counts <- acl[["RNA"]]@data

rownames(counts)[1]
colnames(counts)[1]
unique(acl$Sample)


# Demultiplex Data and export as 10X files --------------------------------

setwd("~/Desktop/UCSF/ACL/Data/Demultiplexed_Samples/")

acl.list <- Seurat::SplitObject(acl, split.by = "Sample")
sample.names <- unique(acl$Sample)

head(sample.names)
tail(sample.names)

demultiplex_convert_to_10x <- function(obj, samples) {
  if(class(acl) != "Seurat") {
    message("WARNING: this rds file does not contain a Seurat object! STOP RUNNING THIS SCRIPT")
    message("Check the data type by running:")
    message("class(acl)")
    stop()
  }
  if(!dir.exists(file.path(getwd(), "demultiplexed"))) {
    dir.create(file.path(getwd(), "demultiplexed"))
  } else {
    print("WARNING! A demultiplexed directory already exists")
    return()
  }
  for (i in 1:length(samples)) {
    print(paste0("Converting sample ", samples[i]))
    obj.sub <- obj[[samples[i]]]
    DropletUtils::write10xCounts(path = paste0(getwd(),"/demultiplexed/",samples[i]), x = obj.sub[["RNA"]]@data, type = "sparse", version="3")
  }
}


demultiplex_convert_to_10x(obj = acl.list, 
                           samples = sample.names)














