suppressMessages(
  {
    library(dplyr)
    library(Seurat)
    library(patchwork)
    library(tibble)
    library(DoubletFinder)
    library(scCustomize)
    library(ggplot2)
    #library(loupeR)
    library(qs)
  })

torn_acl <- readRDS("/Volumes/T7 Shield/UCSF/ACL/Objects/ACL_down_integrated.rds")
torn_acl$Gender <- factor(torn_acl$Gender)
torn_acl$group <- "torn"
torn_acl$Sample <- factor(torn_acl$Sample)

torn_acl <- SCTransform(torn_acl)
# These are now standard steps in the Seurat workflow for visualization and clustering
torn_acl <- RunPCA(torn_acl, verbose = FALSE)
torn_acl <- RunUMAP(torn_acl, dims = 1:30, verbose = FALSE)

torn_acl <- FindNeighbors(torn_acl, dims = 1:30, verbose = FALSE)
torn_acl <- FindClusters(torn_acl, verbose = FALSE)
DimPlot(torn_acl, label = TRUE)



saveRDS(torn_acl,"/Volumes/T7 Shield/UCSF/ACL/Objects/ACL_sct_down_integrated.rds")



# L1 ----------------------------------------------------------------------
l1_data <- Read10X("D:/10X_CLLI/txg-macos-v3.0.0/healthy_acl/HRR865576/filtered_feature_bc_matrix/")
l1 <- CreateSeuratObject(counts = l1_data)

obj <- l1
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj)
l1 <- obj

l1$Gender <- "Female"
l1$orig.ident <- "HRR865576"
l1$Sample <- "L1"
l1$group <- "healthy"
#saveRDS(l1, "D:/10X_CLLI/txg-macos-v3.0.0/healthy_acl/HRR865576/l1.rds")


# L3 ----------------------------------------------------------------------
l3_data <- Read10X("D:/10X_CLLI/txg-macos-v3.0.0/healthy_acl/HRR865577/filtered_feature_bc_matrix/")
l3 <- CreateSeuratObject(counts = l3_data)

obj <- l3
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj)
l3 <- obj

l3$Gender <- "Male"
l3$orig.ident <- "HRR865577"
l3$Sample <- "L3"
l3$group <- "healthy"
#saveRDS(l3, "/Volumes/T7 Shield/10X_CLLI/txg-macos-v3.0.0/healthy_acl/HRR865577/l3.rds")


# L4 ----------------------------------------------------------------------
l4_data <- Read10X("D:/10X_CLLI/txg-macos-v3.0.0/healthy_acl/HRR865578/filtered_feature_bc_matrix/")
l4 <- CreateSeuratObject(counts = l4_data)
obj <- l4
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj)
l4 <- obj

l4$Gender <- "Female"
l4$orig.ident <- "HRR865578"
l4$Sample <- "L4"
l4$group <- "healthy"
#saveRDS(l4, "/Volumes/T7 Shield/10X_CLLI/txg-macos-v3.0.0/healthy_acl/HRR865578/l4.rds")

# L10 ---------------------------------------------------------------------
l10_data <- Read10X("D:/10X_CLLI/txg-macos-v3.0.0/healthy_acl/HRR865579/filtered_feature_bc_matrix/")
l10 <- CreateSeuratObject(counts = l10_data)
obj <- l10
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj)
l10 <- obj
l10$Gender <- "Male"
l10$orig.ident <- "HRR865579"
l10$Sample <- "L10"
l10$group <- "healthy"
#saveRDS(l10, "/Volumes/T7 Shield/10X_CLLI/txg-macos-v3.0.0/healthy_acl/HRR865579/l10.rds")




# Read Data ---------------------------------------------------------------
torn_acl <- readRDS("D:/UCSF/ACL/Objects/ACL_sct_down_integrated.rds")
#torn_acl$Gender <- factor(torn_acl$Gender)
#torn_acl$group <- "torn"
#torn_acl$Sample <- factor(torn_acl$Sample)

#l1 <- readRDS("D:/10X_CLLI/txg-macos-v3.0.0/healthy_acl/HRR865576/l1.rds")
#l3 <- readRDS("D:/10X_CLLI/txg-macos-v3.0.0/healthy_acl/HRR865577/l3.rds")
#l4 <- readRDS("D:/10X_CLLI/txg-macos-v3.0.0/healthy_acl/HRR865578/l4.rds")
#l10 <- readRDS("D:/10X_CLLI/txg-macos-v3.0.0/healthy_acl/HRR865579/l10.rds")

# Merge/ Integrate --------------------------------------------------------
DefaultAssay(torn_acl) <- "RNA"
torn_acl <- JoinLayers(torn_acl)
torn_acl <- NormalizeData(torn_acl)
torn_acl <- FindVariableFeatures(torn_acl)
torn_acl <- ScaleData(torn_acl)
torn_acl <- RunPCA(torn_acl)
torn_acl <- FindNeighbors(torn_acl, dims = 1:30)
torn_acl <- FindClusters(torn_acl)

acl <- merge(torn_acl, c(l1,l3,l4,l10),add.cell.ids = c("torn",'l1','l3','l4','l10'),merge.data = F)
rm(l1)
rm(l1_data)
rm(l3)
rm(l3_data)
rm(l4)
rm(l4_data)
rm(l10)
rm(l10_data)
rm(obj)
rm(torn_acl)
gc()

DefaultAssay(acl) <- "RNA"
acl <- JoinLayers(acl, assay = "RNA")
acl[["RNA"]] <- split(acl[["RNA"]], f = acl$group)

acl <- NormalizeData(acl)
acl <- FindVariableFeatures(acl)
acl <- ScaleData(acl)
acl <- RunPCA(acl)
acl <- FindNeighbors(acl, dims = 1:30)
acl <- FindClusters(acl)

acl <- IntegrateLayers(object = acl, method = CCAIntegration, orig.reduction = "pca",
                                new.reduction = "integrated.cca", verbose = FALSE)
acl <- FindNeighbors(acl, reduction = "integrated.cca", dims = 1:30)
acl <- FindClusters(acl)
acl <- RunUMAP(acl, reduction = "integrated.cca", dims = 1:30)
DimPlot(acl, label = T)

#saveRDS(acl, "D:/UCSF/ACL/Objects/ACL_integrated_torn_healthy.rds")

#acl <- readRDS("D:/UCSF/ACL/Objects/ACL_integrated_torn_healthy.rds")


# all_cells - label transfer ----------------------------------------------------------
Idents(acl) <- "group"
acl.ref <- subset(acl, idents = "torn")
acl.query <- acl #subset(acl, idents = "healthy")


acl.query <- NormalizeData(acl.query)
acl.anchors <- FindTransferAnchors(reference = acl.ref, query = acl.query, dims = 1:30,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = acl.anchors, refdata = acl.ref$celltype, dims = 1:30)
acl.query <- AddMetaData(acl.query, metadata = predictions)


acl.query$prediction.match <- acl.query$predicted.id == acl.query$celltype
table(acl.query$prediction.match)
table(acl.query$predicted.id)

acl.ref <- RunUMAP(acl.ref, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
acl.query <- MapQuery(anchorset = acl.anchors, reference = acl.ref, query = acl.query,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(acl.ref, reduction = "umap", group.by = "celltype", cols = c("#f09234","#54ade0","#7f4d9f","#bc3930","#458a40")) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(acl.query, reduction = "ref.umap", group.by = "predicted.celltype", cols = c("#f09234","#54ade0","#7f4d9f","#bc3930","#458a40")) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

DimPlot(acl.query, group.by = "predicted.id", split.by = "group",cols = c("#f09234","#54ade0","#7f4d9f","#bc3930","#458a40"))

FeaturePlot(acl.query, feature = c("PDGFRA","TPPP3","PECAM1","CD34","RGS5","PTPRC","MRC1","CD3E"), cols = c("#b0e0e6","#8b0000"), ncol = 4, split.by = "group")


# Fibroblast - label transffer --------------------------------------------------------------
torn_fibro <- readRDS("D:/UCSF/ACL/Objects/acl_fibro.rds")
torn_fibro$group <- "torn"

Idents(acl) <- "seurat_clusters"
fibro <- subset(acl, idents = c("4",'11','5','12','10','1','8','0','20'))
Idents(fibro) <- "group"
healthy_fibro <- subset(fibro, idents = "healthy")

fibro <- merge(torn_fibro, y = healthy_fibro, add.cell.ids = c("torn","healthy"))
DefaultAssay(fibro) <- "RNA"
fibro <- JoinLayers(fibro, assay = "RNA")
fibro[["RNA"]] <- split(fibro[["RNA"]], f = fibro$group)

fibro <- NormalizeData(fibro)
fibro <- FindVariableFeatures(fibro)
fibro <- ScaleData(fibro)
fibro <- RunPCA(fibro)
fibro <- FindNeighbors(fibro, dims = 1:30)
fibro <- FindClusters(fibro)

fibro <- IntegrateLayers(object = fibro, method = CCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.cca", verbose = FALSE)
fibro <- FindNeighbors(fibro, reduction = "integrated.cca", dims = 1:30)
fibro <- FindClusters(fibro)
fibro <- RunUMAP(fibro, reduction = "integrated.cca", dims = 1:30)
DimPlot(fibro, label = T)



Idents(fibro) <- "group"
fibro.ref <- subset(fibro, idents = "torn")
fibro.query <- fibro #subset(fibro, idents = "healthy")

fibro.query <- NormalizeData(fibro.query)
fibro.anchors <- FindTransferAnchors(reference = fibro.ref, query = fibro.query, dims = 1:30,
                                   reference.reduction = "pca")
predictions <- TransferData(anchorset = fibro.anchors, refdata = fibro.ref$fib_cells, dims = 1:30)
fibro.query <- AddMetaData(fibro.query, metadata = predictions)


fibro.query$prediction.match <- fibro.query$predicted.id == fibro.query$fib_cells
table(fibro.query$prediction.match)
table(fibro.query$predicted.id)

fibro.ref <- RunUMAP(fibro.ref, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
fibro.query <- MapQuery(anchorset = fibro.anchors, reference = fibro.ref, query = fibro.query,
                      refdata = list(celltype = "fib_cells"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(fibro.ref, reduction = "umap", group.by = "fib_cells", cols = c("salmon1","plum","lightblue","lightpink","lightgreen")) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(fibro.query, reduction = "ref.umap", group.by = "predicted.celltype", cols = c("salmon1","plum","lightblue","lightpink","lightgreen")) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

DimPlot(fibro.query, group.by = "predicted.id", split.by = "group",cols = c("salmon1","plum","lightblue","lightpink","lightgreen"))

FeaturePlot(fibro.query, feature = c('ACAN','CD55','CHI3L2','NOTCH3','TPPP3'), cols = c("#b0e0e6","#8b0000"), ncol = 3, split.by = "group")


# Male --------------------------------------------------------------------
# all_cells - label transfer ----------------------------------------------------------
Idents(acl) <- "Gender"
male_acl <- subset(acl, idents = "Male")
Idents(male_acl) <- "group"
male_acl.ref <- subset(male_acl, idents = "torn")
male_acl.query <- male_acl #subset(male_acl, idents = "healthy")


male_acl.query <- NormalizeData(male_acl.query)
male_acl.anchors <- FindTransferAnchors(reference = male_acl.ref, query = male_acl.query, dims = 1:30,
                                   reference.reduction = "pca")
predictions <- TransferData(anchorset = male_acl.anchors, refdata = male_acl.ref$celltype, dims = 1:30)
male_acl.query <- AddMetaData(male_acl.query, metadata = predictions)


male_acl.query$prediction.match <- male_acl.query$predicted.id == male_acl.query$celltype
table(male_acl.query$prediction.match)
table(male_acl.query$predicted.id)

male_acl.ref <- RunUMAP(male_acl.ref, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
male_acl.query <- MapQuery(anchorset = male_acl.anchors, reference = male_acl.ref, query = male_acl.query,
                      refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(male_acl.ref, reduction = "umap", group.by = "celltype", cols = c("#f09234","#54ade0","#7f4d9f","#bc3930","#458a40")) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(male_acl.query, reduction = "ref.umap", group.by = "predicted.celltype", cols = c("#f09234","#54ade0","#7f4d9f","#bc3930","#458a40")) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

DimPlot(male_acl.query, group.by = "predicted.id", split.by = "group",cols = c("#f09234","#54ade0","#7f4d9f","#bc3930","#458a40"))

FeaturePlot(male_acl, feature = c("PDGFRA","TPPP3","PECAM1","CD34","RGS5","PTPRC","MRC1","CD3E"), cols = c("#b0e0e6","#8b0000"), ncol = 4,split.by = "group")

# Fibroblast - label transffer --------------------------------------------------------------
Idents(fibro) <- "Gender"
male_fibro <- subset(fibro, idents = "Male")
Idents(male_fibro) <- "group"
male_fibro.ref <- subset(male_fibro, idents = "torn")
male_fibro.query <- male_fibro #subset(male_fibro, idents = "healthy")

male_fibro.query <- NormalizeData(male_fibro.query)
male_fibro.anchors <- FindTransferAnchors(reference = male_fibro.ref, query = male_fibro.query, dims = 1:30,
                                     reference.reduction = "pca")
predictions <- TransferData(anchorset = male_fibro.anchors, refdata = male_fibro.ref$fib_cells, dims = 1:30)
male_fibro.query <- AddMetaData(male_fibro.query, metadata = predictions)


male_fibro.query$prediction.match <- male_fibro.query$predicted.id == male_fibro.query$fib_cells
table(male_fibro.query$prediction.match)
table(male_fibro.query$predicted.id)

male_fibro.ref <- RunUMAP(male_fibro.ref, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
male_fibro.query <- MapQuery(anchorset = male_fibro.anchors, reference = male_fibro.ref, query = male_fibro.query,
                        refdata = list(celltype = "fib_cells"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(male_fibro.ref, reduction = "umap", group.by = "fib_cells", cols = c("salmon1","plum","lightblue","lightpink","lightgreen")) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(male_fibro.query, reduction = "ref.umap", group.by = "predicted.celltype", cols = c("salmon1","plum","lightblue","lightpink","lightgreen")) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

DimPlot(male_fibro.query, group.by = "predicted.id", split.by = "group",cols = c("salmon1","plum","lightblue","lightpink","lightgreen"))

FeaturePlot(male_fibro, feature = c('ACAN','CD55','CHI3L2','NOTCH3','TPPP3'), cols = c("#b0e0e6","#8b0000"), ncol = 3, split.by = "group")

# FEMALE --------------------------------------------------------------------
# all_cells - label transfer ----------------------------------------------------------
Idents(acl) <- "Gender"
female_acl <- subset(acl, idents = "Female")
Idents(female_acl) <- "group"
female_acl.ref <- subset(female_acl, idents = "torn")
female_acl.query <- female_acl #subset(female_acl, idents = "healthy")


female_acl.query <- NormalizeData(female_acl.query)
female_acl.anchors <- FindTransferAnchors(reference = female_acl.ref, query = female_acl.query, dims = 1:30,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = female_acl.anchors, refdata = female_acl.ref$celltype, dims = 1:30)
female_acl.query <- AddMetaData(female_acl.query, metadata = predictions)


female_acl.query$prediction.match <- female_acl.query$predicted.id == female_acl.query$celltype
table(female_acl.query$prediction.match)
table(female_acl.query$predicted.id)

female_acl.ref <- RunUMAP(female_acl.ref, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
female_acl.query <- MapQuery(anchorset = female_acl.anchors, reference = female_acl.ref, query = female_acl.query,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(female_acl.ref, reduction = "umap", group.by = "celltype", cols = c("#f09234","#54ade0","#7f4d9f","#bc3930","#458a40")) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(female_acl.query, reduction = "ref.umap", group.by = "predicted.celltype", cols = c("#f09234","#54ade0","#7f4d9f","#bc3930","#458a40")) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

DimPlot(female_acl.query, group.by = "predicted.id", split.by = "group",cols = c("#f09234","#54ade0","#7f4d9f","#bc3930","#458a40"))

FeaturePlot(female_acl, feature = c("PDGFRA","TPPP3","PECAM1","CD34","RGS5","PTPRC","MRC1","CD3E"), cols = c("#b0e0e6","#8b0000"), ncol = 4, split.by = "group")

# Fibroblast - label transffer --------------------------------------------------------------
Idents(fibro) <- "Gender"
female_fibro <- subset(fibro, idents = "Female")
Idents(female_fibro) <- "group"
female_fibro.ref <- subset(female_fibro, idents = "torn")
female_fibro.query <- female_fibro #subset(female_fibro, idents = "healthy")

female_fibro.query <- NormalizeData(female_fibro.query)
female_fibro.anchors <- FindTransferAnchors(reference = female_fibro.ref, query = female_fibro.query, dims = 1:30,
                                          reference.reduction = "pca")
predictions <- TransferData(anchorset = female_fibro.anchors, refdata = female_fibro.ref$fib_cells, dims = 1:30)
female_fibro.query <- AddMetaData(female_fibro.query, metadata = predictions)


female_fibro.query$prediction.match <- female_fibro.query$predicted.id == female_fibro.query$fib_cells
table(female_fibro.query$prediction.match)
table(female_fibro.query$predicted.id)

female_fibro.ref <- RunUMAP(female_fibro.ref, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
female_fibro.query <- MapQuery(anchorset = female_fibro.anchors, reference = female_fibro.ref, query = female_fibro.query,
                             refdata = list(celltype = "fib_cells"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(female_fibro.ref, reduction = "umap", group.by = "fib_cells", cols = c("salmon1","plum","lightblue","lightpink","lightgreen")) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(female_fibro.query, reduction = "ref.umap", group.by = "predicted.celltype", cols = c("salmon1","plum","lightblue","lightpink","lightgreen")) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

DimPlot(female_fibro.query, group.by = "predicted.id", split.by = "group",cols = c("salmon1","plum","lightblue","lightpink","lightgreen"))

FeaturePlot(female_fibro, feature = c('ACAN','CD55','CHI3L2','NOTCH3','TPPP3'), cols = c("#b0e0e6","#8b0000"), ncol = 3, split.by = "group")
