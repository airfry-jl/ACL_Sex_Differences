# Library -----------------------------------------------------------------
library(data.table)
library(dplyr)
library(Seurat)
library(patchwork)
library(tibble)
library(CellChat)
library(plyr)
library(NMF)
library(ggalluvial)
library(devtools)
library(sctransform)
library(Matrix)
library(fields)
library(KernSmooth)
library(ROCR)
library(parallel)
library(DoubletFinder)

# Load Data ---------------------------------------------------------------
setwd("/path/to/wd")
acl <- readRDS(file = "acl_firstfive.rds")
DimPlot(acl, reduction = "umap")
freemuxlet <- fread("//freemuxlet.clust1.samples")

# Recluster + Add Freemux ---------------------------------------------------------------
#Add sample information (0-4)
df <- acl@meta.data
d <- tibble::rownames_to_column(df, "BARCODE")
f <- subset(freemuxlet, select = c("BARCODE","SNG.BEST.GUESS"))
colnames(f)[2] <-'Sample'
d <- left_join(d,f,by = "BARCODE")
d_col <- d[,-1]
rownames(d_col) <- d[,1]
d_col$Sample <- as.factor(d_col$Sample)

acl <- AddMetaData(object = acl, metadata = d_col$Sample, col.name = "Sample")

#test that the data is correct
meta <- acl@meta.data
all.equal(d_col, meta)

#Recluster
new.cluster.ids <- c("pericyte", "macrophage(M2)", "endothelial", "fibroblast",
                     "fibroblast","endothelial", "endothelial","fibroblast","T-cell", 
                     "fibroblast", "macrophage(M1)", "fibroblast/ FAP", "fibroblast", 
                     "fibroblast", "B-cells","fibroblast")
names(new.cluster.ids) <- levels(acl)
acl <- RenameIdents(acl, new.cluster.ids)
acl <- AddMetaData(object = acl, acl@active.ident, col.name = "cluster_id")

DimPlot(acl, reduction = "umap")

DimPlot(acl, reduction = "umap", group.by=c("ident","Sample"), split.by = "Sample")

# Male or Female ----------------------------------------------------------
genes.file = "~/Desktop/UCSF/ACL/genes.table.csv"

if (!file.exists(genes.file)) {
  suppressMessages(require(biomaRt))
  
  # initialize connection to mart, may take some time if the sites are
  # unresponsive.
  mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  
  # fetch chromosome info plus some other annotations
  genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                                   "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart,
                                    useCache = F))
  
  if (!dir.exists("~/Desktop/UCSF/ACL/results")) {
    dir.create("~/Desktop/UCSF/ACL/results")
  }
  if (is.data.frame(genes.table)) {
    write.csv(genes.table, file = genes.file)
  }
  
  if (!file.exists(genes.file)) {
    download.file("https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/master/labs/misc/genes.table.csv",
                  destfile = "~/Desktop/UCSF/ACL/genes.table.csv")
    genes.table = read.csv(genes.file)
  }
  
} else {
  genes.table = read.csv(genes.file)
}

genes.table <- genes.table[genes.table$external_gene_name %in% rownames(acl),
]

chrY.gene = genes.table$external_gene_name[genes.table$chromosome_name == "Y"]

acl$pct_chrY = colSums(acl@assays$RNA@counts[chrY.gene, ])/colSums(acl@assays$RNA@counts)


acl <- SetIdent(acl, value = acl@meta.data$Sample)


FeatureScatter(acl, feature1 = "XIST", feature2 = "pct_chrY")
VlnPlot(acl, features = c("XIST", "pct_chrY"))

gender.ids <- c("Male","Female","Male","Female","Male")
names(gender.ids) <- levels(acl)
acl <- RenameIdents(acl, gender.ids)
acl <- AddMetaData(object = acl, acl@active.ident, col.name = "Gender")


# Doublets ----------------------------------------------------------------
#acl <- SetIdent(acl, value = acl@meta.data$orig.ident)
#acl <- SCTransform(acl)
#acl <- RunPCA(acl)

#acl <- FindNeighbors(acl, dims = 1:15)
#acl <- FindClusters(acl, resolution = 0.8)
#acl <- RunUMAP(acl, dims = 1:15)

# pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_acl <- paramSweep_v3(acl, PCs = 1:15, sct = FALSE)
sweep.stats_acl <- summarizeSweep(sweep.res.list_acl, GT = FALSE)
bcmvn_acl <- find.pK(sweep.stats_acl)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- acl@meta.data$cluster_id
homotypic.prop <- modelHomotypic(annotations)           ## 
nExp_poi <- round(0.031*nrow(acl@meta.data))  ## Assuming 3.1% doublet formation rate with 4000 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
acl <- doubletFinder_v3(acl, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

#Remove Doublet
acl <- SetIdent(acl, value = acl@meta.data$DF.classifications_0.25_0.09_294)
DimPlot(acl, reduction = "umap")

doub <- as.factor(acl@meta.data$DF.classifications_0.25_0.09_294)
acl <- AddMetaData(object = acl, metadata = doub, col.name = "Singlet_Doublet")

acl <- subset(acl, subset = Singlet_Doublet == "Singlet")

acl <- SetIdent(acl, value = acl@meta.data$DF.classifications_0.25_0.09_294)
DimPlot(acl, reduction = "umap")


# Plots -------------------------------------------------------------------
#Load Data
setwd(dir = '~/Desktop/UCSF/ACL/Objects/')
acl <- readRDS(file = "~/Desktop/UCSF/ACL/acl_firstfive_gender_df.rds")

#Plot original data
acl <- SetIdent(acl, value = acl@meta.data$cluster_id)

DimPlot(acl, reduction = 'umap')

DimPlot(acl, reduction = "umap", group.by=c("ident","Sample"), split.by = "Sample")

#male v female
acl <- SetIdent(acl, value = acl@meta.data$Sample)
FeatureScatter(acl, feature1 = "XIST", feature2 = "pct_chrY")
VlnPlot(acl, features = c("XIST", "pct_chrY"))

acl <- SetIdent(acl, value = acl@meta.data$cluster_id)
DimPlot(acl, reduction = "umap", group.by=c("ident","Gender"), split.by = "Gender")
DimPlot(acl, reduction = "umap", group.by=c("ident","Gender", "Sample"), split.by = "Sample")




# Save Seurat Objects ---------------------------------------------------------------
saveRDS(acl,file ="~/Desktop/UCSF/ACL/acl_firstfive_gender_df.rds")

setwd(dir = '~/Desktop/UCSF/ACL/Objects/')
acl <- readRDS(file = "ACLmerge.rds")
# CellChat ----------------------------------------------------------------
male <- subset(acl, subset = Gender == "Male")
female <- subset(acl, subset = Gender == "Female")

# Male --------------------------------------------------------------------
male <- SetIdent(male, value = male@meta.data$celltype)

data.input_m <- GetAssayData(male, assay = "RNA", slot = "data") # normalized data matrix
labels_m <- Idents(male)
meta_m <- data.frame(group = labels_m, row.names = names(labels_m)) # create a dataframe of the cell labels
names(meta_m)[1] <- "labels"

cellchat_m <- createCellChat(object = data.input_m, meta = meta_m, group.by = "labels")

cellchat_m <- addMeta(cellchat_m, meta = meta_m, meta.name = "labels")
cellchat_m <- setIdent(cellchat_m, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_m@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_m@idents))

CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellchat_m@DB <- CellChatDB.use

cellchat_m <- subsetData(cellchat_m) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4)

cellchat_m <- identifyOverExpressedGenes(cellchat_m)
cellchat_m <- identifyOverExpressedInteractions(cellchat_m)  #LR

cellchat_m <- computeCommunProb(cellchat_m) #net
cellchat_m <- filterCommunication(cellchat_m, min.cells = 10)


cellchat_m <- computeCommunProbPathway(cellchat_m) #netP
cellchat_m <- aggregateNet(cellchat_m)


# Female ------------------------------------------------------------------

female <- SetIdent(female, value = female@meta.data$celltype)

data.input_f <- GetAssayData(female, assay = "RNA", slot = "data") # normalized data matrix
labels_f <- Idents(female)
meta_f <- data.frame(group = labels_f, row.names = names(labels_f)) # create a dataframe of the cell labels
names(meta_f)[1] <- "labels"

cellchat_f <- createCellChat(object = data.input_f, meta = meta_f, group.by = "labels")

cellchat_f <- addMeta(cellchat_f, meta = meta_f, meta.name = "labels")
cellchat_f <- setIdent(cellchat_f, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_f@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_f@idents))

CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellchat_f@DB <- CellChatDB.use

cellchat_f <- subsetData(cellchat_f) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4)

cellchat_f <- identifyOverExpressedGenes(cellchat_f)
cellchat_f <- identifyOverExpressedInteractions(cellchat_f)  #LR

cellchat_f <- computeCommunProb(cellchat_f) #net
cellchat_f <- filterCommunication(cellchat_f, min.cells = 10)


cellchat_f <- computeCommunProbPathway(cellchat_f) #netP
cellchat_f <- aggregateNet(cellchat_f)


# Incoming/ Outgoing ------------------------------------------------------
#Male
cellchat_m<- netAnalysis_computeCentrality(cellchat_m, slot.name = "netP")
selectK(cellchat_m, pattern = "outgoing")
nPatterns <- 3 
cellchat_m <- identifyCommunicationPatterns(cellchat_m, pattern = "outgoing", k = nPatterns)

selectK(cellchat_m, pattern = "incoming")
nPatterns <- 3
cellchat_m <- identifyCommunicationPatterns(cellchat_m, pattern = "incoming", k = nPatterns)

#Female
cellchat_f<- netAnalysis_computeCentrality(cellchat_f, slot.name = "netP")
selectK(cellchat_f, pattern = "outgoing")
nPatterns <- 3
cellchat_f <- identifyCommunicationPatterns(cellchat_f, pattern = "outgoing", k = nPatterns)

selectK(cellchat_f, pattern = "incoming")
nPatterns <- 3
cellchat_f <- identifyCommunicationPatterns(cellchat_f, pattern = "incoming", k = nPatterns)

# Save Cellchat Objects --------------------------------------------------------------
saveRDS(cellchat_m, "cellchat_m.rds")
saveRDS(cellchat_f, "cellchat_f.rds")

# Load Files ------------------------------------------------------------
cellchat_m <- readRDS("cellchat_m.rds")
cellchat_f <- readRDS("cellchat_f.rds")

updateCellChat(cellchat_m)
updateCellChat(cellchat_f)

object.list <- list(`male` = cellchat_m, 
                    `female` = cellchat_f)
cellchat <- mergeCellChat(object.list,add.names=names(object.list))

# Circle Plots ------------------------------------------------------------
all_male <- cellchat@netP$male$pathways
all_female <- cellchat@netP$female$pathways


all_male_df <- data.frame(all_male)
colnames(all_male_df) <- "Path"
all_female_df <- data.frame(all_female)
colnames(all_female_df) <- "Path"

in_both <- inner_join(all_male_df,all_female_df)
in_male <- anti_join(all_male_df, all_female_df)
in_female <- anti_join(all_female_df, all_male_df)

both <- all_male
both <- in_both$Path
male <- all_male
male <- in_male$Path
female <- all_female
female <- in_female$Path

pdf("Circle_Both.pdf")
for (j in 1:length(both)){
  pathway.show <- both[j]
  par(mfrow = c(1,2), xpd=TRUE)
  weight.max <- getMaxWeight(object.list[1:2], slot.name = c("netP"), attribute = pathway.show)
  for (i in 1:length(object.list)) {
      netVisual_aggregate(object.list[[i]], 
                          vertex.receiver = vertex.receiver, 
                          signaling = pathway.show, 
                          layout = "circle", 
                          vertex.weight = NULL, 
                          #edge.weight.max = weight.max[1], 
                          signaling.name = paste(pathway.show, 
                                                 names(object.list)[i]), 
                          arrow.size = 0.3)
  
  }
}
dev.off()

pdf("Circle_Male.pdf")
for (j in 1:length(male)){
  pathway.show <- male[j]
  par(mfrow = c(1,1), xpd=TRUE)
  weight.max <- getMaxWeight(object.list[1], slot.name = c("netP"), attribute = pathway.show)
  netVisual_aggregate(object.list[[1]], 
                        vertex.receiver = vertex.receiver, 
                        signaling = pathway.show, 
                        layout = "circle", 
                        vertex.weight = NULL, 
                        #edge.weight.max = weight.max[1], 
                        signaling.name = paste(pathway.show, 
                                               names(object.list)[1]), 
                        arrow.size = 0.3)
    
}
dev.off()

pdf("Circle_Female.pdf")
for (j in 1:length(female)){
  pathway.show <- female[j]
  par(mfrow = c(1,1), xpd=TRUE)
  weight.max <- getMaxWeight(object.list[2], slot.name = c("netP"), attribute = pathway.show)
  netVisual_aggregate(object.list[[2]], 
                      vertex.receiver = vertex.receiver, 
                      signaling = pathway.show, 
                      layout = "circle", 
                      vertex.weight = NULL, 
                      #edge.weight.max = weight.max[1], 
                      signaling.name = paste(pathway.show, 
                                             names(object.list)[2]), 
                      arrow.size = 0.3)
  
}
dev.off()

pdf("Circle_I.pdf")
path <- c("OSM", "CALCR", "NECTIN", "VCAM", "COLLAGEN", "PDGF", "TWEAK")
for (j in 1:length(path)){
  pathway.show <- path[j]
  par(mfrow = c(1,2), xpd=TRUE)
  weight.max <- getMaxWeight(object.list[1:2], slot.name = c("netP"), attribute = pathway.show)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], 
                        vertex.receiver = vertex.receiver, 
                        signaling = pathway.show, 
                        layout = "circle", 
                        vertex.weight = NULL, 
                        #edge.weight.max = weight.max[1], 
                        signaling.name = paste(pathway.show, 
                                               names(object.list)[i]), 
                        arrow.size = 0.3)
    
  }
}

dev.off()





#When signaling for certain groups does not show up 
netVisual_aggregate(cellchat_pc, signaling = pathways.show, signaling.name = paste(pathways.show, names(object.list)[1]),layout = "circle")

netVisual_aggregate(cellchat_pd, signaling = pathways.show,signaling.name = paste(pathways.show, names(object.list)[2]), layout = "circle")

netVisual_aggregate(cellchat_fc, signaling = pathways.show, signaling.name = paste(pathways.show, names(object.list)[3]),layout = "circle")

netVisual_aggregate(cellchat_fd, signaling = pathways.show,signaling.name = paste(pathways.show, names(object.list)[4]), layout = "circle")


# LR Plots ----------------------------------------------------------------------
pathways <- c("TGFb", "APP", "THY1", "TENASCIN", "VISFATIN", "SEMA3", "HGF", "BMP")
pathways.show <- pathways

path <- "~/Desktop/UCSF/Steven/Cellchat/Paper/General/LR/"
for(i in 1:8){
  tryCatch({
    pdf(paste0(path,"general_LR_",pathways.show[i],".pdf"))
    print(netAnalysis_contribution(cellchat_pc,title = paste("Contribution of each L-R pair in", pathways.show[i], "Signaling Pathway, Partial Cuff"), signaling = pathways.show[i]))
    print(netAnalysis_contribution(cellchat_fc,title = paste("Contribution of each L-R pair in", pathways.show[i], "Signaling Pathway, Full Cuff"), signaling = pathways.show[i]))
    print(netAnalysis_contribution(cellchat_pd,title = paste("Contribution of each L-R pair in", pathways.show[i], "Signaling Pathway, Partial Delt"), signaling = pathways.show[i]))
    print(netAnalysis_contribution(cellchat_fd,title = paste("Contribution of each L-R pair in", pathways.show[i], "Signaling Pathway, Full Delt"), signaling = pathways.show[i]))
    dev.off()
  }, error=function(e){})
}

#In case some objects do not have significant communication
pathways.show <- "HGF"
pdf(paste0("~/Desktop/UCSF/Steven/Cellchat/Paper/General/LR/general_LR_",pathways.show,".pdf"))
netAnalysis_contribution(cellchat_pc,title = paste("Contribution of each L-R pair in", pathways.show, "Signaling Pathway, Partial Cuff"), signaling = pathways.show)
netAnalysis_contribution(cellchat_fc,title = paste("Contribution of each L-R pair in", pathways.show, "Signaling Pathway, Full Cuff"), signaling = pathways.show)
netAnalysis_contribution(cellchat_pd,title = paste("Contribution of each L-R pair in", pathways.show, "Signaling Pathway, Partial Delt"), signaling = pathways.show)
netAnalysis_contribution(cellchat_fd,title = paste("Contribution of each L-R pair in", pathways.show, "Signaling Pathway, Full Delt"), signaling = pathways.show)
dev.off()



# HeatMap -----------------------------------------------------------------
cellchat_pc <- netAnalysis_computeCentrality(cellchat_pc, slot.name = "netP") 
cellchat_pd <- netAnalysis_computeCentrality(cellchat_pd, slot.name = "netP")
cellchat_fc <- netAnalysis_computeCentrality(cellchat_fc, slot.name = "netP")
cellchat_fd <- netAnalysis_computeCentrality(cellchat_fd, slot.name = "netP")

pathways.show <- pathways

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Heatmap/general_heatmap_pc.pdf")

#pc
for(i in 1:8){
  netAnalysis_signalingRole_network(cellchat_pc, 
                                    signaling = pathways.show[i], 
                                    width = 8, 
                                    height = 2.5, 
                                    font.size = 10)
}
dev.off()

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Heatmap/general_heatmap_fc.pdf")
#fc
for(i in 1:8){
  netAnalysis_signalingRole_network(cellchat_fc, 
                                    signaling = pathways.show[i], 
                                    width = 8, 
                                    height = 2.5, 
                                    font.size = 10)
}
dev.off()

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Heatmap/general_heatmap_pd.pdf")
#pd
for(i in 1:8){
  tryCatch({
    netAnalysis_signalingRole_network(cellchat_pd, 
                                      signaling = pathways.show[i], 
                                      width = 8, 
                                      height = 2.5, 
                                      font.size = 10)
  }, error=function(e){})
}
dev.off()

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Heatmap/general_heatmap_fd.pdf")
#fd
for(i in 1:8){
  netAnalysis_signalingRole_network(cellchat_fd, 
                                    signaling = pathways.show[i], 
                                    width = 8, 
                                    height = 2.5, 
                                    font.size = 10)
}
dev.off()

# Fold Change -------------------------------------------------------------
pathways.show <- c("HGF", "THY1", "SEMA3","BMP")
#PC vs FC
gg1<- netVisual_bubble(cellchat, 
                       sources.use = c(1:12), 
                       targets.use = c(1:12),  
                       comparison = c(1,3),
                       title.name = "Increased signaling in Full Cuff", 
                       angle.x = 45, 
                       remove.isolate = T,
                       max.dataset = 3,
                       return.data = T)


gg2<- netVisual_bubble(cellchat, 
                       sources.use = c(1:12), 
                       targets.use = c(1:12),  
                       comparison = c(1,3),
                       title.name = "Decreased signaling in Full Cuff", 
                       angle.x = 45, 
                       remove.isolate = T,
                       max.dataset = 1,
                       return.data = T)

write.csv(gg1$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/General/FC/gen_Inc_FCvPC.csv")
write.csv(gg2$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/General/FC/gen_Dec_FCvPC.csv")

#dot plot
pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/UpDown/fc_pc.pdf")
pathways <- pathways.show[4]



netVisual_bubble(cellchat, 
                 sources.use = c(1:7), 
                 targets.use = c(1:7),
                 signaling = pathways,
                 comparison = c(1,3),
                 color.heatmap = c("viridis"),
                 title.name = paste0("Increased ",pathways, " signaling in Full Cuff vs Partial Cuff"), 
                 angle.x = 45, 
                 remove.isolate = T,
                 max.dataset = 3) +
  scale_colour_gradient2(high = "#8b0000",                                  
                         low = "#B6D0E2",
                         midpoint = 0.16)

netVisual_bubble(cellchat, 
                 sources.use = c(1:7), 
                 targets.use = c(1:7),  
                 signaling = pathways,
                 comparison = c(1,3),
                 title.name = paste0("Decreased ",pathways, " signaling in Full Cuff vs Partial Cuff"), 
                 angle.x = 45, 
                 remove.isolate = T,
                 max.dataset = 1) +
  scale_colour_gradient2(high = "#8b0000",                                  
                         low = "#B6D0E2",
                         midpoint = 0.137)
dev.off()



#PC vs PD
gg3<- netVisual_bubble(cellchat, 
                       sources.use = c(1:12), 
                       targets.use = c(1:12),  
                       comparison = c(1,2),
                       title.name = "Increased signaling in Partial Cuff", 
                       angle.x = 45, 
                       remove.isolate = T,
                       max.dataset = 1,
                       return.data = T)


gg4<- netVisual_bubble(cellchat, 
                       sources.use = c(1:12), 
                       targets.use = c(1:12),  
                       comparison = c(1,2),
                       title.name = "Decreased signaling in Partial Cuff", 
                       angle.x = 45, 
                       remove.isolate = T,
                       max.dataset = 2,
                       return.data = T)

write.csv(gg3$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/General/FC/gen_Inc_PCvPD.csv")
write.csv(gg4$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/General/FC/gen_Dec_PCvPD.csv")

#FC vs FD
gg5<- netVisual_bubble(cellchat, 
                       sources.use = c(1:12), 
                       targets.use = c(1:12),  
                       comparison = c(3,4),
                       title.name = "Increased signaling in Full Cuff", 
                       angle.x = 45, 
                       remove.isolate = T,
                       max.dataset = 3,
                       return.data = T)


gg6<- netVisual_bubble(cellchat, 
                       sources.use = c(1:12), 
                       targets.use = c(1:12),  
                       comparison = c(3,4),
                       title.name = "Decreased signaling in Full Cuff", 
                       angle.x = 45, 
                       remove.isolate = T,
                       max.dataset = 4,
                       return.data = T)

write.csv(gg5$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/General/FC/gen_Inc_FCvFD.csv")
write.csv(gg6$communication, file = "~/Desktop/UCSF/Steven/Cellchat/Paper/General/FC/gen_Dec_FCvFD.csv")

#dot plot
pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/UpDown/fc_fd.pdf")

pathways <- pathways.show[1]

netVisual_bubble(cellchat, 
                 sources.use = c(1:7), 
                 targets.use = c(1:7),
                 signaling = pathways,
                 comparison = c(3,4),
                 title.name = paste0("Increased ",pathways, " signaling in Full Cuff vs Full Delt"), 
                 angle.x = 45, 
                 remove.isolate = T,
                 max.dataset = 3) +
  scale_colour_gradient2(high = "#8b0000",                                  
                         low = "#B6D0E2",
                         midpoint = 0.17)

netVisual_bubble(cellchat, 
                 sources.use = c(1:7), 
                 targets.use = c(1:7),  
                 signaling = pathways,
                 comparison = c(3,4),
                 title.name = paste0("Decreased ",pathways, " signaling in Full Cuff vs Full Delt"), 
                 angle.x = 45, 
                 remove.isolate = T,
                 max.dataset = 4) +
  scale_colour_gradient2(high = "#8b0000",                                  
                         low = "#B6D0E2",
                         midpoint = 0.17)
dev.off()



# Chord -------------------------------------------------------------------

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Chord/general_chord.pdf")
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]]@net$count, net = object.list[[i]]@net$count, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

pdf("~/Desktop/UCSF/Steven/Cellchat/Paper/General/Chord/general_chord_weight.pdf")
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]]@net$weight, net = object.list[[i]]@net$weight, title.name = paste0("Interaction weights/strength - ", names(object.list)[i]))
}
dev.off()
