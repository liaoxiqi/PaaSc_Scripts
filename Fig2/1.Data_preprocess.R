library(Seurat)
library(SeuratData)
library(tidyverse)
library(scRNAseq)

#--- Bmcite ---#
InstallData("bmcite")
data("bmcite")
bmcite.meta <- bmcite@meta.data
bmcite <- CreateSeuratObject(bmcite@assays$RNA@counts, min.cells = 10)
bmcite@meta.data <- bmcite.meta
bmcite$celltype.l2 <- recode(bmcite$celltype.l2, "CD4 Memory" = "CD4 T", "CD4 Naive" = "CD4 T", 
                             "CD8 Effector_1" = "CD8 T", "CD8 Effector_2" = "CD8 T", "CD8 Memory_1" = "CD8 T", 
                             "CD8 Memory_2" = "CD8 T", "CD8 Naive" = "CD8 T", "cDC2" = "cDC", "Memory B" = "B", 
                             "Naive B" = "B")
Idents(bmcite) <- "celltype.l2"
bmcite <- subset(bmcite, idents = c("CD4 T", 
                                    "B",  
                                    "CD14 Mono", 
                                    "NK",  "CD8 T", 
                                    "cDC", "GMP", "HSC", "pDC",
                                    "CD16 Mono"))
bmcite <- subset(bmcite, downsample = 1500)
bmcite@assays$ADT <- NULL
saveRDS(bmcite, "bmcite.RDS")

#--- liver ---#
liver.imm <- ZhaoImmuneLiverData()
liver.imm <- as.Seurat(liver.imm, data = NULL)

Idents(liver.imm) <- "sample"
liver.immune <- subset(liver.imm, idents = "donor 1 liver")
liver.immune.meta <- liver.immune@meta.data
matx <- liver.immune@assays$originalexp@counts
rownames(matx) <- liver.imm@assays$originalexp@meta.features$Symbol
liver.immune <- CreateSeuratObject(matx, min.cells = 10)
liver.immune@meta.data <- liver.immune.meta
liver.immune$cell.type <- as.character(liver.immune$fine)
Idents(liver.immune) <- "cell.type"
liver.immune <- liver.immune %>%
  NormalizeData(verbose = F) %>%
  FindVariableFeatures(verbose = F) %>%
  ScaleData(verbose = F) %>%
  RunPCA(verbose = F) %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:10, verbose = F)
liver.immune$cell.type <- recode(liver.immune$cell.type, "CCR7+ CD4" = "CD4", "CCR7+ CD8" = "CD8","CRCR3+ CD8" = "CD8","CRCR3+ CD4" = "CD4",
                                 "CX3CR1+ CD8" = "CD8", "CXCR6+ CD8" = "CD8","Naive B" = "Naive B", "cycling ASC" = "Plasma", "noncycline ASC" = "Plasma",
                                 "CD11c+ mem B" = "mem B","Classical mem B" = "mem B","CX3CR1+ NK" = "NK","CXCR6+ NK" = "NK", "Treg" = "Treg",
                                 "CD14+ Mo" = "CD14+ Mo", "TCL1A+ naive B" = "Naive B") 
saveRDS(liver.immune, "liver.immune.RDS")

#--- spleen ---#
spleen.immune <- subset(liver.imm, idents = "donor 1 spleen")

spleen.immune.meta <- spleen.immune@meta.data
matx <- spleen.immune@assays$originalexp@counts
rownames(matx) <- liver.imm@assays$originalexp@meta.features$Symbol
spleen.immune <- CreateSeuratObject(matx, min.cells = 10)
spleen.immune@meta.data <- spleen.immune.meta

spleen.immune$cell.type <- as.character(spleen.immune$fine)
Idents(spleen.immune) <- "cell.type"
spleen.immune <- spleen.immune %>%
                  NormalizeData(verbose = F) %>%
                  FindVariableFeatures(verbose = F) %>%
                  ScaleData(verbose = F) %>%
                  RunPCA(verbose = F) %>%
                  FindNeighbors() %>%
                  RunUMAP(dims = 1:10, verbose = F)

spleen.immune$cell.type <- recode(spleen.immune$cell.type, "CCR7+ CD4" = "CD4", "CCR7+ CD8" = "CD8","CRCR3+ CD8" = "CD8","CRCR3+ CD4" = "CD4",
                                 "CX3CR1+ CD8" = "CD8", "CXCR6+ CD8" = "CD8","Naive B" = "Naive B", "cycling ASC" = "Plasma", "noncycline ASC" = "Plasma",
                                 "CD11c+ mem B" = "mem B","Classical mem B" = "mem B","CX3CR1+ NK" = "NK","CXCR6+ NK" = "NK", "Treg" = "Treg",
                                 "CD14+ Mo" = "CD14+ Mo", "CD16+ Mo" = "CD16+ Mo", "TCL1A+ naive B" = "Naive B", "CD1c+ DCs" = "DC") 

spleen.immune <- readRDS("spleen.immune.RDS")


#--- Hcortex ---#
hcortex.sce <- NowakowskiCortexData()
hcortex <- CreateSeuratObject(hcortex.sce@assays@data$tpm, min.cells = 10)
hcortex$cell.type <- hcortex.sce$WGCNAcluster
hcortex$cell.type <- recode(hcortex$cell.type, "EN-PFC1" = "EN", "EN-PFC2" = "EN","EN-PFC3" = "EN","EN-V1-1" = "EN",
                            "EN-V1-2" = "EN","EN-V1-3" = "EN","IN-CTX-CGE1" = "IN", "IN-CTX-CGE2" = "IN", "IN-CTX-MGE1" = "IN", 
                            "IN-CTX-MGE2" = "IN", "IN-STR" = "IN", "nIN1" = "IN", "nIN2" = "IN", "nIN3" = "IN", "nIN4" = "IN", "nIN5" = "IN", 
                            "IPC-div1" = "IPC", "IPC-div2" = "IPC", "IPC-nEN1" = "IPC", 
                            "IPC-nEN2" = "IPC", "IPC-nEN3" = "IPC", "MGE-IPC1" = "IPC", "MGE-IPC2" = "IPC", "MGE-IPC3" = "IPC", 
                            "MGE-RG1" = "RG", "MGE-RG2" = "RG", "RG-div1" = "RG", "RG-div2" = "RG","RG-early" = "RG", "tRG" = "RG",
                            "vRG" = "RG","oRG" = "RG")
saveRDS(hcortex, "hcortex.RDS")

#--- pancrea ---#
pancrea <- BaronPancreasData('human')
pancrea <- as.Seurat(pancrea, data = NULL)

pancrea.meta <- pancrea@meta.data
pancrea <- CreateSeuratObject(pancrea@assays$originalexp@counts, min.cells = 10)
pancrea@meta.data <- pancrea.meta
Idents(pancrea) <- "donor"
pancrea <- subset(pancrea, idents = "GSM2230759")
pancrea$cell.type <- as.character(pancrea$label)
Idents(pancrea) <- "cell.type"

pancrea <- pancrea %>%
  NormalizeData(verbose = F) %>%
  FindVariableFeatures(verbose = F) %>%
  ScaleData(verbose = F) %>%
  RunPCA(verbose = F) %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:10, verbose = F)

saveRDS(pancrea, "pancrea.RDS")

