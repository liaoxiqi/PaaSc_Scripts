library(BiocParallel)
library(GSVA)
library(Seurat)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(gsdensity)
library(CelliD)
library(VAM)
library(DescTools)
library(GSEABase)
library(AUCell)
library(nnet)
library(reshape2)
library(RColorBrewer)
library(PaaSc)





###### scGSEA
library(Matrix)
library(fgsea)
library(gficf)

#' Calculate scGSEA activity scores
#'
#' @param glist Gene set list
#' @param pos.id Pathway identifier
#' @param tmp NMF result object containing w and h matrices
#' @param data Main data object containing rawCounts
#' @param sr Seurat object containing cell type metadata
#' @param nsim Number of permutations (default=10000)
#' @param fdr.th FDR threshold (default=0.05)
#' @param minSize Minimum gene set size (default=5)
#' @param maxSize Maximum gene set size (default=Inf)
#' @param gp GSEA parameter (default=0)
#' 
#' @return A data.frame containing scGSEA scores and cell types

calculate_scgsea <- function(glist, pos.id, tmp, data, sr, 
                             nsim = 10000, fdr.th = 0.05,
                             minSize = 5, maxSize = Inf, gp = 0) {
  data$scgsea <- list(
    nmf.w = Matrix::Matrix(tmp$w, sparse = TRUE),
    nmf.h = t(Matrix::Matrix(tmp$h, sparse = TRUE)),
    pathways = setNames(list(glist), pos.id),
    es = NULL, nes = NULL, pval = NULL, fdr = NULL,
    x = NULL, stat = NULL
  )
  
 
  n_pathways <- length(data$scgsea$pathways)
  n_components <- ncol(data$scgsea$nmf.w)
  
  data$scgsea$es <- data$scgsea$nes <- 
  data$scgsea$pval <- data$scgsea$fdr <- 
  Matrix::Matrix(0, n_pathways, n_components)
  
  rownames(data$scgsea$es) <- rownames(data$scgsea$nes) <-
  rownames(data$scgsea$pval) <- rownames(data$scgsea$fdr) <- pos.id
  
  nt_fgsea <- ceiling(n_pathways / 100)
  bpparameters <- BiocParallel::SnowParam(nt_fgsea)
  rownames(data$scgsea$nmf.w) <- rownames(data$rawCounts)
  
  pb <- utils::txtProgressBar(0, n_components, style = 3)
  for(i in seq_len(n_components)) {
    stats <- as.matrix(data$scgsea$nmf.w)[, i]
    df <- fgsea::fgseaMultilevel(
      pathways = data$scgsea$pathways,
      stats = stats,
      nPermSimple = nsim,
      gseaParam = gp,
      BPPARAM = bpparameters,
      minSize = minSize,
      maxSize = maxSize
    )[, 1:7]
    
    data$scgsea$es[df$pathway, i] <- df$ES
    data$scgsea$nes[df$pathway, i] <- df$NES
    data$scgsea$pval[df$pathway, i] <- df$pval
    data$scgsea$fdr[df$pathway, i] <- df$padj
    
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  handle_nas <- function(x) replace(x, is.na(x), 0)
  data$scgsea$nes <- handle_nas(data$scgsea$nes)
  data$scgsea$pval <- handle_nas(data$scgsea$pval)
  data$scgsea$fdr <- handle_nas(data$scgsea$fdr)
  
  data$scgsea$x <- data$scgsea$nes
  data$scgsea$x[data$scgsea$x < 0 | data$scgsea$fdr >= fdr.th] <- 0
  data$scgsea$x <- Matrix::Matrix(
  data$scgsea$nmf.h %*% t(data$scgsea$x), 
    sparse = TRUE
  )
  
  oo <- as.data.frame(data$scgsea$x)
  colnames(oo) <- "scGSEA"
  oo$cell.type <- sr@meta.data$cell.type
  
  return(oo)
}

get.activity.list <- function(glist, pos.id, object, el, cells, cells_rankings){
  lst <- list()
  print(glist)

  scgsea.res <- calculate_scgsea(
    glist = glist,
    pos.id = pos.id,
    tmp = tmp,         
    data = data,       
    sr = sr,           
    nsim = 10000,
    fdr.th = 0.05
  )
  
  lst[['scGSEA']] <- scgsea.res


  lst1 <- list()
  for (method in names(lst)){
    lst1[[method]] <- calc.auc(input = lst[[method]], method.id = method, pos.id = pos.id)
  }
  return(list(score.list = lst, auc.list = lst1))
}


