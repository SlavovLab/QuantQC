#' Computes PCA
#'
#' This function takes a QQC object and computes PCA on the unimputed protein level data taking the eigen
#' values of a correlation matrix computed off pairwise observations.
#'
#' @param QQC a QuantQC object
#' @return A \code{QQC object} with the reductions slot containing a PCA data.frame in list
#' @examples
#' ComputePCA(TestSamples)
#' @export
ComputePCA <- function(QQC,imputed = T){
  if(imputed == T){
    sc.batch_cor <- QQC@matricies@protein.imputed
  }else{
    sc.batch_cor <- QQC@matricies@protein
  }

  cellenONE_meta <- QQC@meta.data

  # Correlation matrix for PCA
  cor_mat <- cor(sc.batch_cor,use = 'pairwise.complete.obs')

  # Eigen values of correlation matrix for PCA
  sc.pca <- eigen(cor_mat)
  scx<-as.data.frame(sc.pca$vectors)
  colnames(scx)<-paste0("PC",1:ncol(scx))

  # Calculate and plot % variance of PCs
  percent_var <- sc.pca$values/sum(sc.pca$values)*100
  plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")



  # Map meta data for plotting PCAs
  scx$ID <-colnames(sc.batch_cor)
  scx <- scx %>% left_join(cellenONE_meta,by = c('ID'))

  QQC@reductions <- list(PCA = scx)

  QQC@misc[['pct_var']] <- round(percent_var,2)

  return(QQC)

}

#' Plots PCA
#'
#' This function takes a QQC object and plots PC1 and PC2 as a scatter plot. It can be colored by different
#' options like sort condition, or various potential batch effect sorces
#'
#' @param QQC a QuantQC object
#' @param by a string that specifies how to color code PCA
#' @return A \code{QQC object} with the reductions slot containing a PCA data.frame in list
#' @examples
#' ComputePCA(TestSamples)
#' @export
PlotPCA <- function(QQC, by = 'Condition'){
  PCA_plot <- QQC@reductions[['PCA']]

  pct_var <- QQC@misc[['pct_var']]

  if(by == 'Condition'){
    pca_plot <- ggplot(PCA_plot, aes(x = PC1, y = PC2, color = sample)) + geom_point()+
      dot_plot + ggtitle('PCA by Condition') + xlab(paste0('PC1(',pct_var[1],'%)'))+
      ylab(paste0('PC2(',pct_var[2],'%)'))
  }
  if(by == 'Total protein'){
    pca_plot <- ggplot(PCA_plot, aes(x = PC1, y = PC2, color = prot_total)) + geom_point()+ ggtitle('PCA by Total Cell Intensity')+
      dot_plot + xlab(paste0('PC1(',pct_var[1],'%)'))+
      ylab(paste0('PC2(',pct_var[2],'%)'))+
      scale_color_gradient2(midpoint = median(PCA_plot$prot_total), low = 'blue',mid = 'white', high = 'red')
  }

  if(by == 'Label'){
    pca_plot <- ggplot(PCA_plot, aes(x = PC1, y = PC2, color = label)) + geom_point()+
      dot_plot + ggtitle('PCA by Label')+ xlab(paste0('PC1(',pct_var[1],'%)'))+
      ylab(paste0('PC1(',pct_var[2],'%)'))
  }

  if(by == 'Run order'){
    pca_plot <- ggplot(PCA_plot, aes(x = PC1, y = PC2, color = Order)) + geom_point() + ggtitle('PCA by Run Order')+
      dot_plot + xlab(paste0('PC1(',pct_var[1],'%)'))+
      ylab(paste0('PC2(',pct_var[2],'%)'))+
      scale_color_gradient2(midpoint = median(PCA_plot$Order,na.rm = T), low = 'blue',mid = 'white', high = 'red')
  }

  pca_plot

}




#' Computes PCA
#'
#' This function takes a QQC object and computes PCA on the unimputed protein level data taking the eigen
#' values of a correlation matrix computed off pairwise observations.
#'
#' @param QQC a QuantQC object
#' @return A \code{QQC object} with the reductions slot containing a PCA data.frame in list
#' @examples
#' ComputePCA(TestSamples)
#' @export
FeaturePCA <- function(QQC, prot = NA, imputed = T){
  PCA_plot <- QQC@reductions[['PCA']]

  if(imputed == T){
    prot_mat <- QQC@matricies@protein.imputed
  }
  if(imputed == F){
    prot_mat <- QQC@matricies@protein
  }

  PCA_plot$protein <- prot_mat[prot,]

  pca_plot <- ggplot(PCA_plot, aes(x = PC1, y = PC2, color = protein)) + geom_point()+
    dot_plot + scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')



  pca_plot

}











#' Computes UMAP
#'
#' This function takes a QQC object and computes UMAP on the imputed protein level data with the UMAP function
#' from the Seurat package
#'
#' @param QQC a QuantQC object
#' @return A \code{QQC object} with the reductions slot containing a UMAP data.frame in list
#' @examples
#' ComputeUMAP(TestSamples)
#' @export
ComputeUMAP <- function(QQC){
  protein_Data <- QQC@matricies@protein.imputed
  scx <- QQC@reductions[['PCA']]

  prot_umap <- CreateSeuratObject(counts = protein_Data, project = "prot_mat")
  prot_umap <- NormalizeData(prot_umap, normalization.method = "LogNormalize", scale.factor = 10000)
  prot_umap@assays$RNA@layers$data <- protein_Data

  all.genes <- rownames(protein_Data)
  prot_umap <- ScaleData(prot_umap, features = all.genes)
  prot_umap@assays$RNA@layers$data <- protein_Data
  prot_umap <- Seurat::RunPCA(prot_umap, features = all.genes)
  prot_umap <- FindNeighbors(prot_umap, dims = 1:6)
  prot_umap <- FindClusters(prot_umap, resolution = 0.5)
  prot_umap <- RunUMAP(prot_umap, dims = 1:6)

  um_plot <- as.data.frame(prot_umap@reductions[["umap"]]@cell.embeddings)

  um_plot$sample <- scx$sample
  um_plot$cluster <- prot_umap@meta.data[["RNA_snn_res.0.5"]]
  um_plot$lab <- scx$label
  um_plot$prot_total <- scx$prot_total
  um_plot$Order <- scx$Order
  um_plot$diameter <- scx$diameter




  # ggplot(um_plot, aes(x = UMAP_1,y = UMAP_2, color = cluster)) + geom_point(size = 3) + theme_classic()+
  #   ggtitle('UMAP colored by Sample')

  QQC@reductions[['UMAP']] <- um_plot

  return(QQC)
}

#' Computes PCA
#'
#' This function takes a QQC object and computes PCA on the unimputed protein level data taking the eigen
#' values of a correlation matrix computed off pairwise observations.
#'
#' @param QQC a QuantQC object
#' @return A \code{QQC object} with the reductions slot containing a PCA data.frame in list
#' @examples
#' ComputePCA(TestSamples)
#' @export
PlotUMAP <- function(QQC, by = 'Cluster'){
  UMAP_plot <- QQC@reductions[['UMAP']]

  if(by == 'Cluster'){
    umap_plot <- ggplot(UMAP_plot, aes(x = umap_1, y = umap_2, color = cluster)) + geom_point()+
      um_plot
  }

  if(by == 'Condition'){
    umap_plot <- ggplot(UMAP_plot, aes(x = umap_1, y = umap_2, color = sample)) + geom_point()+
      um_plot
  }
  if(by == 'Total protein'){
    umap_plot <- ggplot(UMAP_plot, aes(x = umap_1, y = umap_2, color = prot_total)) + geom_point()+
      um_plot + scale_color_gradient2(midpoint = median(UMAP_plot$prot_total), low = 'blue',mid = 'white', high = 'red')
  }

  if(by == 'Label'){
    umap_plot <- ggplot(UMAP_plot, aes(x = umap_1, y = umap_2, color = lab)) + geom_point()+
      um_plot
  }

  if(by == 'Run order'){
    umap_plot <- ggplot(UMAP_plot, aes(x = umap_1, y = umap_2, color = Order)) + geom_point()+
      um_plot+ scale_color_gradient2(midpoint = median(UMAP_plot$Order), low = 'blue',mid = 'white', high = 'red')

  }

  if(by == 'Diameter'){
    umap_plot <- ggplot(UMAP_plot, aes(x = umap_1, y = umap_2, color = diameter)) + geom_point()+
      um_plot+ scale_color_gradient2(midpoint = mean(UMAP_plot$diameter), low = 'blue',mid = 'white', high = 'red')

  }

  umap_plot

}

#' Computes PCA
#'
#' This function takes a QQC object and computes PCA on the unimputed protein level data taking the eigen
#' values of a correlation matrix computed off pairwise observations.
#'
#' @param QQC a QuantQC object
#' @return A \code{QQC object} with the reductions slot containing a PCA data.frame in list
#' @examples
#' ComputePCA(TestSamples)
#' @export
FeatureUMAP <- function(QQC, prot = NA, imputed = T){
  UMAP_plot <- QQC@reductions[['UMAP']]

  if(imputed == T){
    prot_mat <- QQC@matricies@protein.imputed
  }
  if(imputed == F){
    prot_mat <- QQC@matricies@protein
  }

  UMAP_plot$protein <- prot_mat[prot,]

  umap_plot <- ggplot(UMAP_plot, aes(x = umap_1, y = umap_2, color = protein)) + geom_point()+
      dot_plot + scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')



  umap_plot

}







FeatureUMAP_abs <- function(QQC, prot = NA){
  UMAP_plot <- QQC@reductions[['UMAP']]

  prot_mat <- QQC@matricies@protein_abs

  UMAP_plot$protein <- prot_mat[prot,]

  umap_plot <- ggplot(UMAP_plot, aes(x = umap_1, y = umap_2, color = log2(protein))) + geom_point()+
    dot_plot + scale_color_gradient2(midpoint = median(log2(UMAP_plot$protein),na.rm = T), low = 'blue',mid = 'white', high = 'red')


  umap_plot

}



ClustBoxPlot <- function(QQC, prot = NA, imputed = F){
  UMAP_plot <- QQC@reductions[['UMAP']]

  if(imputed == T){
    prot_mat <- QQC@matricies@protein.imputed
  }
  if(imputed == F){
    prot_mat <- QQC@matricies@protein
  }

  UMAP_plot$protein <- prot_mat[prot,]

  box_plot <- ggplot(UMAP_plot, aes(x = cluster, y = protein,color = cluster)) + geom_boxplot()+
    dot_plot



  box_plot

}


