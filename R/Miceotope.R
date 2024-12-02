PepCorMini <- function(peptide_data,peptide_protein_map){

  # Initialized empty matricies to store correlations between peptides
  mat_stor = matrix(data = NA,nrow = 10000,ncol = 3)
  mat_stor_fake = matrix(data = NA,nrow = 10000,ncol = 2)

  # List of all the unique proteins in data
  upl <- unique(peptide_protein_map$Protein)

  # List describing which protein each peptide comes from
  prot_list <- peptide_protein_map$Protein


  peptide_data <- Normalize_reference_vector_log(log2(peptide_data))

  count <- 0

  # Loop over each protein, calculate correlations between peptides mapping to a protein
  for(i in upl){

    # Matrix for a single protein
    mat_p1 <- peptide_data[which(prot_list == i),]
    if(is.null(nrow(mat_p1)) == FALSE ){

      # calculate pairwise observations (how many times peptides are observed in same cells)
      obs_mat <- psych::pairwiseCount(t(mat_p1))

      # calculate correlations between peptides
      cor_mat <- cor(t(mat_p1),use = 'pairwise.complete.obs')

      # Require atleast 5 pairwise observations to plot the correlations
      cor_mat[obs_mat < 4] <- NA

      # Turn correlation matrix to a vector for storage, store correlations and protein identity
      cor_mat <- cor_mat[lower.tri(cor_mat)]
      if(is.na(median(cor_mat,na.rm = T)) == F){
        count = count + 1
        mat_stor[count,1] <- i
        mat_stor[count,2] <- median(cor_mat,na.rm = T)
        mat_stor[count,3] <- nrow(mat_p1)
      }
    }
  }

  mat_stor <- as.data.frame(mat_stor)
  colnames(mat_stor) <- c('Protein','Cor','Obs')
  mat_stor <- mat_stor %>% filter(is.na(Protein) == F)
  mat_stor$Cor <- as.numeric(mat_stor$Cor)
  return(mat_stor)


}


Miceotope_cellXpeptide <- function(QQC,TQVal = 1, chQVal = 1, t = 5){


  sc.data <- QQC@raw_data

  sc.data <- sc.data %>% filter(Channel.Q.Value < chQVal)
  sc.data <- sc.data %>% filter(Translated.Q.Value < TQVal)

  sc.data$plex_new <- NA
  sc.data$plex_new[sc.data$plex == "0" | sc.data$plex == "1"] <- '0'
  sc.data$plex_new[sc.data$plex == "2" | sc.data$plex == "3"] <- '4'
  sc.data$plex <- sc.data$plex_new
  sc.data$plex_new <- NULL

  sc.data$ID <- paste0(sc.data$Well,sc.data$plate,'.',sc.data$plex)

  sc.data$pep_type <- str_sub(sc.data$Stripped.Sequence,-1,-1)

  sc.data <- sc.data %>% filter(Ms1.Area != 0)

  mice_K <- sc.data %>% filter(pep_type == 'K')
  mice_R <- sc.data %>% filter(pep_type == 'R')

  mice_K <- mice_K %>% filter(Run %in% mice_R$Run)
  mice_R <- mice_R %>% filter(Run %in% mice_K$Run)

  mice_R$seqRun <- paste0(mice_R$seqcharge,mice_R$ID)
  mice_R <- mice_R %>% distinct(seqRun,.keep_all = T)

  mice_R$seqRun <- NULL




  mice_R_ <- reshape2::dcast(mice_R,Protein.Group+seqcharge ~ ID, value.var = 'Ms1.Area')

  n_occur <- data.frame(table(mice_R_$seqcharge))
  n_occur <- n_occur[n_occur$Freq > 1,]
  if(nrow(n_occur) > 0){
    mice_R_ <- mice_R_ %>% filter(seqcharge %in% n_occur$Var1)
    mice_R_ <- mice_R_ %>% distinct(seqcharge,.keep_all = T)

    for(i in 1:nrow(mice_R_)){
      mice_R$Protein.Group[mice_R$seqcharge == mice_R_$seqcharge[i]] <- mice_R_$Protein.Group[i]
    }
  }
  mice_R <- reshape2::dcast(mice_R,Protein.Group+seqcharge ~ ID, value.var = 'Ms1.Area')


  Prot_pep_mapR <- as.data.frame(cbind(mice_R$Protein.Group,mice_R$seqcharge))
  colnames(Prot_pep_mapR) <- c('Protein','seqcharge')

  rownames(mice_R) <- mice_R$seqcharge
  mice_R$seqcharge <- NULL
  mice_R$Protein.Group <- NULL

  mice_K$Iso <- str_sub(mice_K$Precursor.Id,-3,-3)
  mice_K$Iso[mice_K$Iso == "0" | mice_K$Iso == "2"] <- 'L'
  mice_K$Iso[mice_K$Iso == "1" | mice_K$Iso == "3"] <- 'H'

  mice_K_count <- mice_K %>% group_by(ID,seqcharge) %>% summarise(numbObs = sum(is.na(Ms1.Area)==F))
  mice_K_count <- mice_K_count %>% filter(numbObs == 2)


  mice_K <- mice_K %>% filter(seqcharge %in% mice_K_count$seqcharge)

  mice_K_H <- mice_K %>% filter(Iso == 'H')

  mice_K_H_ <- reshape2::dcast(mice_K_H,Protein.Group+seqcharge ~ ID, value.var = 'Ms1.Area')

  n_occur <- data.frame(table(mice_K_H_$seqcharge))
  n_occur <- n_occur[n_occur$Freq > 1,]

  if(nrow(n_occur) > 0){
    mice_K_H_ <- mice_K_H_ %>% filter(seqcharge %in% n_occur$Var1)
    mice_K_H_ <- mice_K_H_ %>% distinct(seqcharge,.keep_all = T)

    for(i in 1:nrow(mice_K_H_)){
      mice_K_H$Protein.Group[mice_K_H$seqcharge == mice_K_H_$seqcharge[i]] <- mice_K_H_$Protein.Group[i]
    }
  }


  mice_K_H <- reshape2::dcast(mice_K_H,Protein.Group+seqcharge ~ ID, value.var = 'Ms1.Area')



  mice_K_L <- mice_K %>% filter(Iso == 'L')
  mice_K_L_ <- reshape2::dcast(mice_K_L,Protein.Group+seqcharge ~ ID, value.var = 'Ms1.Area')

  n_occur <- data.frame(table(mice_K_L_$seqcharge))
  n_occur <- n_occur[n_occur$Freq > 1,]
  if(nrow(n_occur) > 0){
    mice_K_L_ <- mice_K_L_ %>% filter(seqcharge %in% n_occur$Var1)
    mice_K_L_ <- mice_K_L_ %>% distinct(seqcharge,.keep_all = T)

    for(i in 1:nrow(mice_K_L_)){
      mice_K_L$Protein.Group[mice_K_L$seqcharge == mice_K_L_$seqcharge[i]] <- mice_K_L_$Protein.Group[i]
    }

  }

  mice_K_L <- reshape2::dcast(mice_K_L,Protein.Group+seqcharge ~ ID, value.var = 'Ms1.Area')


  Prot_pep_mapK <- as.data.frame(cbind(mice_K_L$Protein.Group,mice_K_L$seqcharge))
  colnames(Prot_pep_mapK) <- c('Protein','seqcharge')

  rownames(mice_K_L) <- mice_K_L$seqcharge
  mice_K_L$seqcharge <- NULL
  mice_K_L$Protein.Group <- NULL
  mice_K_L <- as.matrix(mice_K_L)

  rownames(mice_K_H) <- mice_K_H$seqcharge
  mice_K_H$seqcharge <- NULL
  mice_K_H$Protein.Group <- NULL
  mice_K_H <- as.matrix(mice_K_H)

  sect_col <- intersect(colnames(mice_K_H),colnames(mice_K_L))
  #sect_row <- intersect(rownames(mice_K_H),rownames(mice_K_L))
  mice_K_H <- mice_K_H[,sect_col]
  mice_K_L <- mice_K_L[,sect_col]
  #Prot_pep_mapK <- Prot_pep_mapK[sect_row,]
  mice_R <- mice_R[,sect_col]

  mice_K_H_raw = mice_K_H
  mice_K_L_raw = mice_K_L

  mice_K_all <- mice_K_H+mice_K_L

  mice_all_pep <- as.matrix(rbind(mice_R,mice_K_all))

  std_matricies <- new('matricies_DIA',peptide = mice_all_pep,peptide_protein_map = rbind(Prot_pep_mapR,Prot_pep_mapK))

  QQC@matricies <- std_matricies


  mice_K_H_ov_L <- mice_K_H/mice_K_L

  mice_alpha <- log(mice_K_H/mice_K_L+1)/t

  mice_K_all <- QuantQC::normalize(mice_K_all)
  mice_K_H <- mice_K_H/(mice_K_H+mice_K_L)
  mice_K_L <- 1 - mice_K_H
  mice_K_H <- mice_K_all*mice_K_H
  mice_K_L <- mice_K_all*mice_K_L

  mice_beta <- mice_K_H*mice_alpha/(1-exp(-mice_alpha*t)) # Size adjusted translation rate

  mice_beta <- QuantQC::normalize(mice_beta)

  miceotope_matricies <- new('matricies_Miceotopes',Raw_H = as.matrix(mice_K_H_raw),Raw_L = as.matrix(mice_K_L_raw), HovL_pep = as.matrix(mice_K_H_ov_L),Beta_pep = as.matrix(mice_beta),Alpha_pep = as.matrix(mice_alpha), peptide_protein_map=Prot_pep_mapK)


  QQC@miceotopes <- miceotope_matricies


  return(QQC)

}


MicePepCorPlot <- function(QQC){

  Prot_pep_mapK <- QQC@miceotopes@peptide_protein_map

  mice_beta <- QQC@miceotopes@Beta_pep
  mice_alpha <- QQC@miceotopes@Alpha_pep

  PepCorAlpha <- PepCorMini(mice_alpha,Prot_pep_mapK)
  PepCorBeta  <- PepCorMini(mice_beta,Prot_pep_mapK)

  Beta_plot <- ggplot(PepCorBeta,aes(x = Cor)) + geom_histogram() +
    xlab('Correlations') + ylab('# proteins')+
    ggtitle('Beta value')
  Alpha_plot <- ggplot(PepCorAlpha,aes(x = Cor)) + geom_histogram() +
    xlab('Correlations') + ylab('# proteins') +
    ggtitle('Alpha value')

  return(Beta_plot+Alpha_plot)

}


Miceotope_protein_collapse <- function(QQC){

  meta <- QQC@meta.data

  good_cells <- colnames(QQC@matricies@peptide)

  mice_K_H_ov_L <-  QQC@miceotopes@HovL_pep[,good_cells]
  mice_beta <-QQC@miceotopes@Beta_pep[,good_cells]
  mice_alpha <- QQC@miceotopes@Alpha_pep[,good_cells]

  Prot_pep_mapK <- QQC@miceotopes@peptide_protein_map

  mice_alpha_prot <- as.data.frame(mice_alpha)
  mice_alpha_prot$prot <- Prot_pep_mapK$Protein
  mice_alpha_prot <- melt(mice_alpha_prot, ids = c('prot'))
  mice_beta_prot <- as.data.frame(mice_beta)
  mice_beta_prot$prot <- Prot_pep_mapK$Protein
  mice_beta_prot <- melt(mice_beta_prot, ids = c('prot'))
  mice_K_H_ov_L_prot <- as.data.frame(mice_K_H_ov_L)
  mice_K_H_ov_L_prot$prot <- Prot_pep_mapK$Protein
  mice_K_H_ov_L_prot <- melt(mice_K_H_ov_L_prot, ids = c('prot'))

  mice_beta_prot <- mice_beta_prot %>% group_by(variable,prot) %>% summarise(value = median(value,na.rm=T))
  mice_alpha_prot <- mice_alpha_prot %>% group_by(variable,prot) %>% summarise(value = median(value,na.rm=T))
  mice_K_H_ov_L_prot <- mice_K_H_ov_L_prot %>% group_by(variable,prot) %>% summarise(value = median(value,na.rm=T))


  mice_beta_prot <- reshape2::dcast(mice_beta_prot,prot ~ variable, value.var = 'value')
  mice_alpha_prot <- reshape2::dcast(mice_alpha_prot,prot ~ variable, value.var = 'value')
  mice_K_H_ov_L_prot <- reshape2::dcast(mice_K_H_ov_L_prot,prot ~ variable, value.var = 'value')

  rownames(mice_K_H_ov_L_prot) <- mice_K_H_ov_L_prot$prot
  mice_K_H_ov_L_prot$prot <- NULL

  rownames(mice_beta_prot) <- mice_beta_prot$prot
  mice_beta_prot$prot <- NULL

  mice_beta_prot <- log2(mice_beta_prot)

  rownames(mice_alpha_prot) <- mice_alpha_prot$prot
  mice_alpha_prot$prot <- NULL


  QQC@miceotopes@HovL_prot <- as.matrix(mice_K_H_ov_L_prot)
  QQC@miceotopes@Beta_prot <- as.matrix(mice_beta_prot)
  QQC@miceotopes@Alpha_prot <- as.matrix(mice_alpha_prot)


  Turnover_all <- colMedians(as.matrix(mice_K_H_ov_L_prot),na.rm = T)
  Turnover_all_df <- as.data.frame(Turnover_all)
  Turnover_all_df$ID <- colnames(as.matrix(mice_K_H_ov_L_prot))

  meta <- meta %>% left_join(Turnover_all_df, by = c('ID'))


  QQC@meta.data <- meta


  return(QQC)


}


Mice_DimPlot_turnover <-function(QQC, reuct = 'PCA', by = 'Total'){

  meta <- QQC@meta.data
  meta <- meta %>% dplyr::select(ID,Turnover_all)


  if(by == 'Total'){

    if(reuct == 'PCA'){
      reduct <-QQC@reductions[['PCA']]
      reduct$ID <- rownames(reduct)
      reduct <- reduct %>% left_join(meta, by = c('ID'))
      dim_plot <- ggplot(reduct,aes(x = PC1, y = PC2,color = Turnover_all)) + geom_point()+
        um_plot + scale_color_gradient2(midpoint = median(reduct$Turnover_all), low = 'blue',mid = 'white', high = 'red')


    }

    if(reuct == 'UMAP'){

      reduct <- QQC@reductions[['UMAP']]
      reduct$ID <- rownames(reduct)
      reduct <- reduct %>% left_join(meta, by = c('ID'))
      dim_plot <- ggplot(reduct,aes(x = umap_1, y = umap_2,color = Turnover_all)) + geom_point()+
        um_plot + scale_color_gradient2(midpoint = median(reduct$Turnover_all), low = 'blue',mid = 'white', high = 'red')


    }

  }else{

    mice_K_H_ov_L_prot <- QQC@miceotopes@HovL_prot
    Prot_turnover <- mice_K_H_ov_L_prot[by,]
    Prot_turnover[Prot_turnover==Inf] <- NA
    Prot_turnover_df <- as.data.frame(Prot_turnover)
    colnames(Prot_turnover_df) <- 'Prot'

    Prot_turnover_df$ID <- colnames(mice_K_H_ov_L_prot)

    if(reuct == 'PCA'){
      reduct <-QQC@reductions[['PCA']]
      reduct <- reduct %>% left_join(Prot_turnover_df, by = c('ID'))
      dim_plot <- ggplot(reduct,aes(x = PC1, y = PC2,color = log2(Prot))) + geom_point()+ggtitle(by)+
        um_plot + scale_color_gradient2(midpoint = median(log2(reduct$Prot),na.rm = T), low = 'blue',mid = 'white', high = 'red')
    }

    if(reuct == 'UMAP'){

      reduct <-QQC@reductions[['UMAP']]
      reduct$ID <- rownames(reduct)
      reduct <- reduct %>% left_join(Prot_turnover_df, by = c('ID'))
      dim_plot <- ggplot(reduct,aes(x = umap_1, y = umap_2,color = log2(Prot)))+ geom_point()+ggtitle(by)+
        um_plot + scale_color_gradient2(midpoint = median(log2(reduct$Prot),na.rm = T), low = 'blue',mid = 'white', high = 'red')
    }
  }

  return(dim_plot)


}


