




### Peptide correlations


#' Add two numbers.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
SharedPeptideCor <- function(QQC, res = 'sc'){
  peptide_data <- QQC@matricies@peptide
  protein_dat <- QQC@matricies@protein
  peptide_protein_map <- QQC@matricies@peptide_protein_map



  #peptide_data <- Normalize_reference_vector(peptide_data,log = T)

  if(res == 'clust'){

    clusts <- QQC@reductions$UMAP$cluster

    peptide_data_clust <- matrix(data = NA, nrow = nrow(peptide_data),ncol = length(unique(clusts)))
    count <- 0
    for(i in unique(clusts)){
      count <- count + 1
      peptide_data_clust[,count] <- rowMeans(peptide_data[,clusts == i],na.rm = T)

    }

    peptide_data <- peptide_data_clust

  }

  # Initialized empty matricies to store correlations between peptides
  mat_stor = matrix(data = NA,nrow = 10000,ncol = 3)
  mat_stor_fake = matrix(data = NA,nrow = 10000,ncol = 2)

  # List of all the unique proteins in data
  upl <- unique(peptide_protein_map$Protein)

  # List describing which protein each peptide comes from
  prot_list <- peptide_protein_map$Protein


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


  # Shuffle around the protein list so its correlating peptides from different proteins to eachother
  # This should form our "NULL" distribution
  prot_list_mixed <- sample(peptide_protein_map$Protein)

  # Same loop as above but on missmatched peptides
  for(i in upl){
    mat_p1 <- peptide_data[which(prot_list_mixed == i),]
    if(is.null(nrow(mat_p1)) == FALSE ){
      obs_mat <- psych::pairwiseCount(t(mat_p1))
      cor_mat <- cor(t(mat_p1),use = 'pairwise.complete.obs')
      cor_mat[obs_mat < 4] <- NA
      cor_mat <- cor_mat[lower.tri(cor_mat)]
      if(is.na(median(cor_mat,na.rm = T)) == F){
        count = count + 1
        mat_stor_fake[count,1] <- i
        mat_stor_fake[count,2] <- median(cor_mat,na.rm = T)
      }
    }
  }

  # Format storage matrix
  mat_stor <- as.data.frame(mat_stor)
  colnames(mat_stor) <- c('Protein','Cor','Obs')
  mat_stor <- mat_stor %>% filter(is.na(Cor) == F)
  mat_stor$Cor <- as.numeric(mat_stor$Cor)

  # Store a protein matrix with of protein quant for proteins with multiple peptides
  mat_stor <- mat_stor %>% filter(Protein %in% rownames(protein_dat))
  rownames(mat_stor) <- mat_stor$Protein
  protein_mat_mult_pep <- protein_dat[mat_stor$Protein,]
  pep_cor <- mat_stor[rownames(protein_mat_mult_pep),]

  # Binning protein abundance by fold change for faceting correlations
  # Peptides from proteins that change more should correlate better (more signal = better signal/noise)

  pep_cor$FC <- rowMeans(abs(protein_mat_mult_pep),na.rm = T)
  pep_cor$FC[pep_cor$FC <= .4] <- .4
  pep_cor$FC[pep_cor$FC > .4 & pep_cor$FC <= .8] <- .8
  pep_cor$FC[pep_cor$FC > .8 & pep_cor$FC <= 1.2] <- 1.2
  pep_cor$FC[pep_cor$FC > 1.2 & pep_cor$FC <= 2] <- 2
  pep_cor$FC[pep_cor$FC > 2] <- 3
  pep_cor$FC <- as.character(pep_cor$FC)


  QQC@pep.cor <- list(pep_cor,median(as.numeric(mat_stor_fake[,2]),na.rm = T))

  return(QQC)

}

#' Add two numbers.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
PlotPepCor <- function(QQC,type = 'box'){

  pep_cor <- QQC@pep.cor[[1]]
  null_dist <- QQC@pep.cor[[2]]

  if(type == 'box'){

    plot_ <- ggplot(pep_cor, aes(y = Cor, x = FC)) + geom_boxplot(color="black", fill = 'gray') + xlab('Mean abs(protein fold change)') +
      ylab('Cor.; peptides mapping to protein')+
      stat_summary(fun.data=f, geom="text", vjust=-0.5, col="blue")+
      geom_hline(yintercept = null_dist, col = "red")+dot_plot

    return(plot_)

  }
  if(type == 'hist'){

    plot_ <- ggplot(pep_cor, aes(x = Cor)) + geom_histogram(color="black", fill = 'gray') + ylab('# Proteins') +
      xlab('Correlation; peptides mapping to protein')+
      geom_vline(xintercept = null_dist, col = "red")+dot_plot

    return(plot_)


  }






}






#####
# Counting Peptide and protein numbers
####

#' Add two numbers.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
inSet_completness <- function(Raw_data, cellenOne_meta){
  count = 0
  for(i in unique(paste0(cellenONE_meta$injectWell,cellenONE_meta$plate))){

    Data_inset <- as.data.frame(Raw_cell_mat)[,which(grepl(i,colnames(Raw_cell_mat)))]

    Data_inset <- Data_inset[which(is.na(rowMeans(Data_inset,na.rm = T))==F),]

    One_set <- Data_completeness(Data_inset)

    if(count == 0){
      Cell_complete <- One_set[[1]]
      #Protein_complete <- One_set[[2]]
    }else{
      Cell_complete <- rbind(Cell_complete,One_set[[1]])
      #Protein_complete <- rbind(Cell_complete,One_set[[2]])
    }

    count = count+1
  }

  return(Cell_complete)
}



### General Data Analysis functions








#' Add two numbers.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
PlotProtAndPep <- function(QQC){
  # count peptide numbers
  numb_pep <- colSums(is.na(QQC@matricies@peptide)==F)
  numb_pep <- as.data.frame(numb_pep)
  colnames(numb_pep) <- 'Number_precursors'




  # count protein numbers
  numb_prot <- colSums(is.na(QQC@matricies@protein)==F)
  numb_prot <- as.data.frame(numb_prot)
  colnames(numb_prot) <- 'Number_proteins'



  # Plot peptide and protein numbers
  if(QQC@ms_type == 'DDA'){
    numb_pep$ID <- colnames(QQC@matricies@peptide)
    numb_prot$ID <- colnames(QQC@matricies@protein)
    numb_pep <- numb_pep %>% left_join(QQC@meta.data, by = c('ID'))
    numb_prot <- numb_prot %>% left_join(QQC@meta.data, by = c('ID'))

    #pep_number <- ggplot(numb_pep, aes(y = Number_precursors, x = sample)) + geom_boxplot() + ggtitle('# precursors per sample') + rremove('legend')+ylab('# Peptides')+dot_plot + ylim(c((min(numb_pep$Number_precursors)-300)),(max(numb_pep$Number_precursors)+300))

    #prot_number<- ggplot(numb_prot, aes(y = Number_proteins, x = sample)) + geom_boxplot() + ggtitle('# proteins per sample') + ylab('# Proteins')+dot_plot + ylim(c((min(numb_prot$Number_proteins)-200)),(max(numb_prot$Number_proteins)+200))



    pep_number <- ggplot(numb_pep, aes(x = Number_precursors)) + geom_histogram(bins = 15) + ggtitle('# precursors per sample') + rremove('legend')+ylab('# of single cells')+dot_plot + xlim(c((min(numb_pep$Number_precursors)-300)),(max(numb_pep$Number_precursors)+300))

    prot_number<- ggplot(numb_prot, aes(x = Number_proteins)) + geom_histogram(bins = 15) + ggtitle('# proteins per sample') + ylab('# of single cells')+dot_plot + xlim(c((min(numb_prot$Number_proteins)-200)),(max(numb_prot$Number_proteins)+200))


    return(pep_number+prot_number)
  }

  if(QQC@ms_type == 'DIA' | QQC@ms_type == 'DIA_C'){

    numb_pep_NF <- colSums(QQC@matricies@peptide_mask==T)
    numb_pep_NF <- as.data.frame(numb_pep_NF)
    colnames(numb_pep_NF) <- 'Number_precursors'
    numb_pep_NF$Filter <- 'No Ch Qvalue Filter'
    numb_pep$Filter <- paste0(QQC@misc[['ChQ']], ' Channel Q Value')

    numb_pep <- rbind(numb_pep,numb_pep_NF)

    # # count protein numbers
    numb_prot_NF <- colSums(QQC@matricies@protein_mask==T)
    numb_prot_NF <- as.data.frame(numb_prot_NF)
    colnames(numb_prot_NF) <- 'Number_proteins'
    numb_prot_NF$Filter <- 'No Ch Qvalue Filter'
    numb_prot$Filter <- paste0(QQC@misc[['ChQ']], ' Channel Q Value')

    numb_prot <- rbind(numb_prot,numb_prot_NF)

    pep_number <- ggplot(numb_pep, aes(x = Number_precursors, fill = Filter)) + geom_histogram(bins = 30,position = 'identity',alpha = .5) + ggtitle('# precursors per sample') + rremove('legend')+ylab('# of single cells')+dot_plot

    prot_number<- ggplot(numb_prot, aes(x = Number_proteins,fill = Filter)) + geom_histogram(bins = 30,position = 'identity',alpha = .5) + ggtitle('# proteins per sample') + ylab('# of single cells')+dot_plot


    return(pep_number+prot_number)
  }

}

#' Add two numbers.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
PlotDataComplete <- function(QQC){

  data <- QQC@matricies@protein

  missingness_cell_filt <- colSums(is.na(data) == F)/nrow(data)
  missingness_cell_mat <- reshape2::melt(missingness_cell_filt)

  missingness_prot_filt <- rowSums(is.na(data)==F)/ncol(data)
  missingness_prot_mat <- reshape2::melt(missingness_prot_filt)

  if(QQC@ms_type == 'DDA'){
    mp <- ggplot(missingness_prot_mat, aes(x = value)) +
      geom_histogram(bins = 20,position = 'identity',alpha = .5) + ggtitle(paste0('Protein completness, ', nrow(data) ,' proteins'))+rremove('legend') +ylab('# of proteins')+xlab('fraction values present')+dot_plot

    mc <- ggplot(missingness_cell_mat, aes(x = value)) +
      geom_histogram(bins = 20,position = 'identity',alpha = .5) + ggtitle('Cell completness') + ylab('# of single cells')+
      xlab('fraction values present')+dot_plot


    return(mp+mc)
  }


  if(QQC@ms_type == 'DIA' | QQC@ms_type == 'DIA_C'){

    missingness_cell_NF <- colSums(QQC@matricies@protein_mask)/nrow(data)
    missingness_cell_mat_NF <- reshape2::melt(missingness_cell_NF)
    colnames(missingness_cell_mat_NF) <- 'Cell_miss'
    missingness_cell_mat_NF$Filter <- 'No Ch Qvalue Filter'
    colnames(missingness_cell_mat) <- 'Cell_miss'
    missingness_cell_mat$Filter <- paste0(QQC@misc[['ChQ']], ' Channel Q Value')


    missingness_cell_mat <- rbind(missingness_cell_mat,missingness_cell_mat_NF)

    # count protein numbers
    missingness_prot_NF <- rowSums(QQC@matricies@protein_mask)/ncol(data)
    missingness_prot_mat_NF <- reshape2::melt(missingness_prot_NF)
    colnames(missingness_prot_mat_NF) <- 'Miss_proteins'
    missingness_prot_mat_NF$Filter <- 'No Ch Qvalue Filter'
    colnames(missingness_prot_mat) <- 'Miss_proteins'
    missingness_prot_mat$Filter <- paste0(QQC@misc[['ChQ']], ' Channel Q Value')

    missingness_prot_mat <- rbind(missingness_prot_mat,missingness_prot_mat_NF)

    mp <- ggplot(missingness_cell_mat, aes(x = Cell_miss, fill = Filter)) +
      geom_histogram(bins = 20,position = 'identity',alpha = .5) + ggtitle('Cell completness') + ylab('# of single cells')+
      xlab('fraction values present')+dot_plot
    mc <-  ggplot(missingness_prot_mat, aes(x = Miss_proteins,fill = Filter)) +
      geom_histogram(bins = 20,position = 'identity',alpha = .5) + ggtitle(paste0('Protein completness, ', nrow(data) ,' proteins'))+rremove('legend') +ylab('# of proteins')+xlab('fraction values present')+dot_plot



    return(mp+mc)
  }








}





#' Add two numbers.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
PlotSCtoCarrierRatio <- function(QQC){

  do_plot <- "No carrier used"

  if(QQC@ms_type == "DDA"){

    sc.data <- QQC@raw_data
    good_cells <- colnames(QQC@matricies@peptide)

    cellenOne_meta <- QQC@meta.data
    good_cells_p_negs <- c(good_cells,cellenOne_meta$ID[cellenOne_meta$sample == 'neg'])

    ri.index<-which(colnames(sc.data)%in%paste0("Reporter.intensity.",1:18))

    sc.data[sc.data==0] <- NA
    sc.data[, ri.index] <- sc.data[, ri.index] / sc.data[, ri.index[1]]


    sc.data <- sc.data[,c('Raw.file','Well','plate',paste0("Reporter.intensity.",5:18))]
    sc.data <- reshape2::melt(sc.data, id = c('Raw.file','Well','plate'))

    sc.data$ID <- paste0(sc.data$Well,sc.data$plate,sc.data$variable)
    sc.data <- sc.data %>% filter(ID %in% good_cells_p_negs)
    sc.data <- sc.data %>% group_by(ID) %>% dplyr::summarise(value = median(value,na.rm = T))

    sc.data <- sc.data %>% dplyr::left_join(cellenOne_meta, by = c('ID'))

    sc.data$value <- log10(sc.data$value)#1/sc.data$value

    do_plot <- ggplot(sc.data, aes(x = sample, y = value)) + geom_boxplot() + ylab('Single cells / Carrier') +dot_plot

  }

  if(QQC@ms_type == "DIA_C"){

    cellenOne_meta <- QQC@meta.data
    Filt_ratio <- QQC@matricies@peptide #%>% dplyr::select(good_cells$variable)
    Filt_ratio <- reshape2::melt(Filt_ratio)
    Filt_ratio_cell <- Filt_ratio %>% group_by(Var2) %>% dplyr::summarise(value_new = median(value,na.rm = T))
    Filt_ratio_cell <- Filt_ratio_cell %>% left_join(cellenOne_meta,by = c('Var2'='ID'))


    Filt_rat_plot <- ggplot(Filt_ratio_cell, aes(y = 1/(value_new), x =  sample)) + geom_boxplot() +
      ylab('# single cells in carrier') + ggtitle( paste0('Carrier is median ', round(median(1/Filt_ratio$value,na.rm = T),2),' single cells'))

    do_plot <- Filt_rat_plot


  }

  do_plot

}


#' Add two numbers.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
PlotDigestEff <- function(QQC){
  do_plot <- "No carrier used"

  if(QQC@ms_type == "DDA"){


    sc.data <- QQC@matricies@peptide
    peps <- QQC@matricies@peptide_protein_map$seqcharge
    peps <- str_sub(peps,1,-4)

    vect_MC <- ContainsMissedCleaved(peps)

    for(i in 1:ncol(sc.data)){
      sc.data[,i] <- sc.data[,i]/median(sc.data[,i],na.rm=T)

    }

    sc.data.reg <- sc.data[vect_MC == 0,]
    sc.data.reg <- reshape2::melt(sc.data.reg)
    sc.data.reg$cond <- 'Regular'

    sc.data.MC <- sc.data[vect_MC != 0,]
    sc.data.MC <- reshape2::melt(sc.data.MC)
    sc.data.MC$cond <- 'Missed Cleaved'

    sc.data <- rbind(sc.data.reg,sc.data.MC)

    sc.data <- sc.data %>% filter(is.na(value)==F)

    dig_eff <- round(median(sc.data.reg$value,na.rm = T)/median(sc.data.MC$value,na.rm = T),2)

    do_plot <- ggplot(sc.data,aes(y = log2(value), x = cond))+
      geom_boxplot() + ggtitle(paste0(dig_eff,'X digest efficiency of carrier'))+
      ylab('log2(Intensity vs Carrier)')+xlab('')+dot_plot

  }

  if(QQC@ms_type == "DIA_C"){


    sc.data <- QQC@matricies@peptide
    peps <- QQC@matricies@peptide_protein_map$seqcharge
    peps <- str_sub(peps,1,-3)



    vect_MC <- ContainsMissedCleaved(peps)

    for(i in 1:ncol(sc.data)){
      sc.data[,i] <- sc.data[,i]/median(sc.data[,i],na.rm=T)

    }

    sc.data.reg <- sc.data[vect_MC == 0,]
    sc.data.reg <- reshape2::melt(sc.data.reg)
    sc.data.reg$cond <- 'Reg'

    sc.data.MC <- sc.data[vect_MC != 0,]
    sc.data.MC <- reshape2::melt(sc.data.MC)
    sc.data.MC$cond <- 'MC'

    sc.data <- rbind(sc.data.reg,sc.data.MC)

    sc.data <- sc.data %>% filter(is.na(value)==F)

    dig_eff <- round(median(sc.data.reg$value,na.rm = T)/median(sc.data.MC$value,na.rm = T),2)

    do_plot <- ggplot(sc.data,aes(y = log2(value), x = cond))+
      geom_boxplot() + ggtitle(paste0(dig_eff,'X digest efficiency of carrier'))+
      ylab('log2(Intensity vs Carrier)')+xlab('')+dot_plot

  }

  do_plot
}




#' Add two numbers.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
PlotCellSizeVsIntensity <- function(QQC, type = 'sample'){
  meta <- QQC@meta.data

  good_cells <- colnames(QQC@matricies@peptide)
  meta <- meta %>% dplyr::filter(ID %in% good_cells)

  title_text <- paste0('Correlation = ', round(cor(meta$diameter,meta$prot_total,use = 'pairwise.complete.obs'),2))

  if(type == 'sample'){
    plot_ <-    ggplot(meta, aes(x = log2((diameter/2)^3),y = log2(10^prot_total),color = sample)) +
      geom_point() + ggtitle(title_text)+ ylab('log2(Sum cell intensity)')+ xlab('log2(Vol.) cubic uM')+
      dot_plot

    return(plot_)
  }

  if(type == 'Run order'){


    reg_ord <- coef(lm(prot_total~Order,data = meta))
    meta$prot_total2 <- meta$prot_total - meta$Order*reg_ord[2]

    plot_ <-    ggplot(meta, aes(x = log2((diameter/2)^3),y = log2(10^prot_total),color = Order)) +
      geom_point() + ggtitle(title_text)+ ylab('log2(Sum cell intensity)')+ xlab('log2(Vol.) cubic uM')+
      scale_color_gradient2(midpoint = median(meta$Order,na.rm = T), low = 'blue',mid = 'white', high = 'red')+
      dot_plot

    title_text2 <- paste0('adj cor = ', round(cor((meta$diameter/2)^3,10^meta$prot_total2,use = 'pairwise.complete.obs'),2))

    plot_2 <-    ggplot(meta, aes(x = log2((diameter/2)^3),y = log2(10^prot_total2),color = Order)) +
      geom_point() + ggtitle(title_text2)+ylab('log2(Sum cell intens)')+ xlab('log2(Vol.) cubic uM')+
      scale_color_gradient2(midpoint = median(meta$Order,na.rm = T), low = 'blue',mid = 'white', high = 'red')+
      dot_plot

    return(plot_+plot_2)
  }


}




#' Add two numbers.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export
ProteinClustConsistency <- function(QQC, prot = NA, type = 'line',fasta_path = NA){

  if(is.na(fasta_path) == T){
    if(QQC@misc[['Species']] == 'Mouse'){
      fasta_seq <- read.fasta(system.file("extdata", "Mouse.fasta", package = "QuantQC"),seqtype =  "AA",as.string = T)
    }
    if(QQC@misc[['Species']] == 'Human'){
      fasta_seq <- read.fasta(system.file("extdata", "Human.fasta", package = "QuantQC"),seqtype =  "AA",as.string = T)
    }
  }

  if(is.na(fasta_path) == F){
    fasta_seq <- read.fasta(fasta_path,seqtype =  "AA",as.string = T)

  }

  names(fasta_seq) <- sapply(names(fasta_seq), extract_accession)

  fasta_seq_prot <- fasta_seq[[prot]][[1]]

  if(is.null(fasta_seq_prot)==T){
    return('Protein not found in fasta')
  }

  protein_sequence <- data.frame(
    sequence = 1:nchar(fasta_seq_prot),
    color = 'black'
  )

  if(is.null(QQC@reductions[['UMAP']])==F){



    if(QQC@ms_type =='DDA'){
      raw_pep_intense <- QQC@raw_data
      raw_pep_intense <- raw_pep_intense %>% filter(Leading.razor.protein == prot)
      raw_pep_intense <- raw_pep_intense %>% group_by(seqcharge) %>% summarise(pep_raw = median(Reporter.intensity.1,na.rm=T))
      raw_pep_intense <- raw_pep_intense[order(-raw_pep_intense$pep_raw),]

    }
    if(QQC@ms_type =='DIA' | QQC@ms_type =='DIA_C'){
      raw_pep_intense <- QQC@raw_data
      raw_pep_intense <- raw_pep_intense %>% filter(Protein.Group == prot)
      raw_pep_intense <- raw_pep_intense %>% group_by(seqcharge) %>% summarise(pep_raw = median(Ms1.Area,na.rm=T))
      raw_pep_intense <- raw_pep_intense[order(-raw_pep_intense$pep_raw),]

    }





    # Get 4 best peptides


    clusters <- QQC@reductions[['UMAP']]
    prot_map <- QQC@matricies@peptide_protein_map %>% filter(Protein == prot)

    raw_pep_intense_prot <- raw_pep_intense %>% filter(seqcharge %in% prot_map$seqcharge)
    if(nrow(raw_pep_intense_prot)>4){
      raw_pep_intense_prot <- raw_pep_intense_prot[1:4,]
    }
    if(nrow(raw_pep_intense_prot)==1){
      return('Only 1 peptide')
    }

    prot_map <- prot_map %>% filter(seqcharge %in% raw_pep_intense_prot$seqcharge)
    pep_mat <- (Normalize_reference_vector(QQC@matricies@peptide,log = T))
    rownames(pep_mat) <- QQC@matricies@peptide_protein_map$seqcharge

    pep_mat <- pep_mat[prot_map$seqcharge,]
    if(nrow(pep_mat) > 2){
      pep_mat_cor <- cor(t(pep_mat),use = 'pairwise.complete.obs')
      pep_mat_cor <- round(median(pep_mat_cor[lower.tri(pep_mat_cor)],na.rm = T),3)

    }
    if(nrow(pep_mat) == 2){
      pep_mat_cor <- round(cor(pep_mat,use = 'pairwise.complete.obs'))
    }

    pep_mat <- reshape2::melt(pep_mat)
    clusters$ID <- rownames(clusters)
    pep_mat <- pep_mat %>%  dplyr::left_join(clusters, by = c('Var2' = 'ID'))

    pep_mat_values <-  pep_mat %>% dplyr::group_by(Var1,cluster) %>%
      dplyr::summarise(med_abs = median(value,na.rm = T))

    pep_mat_SD<- pep_mat %>%  dplyr::group_by(Var1,cluster) %>%
      dplyr::summarise(sd_error = sd(value,na.rm = T)/sum(is.na(value)==F))

    pep_mat_count <- pep_mat %>%  dplyr::group_by(Var1,cluster) %>%
      dplyr::summarise(numb_data = sum(is.na(value)==F)/(sum(is.na(value)==F)+sum(is.na(value)==T)))

    prot_values <- pep_mat %>%  dplyr::group_by(cluster) %>%
      dplyr::summarise(med_abs = median(value,na.rm = T))



    for(i in 1:nrow(raw_pep_intense_prot)){
      prot_values_pep1 <- prot_values
      prot_values_pep1$Var1 <- unique(pep_mat_values$Var1)[i]
      if(i == 1){
        prot_values <- prot_values_pep1
      }else{
        prot_values <- rbind(prot_values,prot_values_pep1)
      }

    }

    prot_values$sd_error <- 0
    prot_values$numb_dp <- 0
    prot_values$numb_dp_z <- 0
    prot_values$PctCellsMeasured <- NA

    pep_mat_values$sd_error <- pep_mat_SD$sd_error
    pep_mat_values$numb_dp <- pep_mat_count$numb_data
    pep_mat_values$numb_dp_z <- scale(pep_mat_values$numb_dp)[,1]

    pep_mat_values$Var3 <- pep_mat_values$Var1
    prot_values$Var3 <- ' Protein'


    # Transform z-scores into quantiles
    pep_mat_values$PctCellsMeasured <- cut(pep_mat_values$numb_dp, breaks = 3,
                                 labels = c( " Less 33%", "Greater 33%", "Greater 66%"))



    all <- rbind(pep_mat_values,prot_values)

    #pep_mat_values_test <- pep_mat_values

    #pep_mat_values_test$prote <- prot_values$med_abs


    all$cluster <- (as.character(all$cluster))

    if(type == 'line'){
      if(length(all$Var1) > 0){

        create_peptide_vector

        plot_ <- ggplot(all, aes(x = cluster, y = med_abs, color = Var3, group = Var3)) +
          geom_line( aes(color = Var3)) +
          geom_point(aes(size = PctCellsMeasured,color = Var3)) +
          scale_size_manual(values = c(" Less 33%" = 2, "Greater 33%" = 4, "Greater 66%" = 6)) +
          geom_errorbar(aes(ymin = med_abs - sd_error, ymax = med_abs + sd_error, color = Var3), width = 0.2) +
          labs(x = "X-axis", y = "Y-axis") +
          ggtitle(paste0("Individual peptide vs all averaged peptides across clusters")) +
          theme_bw() + scale_colour_manual(values = c("black","red", "blue", "green","purple"))+
          xlab('Clusters') +
          ggplot2::facet_wrap(~Var1, scales = "free_y") +
          theme(strip.background = element_blank())+ylab('Log2(FoldChange)')+guides(color = "none")

        col_vect <- c("red", "blue", "green","purple")
        count = 0
        for(i in unique(all$Var1)){
          count = count+1
          peptide <- gsub('_','',i)
          peptide <- gsub('[[:digit:]]+', '', peptide)

          protein_sequence$color[create_peptide_vector(fasta_seq_prot, peptide) == 1] <- col_vect[count]
        }

        prot_for_plot <- ggplot(protein_sequence, aes(xmin = sequence - 0.5, xmax = sequence + 0.5, ymin = -0.2, ymax = 0.2, fill = color)) +
          geom_rect() +
          coord_cartesian(ylim = c(-0.5, 0.5)) +
          scale_fill_identity() +
          theme_void()+ggtitle('Peptide locations across protein')

        return((prot_for_plot/plot_)+plot_layout(heights = c(1,10)) + plot_annotation(paste0(prot,' Median cor. = ',pep_mat_cor)))
      }
    }

    if(type == 'bubble'){

      all$var4 <- NA
      all$var4[all$Var3 == 'Protein'] <- 'Protein'
      all$var4[all$Var3 != 'Protein'] <- 'Peptide'

      all$FractionCellsExpress[all$Var3 == 'Protein'] <- 'Greater 66%'
      all$FractionCellsExpress[all$FractionCellsExpress == 'Less 33%'] <- ' Less 33%'

      plot_ <- ggplot(all, aes(x = cluster, y = Var3, fill = med_abs)) +
        geom_point(aes(size = FractionCellsExpress), shape = 21) + theme_bw() +
        scale_size_manual(values = c(" Less 33%" = 3, "Greater 33%" = 5, "Greater 66%" = 8))+
        scale_fill_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')+
        facet_grid(vars(var4),scales = "free_y", space = "free_y") + ylab('')

      return(plot_)



    }



  }


}



#' Add two numbers.
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(2, 3)
#' @export

PlotMS1vMS2 <- function(QQC){

  peps <- QQC@matricies@peptide
  raw_data <- QQC@raw_data
  raw_data <- raw_data %>% filter(Channel.Q.Value < QQC@misc[['ChQ']])

  rownames(peps) <- QQC@matricies@peptide_protein_map$seqcharge

  raw_data <- raw_data %>% filter(ID %in% colnames(peps))

  MS2_mat <- reshape2::dcast(raw_data, seqcharge ~ ID, value.var = 'Precursor.Quantity')

  rownames(MS2_mat)<- MS2_mat$seqcharge
  MS2_mat$seqcharge <- NULL
  MS2_mat <- as.matrix(MS2_mat)

  sect_pep <- intersect(rownames(MS2_mat),  rownames(peps) )
  sect_col <- intersect(colnames(MS2_mat),  colnames(peps) )

  MS2_mat <- MS2_mat[sect_pep,sect_col]
  peps <- peps[sect_pep,sect_col]

  MS2_mat <- Normalize_reference_vector(MS2_mat, log = T)
  peps <- Normalize_reference_vector(peps, log = T)

  cors <- c()
  pep <- c()
  for(i in 1:nrow(MS2_mat)){
    same <- psych::pairwiseCount(MS2_mat[i,],peps[i,])
    if(same > 10){
      pep <- c(pep, rownames(peps)[i])
      cors <- c(cors, cor(MS2_mat[i,],peps[i,], use = 'pairwise.complete.obs'))
    }

  }


  df <- as.data.frame(cors)
  df$FC <- rowMeans(abs(peps[pep,]),na.rm = T)
  df$FC[df$FC < .4] <- .4
  df$FC[df$FC > .4 & df$FC < .8] <- .8
  df$FC[df$FC > .8 & df$FC < 1.2] <- 1.2
  df$FC[df$FC > 1.2 & df$FC < 1.6] <- 1.6
  df$FC[df$FC > 1.6] <- 2

  df$FC <- as.character(df$FC)


  ggplot(df, aes(x = FC, y = cors))+ geom_boxplot(color="black", fill = 'gray') +
    ggtitle('Ms1, Ms2 Peptide Correlation') + ylab('Correlations') +
    xlab('Average abs(log2 Fold Change)')+ theme_linedraw()+
    stat_summary(fun.data=f, geom="text", vjust=-0.5, col="blue")





}





ImputationComparison <- function(QQC, cluster = 1){
  cluster = 1
  clusters <- QQC@reductions[['UMAP']]
  clusters$ID <- rownames(clusters)



  prot_imp <- QQC@protein.imputed
  prot_noimp <- QQC@protein


  prot_imp <- reshape2::melt(prot_imp)
  prot_imp <- prot_imp %>% left_join(clusters, by = c('Var2' = 'ID'))
  prot_imp <- prot_imp %>% group_by(cluster,Var1) %>% summarise(prot_score = median(value,na.rm=T))
  prot_imp <- prot_imp %>% filter(cluster == 2)

  prot_noimp <- reshape2::melt(prot_noimp)
  prot_noimp <- prot_noimp %>% left_join(clusters, by = c('Var2' = 'ID'))
  prot_noimp <- prot_noimp %>% group_by(cluster,Var1) %>% summarise(prot_score = median(value,na.rm=T))
  prot_noimp <- prot_noimp %>% filter(cluster == 2)


  prot_noimp$imp <- prot_imp$prot_score

  ggplot(prot_noimp, aes(x = prot_score, y = imp)) + geom_point()+ theme_bw()+
    xlab('No Imputation') + ylab('With Imputation') +
    ggtitle(paste0('Cluster ',cluster, ' Protein Averages, Cor = ',
                   round(cor(prot_noimp$prot_score,prot_imp$prot_score, use = 'pairwise.complete.obs'),3)))






}


