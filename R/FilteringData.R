

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
EvaluateNegativeControls <- function(QQC,CV_thresh){

  if(nrow(QQC@meta.data)==0){
    return("Must map sample identites to data first")
  }


  if(QQC@ms_type == 'DDA'){
    QQC <- EvaluateNegativeControls_DDA(QQC,CV_thresh)
  }

  if(QQC@ms_type == 'DIA' | QQC@ms_type ==  'DIA_C'){
    QQC <- EvaluateNegativeControls_DIA(QQC)
  }


  return(QQC)


}




CVs <- function(QQC,thresh){
  cell_id <- QQC@meta.data
  # Normalize peptide data
  mat_norm <- Normalize_reference_vector(QQC@matricies@peptide)

  mat_norm <- as.data.frame(mat_norm)
  mat_norm$Protein <- QQC@matricies@peptide_protein_map$Protein
  mat_norm$pep <- QQC@matricies@peptide_protein_map$seqcharge

  # convert to data.table for fast computation
  mat_norm <- as.data.table(mat_norm)

  # Melt data frame
  mat_norm.melt <- data.table::melt(mat_norm, id = c('Protein','pep'))


  # Call the function on the data.table
  CV_mat <- fast_cv(mat_norm.melt)

  # Revert back to DF
  CV_mat <- as.data.frame(CV_mat)
  CV_mat$value <- NULL

  # count number for protein with multiple peptides in each cell
  CV_mat_count <- CV_mat %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(counter = sum(is.na(cvq)==F))


  # only look at cells with atleast 20 proteins with multiple peptides
  CV_mat_count <- CV_mat_count %>% filter(counter > 20)


  # Take median protein CV per cell
  CV_mat <- CV_mat %>%
    dplyr::group_by( variable) %>%
    dplyr::summarise(cvq = median(cvq,na.rm = T))

  CV_mat <- CV_mat %>% filter(variable %in% CV_mat_count$variable)


  # store values from cells, not negative ctrls
  pos <- cell_id %>% filter(sample != 'neg')

  CV_mat$value <- NA
  CV_mat$value[CV_mat$variable %in% pos$ID] <- 'cell'
  CV_mat$value[!CV_mat$variable %in% pos$ID] <- 'neg'

  CV_mat_pos <- CV_mat %>% filter(value == 'cell')
  CV_mat_neg <- CV_mat %>% filter(value == 'neg')

  return(CV_mat)

}



Count_peptides_per_cell <- function(sc.data,cellenONE_meta,good_cells = NULL){




  # Get IDs for negative controls
  negative_IDs <- cellenONE_meta$ID[cellenONE_meta$sample == 'neg']


  # If negative control has 0 peptides detected, it will be left out of sc.data,
  # find negative controls with 0 peptides detected so 0s can be included on plot

  numb_neg_controls <- length(intersect(negative_IDs,colnames(sc.data)))


  #If negative controls are not all 0 peptides and there are more than 2, count number peptides for non 0 negative controls
  if(numb_neg_controls > 1){

    # Data matrix for negative controls

    neg_mat <- sc.data[,negative_IDs]

    # Count number peptides in each negative control
    sum_neg <- colSums(is.na(neg_mat)==F)
    sum_neg_int <- colSums(neg_mat,na.rm = T)
    # Make data frame for plotting
    neg_df <- as.data.frame(sum_neg)
    colnames(neg_df) <- 'Number_precursors'
    neg_df$intense <- sum_neg_int
    neg_df$type <- 'negative ctrl'



    # # Add in 0s for any negative controls with 0 peptides
    # if(zero_peptide_negs != 0 ){
    #   for(i in 1:zero_peptide_negs){
    #     neg_df[nrow(neg_df) + 1,] = c(0,0,"negative ctrl")
    #   }
    # }


  }else if(numb_neg_controls == 1){
    # If the negative controls all have 0 peptides measured
    neg_vect <- sc.data[,negative_IDs]

    # Make data frame for plotting
    neg_df <- as.data.frame(matrix(data = 0,ncol = 1,nrow = length(negative_IDs)))
    colnames(neg_df) <- 'Number_precursors'
    neg_df$intense <- sum(is.na(neg_vect)==F)
    neg_df$type <- 'negative ctrl'

  }else{
    # Make data frame for plotting
    neg_df <- as.data.frame(matrix(data = 0,ncol = 1,nrow = 1))
    colnames(neg_df) <- 'Number_precursors'
    neg_df$intense <- 0
    neg_df$type <- 'negative ctrl'
  }



  if(is.null(good_cells)==F){
    # Filter for single cells that passed CV test (for pSCoPE data ONLY)
    sum_cell <- sc.data[,colnames(sc.data) %in% good_cells$variable]

  }else{

    # for DIA its all non zero single cells
    good_cells <- cellenONE_meta %>% filter(sample != 'neg')
    sum_cell <- sc.data[,colnames(sc.data) %in% (good_cells$ID)]
  }

  # Sum number peptides for all real cells
  sum_cell_id <- colSums(is.na(sum_cell)==F)
  sum_cell_intense <- colSums(sum_cell,na.rm = T)

  # Make cell count data frame for plotting
  pos_df <- as.data.frame(sum_cell_id)
  colnames(pos_df) <- 'Number_precursors'
  pos_df$intense <-  sum_cell_intense
  pos_df$type <- 'single cells'

  # Combind positive and negative control data frames
  neg_vs_pos_DF <- rbind(pos_df,neg_df)

  return(neg_vs_pos_DF)
}



EvaluateNegativeControls_DDA <- function(QQC,CV_thresh){


  # Compute CVs of cells and negative controls, function outputs CV plot and list of cells with CVs
  CVm <- CVs(QQC,CV_thresh)


  # Get IDs of cells with median protein CVs (good_cells is a df with cell ID and CV)
  good_cells <- CVm %>% filter(cvq < CV_thresh)


  if(length(good_cells) < 3){

    return('less than three good cells, try increasing CV filter')

  }


  # Count the number of peptides in negative controls and good single cells
  Peptide_counts_by_sample <- Count_peptides_per_cell(QQC@matricies@peptide,QQC@meta.data,good_cells)


  # Plot distributions
  Numb_data_points <- ggplot(Peptide_counts_by_sample, aes(x = Number_precursors, fill = type)) + geom_histogram(position = 'identity', alpha = .5) + ggtitle(paste0('# precursors per sample')) + ylab('# of samples')+dot_plot



  # Filter for only good cells
  #cols_to_keep <- c('Protein.Group','seqcharge',as.character(good_cells$variable))

  Peptide_counts_by_sample$variable <- rownames(Peptide_counts_by_sample)
  neg_meta <-  CVm %>% left_join(Peptide_counts_by_sample, by = c('variable'))

  QQC@neg_ctrl.info <- neg_meta

  return(QQC)


}



EvaluateNegativeControls_DIA <- function(QQC){


  # Count the number of peptides in negative controls and good single cells
  Peptide_counts_by_sample <- Count_peptides_per_cell(QQC@matricies@peptide,QQC@meta.data)

  #Ref_norm_data_filtered <- Ref_norm_data[,colnames(Ref_norm_data) %in% cols_to_keep]
  Peptide_counts_by_sample$variable <- rownames(Peptide_counts_by_sample)

  QQC@neg_ctrl.info <- Peptide_counts_by_sample

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
PlotNegCtrl <- function(QQC,CV_thresh){


  plot_data <- QQC@neg_ctrl.info

  if(QQC@ms_type == 'DDA'){

    #peps <- ggplot(plot_data, aes(x = Number_precursors, fill = type)) + geom_histogram(position = 'identity', alpha = .5) + ggtitle(paste0('precursors per cell')) + ylab('# of samples')+dot_plot

    peps <- ggplot(plot_data, aes( x = log10(intense), fill = type)) + geom_histogram(position = 'identity',alpha = .5)  + ylab('# of samples')+dot_plot+
      xlab('sum(log10(Intensity))')

    CV_mat_pos <- plot_data %>% filter(value == 'cell')
    CV_mat_neg <- plot_data %>% filter(value == 'neg')

    cvs <- ggplot(data=plot_data, aes(x=cvq,fill=value)) + geom_density( alpha=0.5,adjust=1.5) + dot_plot+
      scale_fill_manual(values=my_col3) +
      xlab("CV, peptides from same protein") + ylab("Fraction of cells") + rremove("y.ticks") + rremove("y.text") +
      coord_cartesian(xlim=c(.1,.65))+
      annotate("text", x=0.2, y= 14, label=paste0(sum(CV_mat_pos$cvq < CV_thresh)," cells"), size=10, color=my_col3[c(1)])+
      annotate("text", x=0.64, y= 12, label=paste0(sum(CV_mat_neg$cvq > CV_thresh,na.rm = T)," Ctr -"), size=10, color=my_col3[c(2)])+
      annotate("text", x=0.63, y= 14, label=paste0(sum(CV_mat_pos$cvq > CV_thresh)," cells"), size=10, color=my_col3[c(1)])+
      annotate("text", x=0.2, y= 12, label=paste0((sum(CV_mat_neg$cvq < CV_thresh,na.rm = T)-1)," Ctr -"), size=10, color=my_col3[c(2)])+
      ggtitle('Min 3 proteins w/ many peps')+
      rremove("legend") + geom_vline(xintercept=CV_thresh, lty=2, size=2, color="gray50") + theme(plot.margin = margin(1, 1, 0, 1, "cm"))


    return(peps+cvs)
  }


  if(QQC@ms_type == 'DIA' | QQC@ms_type ==  'DIA_C'){
    peps <- ggplot(plot_data, aes( x = log10(intense), fill = type)) + geom_histogram(position = 'identity',alpha = .5) + ggtitle(paste0('Negative ctrl Vs Single cells')) + ylab('# of samples')+dot_plot+
      xlab('log10(Intensity)')

    return(peps)
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
FilterBadCells <- function(QQC, CV_thresh = NA, min_intens = NA){

  neg_filter <- QQC@neg_ctrl.info
  peptide_data <- QQC@matricies@peptide


  if(QQC@ms_type == 'DIA' | QQC@ms_type == 'DIA_C'){
    peptide_mask <- QQC@matricies@peptide_mask
    neg_filter <- neg_filter %>% dplyr::filter(type != 'negative ctrl')
    neg_filter <- neg_filter %>% dplyr::filter(log10(intense) > min_intens)
    peptide_mask <- peptide_mask[,colnames(peptide_mask) %in% neg_filter$variable]
    QQC@matricies@peptide_mask <- peptide_mask

  }

  if(QQC@ms_type == 'DDA'){
    neg_filter <- neg_filter %>% dplyr::filter(value != 'neg')
    if(is.na(min_intens)==F){
      neg_filter <- neg_filter %>% dplyr::filter(log10(intense) > min_intens)
    }
    neg_filter <- neg_filter %>% dplyr::filter(cvq < CV_thresh)


  }

  peptide_data <- peptide_data[,colnames(peptide_data) %in% neg_filter$variable]

  QQC@matricies@peptide <- peptide_data


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
filt.mat.cr<-function(mat, pct.r,pct.c){

  kc<-c()
  for(k in 1:ncol(mat)){

    pct.na<-length(which(is.na(mat[,k]))) / length(mat[,k])
    if(pct.na <= pct.c){ kc<-c(kc,k)}
    #print(pct.na)


  }

  mat<-mat[,kc]

  kr<-c()
  for(k in 1:nrow(mat)){

    pct.na<-length(which(is.na(mat[k,]))) / length(mat[k,])
    if(pct.na <= pct.r){ kr<-c(kr,k)}
    #print(pct.na)


  }

  mat<-mat[kr,]



  return(mat)

}
