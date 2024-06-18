#####
# TMT Processing
#####


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
cellXpeptide <- function(QQC, TQVal = 1, chQVal = 1){

  if(QQC@ms_type == 'DIA' |QQC@ms_type == 'DIA_C' ){
    QQC <- cellXpeptide_DIA(QQC,TQVal, chQVal)
  }

  if(QQC@ms_type == 'DDA'){
    QQC <- TMT_Reference_channel_norm(QQC)
  }

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
TMT_Reference_channel_norm <- function(QQC){

  plex <- QQC@misc[['plex']]

  sc.data <- QQC@raw_data

  sc.data <- sc.data %>% filter(Reporter.intensity.2 != 0)

  if(plex == 14){
    ri.index<-which(colnames(sc.data)%in%paste0("Reporter.intensity.",2:18))
  }
  if(plex == 29){
    ri.index<-which(colnames(sc.data)%in%paste0("Reporter.intensity.",2:32))
  }



  sc.data[, ri.index] <- sc.data[, ri.index] / sc.data[, ri.index[1]]

  sc.data <- as.data.table(sc.data)

  if(plex == 14){
    sc.data <- sc.data[,c('seqcharge','Leading.razor.protein','Raw.file','Well','plate',paste0("Reporter.intensity.",5:18))]
  }
  if(plex == 29){
    sc.data <- sc.data[,c('seqcharge','Leading.razor.protein','Raw.file','Well','plate',paste0("Reporter.intensity.",4:32))]
  }

  sc.data <- data.table::melt(sc.data, id = c('seqcharge','Leading.razor.protein','Raw.file','Well','plate'))

  if(plex == 14){
    sc.data <- sc.data[sc.data$value < 2.5,]
  }
  if(plex == 29){
    #sc.data <- sc.data[sc.data$value < 10,]
  }

  sc.data$ID <- paste0(sc.data$Well,sc.data$plate,sc.data$variable)

  sc.data <- data.table::dcast(sc.data,Leading.razor.protein+seqcharge~ID,value.var = 'value')


  # Add in 0 peptides for negative controls that were totally filtered out
  #sc.data[sc.data==0] <- NA

  sc.data <- as.data.frame(sc.data)

  sc.data <- sc.data %>% distinct(seqcharge,.keep_all = T)


  prot_pep_map <- as.data.frame(cbind(sc.data$Leading.razor.protein,sc.data$seqcharge))
  colnames(prot_pep_map) <- c('Protein','seqcharge')

  data_matricies <- new('matricies_DDA',peptide = as.matrix(sc.data[,3:ncol(sc.data)]),peptide_protein_map = prot_pep_map)

  QQC@matricies <- data_matricies

  return(QQC)

}



#####
# mTRAQ Processing
#####

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
cellXpeptide_DIA <- function(QQC,TQVal, chQVal){

  Raw_data <- QQC@raw_data
  plex <- QQC@misc[['plex']]
  type <- QQC@ms_type



  if(plex == 2 & type == 'DIA_C'){
    plex_used <- c(0,4,8)

  }else if(plex == 2 & type == 'DIA'){
    plex_used <- c(0,4)
  }else if(plex == 3){
    plex_used <- c(0,4,8)
  }else{
    return('plex not valid')
  }


  Raw_data <- Raw_data %>% filter(plex %in% plex_used)
  Raw_data_filt <- Raw_data %>% filter(Channel.Q.Value < chQVal)
  Raw_data_filt <- Raw_data_filt %>% filter(Translated.Q.Value < TQVal)

  Raw_data_lim_filt <- Raw_data_filt %>% dplyr::select(Protein.Group,seqcharge,Ms1.Area,File.Name)
  Raw_data_lim.d_filt <- reshape2::dcast(Raw_data_lim_filt,Protein.Group+seqcharge~File.Name,value.var = 'Ms1.Area')
  Raw_data_lim.d_filt[Raw_data_lim.d_filt==0] <- NA

  Raw_data_lim_NF <- Raw_data %>% dplyr::select(Protein.Group,seqcharge,Ms1.Area,File.Name)
  Raw_data_lim.d_NF <- reshape2::dcast(Raw_data_lim_NF,Protein.Group+seqcharge~File.Name,value.var = 'Ms1.Area')

  Raw_data_lim.d_NF <- Raw_data_lim.d_NF %>% filter(seqcharge %in% Raw_data_lim.d_filt$seqcharge)



  # Normalize data by carrier
  ## This code also removes all sets where a single cell is larger in mean intensity than the carrier
  if(type == 'DIA_C'){

    Raw_data_lim.d_filt <- DIA_carrier_norm(Raw_data_lim.d_filt,8,plex_used)

    Raw_data_lim.d_NF <- Raw_data_lim.d_NF %>% dplyr::select(colnames(Raw_data_lim.d_filt))


  }



  peptide_protein_map <- Raw_data_lim.d_filt %>% dplyr::select(seqcharge,Protein.Group)

  colnames(peptide_protein_map) <- c('seqcharge','Protein')

  Raw_data_lim.d_filt <- as.matrix(Raw_data_lim.d_filt[,3:ncol(Raw_data_lim.d_filt)])
  Raw_data_lim.d_NF <- as.matrix(Raw_data_lim.d_NF[,3:ncol(Raw_data_lim.d_NF)])
  Raw_data_lim.d_NF[Raw_data_lim.d_NF==0] <- NA
  pep_mask <- is.na(Raw_data_lim.d_NF)==F


  data_matricies <- new('matricies_DIA',peptide = Raw_data_lim.d_filt,peptide_mask = pep_mask,peptide_protein_map = peptide_protein_map)


  QQC@matricies <- data_matricies
  QQC@misc[['ChQ']] <- chQVal
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
remove_sets_carrier_lessthan_SCs <- function(df,files){
  files_remove <- c()
  #df[is.nan(df)] <- NA
  df[df ==0] <- NA
  df[df == Inf] <- NA
  df[df == -Inf] <- NA


  for(i in unique(files)){
    if(length(intersect(colnames(df),(paste0(i,'8')))) > 0){
      count = FALSE
      count1 = FALSE
      count2 = FALSE
      if(length(intersect(colnames(df),(paste0(i,'0')))) > 0){
        count1 = TRUE
        dif1 <- median(as.matrix((df %>% select(paste0(i,'8')))/(df %>% select(paste0(i,'0'))))[,1],na.rm = T)
        if(dif1 > 3){
          count = TRUE
          count1 = TRUE
        }
      }
      if(length(intersect(colnames(df),(paste0(i,'4')))) > 0){
        count2 = TRUE
        dif2 <- median(as.matrix((df %>% select(paste0(i,'8')))/(df %>% select(paste0(i,'4'))))[,1],na.rm = T)
        if(dif2 > 3){
          count = TRUE

        }else{
          count = FALSE
        }
      }

      if(count == FALSE){
        if(count1 == TRUE){
          df <- df %>% select(-paste0(i,'0'))
        }
        if(count2 == TRUE){
          df <- df %>% select(-paste0(i,'4'))
        }

        df <- df %>% select(-paste0(i,'8'))
        files_remove <- c(files_remove,i)
      }

    }else{
      df <- df %>% select(-paste0(i,'4'))
      df <- df %>% select(-paste0(i,'0'))
      files_remove <- c(files_remove,i)
    }
  }
  files <- files[!files %in% files_remove]
  return(list(files,df))

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
DIA_carrier_norm <- function(Raw_data_lim.d,carrier_CH,plex_used){

  Carrier_list <- colnames(Raw_data_lim.d)[str_sub(colnames(Raw_data_lim.d), start= -1) == carrier_CH]

  Remove_cells <- c()
  for(i in Carrier_list){
    Flag = FALSE
    i <- str_sub(i, 1, -2)
    for(j in 1:length(plex_used)){
      if(paste0(i,plex_used[j]) %in% colnames(Raw_data_lim.d)){
        if(median(as.matrix((Raw_data_lim.d %>% dplyr::select(paste0(i,plex_used[j])))/(Raw_data_lim.d %>% dplyr::select(paste0(i,carrier_CH)))),na.rm = T) > 1){
          Flag = TRUE
        }
      }
    }
    if(Flag == T){
      Remove_cells <- c(Remove_cells,i)
    }
  }

  non_carrier <- colnames(Raw_data_lim.d)[!colnames(Raw_data_lim.d) %in% Carrier_list]
  for(i in Remove_cells){
    Carrier_list <- Carrier_list[Carrier_list != paste0(i,carrier_CH)]
  }
  for(i in 1:length(plex_used)){
    non_carrier <-non_carrier[non_carrier != paste0(Remove_cells,plex_used[i])]
  }
  all_cols <- c('Protein.Group','seqcharge',non_carrier)
  Carrier_list <- str_sub(Carrier_list, 1, -2)



  for(i in Carrier_list){
    for(j in 1:length(plex_used)){
      if(paste0(i,plex_used[j]) %in% colnames(Raw_data_lim.d)){
        Raw_data_lim.d[,paste0(i,plex_used[j])] <- (Raw_data_lim.d %>% dplyr::select(paste0(i,plex_used[j])))/(Raw_data_lim.d %>% dplyr::select(paste0(i,carrier_CH)))
      }
    }
  }


  Raw_data_lim.d <- Raw_data_lim.d %>% dplyr::select(intersect(all_cols,colnames(Raw_data_lim.d)))

  Raw_data_lim.d[Raw_data_lim.d == Inf] <- NA
  Raw_data_lim.d[Raw_data_lim.d == -Inf] <- NA
  Raw_data_lim.d[Raw_data_lim.d ==0] <- NA

  return(Raw_data_lim.d)
}


#####
#General Proc
#####

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
KNN_impute<-function(QQC, k = 3){


  sc.data <- QQC@matricies@protein



  # Create a copy of the data, NA values to be filled in later
  sc.data.imp<-sc.data

  # Calculate similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(sc.data)) )
  #dist.mat<- 1-as.matrix(cor((dat), use="pairwise.complete.obs"))

  #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))

  # Column names of the similarity matrix, same as data matrix
  cnames<-colnames(dist.mat)

  # For each column in the data...
  for(X in cnames){

    # Find the distances of all other columns to that column
    distances<-dist.mat[, X]

    # Reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered<-distances[order(distances, decreasing = F)]

    # Reorder the data matrix columns, smallest distance to largest from the column of interest
    # Obviously, first column will be the column of interest, column X
    dat.reordered<-sc.data[ , names(distances.ordered ) ]

    # Take the values in the column of interest
    vec<-sc.data[, X]

    # Which entries are missing and need to be imputed...
    na.index<-which( is.na(vec) )

    # For each of the missing entries (rows) in column X...
    for(i in na.index){

      # Find the most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )

      #print(length(closest.columns))

      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){

        # Replace NA in column X with the mean the same row in k of the most similar columns
        vec[i]<-mean( sc.data[ i, closest.columns[1:k] ] )

      }


      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){

        # Replace NA in column X with the mean the same row in all of the most similar columns
        vec[i]<-mean( sc.data[ i, closest.columns ] )

      }


    }

    # Populate a the matrix with the new, imputed values
    sc.data.imp[,X]<-vec

  }

  # Normalize imputed data
  sc.data.imp <- Normalize_reference_vector_log(sc.data.imp)



  QQC@matricies@protein.imputed <- sc.data.imp



  return(QQC)

}



MinValue_impute <- function(QQC){
  prot_dat <- QQC@matricies@protein

  for(i in 1:nrow(prot_dat)){
    prot_vect <- prot_dat[i,]
    prot_vect[is.na(prot_vect) == T] <- min(prot_vect,na.rm = T)
    prot_dat[i,] <- prot_vect

  }


  QQC@matricies@protein.imputed <- (prot_dat)

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
Normalize_reference_vector <- function(dat,log = F){
  dat <- as.matrix(dat)
  dat[dat==0] <- NA

  refVec <- matrixStats::rowMedians(x = dat, cols = 1:ncol(dat), na.rm = T)
  for(k in 1:ncol(dat)){
    dat[,k]<-dat[,k] * median(refVec/dat[,k], na.rm = T)
  }
  for(k in 1:nrow(dat)){
    dat[k,]<-dat[k,]/mean(dat[k,], na.rm = T)
  }
  if(log == T){
    dat <- log2(dat)
  }
  return(dat)
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
Normalize_reference_vector_log <- function(dat){
  dat <- as.matrix(dat)

  refVec <- matrixStats::rowMedians(x = dat, cols = 1:ncol(dat), na.rm = T)
  for(k in 1:ncol(dat)){
    dat[,k]<-dat[,k] + median(refVec - dat[,k], na.rm = T)
  }
  for(k in 1:nrow(dat)){
    dat[k,]<-dat[k,]-mean(dat[k,], na.rm = T)
  }
  return(dat)
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
CollapseToProtein <- function(QQC, opt, norm = 'ref'){

  sc.data <- QQC@matricies@peptide


  ## Store absolute abundances
  Abs_peptide_data <- as.data.table(cbind(QQC@matricies@peptide_protein_map,sc.data))
  Abs_peptide_data <- data.table::melt(Abs_peptide_data,id.vars = c('Protein','seqcharge'))
  Abs_peptide_data$seqcharge <- NULL
  Abs_peptide_data <- Abs_peptide_data[, lapply(.SD, median,na.rm = TRUE), by = c('Protein','variable')]

  # Create Protein x Cell matrix

  Abs_peptide_data <- data.table::dcast(Abs_peptide_data, Protein ~ variable, value.var = 'value')

  Abs_peptide_data <- as.data.frame(Abs_peptide_data)
  rownames(Abs_peptide_data) <- Abs_peptide_data$Protein
  Abs_peptide_data$Protein <- NULL



  QQC@matricies@protein_abs <- as.matrix(Abs_peptide_data)


  # This function colapses peptide level data to the protein level
  # There are different ways to collapse peptides mapping from the same
  # protein to a single data point. The simplest way is to take the median
  # of the peptide values after normalizing. There are also more sophisticated
  # ways that rely on the raw peptide intensities like MaxLFQ

  if(opt == 1){

    # Remove peptide protein name columns
    #Normalize_peptide_data <- as.matrix(sc.data)
    #
    # # Normalize peptide data for cell size and then to relative abundances and log transform
    # if(norm == 'std'){
    #   Normalize_peptide_data <- normalize(Normalize_peptide_data,log = T)
    # }
    # if(norm == 'ref'){
    #   Normalize_peptide_data <- Normalize_reference_vector(Normalize_peptide_data,log = T)
    # }



    # Remove unwanted values
    #Normalize_peptide_data[Normalize_peptide_data == Inf] <- NA
    #Normalize_peptide_data[Normalize_peptide_data == -Inf] <- NA

    if(QQC@ms_type == 'DDA'){

      Normalize_peptide_data <- Normalize_reference_vector(sc.data,log = T)
      QQC@matricies@peptide <-  Normalize_peptide_data

    }

    if(QQC@ms_type != 'DDA'){

      QQC <- LC_BatchCorrect(QQC)
      Normalize_peptide_data <- QQC@matricies@peptide

    }



    #Re-Join data
    Normalize_peptide_data <- as.data.table(cbind(QQC@matricies@peptide_protein_map,Normalize_peptide_data))

    # Remove peptides observed less than 10 times
    Normalize_peptide_data <- Normalize_peptide_data[rowSums(is.na(Normalize_peptide_data)==F) > 3,]


    # Collapse peptide levels to median protein level, first melt, then collapse then expand back out
    Normalize_peptide_data <- data.table::melt(Normalize_peptide_data,id.vars = c('Protein','seqcharge'))
    Normalize_peptide_data$seqcharge <- NULL
    Normalize_protein_data <- Normalize_peptide_data[, lapply(.SD, median,na.rm = TRUE), by = c('Protein','variable')]

    # Create Protein x Cell matrix

    Normalize_protein_data <- data.table::dcast(Normalize_protein_data, Protein ~ variable, value.var = 'value')

    Normalize_protein_data <- as.data.frame(Normalize_protein_data)
    rownames(Normalize_protein_data) <- Normalize_protein_data$Protein
    Normalize_protein_data$Protein <- NULL

    # Re-column and row normalize:
    if(norm == 'std'){
      Normalize_protein_data <- normalize_log(Normalize_protein_data)
    }
    if(norm == 'ref'){
      Normalize_protein_data <- Normalize_reference_vector_log(Normalize_protein_data)
    }


    QQC@matricies@protein <- Normalize_protein_data


    if(QQC@ms_type == 'DIA' | QQC@ms_type == 'DIA_C'){

      prot_NF <- QQC@matricies@peptide_mask
      prot_NF <- cbind(prot_NF,QQC@matricies@peptide_protein_map)
      prot_NF$seqcharge <- NULL
      prot_NF <- reshape2::melt(prot_NF,id.var = 'Protein')
      prot_NF <- prot_NF %>% dplyr::group_by(Protein,variable) %>% dplyr::summarise(numb = sum(value)>0)

      prot_NF <- reshape2::dcast(prot_NF,Protein~variable,value.var = 'numb')
      prot_NF$Protein <- NULL
      prot_NF <- as.matrix(prot_NF)
      QQC@matricies@protein_mask <- prot_NF

    }

    return(QQC)

  }




  if(opt == 2){
    library(diann)
    #Max LFQ protein level

    sc.data <- QQC@raw_data




    protein_mat <- diann_maxlfq(sc.data,sample.header = "File.Name",group.header = "Protein.Group",id.header = "seqcharge",quantity.header = "Ms1.Area")

    #Normalize protein level data and log transform
    if(norm == 'std'){
      Normalize_peptide_data <- normalize(protein_mat,log = T)
    }
    if(norm == 'ref'){
      Normalize_peptide_data <- Normalize_reference_vector(protein_mat,log = T)
    }
    QQC@matricies@protein <- Normalize_protein_data



    if(QQC@ms_type == 'DIA' | QQC@ms_type == 'DIA_C'){

      prot_NF <- QQC@matricies@peptide_mask
      prot_NF <- cbind(prot_NF,QQC@matricies@peptide_protein_map)
      prot_NF$seqcharge <- NULL
      prot_NF <- reshape2::melt(prot_NF,id.var = 'Protein')
      prot_NF <- prot_NF %>% dplyr::group_by(Protein,variable) %>% dplyr::summarise(numb = sum(value)>0)

      prot_NF <- reshape2::dcast(prot_NF,Protein~variable,value.var = 'numb')
      prot_NF$Protein <- NULL
      prot_NF <- as.matrix(prot_NF)
      QQC@matricies@protein_mask <- prot_NF

    }



    return(QQC)
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
BatchCorrect <- function(QQC, labels = T, run = T, batch = F, norm = 'ref'){

  if(run == T & batch == T){
    return("cant correct on run and batch")
  }

  cellenONE_meta <- QQC@meta.data
  protein_mat_imputed <- QQC@matricies@protein.imputed
  protein_mat <- QQC@matricies@protein

  # Get meta data for batch correction
  batch_label  <- cellenONE_meta %>% dplyr::filter(ID %in% colnames(protein_mat_imputed))
  batch_label <- batch_label[order(match(batch_label$ID,colnames(protein_mat_imputed))),]



  # Perform batch corrections, possible sources label bias, Every LC/MS runs or groups of LC/MS runs
  #sc.batch_cor <- ComBat(protein_mat_imputed, batch=factor(batch_label$label))

  if(labels == T & run == T){
    sc.batch_cor <- limma::removeBatchEffect(protein_mat_imputed,batch = batch_label$injectWell, batch2 = batch_label$label)
  }
  if(labels == T & batch == T){
    sc.batch_cor <- limma::removeBatchEffect(protein_mat_imputed,batch = batch_label$label, batch2 = batch_label$LCMS_Batch)
  }
  if(labels == T & batch == F & run == F){
    sc.batch_cor <- sva::ComBat(protein_mat_imputed,batch = batch_label$label)
  }
  if(labels == F & batch == T & run == F){
    sc.batch_cor <- limma::removeBatchEffect(protein_mat_imputed,batch = batch_label$LCMS_Batch)
  }
  if(labels == F & batch == F & run == T){
    sc.batch_cor <- limma::removeBatchEffect(protein_mat_imputed,batch = batch_label$injectWell)
  }



  # Re normalize and NA out imputed values
  if(norm == 'std'){
    sc.batch_cor <- normalize_log(sc.batch_cor)
  }
  if(norm == 'ref'){
    sc.batch_cor <- Normalize_reference_vector_log(sc.batch_cor)
  }


  # Store unimputed matrix
  sc.batch_cor_noimp <- sc.batch_cor
  sc.batch_cor_noimp[is.na(protein_mat)==T] <- NA

  QQC@matricies@protein.imputed <- sc.batch_cor
  QQC@matricies@protein <- sc.batch_cor_noimp


  return(QQC)

}




LC_BatchCorrect <- function(QQC){

  pep_norm <- QQC@matricies@peptide
  pep_norm[pep_norm==0] <- NA

  pep_norm <- QuantQC::Normalize_reference_vector(pep_norm,log = T)

  order_vect <- QQC@meta.data %>% filter(ID %in% colnames(QQC@matricies@peptide))

  if(sum(order_vect$ID == colnames(pep_norm)) != ncol(pep_norm)){
    return('something went wrong')
  }

  for(i in 1:nrow(pep_norm)){

    set_df <- floor(sum(is.na(pep_norm[i,])==F)/11)
    if(set_df > 20){
      set_df <- 20
    }
    if(set_df == 0){
      set_df <- 2
    }
    if(set_df == 1){
      set_df <- 2
    }

    if(sum(is.na(pep_norm[i,])==F) > 29){


      df_spline <- as.data.frame(order_vect$Order)
      colnames(df_spline) <- 'order'
      df_spline$data <- pep_norm[i,]
      df_spline$holder <- 1:ncol(pep_norm)
      df_spline <- df_spline %>% filter(is.na(data)==F)

      if(QQC@ms_type == 'DDA'){
        set_df <- length(unique(df_spline$order))
        if(set_df > 20){
          set_df <- 20
        }
        if(set_df == 3){
          df_spline$order[1] <- df_spline$order[1] +.001
        }


      }



      smooth = stats::smooth.spline(df_spline$order,df_spline$data,df = set_df)
      predicted_y <- predict(smooth,  df_spline$order)

      RSS <- sum((df_spline$data - predicted_y$y)^2)

      # Calculate total sum of squares (TSS)
      TSS <- sum((df_spline$data - mean(predicted_y$y))^2)

      # Calculate R-squared (R²)
      R_squared <- 1 - (RSS / TSS)

      if(R_squared > .1){

        x = predicted_y$x
        y = predicted_y$y
        df <- as.data.frame(cbind(x,y))
        df <- df[order(df$x),]

        df_spline <- df_spline[order(df_spline$order),]
        df_spline$data = df_spline$data-df$y

        pep_norm[i,df_spline$holder] <- df_spline$data

      }
    }

  }

  QQC@matricies@peptide <- pep_norm

  return(QQC)

}



# BatchCorrectPeptide <- function(QQC){
#
#   cellenONE_meta <- QQC@meta.data
#
#   peptide_mat <- QQC@peptide
#   peptide_mat <- Normalize_reference_vector(peptide_mat, log = T)
#   QQC@peptide <- peptide_mat
#
#   QQC <- KNN_impute(QQC, k = 3, dat = 'peptide')
#   peptide_mat_imputed <- QQC@peptide.imputed
#
#   # Get meta data for batch correction
#   batch_label  <- cellenONE_meta %>% dplyr::filter(ID %in% colnames(peptide_mat_imputed))
#   batch_label <- batch_label[order(match(batch_label$ID,colnames(peptide_mat_imputed))),]
#
#
#
#   # Perform batch corrections, possible sources label bias, Every LC/MS runs or groups of LC/MS runs
#   #sc.batch_cor <- ComBat(protein_mat_imputed, batch=factor(batch_label$label))
#   sc.batch_cor <- limma::removeBatchEffect(peptide_mat_imputed,batch = batch_label$injectWell, batch2 = batch_label$label)
#
#   # Re normalize and NA out imputed values
#   sc.batch_cor <- Normalize_reference_vector_log(sc.batch_cor)
#
#   # Store unimputed matrix
#   sc.batch_cor_noimp <- sc.batch_cor
#   sc.batch_cor_noimp[is.na(peptide_mat)==T] <- NA
#
#   QQC@peptide.imputed <- sc.batch_cor
#   QQC@peptide <- sc.batch_cor_noimp
#
#
#   return(QQC)
#
# }
