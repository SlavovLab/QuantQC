Proc_fasta <- function(path){
  convert_mouse <- read.fasta(path,set.attributes = T,whole.header = T)
  convert_mouse <- names(convert_mouse)
  parse_row<-grep("GN=",convert_mouse, fixed=T)
  split_prot<-str_split(convert_mouse[parse_row], pattern = fixed("GN="))
  gene<-unlist(split_prot)[seq(2,2*length(split_prot),2)]
  prot <- unlist(split_prot)[seq(1,2*length(split_prot),2)]
  prot_parse <- grep("|",prot, fixed=T)
  gene_parse <- grep(" ",gene, fixed=T)
  split_gene<-str_split(gene[parse_row], pattern = fixed(" "))
  split_gene<-unlist(split_gene)[seq(1,3*length(split_gene),3)]
  split_prot<-str_split(prot[parse_row], pattern = fixed("|"))
  split_prot<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
  convert_mouse  <- as.data.frame(cbind(split_prot,split_gene))

  return(convert_mouse)
}




# FindPermeableCells <- function(mat,species = 'Human'){
#
#   mat <- mat/rowSds(mat,na.rm = T)
#
#   training_data <- system.file("extdata", "PermeableClassification.csv", package = "QuantQC")
#   training_data <- read.csv(training_data,row.names = 1)
#
#   training_data_dat <- training_data
#   training_data_dat$train <- NULL
#   training_data_dat <- as.matrix(training_data_dat)
#
#
#   if(!species %in% c('Human','Mouse')){
#     return("species not supported, options are 1. Human 2. Mouse")
#   }
#
#   if(species == 'Human'){
#
#     Human <- system.file("extdata", "Human.fasta", package = "QuantQC")
#     Human <- Proc_fasta(Human)
#
#     Human <- Human %>% filter(split_gene %in% colnames(training_data_dat))
#     Human <- Human %>% filter(split_prot %in% rownames(mat))
#
#     training_data_dat <- training_data_dat[,Human$split_gene]
#
#     mat <- mat[Human$split_prot,]
#     rownames(mat) <- Human$split_gene
#
#     mat <- mat[colnames(training_data_dat),]
#
#     mat <- t(mat)
#     mat <- mat/colSds(mat,na.rm = T)
#
#
#
#   }
#
#   if(species == 'Mouse'){
#
#     Mouse <- system.file("extdata", "Mouse.fasta", package = "QuantQC")
#     Mouse <- Proc_fasta(Mouse)
#     Mouse$split_gene <- toupper(Mouse$split_gene)
#     Mouse <- Mouse %>% filter(split_gene %in% colnames(training_data_dat))
#     Mouse <- Mouse %>% filter(split_prot %in% rownames(mat))
#
#     training_data_dat <- training_data_dat[,Mouse$split_gene]
#
#     mat <- mat[Mouse$split_prot,]
#     rownames(mat) <- Mouse$split_gene
#
#     mat <- mat[colnames(training_data_dat),]
#
#     mat <- t(mat)
#
#     mat <- mat/colSds(mat,na.rm = T)
#
#
#   }
#
#
#   train_matrix <- xgb.DMatrix(data = training_data_dat, label = training_data$train)
#
#   # Set parameters for the model
#   params <- list(objective = "binary:logistic", eval_metric = "logloss")
#
#   # Train the model
#   model <- xgb.train(params, train_matrix, nrounds = 100)
#
#
#
#
#
#   predict_mat <- xgb.DMatrix(data = mat)
#
#   # Make predictions
#   predictions <- predict(model, newdata = predict_mat)
#
#
#
#
#   return(predictions)
#
#
#
# }




#' gives probability cell is permeable
#'
#' This function takes a matrix as an object with cell IDs on columns and uniprot IDs on rownames
#'
#' @param mat an expression matrix (can have NAs) proteins rows cellID columns
#' @param species Human or Mouse
#' @return A vector of probabilities length number columns of matrix
#' @examples
#' FindPermeableCells(protein_matrix, species = 'Mouse')
#' @export
FindPermeableCells <- function(mat, species = 'Human'){

  mat <- mat/rowSds(mat, na.rm = TRUE)

  training_data <- system.file("extdata", "PermeableClassification.csv", package = "QuantQC")
  training_data <- read.csv(training_data, row.names = 1)

  training_data_dat <- training_data
  training_data_dat$train <- NULL
  training_data_dat <- as.matrix(training_data_dat)

  if (!species %in% c('Human', 'Mouse')) {
    stop("Species not supported. Options are: 'Human', 'Mouse'")
  }

  if (species == 'Human') {
    Human <- system.file("extdata", "Human.fasta", package = "QuantQC")
    Human <- Proc_fasta(Human)

    Human <- Human %>% filter(split_gene %in% colnames(training_data_dat))
    Human <- Human %>% filter(split_prot %in% rownames(mat))

    training_data_dat <- training_data_dat[, Human$split_gene]

    mat <- mat[Human$split_prot, ]
    rownames(mat) <- Human$split_gene

    mat <- mat[colnames(training_data_dat), ]
    mat <- t(mat)
    mat <- mat/colSds(mat, na.rm = TRUE)
  }

  if (species == 'Mouse') {
    Mouse <- system.file("extdata", "Mouse.fasta", package = "QuantQC")
    Mouse <- Proc_fasta(Mouse)
    Mouse$split_gene <- toupper(Mouse$split_gene)
    Mouse <- Mouse %>% filter(split_gene %in% colnames(training_data_dat))
    Mouse <- Mouse %>% filter(split_prot %in% rownames(mat))

    training_data_dat <- training_data_dat[, Mouse$split_gene]

    mat <- mat[Mouse$split_prot, ]
    rownames(mat) <- Mouse$split_gene

    mat <- mat[colnames(training_data_dat), ]
    mat <- t(mat)
    mat <- mat/colSds(mat, na.rm = TRUE)
  }

  # Check if xgboost is available
  if (requireNamespace("xgboost", quietly = TRUE)) {
    train_matrix <- xgboost::xgb.DMatrix(data = training_data_dat, label = training_data$train)
    params <- list(objective = "binary:logistic", eval_metric = "logloss")
    model <- xgboost::xgb.train(params, train_matrix, nrounds = 100)

    predict_mat <- xgboost::xgb.DMatrix(data = mat)
    predictions <- predict(model, newdata = predict_mat)
    return(predictions)
  } else {
    stop("The 'xgboost' package is required to use this function.")
  }
}


