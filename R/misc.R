

#formatting

my_col3 <- c("purple2","black")

dot_plot <-  theme_bw()+theme(plot.title = element_text(hjust = .5,size = 22),
                              axis.title.x = element_text(size = 18),
                              axis.title.y = element_text(size = 18),
                              axis.text.x = element_text(size = 16),
                              axis.text.y = element_text(size = 16),
                              legend.text = element_text(size = 18))

um_plot <-  theme_classic()+theme(plot.title = element_text(hjust = .5,size = 22),
                              axis.title.x = element_text(size = 18),
                              axis.title.y = element_text(size = 18),
                              axis.text.x = element_text(size = 12),
                              axis.text.y = element_text(size = 12))



# Other small less complex functions
extract_accession <- function(input_string) {
  result <- str_extract(input_string, "\\|([^|]+)\\|")
  result <- substr(result, 2, nchar(result) - 1)
  return(result)
}

inSet_norm <- function(Raw_data, cellenOne_meta){
  count = 0

  for(i in unique(paste0(cellenONE_meta$injectWell))){

    set <- Raw_data[,which(grepl(i,colnames(Raw_data)))]

    Raw_data[,which(grepl(i,colnames(Raw_data)))] <- set - rowMeans(set,na.rm = T)

  }

  return(Raw_data)
}


f <- function(y) {c(label=length(y), y=median(y))}


cv<-function(x){

  sd(x, na.rm=T) / mean(x, na.rm=T)

}


fast_cv <- function(dt) {
  summary <- dt[!is.na(value), .(sd_value = sd(value), mean_value = mean(value)), by = .(variable, Protein)]
  dt[!is.na(value), cvq := summary[.SD, on = .(variable, Protein), sd_value / mean_value]]
}


normalize <- function(evnew,log = F){
  evnew <- as.matrix(evnew)
  evnew[evnew==0] <- NA
  for(i in 1:ncol(evnew)){
    evnew[,i] <- evnew[,i]/median(evnew[,i],na.rm = T)
  }
  for(i in 1:nrow(evnew)){
    evnew[i,] <- evnew[i,]/mean(evnew[i,],na.rm = T)
  }
  if(log == T){
    evnew <- log2(evnew)
  }

  return(evnew)
}

normalize_log <- function(evnew,log = F){
  evnew <- as.matrix(evnew)
  evnew[evnew==0] <- NA
  for(i in 1:ncol(evnew)){
    evnew[,i] <- evnew[,i] - median(evnew[,i],na.rm = T)
  }
  for(i in 1:nrow(evnew)){
    evnew[i,] <- evnew[i,] - mean(evnew[i,],na.rm = T)
  }

  return(evnew)
}

create_peptide_vector <- function(protein_sequence, peptide) {
  positions <- regexpr(peptide, protein_sequence)
  binary_vector <- numeric(nchar(protein_sequence))

  if (attr(positions, "match.length") > 0) {
    peptide_start <- positions
    peptide_end <- positions + attr(positions, "match.length") - 1
    binary_vector[peptide_start:peptide_end] <- 1
  }

  return(binary_vector)
}



ContainsMissedCleaved <- function(pep_list){

  pep_list <- str_sub(pep_list,1,-2)

  pep_list_K <- grepl("K", pep_list, fixed = TRUE)
  pep_list_R <- grepl("R", pep_list, fixed = TRUE)

  pep_list <- pep_list_K+pep_list_R

  return(pep_list)

}
