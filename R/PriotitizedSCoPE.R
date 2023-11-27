### pSCoPE

Gen_Tlist_DIANN <- function(DIA,start_time){
  DIA_filt <- DIA %>%
    dplyr::select(Stripped.Sequence,Modified.Sequence,RT,PEP,Precursor.Charge,Precursor.Mz,Ms1.Area,Genes,Run,Protein.Group)

  DIA_filt$Charge_double <- as.double(DIA_filt$Precursor.Charge)
  DIA_filt$Mass <- ((DIA_filt$Precursor.Mz)*(DIA_filt$Charge_double))-(1.0072826748*((DIA_filt$Charge_double)))
  DIA_filt$MassRound <- round(DIA_filt$Mass,7)

  DIA_filt$SeqCharge <- paste0(DIA_filt$Stripped.Sequence, DIA_filt$Precursor.Charge)
  DIA_filt_single <- DIA_filt %>% dplyr::group_by(SeqCharge) %>% dplyr::top_n(-1,PEP) %>% dplyr::slice(1)

  DIA_filt_single$RT_sec <- round(DIA_filt_single$RT*60)
  DIA_filt_single$ind <- rownames(DIA_filt_single)
  DIA_filt_single <- DIA_filt_single %>% ungroup()

  filtDIA <- DIA_filt_single %>% group_by(SeqCharge) %>% slice(1)
  #unfiltDIA <- DIA_unfilt_single %>% group_by(SeqCharge) %>% slice(1)

  filtDIA$TargBool <- TRUE
  #unfiltDIA$TargBool <- FALSE

  DIA_input <- filtDIA

  DIA_input$Prot <- gsub(";.*","",DIA_input$Protein.Group)
  DIA_input$Prot <- gsub("-.*","",DIA_input$Prot)


  DIA_input$RTBool <- TRUE
  DIA_input$Masses <- "376.27"
  DIA_input$Leading.razor.protein <- DIA_input$Protein.Group
  DIA_input$Leading.razor.protein <- sub(';.*$','', DIA_input$Leading.razor.protein)


  colsToSelect <- c("Stripped.Sequence","Modified.Sequence","RT","Precursor.Charge","Ms1.Area","MassRound","TargBool","RTBool","Masses", "Leading.razor.protein","SeqCharge")

  FinalDF <- DIA_input %>% dplyr::ungroup() %>% dplyr::select(all_of(colsToSelect))

  #FinalDF$newMod <- paste0(FinalDF$PEPStrippedSequence,"0",rownames(FinalDF))

  colnames(FinalDF) <-c("Sequence","ModSeq","Retention.time","Charge","Apex.intensity","Mass","TargBoolean","RTC_Boolean","Fragments.mz","Leading.razor.protein","Modified")


  FinalDF$Retention.time <- round(FinalDF$Retention.time,3)
  FinalDF$Apex.intensity <- round(FinalDF$Apex.intensity,6)

  ###
  # Uncomment this if you're using a delayed acquisition method, otherwise not
  ###
  FinalDF$Retention.time <- FinalDF$Retention.time - start_time




  # Defining intensity tertiles for intensity-based priority tiers
  IntQuant <-quantile(FinalDF$Apex.intensity, c(.33,.66))

  # To define priority based on precursor intensity
  FinalDF$Priority <- ifelse((FinalDF$TargBoolean == TRUE) & (FinalDF$Apex.intensity <= IntQuant[[1]]),1,ifelse((FinalDF$TargBoolean == TRUE) & (FinalDF$Apex.intensity > IntQuant[[1]])&(FinalDF$Apex.intensity <= IntQuant[[2]]),2,ifelse((FinalDF$TargBoolean == TRUE)& (FinalDF$Apex.intensity > IntQuant[[2]]), 3,0)))

  FinalDF$ModSeq <- NULL

  FinalDF_differential <- FinalDF %>% filter(Priority == 3)

  prot <- c()
  numb <- c()

  for(i in unique(FinalDF_differential$Leading.razor.protein)){
    numb <- c(numb,sum(FinalDF_differential$Leading.razor.protein == i))
    prot <- c(prot,i)
  }

  df_lim <- as.data.frame(prot)
  df_lim$numn <- numb

  df_lim <- df_lim %>% filter(numb > 4)

  FinalDF_differential <- FinalDF_differential[order(-FinalDF_differential$Apex.intensity),]



  downgrade <- c()
  for(i in df_lim$prot){
    df_prot_lim <- FinalDF_differential %>% filter(Leading.razor.protein == i)
    downgrade <- c(downgrade,df_prot_lim$Modified[5:nrow(df_prot_lim)])

  }

  FinalDF$Priority[FinalDF$Modified %in% downgrade] <- 1


  return(FinalDF)
}

reformat <- function(FinalDF,scout){
  scout$seqcharge <- paste0(scout$Sequence, scout$Charge)
  scout <- scout %>% filter(PEP < .04)
  scout <- scout %>% filter(PIF > .75)

  FinalDF$Priority[FinalDF$Priority ==3] <- 2
  FinalDF$Priority[FinalDF$Modified %in% scout$seqcharge ] <- 3


  FinalDF_differential <- FinalDF %>% filter(Priority == 3)

  prot <- c()
  numb <- c()

  for(i in unique(FinalDF_differential$Leading.razor.protein)){
    numb <- c(numb,sum(FinalDF_differential$Leading.razor.protein == i))
    prot <- c(prot,i)
  }

  df_lim <- as.data.frame(prot)
  df_lim$numn <- numb

  df_lim <- df_lim %>% filter(numb > 4)

  FinalDF_differential <- FinalDF_differential[order(-FinalDF_differential$Apex.intensity),]



  downgrade <- c()
  for(i in df_lim$prot){
    df_prot_lim <- FinalDF_differential %>% filter(Leading.razor.protein == i)
    downgrade <- c(downgrade,df_prot_lim$Modified[5:nrow(df_prot_lim)])

  }

  FinalDF$Priority[FinalDF$Modified %in% downgrade] <- 1


  return(FinalDF)

}
