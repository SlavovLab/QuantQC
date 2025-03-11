#' Link Manual Annotations to raw data.
#'
#' This function takes a QuantQC (QQC) object as an input and return a QQC object
#' that has sample meta data linked in the meta.data slot. Its only to be used if
#' cellenONE files are not available
#'
#' @param QQC A QuantQC object.
#' @return \code{QQC} object with linked sample meta data
#' @examples
#' link_manual_Raw(TestSamples)
#' @export
link_manual_Raw <- function(QQC){
  linker <- QQC@meta.data


  linker <- melt(linker, ids = c('Run','well','plate'))



  peptide_data <- QQC@matricies@peptide
  # Get list of unique cell IDs
  cellID <- colnames(peptide_data)
  cellID <- as.data.frame(cellID)
  colnames(cellID) <- 'ID'

  cellenOne_data_small <- cellenOne_data %>% dplyr::select(any_of(c('ID','diameter','sample','label','injectWell','plate')))
  cellenOne_data_small <- as.data.frame(cellenOne_data_small)


  cellID <- cellID %>% left_join(cellenOne_data_small,by = c('ID'))

  cellID$sample[is.na(cellID$sample)==T] <- 'neg'

  cellID$prot_total <- log2(colSums(peptide_data[,1:ncol(peptide_data)],na.rm = T))

  QQC@cellenONE.meta <- cellenOne_data

  cellID$WP <- paste0(cellID$plate,cellID$injectWell)
  cellID$plate <- NULL
  linker$WP <- paste0(linker$plate,linker$Well)

  cellID <- cellID %>% left_join(linker, by = c('WP'))
  cellID$WP <- NULL

  QQC@meta.data <- cellID

  return(QQC)

}



#' Link CellenONE metadata
#'
#' This function takes a QuantQC (QQC) object and file paths to cellenONE isolation files
#' as a named list.
#'
#' @param QQC A QuantQC object
#' @param allCells A named list of file paths to cellenONE isolation files
#' @return A QQC object with \code{QQC at meta.data} slot and \code{cellenONE.meta} slot linking raw data to meta data from cellenONE
#' @examples
#' link_cellenONE_Raw(TestSamples,SortData_3_7)
#' @export
link_cellenONE_Raw <- function(QQC,cells_file){

  linker <- QQC@meta.data



  if(QQC@ms_type == 'DDA'){
    cellenOne_data <- analyzeCellenONE_TMT(cells_file,QQC@misc[['plex']])
  }
  if(QQC@ms_type == 'DIA' | QQC@ms_type =='DIA_C'){
    cellenOne_data <- analyzeCellenONE_mTRAQ(cells_file,QQC@misc[['plex']])
  }

  peptide_data <- QQC@matricies@peptide
  # Get list of unique cell IDs
  cellID <- colnames(peptide_data)
  cellID <- as.data.frame(cellID)
  colnames(cellID) <- 'ID'

  if(sum(colnames(cellenOne_data)=='Intensity') == 0){
    cellenOne_data_small <- cellenOne_data %>% dplyr::select(any_of(c('ID','diameter','sample','label','injectWell','plate')))
    cellenOne_data_small <- as.data.frame(cellenOne_data_small)
  }
  if(sum(colnames(cellenOne_data)=='Intensity') == 1){
    cellenOne_data_small <- cellenOne_data %>% dplyr::select(any_of(c('ID','diameter','sample','label','injectWell','plate','Intensity','Stain_Diameter')))
    cellenOne_data_small <- as.data.frame(cellenOne_data_small)
  }

  #cellenOne_data_small <- cellenOne_data_small %>% filter(injectWell %in%  QQC@raw_data$Well)
  #cellenOne_data_small <- cellenOne_data_small %>% filter(plate %in%  QQC@raw_data$plate)


  cellenOne_data_small <- cellenOne_data_small %>% distinct(ID,.keep_all = T)

  cellID <- cellID %>% dplyr::left_join(cellenOne_data_small,by = c('ID'))

  cellID$sample[is.na(cellID$sample)==T] <- 'neg'

  cellID$prot_total <- log2(colSums(peptide_data,na.rm = T))

  QQC@cellenONE.meta <- cellenOne_data

  cellID$WP <- paste0(cellID$plate,cellID$injectWell)
  cellID$plate <- NULL
  linker$WP <- paste0(linker$plate,linker$Well)

  cellID <- cellID %>% left_join(linker, by = c('WP'))
  cellID$WP <- NULL

  QQC@meta.data <- cellID

  return(QQC)

}



analyzeCellenONE_TMT <- function(cells_file,plex){
  #cells_file <- all_cells
  #plex = 32
  # Code to parse cellenONE files and map cell diameters, a mess and not too important,
  # dont feel obligeted to read
  for(i in 1:length(cells_file)){
    df1 <- read.delim(cells_file[[i]])
    df1$condition <- names(cells_file[i])

    if(i == 1){
      df <- df1
    }else{
      df <- rbind(df,df1)
    }

  }

  cells_file <- df


  #file_paths
  if(plex == 14){
    #File paths to pickup/label files

    # 2plex
    labelPath <- system.file("extdata", "14plex_files/Labels.fld", package = "QuantQC")
    pickupPath1 <-  system.file("extdata", "14plex_files/Pickup_mock.fld", package = "QuantQC")

  }
  if(plex == 12){
    #File paths to pickup/label files
    # 2plex
    labelPath <- system.file("extdata", "12plex_files/Labels.fld", package = "QuantQC")
    pickupPath1 <-  system.file("extdata", "12plex_files/Pickup_mock.fld", package = "QuantQC")

  }

  if(plex == 29){

    # 29plex
    labelPath <- system.file("extdata", "29plex_files/Labels.fld", package = "QuantQC")
    pickupPath1 <-  system.file("extdata", "29plex_files/Pickup_mock.fld", package = "QuantQC")
  }

  if(plex == 32){

    # 29plex
    labelPath <- system.file("extdata", "32plex_files/Labels.fld", package = "QuantQC")
    pickupPath1 <-  system.file("extdata", "32plex_files/Pickup_mock.fld", package = "QuantQC")
  }

  cells_file[grepl("Transmission",cells_file$X),]$X <- NA
  cells_file <- cells_file %>% fill(2:7, .direction = "up") %>% drop_na(XPos)

  if(sum(cells_file$X == 'Blue') > 0){
    cells_file[grepl("Blue",cells_file$X),]$X <- NA
    cells_file <- cells_file %>% filter(is.na(X) == F)
  }
  store <- F
  if(sum(cells_file$X == 'Green') > 0){
    store <- T
    get_green <- seq(2, nrow(cells_file), by = 2)
    green <- cells_file[get_green,]
    cells_file[grepl("Green",cells_file$X),]$X <- NA
    cells_file <- cells_file %>% filter(is.na(X) == F)
  }




  #cells_file1 <- cells_file %>% filter(XPos > 50)
  #cells_file1$XPos <- cells_file1$XPos + 1
  #cells_file2<- cells_file %>% filter(XPos < 50)
  #cells_file <- rbind(cells_file2,cells_file1)


  #### Labelling and Pickup Field Files
  ## labelling file

  con_lab <-file(labelPath)
  lines_lab<- readLines(con_lab)
  close(con_lab)
  slines_lab <- strsplit(lines_lab,"\t")
  colCount_lab <- max(unlist(lapply(slines_lab, length)))

  label <- read.table(labelPath, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 22)

  fieldRaw_label <- label[grepl("\\[\\d",label$position),]$position
  fieldBrack_label <- gsub("\\[|\\]", "", fieldRaw_label)
  fieldComma_label <- strsplit(fieldBrack_label, "," )
  fieldNum_label <- unlist(lapply(fieldComma_label, `[[`, 1))
  label$field[grepl("\\[\\d",label$position)] <- fieldNum_label
  label <- label %>% fill(field, .direction = "down")
  label <-  label[!label$well=="",]
  label$field <- as.numeric(label$field) + 1

  labelxyPos <- strsplit(label$position, "\\/")
  label$yPos <- unlist(lapply(labelxyPos, '[[', 1))
  label$xPos <- unlist(lapply(labelxyPos, '[[', 2))



  ## sample pickup file

  con_pickup <-file(pickupPath1)
  lines_pickup<- readLines(con_pickup)
  close(con_pickup)
  slines_pickup <- strsplit(lines_pickup,"\t")
  colCount_pickup <- max(unlist(lapply(slines_pickup, length)))

  pickup <- read.table(pickupPath1, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 22)


  fieldRaw_pickup <- pickup[grepl("\\[\\d",pickup$position),]$position
  fieldBrack_pickup <- gsub("\\[|\\]", "",fieldRaw_pickup)
  fieldComma_pickup <- strsplit(fieldBrack_pickup, "," )
  fieldNum_pickup <- unlist(lapply(fieldComma_pickup, `[[`, 1))
  pickup$field[grepl("\\[\\d",pickup$position)] <- fieldNum_pickup
  pickup <- pickup %>% fill(field, .direction = "down")
  pickup <-  pickup[!pickup$well=="",]
  pickup$field <- as.numeric(pickup$field) + 1

  ## trying to get pickup working in the same way
  # fixing pickup file
  pickup$target <- 1

  pickup$well <- substring(pickup$well, 2)
  pickup$well <- gsub(",","",pickup$well)



  pickupxyPos <- strsplit(pickup$position, "\\/")
  pickup$yPos <- unlist(lapply(pickupxyPos, '[[', 1))
  pickup$xPos <- unlist(lapply(pickupxyPos, '[[', 2))
  pickup$xPos <- as.numeric(pickup$xPos)
  pickup$yPos <- as.numeric(pickup$yPos)

  ### making sure that label "slide" refers to fields or actual slides.
  ## There is a sequence of 108 x,y and there are 108 spots per slide
  # order is y/x positions

  label$yPos <- as.numeric(label$yPos)
  label$xPos <- as.numeric(label$xPos)
  label$well[label$well == '1H1,'] <- '1G25,'
  label$well[label$well == '1H2,'] <- '1G26,'
  label$well[label$well == '1H3,'] <- '1G27,'
  label$well[label$well == '1H4,'] <- '1G28,'
  label$well[label$well == '1H5,'] <- '1G29,'
  if(plex == 32){
    label$well[label$well == '1P1,'] <- '1G25,'
    label$well[label$well == '1P2,'] <- '1G26,'
    label$well[label$well == '1P3,'] <- '1G27,'
    label$well[label$well == '1P4,'] <- '1G28,'
    label$well[label$well == '1P5,'] <- '1G29,'
    label$well[label$well == '1P6,'] <- '1G30,'
    label$well[label$well == '1P7,'] <- '1G31,'
    label$well[label$well == '1P8,'] <- '1G32,'

  }
  label$well <- substring(label$well, 3)
  label$well <- as.numeric(gsub(",","",label$well))
  matchTMTSCP <- paste0("TMT", 1:length(unique(label$well)))



  for (i in 1:length(matchTMTSCP)) {

    label[which(label$well == i),]$well <- matchTMTSCP[i]

  }

  label$well <- substring(label$well, 4)


  label <- label %>% filter(field %in% unique(cells_file$Field))
  pickup <- pickup %>% filter(field %in% unique(cells_file$Field))

  #### Trying to map label to cell
  ###  lets try and keep all three important in one
  cells_file <- transform(cells_file, xyf = paste0(XPos, YPos, Field))
  label <- transform(label, xyf = paste0(xPos, yPos, field))


  ### Isolation and Label merged

  isoLab <- cells_file %>% group_by(Target) %>% left_join(label, by = 'xyf')
  labelCount <- isoLab %>% group_by(well) %>% dplyr::summarize(count=n())


  ### clustering together pickup points with sample points using ANN
  ## super convolutedly/idiotically

  isoLab$ann <- NA
  isoLab$pickupX <- NA
  isoLab$pickupY <- NA


  if(plex == 32){
    isoLab <- isoLab %>% filter(is.na(yPos)==F)
    isoLab$yPos <- as.numeric(isoLab$yPos)
    isoLab$yPos[isoLab$yPos > 50] <- isoLab$yPos[isoLab$yPos > 50] + 10
    isoLab$yPos[isoLab$yPos > 33] <- isoLab$yPos[isoLab$yPos > 33] + 10
    isoLab$yPos[isoLab$yPos > 17] <- isoLab$yPos[isoLab$yPos > 17] + 10

    pickup$yPos[pickup$yPos > 50] <- pickup$yPos[pickup$yPos > 50] + 10
    pickup$yPos[pickup$yPos > 33] <- pickup$yPos[pickup$yPos > 33] + 10
    pickup$yPos[pickup$yPos > 17] <- pickup$yPos[pickup$yPos > 17] + 10

  }

  ann_ <- yaImpute::ann(ref = as.matrix(unique(pickup[, c("xPos","yPos")])),  target = as.matrix(isoLab[ , c("xPos","yPos")]), k=1)
  isoLab$ann <-  ann_$knnIndexDist[,1]

  if(plex == 32){


    isoLab$yPos[isoLab$yPos > 17] <- isoLab$yPos[isoLab$yPos > 17] - 10
    isoLab$yPos[isoLab$yPos > 33] <- isoLab$yPos[isoLab$yPos > 33] - 10
    isoLab$yPos[isoLab$yPos > 50] <- isoLab$yPos[isoLab$yPos > 50] - 10

    pickup$yPos[pickup$yPos > 17] <- pickup$yPos[pickup$yPos > 17] - 10
    pickup$yPos[pickup$yPos > 33] <- pickup$yPos[pickup$yPos > 33] - 10
    pickup$yPos[pickup$yPos > 50] <- pickup$yPos[pickup$yPos > 50] - 10

  }


  isoLab_new <- unique(pickup[, c("xPos","yPos")])


  isoLab$pickupX <- isoLab_new[isoLab$ann,]$xPos
  isoLab$pickupY <- isoLab_new[isoLab$ann,]$yPos
  isoLab_bound <-isoLab



  ### Merge pickup and isoLab
  isoLab_bound <- transform(isoLab_bound, xyft = paste0(pickupX, pickupY, Field, Target))
  pickup <-  transform(pickup, xyft = paste0(xPos, yPos, field, target))


  isoLab_final <- isoLab_bound %>% left_join(pickup, by = 'xyft')
  wellCount <- isoLab_final %>% group_by(well.y) %>% dplyr::summarize(count=n())

  isoLab_final$ImageFile <- str_sub(isoLab_final$ImageFile, end = -2)
  isoLab_final$ImageFile <- str_sub(isoLab_final$ImageFile,12)

  ### Clean up to yield final dataframe
  cellenOne_data <- data.frame(sample = isoLab_final$condition, isoTime = isoLab_final$Time, diameter = isoLab_final$Diameter, elongation = isoLab_final$Elongation, slide = isoLab_final$Target, field = isoLab_final$Field, dropXPos = isoLab_final$XPos, dropYPos = isoLab_final$YPos, label = isoLab_final$well.x, pickupXPos = isoLab_final$pickupX, pickupYPos = isoLab_final$pickupY, injectWell = isoLab_final$well.y, picture = isoLab_final$ImageFile)


  ##*sigh* not done yet

  cellenOne_data$wellAlph <- substring(cellenOne_data$injectWell, 1, 1)
  cellenOne_data$wellNum <- substring(cellenOne_data$injectWell, 2)

  cellenOne_data <- cellenOne_data %>% arrange(wellAlph, as.numeric(wellNum), as.numeric(label))
  #cellenOne_data <- cellenOne_data %>% dplyr::select(dropYPos,dropXPos,sample,diameter,elongation,field,label,injectWell)

  cellenOne_data$pickupXPos_numb <- as.numeric(cellenOne_data$pickupXPos)
  cellenOne_data$pickupYPos_numb  <- as.numeric(cellenOne_data$pickupYPos)


  # Assigning TMT tags to wells where label was picked up out of during prep
  # Each tag was dispensed multiple times so maps to multiple wells of plate
  if(plex == 14){
    cellenOne_data$label <- paste0('Reporter.intensity.' , (as.numeric(cellenOne_data$label) + 4))

  }
  if(plex == 29){

    labs_map <- c('Reporter.intensity.4', 'Reporter.intensity.5','Reporter.intensity.6',
              'Reporter.intensity.7', 'Reporter.intensity.8','Reporter.intensity.9',
              'Reporter.intensity.10', 'Reporter.intensity.11','Reporter.intensity.12',
              'Reporter.intensity.13','Reporter.intensity.14','Reporter.intensity.15',
              'Reporter.intensity.16','Reporter.intensity.17','Reporter.intensity.18',
              'Reporter.intensity.19','Reporter.intensity.20','Reporter.intensity.21',
              'Reporter.intensity.22','Reporter.intensity.23','Reporter.intensity.24',
              'Reporter.intensity.25','Reporter.intensity.26','Reporter.intensity.27',
              'Reporter.intensity.28','Reporter.intensity.29','Reporter.intensity.30',
              'Reporter.intensity.31','Reporter.intensity.32')

    cellenOne_data$label <- as.numeric(cellenOne_data$label)
    cellenOne_data$label_new <- NA
    for(i in 1:29){
      cellenOne_data$label_new[cellenOne_data$label == i] <- labs_map[i]

    }
    cellenOne_data$label <- cellenOne_data$label_new
    cellenOne_data$label_new <- NULL

    #cellenOne_data$label <- paste0('Reporter.intensity.' , (as.numeric(cellenOne_data$label) + 3))



  }

  if(plex == 32){

    labs_map <- c('Reporter.intensity.4', 'Reporter.intensity.5','Reporter.intensity.6',
                  'Reporter.intensity.7', 'Reporter.intensity.8','Reporter.intensity.9',
                  'Reporter.intensity.10', 'Reporter.intensity.11','Reporter.intensity.12',
                  'Reporter.intensity.13','Reporter.intensity.14','Reporter.intensity.15',
                  'Reporter.intensity.16','Reporter.intensity.17','Reporter.intensity.18',
                  'Reporter.intensity.19','Reporter.intensity.20','Reporter.intensity.21',
                  'Reporter.intensity.22','Reporter.intensity.23','Reporter.intensity.24',
                  'Reporter.intensity.25','Reporter.intensity.26','Reporter.intensity.27',
                  'Reporter.intensity.28','Reporter.intensity.29','Reporter.intensity.30',
                  'Reporter.intensity.31','Reporter.intensity.32','Reporter.intensity.33',
                  'Reporter.intensity.34','Reporter.intensity.35')

    cellenOne_data$label <- as.numeric(cellenOne_data$label)
    cellenOne_data$label_new <- NA
    for(i in 1:32){
      cellenOne_data$label_new[cellenOne_data$label == i] <- labs_map[i]

    }
    cellenOne_data$label <- cellenOne_data$label_new
    cellenOne_data$label_new <- NULL

    cellenOne_data <- cellenOne_data %>% filter(is.na(label)==F)

    #cellenOne_data$label <- paste0('Reporter.intensity.' , (as.numeric(cellenOne_data$label) + 3))



  }

  # only one plate TMT
  cellenOne_data$plate <- 1

  # Create cell ID to match MaxQuant report
  cellenOne_data$ID <- paste0(cellenOne_data$injectWell,cellenOne_data$plate,cellenOne_data$label)


  if(store == T){

    cellenOne_data$match_stain <- paste0('F',cellenOne_data$field,'X',cellenOne_data$dropXPos,'Y',cellenOne_data$dropYPos)
    green$match_stain <- paste0('F',green$Field,'X',green$XPos,'Y',green$YPos)
    green <- green %>% dplyr::select(match_stain,Intensity,Diameter)
    colnames(green)[3] <- 'Stain_Diameter'
    cellenOne_data <- cellenOne_data %>% left_join(green, by = c('match_stain'))

  }

  return(cellenOne_data)


}


analyzeCellenONE_mTRAQ <- function(cells_file,plex){

  for(i in 1:length(cells_file)){
    df1 <- read.delim(cells_file[[i]])
    df1$condition <- names(cells_file[i])

    if(i == 1){
      df <- df1
    }else{
      df <- rbind(df,df1)
    }

  }

  df <- df %>% filter(X != 'Green')
  df <- df %>% filter(X != '0')
  cells_file <- df


  # Code to parse cellenONE files and map cell diameters, a mess and not too important,
  # dont feel obligeted to read

  if(plex == 2){
    #File paths to pickup/label files

    # 2plex
    labelPath <- system.file("extdata", "2plex_files/Labels.fld", package = "QuantQC")
    pickupPath1 <- system.file("extdata", "2plex_files/Pickup_1_mock.fld", package = "QuantQC")
    pickupPath2 <- system.file("extdata", "2plex_files/Pickup_2_mock.fld", package = "QuantQC")

  }
  if(plex == 3){

    # 3plex
    labelPath <- system.file("extdata", "3plex_files/Labels.fld", package = "QuantQC")
    pickupPath1 <- system.file("extdata", "3plex_files/Pickup_mock_1.fld", package = "QuantQC")
    pickupPath2 <- system.file("extdata", "3plex_files/Pickup_mock_2.fld", package = "QuantQC")
  }


  cells_file[grepl("Transmission",cells_file$X),]$X <- NA
  cells_file <- cells_file %>% fill(2:7, .direction = "up") %>% drop_na(XPos)

  cells_file <- cells_file %>% filter(YPos < 69)
  #### Labelling and Pickup Field Files
  ## labelling file

  con_lab <-file(labelPath)
  lines_lab<- readLines(con_lab)
  close(con_lab)
  slines_lab <- strsplit(lines_lab,"\t")
  colCount_lab <- max(unlist(lapply(slines_lab, length)))

  label <- read.table(labelPath, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 22)

  fieldRaw_label <- label[grepl("\\[\\d",label$position),]$position
  fieldBrack_label <- gsub("\\[|\\]", "", fieldRaw_label)
  fieldComma_label <- strsplit(fieldBrack_label, "," )
  fieldNum_label <- unlist(lapply(fieldComma_label, `[[`, 1))
  label$field[grepl("\\[\\d",label$position)] <- fieldNum_label
  label <- label %>% fill(field, .direction = "down")
  label <-  label[!label$well=="",]
  label$field <- as.numeric(label$field) + 1

  labelxyPos <- strsplit(label$position, "\\/")
  label$yPos <- unlist(lapply(labelxyPos, '[[', 1))
  label$xPos <- unlist(lapply(labelxyPos, '[[', 2))

  ## sample pickup file

  con_pickup <-file(pickupPath1)
  lines_pickup<- readLines(con_pickup)
  close(con_pickup)
  slines_pickup <- strsplit(lines_pickup,"\t")
  colCount_pickup <- max(unlist(lapply(slines_pickup, length)))

  pickup <- read.table(pickupPath1, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 22)



  con_pickup <-file(pickupPath2)
  lines_pickup<- readLines(con_pickup)
  close(con_pickup)
  slines_pickup <- strsplit(lines_pickup,"\t")
  colCount_pickup <- max(unlist(lapply(slines_pickup, length)))
  pickup2 <- read.table(pickupPath2, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 27)

  pickup <- rbind(pickup,pickup2)

  fieldRaw_pickup <- pickup[grepl("\\[\\d",pickup$position),]$position
  fieldBrack_pickup <- gsub("\\[|\\]", "",fieldRaw_pickup)
  fieldComma_pickup <- strsplit(fieldBrack_pickup, "," )
  fieldNum_pickup <- unlist(lapply(fieldComma_pickup, `[[`, 1))
  pickup$field[grepl("\\[\\d",pickup$position)] <- fieldNum_pickup
  pickup <- pickup %>% fill(field, .direction = "down")
  pickup <-  pickup[!pickup$well=="",]
  pickup$field <- as.numeric(pickup$field) + 1



  pickupxyPos <- strsplit(pickup$position, "\\/")
  pickup$yPos <- unlist(lapply(pickupxyPos, '[[', 1))
  pickup$xPos <- unlist(lapply(pickupxyPos, '[[', 2))


  ### making sure that label "slide" refers to fields or actual slides.
  ## There is a sequence of 108 x,y and there are 108 spots per slide
  # order is y/x positions
  LabelxyPos <- strsplit(label$position, "\\/")
  label$yPos <- unlist(lapply(LabelxyPos, '[[', 1))
  label$xPos <- unlist(lapply(LabelxyPos, '[[', 2))
  label$well <- substring(label$well, 3)
  label$well <- as.numeric(gsub(",","",label$well))
  matchTMTSCP <- paste0("TMT", 1:length(unique(label$well)))

  for (i in 1:length(matchTMTSCP)) {

    label[which(label$well == i),]$well <- matchTMTSCP[i]

  }

  label$well <- substring(label$well, 4)


  label <- label %>% filter(field %in% unique(cells_file$Field))
  pickup <- pickup %>% filter(field %in% unique(cells_file$Field))


  #### Trying to map label to cell
  ###  lets try and keep all three important in one
  cells_file <- transform(cells_file, xyf = paste0(XPos, YPos, Field))
  label <- transform(label, xyf = paste0(xPos, yPos, field))

  ### Isolation and Label merged

  isoLab <- cells_file %>% group_by(Target) %>% left_join(label, by = 'xyf')
  labelCount <- isoLab %>% group_by(well) %>% dplyr::summarize(count=n())

  ## trying to get pickup working in the same way
  # fixing pickup file
  pickup$target <- 1

  pickup$well <- substring(pickup$well, 2)
  pickup$well <- gsub(",","",pickup$well)

  #pickup <- pickup %>% filter(!well %in% No_use$Wells_no )


  slidesUsedForPrep <- unique(isoLab$Target)

  #pickup <- pickup[which(pickup$target %in% slidesUsedForPrep),]

  ### clustering together pickup points with sample points using ANN
  ## super convolutedly/idiotically

  isoLab$ann <- NA
  isoLab$pickupX <- NA
  isoLab$pickupY <- NA

  ann_ <- ann(ref = as.matrix(unique(pickup[, c("xPos","yPos")])),  target = as.matrix(isoLab[ , c("xPos","yPos")]), k=1)


  isoLab$ann <-  ann_$knnIndexDist[,1]


  isoLab_new <- unique(pickup[, c("xPos","yPos")])

  isoLab$pickupX <- isoLab_new[isoLab$ann,]$xPos
  isoLab$pickupY <- isoLab_new[isoLab$ann,]$yPos

  isoLab_bound <-isoLab
  #isoLab_bound <- rbind(isoLab_123, isoLab_4)



  ### Merge pickup and isoLab
  isoLab_bound <- transform(isoLab_bound, xyft = paste0(pickupX, pickupY, Field, Target))
  pickup <-  transform(pickup, xyft = paste0(xPos, yPos, field, target))

  #intersect(pickup$xyft,isoLab_bound$xyft)

  isoLab_final <- isoLab_bound %>% left_join(pickup, by = 'xyft')
  wellCount <- isoLab_final %>% group_by(well.y) %>% dplyr::summarize(count=n())



  ### Clean up to yield final dataframe
  cellenOne_data <- data.frame (sample = isoLab_final$condition, isoTime = isoLab_final$Time, diameter = isoLab_final$Diameter, elongation = isoLab_final$Elongation, slide = isoLab_final$Target, field = isoLab_final$Field, dropXPos = isoLab_final$XPos, dropYPos = isoLab_final$YPos, label = isoLab_final$well.x, pickupXPos = isoLab_final$pickupX, pickupYPos = isoLab_final$pickupY, injectWell = isoLab_final$well.y, picture = isoLab_final$ImageFile)



  ##*sigh* not done yet

  cellenOne_data$wellAlph <- substring(cellenOne_data$injectWell, 1, 1)
  cellenOne_data$wellNum <- substring(cellenOne_data$injectWell, 2)

  cellenOne_data <- cellenOne_data %>% arrange(wellAlph, as.numeric(wellNum), as.numeric(label))
  #cellenOne_data <- cellenOne_data %>% dplyr::select(dropYPos,dropXPos,sample,diameter,elongation,field,label,injectWell)
  cellenOne_data$pickupXPos_numb <- as.numeric(cellenOne_data$pickupXPos)
  cellenOne_data$pickupYPos_numb  <- as.numeric(cellenOne_data$pickupYPos)


  # if using plex/scopeDIA samples are injected into 2 different plates,
  # there are 16 fields and samples from 1-8 go into plate 1, 9-16 go into plate 2
  cellenOne_data$plate <- NA
  cellenOne_data$plate[cellenOne_data$field > 8] <- 2
  cellenOne_data$plate[cellenOne_data$field <= 8] <- 1


  # Assigning mTRAQ tags to wells where label was picked up out of during prep
  # Each tag was dispensed multiple times so maps to multiple wells of plate
  cellenOne_data$tag <- NA

  if(plex == 2){
    cellenOne_data$tag[cellenOne_data$label %in% c('1','3','5','7')] <- '0'
    cellenOne_data$tag[cellenOne_data$label %in% c('2','4','6','8')] <- '4'
  }
  if(plex == 3){
    cellenOne_data$tag[cellenOne_data$label %in% c('1','4','7','10')] <- '0'
    cellenOne_data$tag[cellenOne_data$label %in% c('2','5','8','11')] <- '4'
    cellenOne_data$tag[cellenOne_data$label %in% c('3','6','9','12')] <- '8'
  }

  cellenOne_data$label <- cellenOne_data$tag
  # Create cell ID to match DIA-NN report
  cellenOne_data$ID <- paste0(cellenOne_data$injectWell,cellenOne_data$plate,'.',cellenOne_data$tag)

  return(cellenOne_data)
}





#' Plot cellenONE droplets by cell types
#'
#' This function takes plots the positions of the single cell droplets over glass slide colored
#' by the type of sample sorted. Wells where sample is placed are written in text in center of drop array.
#'
#' @param QQC A QQC object
#' @return A plot
#' @examples
#' PlotSlideLayout_celltype(TestSamples)
#' @export
PlotSlideLayout_celltype <- function(QQC){

  ggplot(QQC@cellenONE.meta) +
    geom_point(aes(x = dropXPos,y = dropYPos,color = sample)) +
    geom_text(aes(x = pickupXPos_numb,y = pickupYPos_numb,label = injectWell,size = 5),hjust= .5, vjust=-.6) +
    facet_wrap(~field,ncol = 4)+
    scale_y_reverse()+
    theme(plot.title = element_text(hjust = .5,size = 22),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          legend.text = element_text(size = 20))

}


#' Plot cellenONE droplets by labels
#'
#' This function takes plots the positions of the single cell droplets over glass slide colored
#' by the label used.
#'
#' @param QQC A QQC object
#' @return A plot
#' @examples
#' PlotSlideLayout_celltype(TestSamples)
#' @export
PlotSlideLayout_label <- function(QQC){

  # print the mTRAQ labels overlayed on the positions of the slide
  ggplot(QQC@cellenONE.meta, aes(x = dropXPos,y = dropYPos,color = label)) +
    geom_point() +scale_y_reverse()+ facet_wrap(~field,ncol = 4)+
    theme(plot.title = element_text(hjust = .5,size = 22),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          legend.text = element_text(size = 20))

}

