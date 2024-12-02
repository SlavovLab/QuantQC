library(QuantQC)
library(ggpointdensity)
library(viridis)
library(ggpubr)
library(dplyr)
library(matrixStats)
library(reshape2)
library(stringr)


# Function used
# Compute total least squares fit
TLS <- function(vect1,vect2){

  vect1[vect1 == -Inf] <- NA
  vect1[vect1 == Inf] <- NA



  vect2[vect2 == -Inf] <- NA
  vect2[vect2 == Inf] <- NA

  int_x <- mean(vect1,na.rm=T)
  int_y <- mean(vect2,na.rm=T)

  vect1 <- vect1-int_x
  vect2 <- vect2-int_y

  mat <- cbind(vect1,vect2)

  mat <- mat[complete.cases(mat),]

  TLS_mat <- svd(mat)$v

  slope <- TLS_mat[1,1]/TLS_mat[2,1]

  int <- c(int_x,int_y)

  return(list(slope,int))

}

# link gene names and uniprot
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



## mTRAQ Labeling Eff test in single cells by pooled DDA Experiment

mTRAQ_Lab <- read.delim('/Volumes/My Passport/LabTest/combined/txt/evidence.txt')
mTRAQ_Lab <- mTRAQ_Lab %>% filter(Oxidation..M. == 0)
mTRAQ_Lab <- mTRAQ_Lab %>% filter(Acetyl..Protein.N.term. == 0)

# nTerm Labeling Eff = 99%
sum(mTRAQ_Lab$Var_mTRAQ.Nter0,mTRAQ_Lab$Var_mTRAQ.Nter4,mTRAQ_Lab$Var_mTRAQ.Nter8)/nrow(mTRAQ_Lab)

# Lyine labeling Eff = 90%
sum(str_sub(mTRAQ_Lab$Sequence,-1) == 'K')/sum(mTRAQ_Lab$Var_mTRAQ.Lys0+mTRAQ_Lab$Var_mTRAQ.Lys4,mTRAQ_Lab$Var_mTRAQ.Lys8)



# single cell proc
## ------------------------------------------------------------------------------------------------------------


devtools::install_github("https://github.com/Andrew-Leduc/QuantQC")
library(QuantQC)


## -------------------------------------------
#plexDIA Analysis

## Data files,  UPDATE FILE PATHS

#searched SC data, find on massive repo
data_path <- "/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Users/appnote/Searched_DDM/report.tsv"


# link raw file name to well plate, find in folder
linker_path <- "~/Desktop/Github/QuantQC/AnalysisFromPaper/plexDIA/linker.csv"

## Read in cell isolation files from CellenONE and assign cell type, find in folder
one <-"~/Desktop/Github/QuantQC/AnalysisFromPaper/plexDIA/Melanoma.xls"
two <- "~/Desktop/Github/QuantQC/AnalysisFromPaper/plexDIA/PDAC.xls"
three <- "~/Desktop/Github/QuantQC/AnalysisFromPaper/plexDIA/Monocyte.xls"

all_cells <- list(Melanoma = one,
                  PDAC = two,
                  Monocyte = three)



# input character vector of proteins you'd like to see visualized in report
#prot_vis_umap <- c('O00487')

output_path <- "~/Desktop/QQC_DIA_Report.html"

# Generate the HTML one line
Gen_QQC_report_DIA(data_path = data_path,
                   linker_path = linker_path,
                   isolation = all_cells,
                   output_path = output_path,
                   plex_exp = 3,
                   carrier_used = F,
                   ChQ = 1)




######------------------------------------------------------------------------
# Or you can run line by line just like seurat


#Generate nPOP object from raw data
AppNote <- DIANN_to_QQC(data_path,linker_path, plex = 3, carrier = F)


# Normalize single cell runs to reference channel,
# filter out data points over twice reference
# Generate turn long table format to Peptide X Cell matrix
AppNote <- cellXpeptide(AppNote, TQVal = 1, chQVal = .05)


## Mapping cellenONE meta data to raw data
AppNote <- link_cellenONE_Raw(AppNote,all_cells)

  #plot exp design on glass slide
  PlotSlideLayout_celltype(AppNote)
  PlotSlideLayout_label(AppNote)



## Calculate statistics to look for decrease in LC/MS performance
# over time
  AppNote <- Calculate_run_order_statistics(AppNote)

  # Plot MS Intensity Statistics
  PlotIntensityDrift(AppNote)

  #Plot Retention time Statistics
  PlotRTDrift(AppNote)




# Test negative controls, i.e. samples with no cell
AppNote <- EvaluateNegativeControls(AppNote)


# Make the classic neg ctrl plots
PlotNegCtrl(AppNote)


# filter bad cells based off above, put in log10 intensity
AppNote <- FilterBadCells(AppNote, min_intens = 7.5)



# Plot cell size vs intensity in MS, options to color code by" "Run order" or "sample"
PlotCellSizeVsIntensity(AppNote, type = 'sample')


# Collapse and give statistics options:
#1 = median relative peptide, 2 = maxLFQ... soon to be directLFQwhich is much faster

AppNote <- CollapseToProtein(AppNote, 1)



# plot protein numbers
PlotProtAndPep(AppNote)



PlotDataComplete(AppNote)



# Compute correlations between peptides mapping to prot
AppNote <- SharedPeptideCor(AppNote)

PlotPepCor(AppNote)

PlotMS1vMS2(AppNote)

# Impute, still old imputation very slow, soon will update with C++ implementation
AppNote <- KNN_impute(AppNote)


## currently does label and LC/MS run but options comming soon
AppNote <- BatchCorrect(AppNote,run = F,labels = T)



AppNote <- ComputePCA(AppNote)

## plot PCA options are "Run order" "Total protein" "Condition" "Label"
PlotPCA(AppNote, by = "Run order")
PlotPCA(AppNote, by = "Condition")

## also assigns louvain clusters
AppNote <- ComputeUMAP(AppNote)

## plots by cluster
PlotUMAP(AppNote, by = 'Condition')

View(AppNote@reductions$UMAP)
## Plots peptide agreement across clusters and shows location of peptide mapped to protein
## sequence

# AppNote@misc[["Species"]] <- "Human" # or "Mouse" supported
# ProteinClustConsistency(AppNote, prot = 'Q02790', type = 'line')


## Color code umap by proteins

# FeatureUMAP(AppNote, prot = 'Q02750', imputed = F)


convert <- Proc_fasta('/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Users/LK/FASTA/swissprot_human_20211005.fasta')
convert <- convert %>% filter(split_prot %in% test)


Membrane <- read.gmt('/Users/andrewleduc/Downloads/PLASMA_MEMBRANE.v2023.2.Hs.gmt')
Membrane <- Membrane %>% filter(gene %in% convert$split_gene)
Kinases <- read.csv('/Users/andrewleduc/Downloads/Kinases.csv')
Kinases <- Kinases %>% filter(V8 %in% test)
TFs <- read.csv('/Users/andrewleduc/Downloads/TFs.csv')
TFs <- TFs %>% filter(HGNC.symbol %in% convert$split_gene)
TFs <- TFs %>% filter(Is.TF. == 'Yes')
UBligases <-read.csv('/Users/andrewleduc/Downloads/Ubiquitin_ligases.csv')
UBligases <- UBligases %>% filter(Swiss.Prot %in% test)

df_prot_type <- as.data.frame(cbind(c(nrow(Membrane), nrow(Kinases), nrow(TFs), nrow(UBligases)),c(('Membrane'), ('Kinases'), ('TFs'), ('UBligases'))))
colnames(df_prot_type) <- c('Number','Type')
df_prot_type$Number <- as.numeric(df_prot_type$Number)
ggplot(df_prot_type, aes(x = Type, y = Number)) + dot_plot + geom_bar(stat = 'identity')+
  ylab("# Unique proteins") + xlab('')

# Bulk to sc correlations
## ------------------------------------------------------------------------------------------------------------


# Proteins with good peptide correlations single cell data
Pos_pep_cor <- as.data.frame(AppNote@pep.cor)
Pos_pep_cor <- Pos_pep_cor %>% filter(Cor > .2)

# Download data from massive repo
LF <- read.delim('/path/report_2.tsv')
LF_mat <- diann::diann_maxlfq(LF,group.header = "Protein.Group",quantity.header = "Ms1.Area")
sect <- intersect(rownames(LF_mat),Pos_pep_cor$Protein )
LF_mat2 <- LF_mat[sect,]
for(i in 1:ncol(LF_mat)){
  LF_mat2[,i]<-LF_mat2[,i]/median(LF_mat2[,i],na.rm = T)
}

colnames(LF_mat2) <- c('Monocyte_rep1','Monocyte_rep2','PDAC_rep1','PDAC_rep2','Melanoma_rep1','Melanoma_rep2')


Monocyte <- (LF_mat2[,1]+LF_mat2[,2])/2
PDAC <- (LF_mat2[,3]+LF_mat2[,4])/2
Melanoma <- (LF_mat2[,5]+LF_mat2[,6])/2


Melanoma_sc <- rowMedians(AppNote@matricies@protein[sect,colnames(AppNote@matricies@protein) %in% AppNote@meta.data$ID[AppNote@meta.data$sample == 'Melanoma']],na.rm = T)
PDAC_sc <- rowMedians(AppNote@matricies@protein[sect,colnames(AppNote@matricies@protein) %in% AppNote@meta.data$ID[AppNote@meta.data$sample == 'PDAC']],na.rm = T)
Monocyte_sc <- rowMedians(AppNote@matricies@protein[sect,colnames(AppNote@matricies@protein) %in% AppNote@meta.data$ID[AppNote@meta.data$sample == 'Monocyte']],na.rm = T)


## Make plots

vect_sc <- Melanoma_sc-PDAC_sc
df_plot <- as.data.frame(vect_sc)
df_plot$vect_bulk <- log2(Melanoma/PDAC)

ggplot(df_plot, aes(x = vect_bulk, y = vect_sc)) + geom_pointdensity() +
  scale_color_viridis() + dot_plot + xlab('Bulk Label Free') + ylab('Single Cell')+
  ggtitle('Log2(Melanoma/Monocyte)') + rremove('legend') +ylim(c(-4,4)) +xlim(c(-4,4))


cor(df_plot$vect_sc,df_plot$vect_bulk,use = 'pairwise.complete.obs', method = 'pearson')
TLS(df_plot$vect_sc,df_plot$vect_bulk)[[1]]



## -----------------------------------------------
#pSCoPE Analysis


#searched SC data, find on massive repo
data_path <- "/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Users/aleduc/TMT29/Protocol_final_data/All_1p_FDR/evidence.txt"

# link raw file name to well plate, find in the github repo /QuantQC/AnalysisFromPaper/pSCoPE/
linker_path <- "/Users/andrewleduc/Desktop/Github/QuantQC/AnalysisFromPaper/pSCoPE/linker.csv"

# Read in cell isolation files from CellenONE and assign cell type, find in folder /QuantQC/AnalysisFromPaper/pSCoPE/
one <-"/Users/andrewleduc/Desktop/Github/QuantQC/AnalysisFromPaper/pSCoPE/Monocyte_isolated.xls"
two <- "/Users/andrewleduc/Desktop/Github/QuantQC/AnalysisFromPaper/pSCoPE/PDAC_isolated.xls"
three <- "/Users/andrewleduc/Desktop/Github/QuantQC/AnalysisFromPaper/pSCoPE/Melanoma_isolated.xls"

all_cells <- list(Monocyte = one,
                  PDAC = two,
                  Melanoma = three)



# input character vector of proteins you'd like to see visualized in report
#prot_vis_umap <- c('O00487')

output_path <- "/Users/andrewleduc/Desktop/DDA_Report.html"

# Generate the HTML one line
Gen_QQC_report_DDA(data_path = data_path,
                   linker_path = linker_path,
                   isolation = all_cells,
                   output_path = output_path,
                   CV_thresh = .36,
                   plex = 29)



######------------------------------------------------------------------------
# Or you can run line by line just like seurat


#Generate nPOP object from raw data
AppNote <- MQ_to_QQC(data_path,linker_path, plex = 29,PIF_in = .2, PEP_in = 1)


# Normalize single cell runs to reference channel,
# filter out data points over twice reference
# Generate turn long table format to Peptide X Cell matrix

AppNote <- TMT_Reference_channel_norm(AppNote)


## Mapping cellenONE meta data to raw data
AppNote <- link_cellenONE_Raw(AppNote,all_cells)

#plot exp design on glass slide
PlotSlideLayout_celltype(AppNote)

PlotSlideLayout_label(AppNote)


## Calculate statistics to look for decrease in LC/MS performance
# over time
AppNote <- Calculate_run_order_statistics(AppNote)

# Plot MS Intensity Statistics
PlotIntensityDrift(AppNote)

#Plot Retention time Statistics
PlotRTDrift(AppNote)


# Test negative controls, i.e. samples with no cell
AppNote <- EvaluateNegativeControls(AppNote)


# Make the classic neg ctrl plots
PlotNegCtrl(AppNote,CV_thresh = .38)


# filter bad cells based off above, put in log10 intensity
AppNote <- FilterBadCells(AppNote, CV_thresh = .4)



# Compute the size relative to the carrier of each cell
PlotSCtoCarrierRatio(AppNote)


# Plot cell size vs intensity in MS, options to color code by" "Run order" or "sample"
PlotCellSizeVsIntensity(AppNote, type = "sample")


# Collapse and give statistics options:
#1 = median relative peptide, 2 = maxLFQ... soon to be directLFQwhich is much faster

AppNote <- CollapseToProtein(AppNote, 1)



# plot protein numbers
PlotProtAndPep(AppNote)



PlotDataComplete(AppNote)



# Compute correlations between peptides mapping to prot
AppNote <- SharedPeptideCor(AppNote)

PlotPepCor(AppNote)


# Impute, still old imputation very slow, soon will update with C++ implementation
AppNote <- KNN_impute(AppNote)


## currently does label and LC/MS run but options comming soon
AppNote <- BatchCorrect(AppNote,run = T,labels = F)


AppNote <- ComputePCA(AppNote,imputed = F)

## plot PCA options are "Run order" "Total protein" "Condition" "Label"
PlotPCA(AppNote, by = "Condition")

## also assigns louvain clusters
AppNote <- ComputeUMAP(AppNote)

## plots by cluster
PlotUMAP(AppNote)
PlotUMAP(AppNote, by = 'Condition')
FeatureUMAP(AppNote, prot = 'P09429') + theme_classic()



# Compare quant to bulk samples

# Download form massive repo
Bulk <- read.delim('/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Users/aleduc/TMT29/Protocol_final_data/Bulk_TMT/evidence.txt')
Bulk <- Bulk %>% filter(Raw.file != "dAL_1k_bulk_T_45k")
Bulk <- Bulk %>% filter( PEP < .07)
Bulk <- Bulk %>% filter( PIF > .7)
Bulk <- Bulk %>% filter( Potential.contaminant != '+')
Bulk <- Bulk %>% filter( Reverse != '+')
Bulk$seqcharge <- paste0(Bulk$Sequence,Bulk$Charge)


# Grab RIs from TMT labels used for bulk data
bulk_samples <- c('Reporter.intensity.6','Reporter.intensity.7','Reporter.intensity.12','Reporter.intensity.13',
  'Reporter.intensity.15','Reporter.intensity.16')

Bulk <- Bulk %>% select(Leading.razor.protein,any_of(bulk_samples))

# Normalize for first for loading and then to relative levels then convert to log
Bulk[,2:7] <- QuantQC::normalize(Bulk[,2:7], log = T)

# Collapse to protein level
Bulk <- reshape2::melt(Bulk, id.vars = 'Leading.razor.protein')
Bulk <- Bulk %>% group_by(Leading.razor.protein,variable) %>% summarise(value = median(value,na.rm=T))
Bulk <- reshape2::dcast(Bulk,Leading.razor.protein~variable, value.var = 'value')

Bulk$Leading.razor.protein <-  gsub(".*\\|([A-Za-z0-9]+)\\|.*", "\\1", Bulk$Leading.razor.protein)

Bulk <- Bulk %>% distinct(Leading.razor.protein,.keep_all = T)
rownames(Bulk) <- Bulk$Leading.razor.protein
Bulk$Leading.razor.protein <- NULL
Bulk <- as.matrix(Bulk)

# Final normalize to relative levels after protein collapse
for(i in 1:nrow(Bulk)){
  Bulk[i,] <- Bulk[i,] - mean(Bulk[i,],na.rm = T)
}


colnames(Bulk) <- c('Monocyte_rep1','PDAC_rep1','Monocyte_rep2','Melanoma_rep1','PDAC_rep2','Melanoma_rep2')
write.csv(Bulk,'/Users/andrewleduc/Desktop/Figshare_Protocol/TMT_1000SPD/bulk.csv')


sect <- intersect(rownames(Bulk), rownames(AppNote@matricies@protein))

Bulk <- Bulk[sect,]





Monocyte <- (Bulk[,1]+Bulk[,3])/2
Melanoma<- (Bulk[,2]+Bulk[,5])/2
PDAC <- (Bulk[,4]+Bulk[,6])/2


Melanoma_sc <- rowMeans(AppNote@matricies@protein[sect,colnames(AppNote@matricies@protein) %in% AppNote@meta.data$ID[AppNote@meta.data$sample == 'Melanoma']],na.rm = T)
PDAC_sc <- rowMeans(AppNote@matricies@protein[sect,colnames(AppNote@matricies@protein) %in% AppNote@meta.data$ID[AppNote@meta.data$sample == 'PDAC']],na.rm = T)
Monocyte_sc <- rowMeans(AppNote@matricies@protein[sect,colnames(AppNote@matricies@protein) %in% AppNote@meta.data$ID[AppNote@meta.data$sample == 'Monocyte']],na.rm = T)


## Make plots

vect_sc <- Monocyte_sc-PDAC_sc
df_plot <- as.data.frame(vect_sc)

df_plot$vect_bulk <- Monocyte - PDAC
df_plot$vect_bulk[df_plot$vect_bulk==Inf] <- NA
df_plot$vect_bulk[df_plot$vect_bulk==-Inf] <- NA

rat1 <- ggplot(df_plot, aes(x = vect_bulk, y = vect_sc)) + geom_pointdensity() +
  scale_color_viridis() + dot_plot + xlab('Bulk TMT') + ylab('Single Cell')+
  ggtitle('Log2(PDAC/Monocyte)') + rremove('legend') +ylim(c(-2.5,2.5)) + theme(plot.title = element_text(size = 14))+
  geom_abline(intercept = 0, slope = TLS(df_plot$vect_sc,df_plot$vect_bulk)[[1]], col = "red")

cor(df_plot$vect_sc,df_plot$vect_bulk,use = 'pairwise.complete.obs', method = 'pearson')
TLS(df_plot$vect_sc,df_plot$vect_bulk)[[1]]

#------------------------------------------------------------------

vect_sc <- Monocyte_sc-Melanoma_sc
df_plot <- as.data.frame(vect_sc)

df_plot$vect_bulk <- Monocyte-Melanoma
df_plot$vect_bulk[df_plot$vect_bulk==Inf] <- NA
df_plot$vect_bulk[df_plot$vect_bulk==-Inf] <- NA

rat2 <- ggplot(df_plot, aes(x = vect_bulk, y = vect_sc)) + geom_pointdensity() +
  scale_color_viridis() + dot_plot + xlab('Bulk TMT') + ylab('Single Cell')+
  ggtitle('Log2(Melanoma/Monocyte)') + rremove('legend') +ylim(c(-2.5,2.5))+
  theme(plot.title = element_text(size = 14))+
  geom_abline(intercept = 0, slope = TLS(df_plot$vect_sc,df_plot$vect_bulk)[[1]], col = "red")

cor(df_plot$vect_sc,df_plot$vect_bulk,use = 'pairwise.complete.obs', method = 'pearson')
TLS(df_plot$vect_sc,df_plot$vect_bulk)[[1]]

#------------------------------------------------------------------

vect_sc <- Melanoma_sc-PDAC_sc
df_plot <- as.data.frame(vect_sc)

df_plot$vect_bulk <- Melanoma-PDAC
df_plot$vect_bulk[df_plot$vect_bulk==Inf] <- NA
df_plot$vect_bulk[df_plot$vect_bulk==-Inf] <- NA

rat3 <- ggplot(df_plot, aes(x = vect_bulk, y = vect_sc)) + geom_pointdensity() +
  scale_color_viridis() + dot_plot + xlab('Bulk TMT') + ylab('Single Cell')+
  ggtitle('Log2(Melanoma/PDAC)') + rremove('legend') +ylim(c(-2.5,2.5)) + theme(plot.title = element_text(size = 14))+
  geom_abline(intercept = 0, slope = TLS(df_plot$vect_sc,df_plot$vect_bulk)[[1]], col = "red")

cor(df_plot$vect_sc,df_plot$vect_bulk,use = 'pairwise.complete.obs', method = 'pearson')
TLS(df_plot$vect_sc,df_plot$vect_bulk)[[1]]


rat1+rat2+rat3

numb_cells = c(1280,1536,3584,3712)
plex = c(' 2-plex',' 3-plex', '18-plex', '32-plex')
df <- as.data.frame(cbind(numb_cells,plex))
df$numb_cells <- as.numeric(numb_cells)
ggplot(df, aes(x = plex,y = numb_cells)) + geom_bar(stat = 'identity',width = .7)+dot_plot+
  ylim(c(0,4200)) + ylab('# of cells / single prep') + xlab('') + ggtitle('nPOP sample prep throughput')






dead <- AppNote@reductions$UMAP %>% filter(cluster %in% c(5,6,7))

alive <- AppNote@reductions$UMAP %>% filter(!cluster %in% c(5,6,7))

dead <- rowMeans(AppNote@matricies@protein[,rownames(alive)],na.rm = T) - rowMeans(AppNote@matricies@protein[,rownames(dead)],na.rm = T)


dead <- as.data.frame(dead)

FeatureUMAP(AppNote, prot = 'P20290') + theme_classic()

dead$prot <- rownames(dead)


Hum <- Proc_fasta('/Users/andrewleduc/Desktop/Github/QuantQC/inst/extdata/Human.fasta')
Mouse <- Proc_fasta('/Users/andrewleduc/Desktop/Github/QuantQC/inst/extdata/Mouse.fasta')
Mouse$split_gene <- toupper(Mouse$split_gene)


Mouse <- Mouse %>% filter(split_gene %in% Hum$split_gene)
Hum <- Hum %>% filter(split_gene %in% Mouse$split_gene)

Mouse <- Mouse[order(Mouse$split_gene),]
Hum <- Hum[order(Hum$split_gene),]

convert <- Mouse %>% left_join(Hum, by = c('split_gene'))
convert <- convert %>% filter(split_prot.y %in% dead$prot)
convert <- convert %>% filter(split_prot.x %in% mat_comp$prot)
convert <- convert %>% distinct(split_prot.y,.keep_all = T)

mat_comp <- mat %>% filter(prot %in% convert$split_prot.x)
dead <- dead %>% filter(prot %in% convert$split_prot.y)

dead <- dead[convert$split_prot.y,]
mat_comp <- mat_comp[convert$split_prot.x,]

dead$dead <- -dead$dead - median(dead$dead)

cor(mat_comp$bas_prot,dead$dead -1)

both_spec <- cbind(mat_comp$bas_prot,dead$dead -1)

colnames(both_spec) <- c('Mouse_Primary_Cells','Human_CellLines')
both_spec <- as.data.frame(both_spec)

ggplot(both_spec, aes(x = Mouse_Primary_Cells,y = Human_CellLines)) + geom_point() + dot_plot +
  ggtitle('Cor = 0.50')






