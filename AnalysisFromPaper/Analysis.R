
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





# single cell proc
## ------------------------------------------------------------------------------------------------------------


devtools::install_github("https://github.com/Andrew-Leduc/QuantQC")
library(QuantQC)


## Data files,  UPDATE FILE PATHS

#searched SC data, find on massive repo
data_path <- "path/report.tsv"

# link raw file name to well plate, find in folder
linker_path <- "path/linker.csv"

## Read in cell isolation files from CellenONE and assign cell type, find in folder
one <-"path/appnote/Monocyte.xls"
two <- "path/PDAC.xls"
three <- "path/Melanoma.xls"

all_cells <- list(Monocyte = one,
                  PDAC = two,
                  Melanoma = three)



# input character vector of proteins you'd like to see visualized in report
#prot_vis_umap <- c('O00487')

output_path <- ".../path/QQC_DIA_Report.html"

# Generate the HTML one line
Gen_QQC_report_DIA(data_path = data_path,
                   linker_path = linker_path,
                   isolation = all_cells,
                   output_path = output_path,
                   plex_exp = 3,
                   carrier_used = F,
                   ChQ = .1)




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
AppNote <- FilterBadCells(AppNote, min_intens = 7)


# Compute the size relative to the carrier of each cell
PlotSCtoCarrierRatio(AppNote)


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


## also assigns louvain clusters
AppNote <- ComputeUMAP(AppNote)

## plots by cluster
PlotUMAP(AppNote, by = 'Condition')


## Plots peptide agreement across clusters and shows location of peptide mapped to protein
## sequence

# AppNote@misc[["Species"]] <- "Human" # or "Mouse" supported
# ProteinClustConsistency(AppNote, prot = 'Q02790', type = 'line')


## Color code umap by proteins

# FeatureUMAP(AppNote, prot = 'Q02750', imputed = F)


# Bulk to sc correlations
## ------------------------------------------------------------------------------------------------------------


library(viridis)
library(ggpointdensity)

# Proteins with good peptide correlations single cell data
Pos_pep_cor <- as.data.frame(AppNote@pep.cor)
Pos_pep_cor <- Pos_pep_cor %>% filter(Cor > .2)



LF <- read.delim('/Users/andrewleduc/Library/CloudStorage/GoogleDrive-leduc.an@husky.neu.edu/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Users/appnote/LF/report_2.tsv')
LF_mat <- diann::diann_maxlfq(LF,group.header = "Protein.Group",quantity.header = "Ms1.Area")
sect <- intersect(rownames(LF_mat),Pos_pep_cor$Protein )
LF_mat2 <- LF_mat[sect,]
for(i in 1:ncol(LF_mat)){
  LF_mat2[,i]<-LF_mat2[,i]/median(LF_mat2[,i],na.rm = T)
}

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

