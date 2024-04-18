################################################################################
###############################      Author : Kyle R. Upton      ###############
################################################################################

### R version 4.3.2 (2023-10-31)

## Attributions: Adapted from StandR Demo scripts

################################################################################
################################### Required R packages ########################
# library(standR)
# library(SpatialExperiment)
# library(limma)
# library(ExperimentHub)
# library(ggalluvial)

Libraries <- c('standR', 'SpatialExperiment', 'limma', 'ExperimentHub',
               'ggalluvial','optparse')
lapply(Libraries, require, character.only = TRUE)
rm(Libraries)

################################################################################
##########################  Read-in arguments for normalisation run ############


# option_list = list(
#   make_option(c("-d", "--datadir"), type="character", default=NULL,
#               help="root directory for data files", metavar="character"),
#   make_option(c("-c", "--countFile"), type="character", default=NULL,
#               help="count data file name", metavar="character"),
#   make_option(c("-s", "--sampleAnnoFile"), type="character", default=NULL,
#               help="sample annotation data file name", metavar="character"),
#   make_option(c("-b", "--biolFactors"), type="character", default=NULL,
#               help="sample annotation data file name", metavar="character"),
#   make_option(c("-f", "--featureAnnoFile"), type="character", default=NULL,
#               help="feature annotation data file name", metavar="character"),
# );
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# # opt = parse_args(opt_parser);


################################################################################
####################  Read-in probe counts, anno and feature files (TSVs). #####

# dataDir = opt$datadir
# setwd(dataDir)
# 
# countFile = opt$countFile
# sampleAnnoFile = opt$sampleAnnoFile
# featureAnnoFile = opt$featureAnnoFile

countFile <- "/Users/upton6/Documents/Nanostring/projects/NS_msWTA/DSP_Output/countFile_SKC.tsv"
sampleAnnoFile <- "/Users/upton6/Documents/Nanostring/projects/NS_msWTA/DSP_Output/sampleAnnoFile_SKC.tsv"
featureAnnoFile <- "/Users/upton6/Documents/Nanostring/projects/NS_msWTA/DSP_Output/featureAnnoFile.tsv"


################################################################################
#####################################  Run through STANDR ######################


# Set up spatial experiment object
spe <- readGeoMx(countFile, sampleAnnoFile, 
                 featureAnnoFile = featureAnnoFile, rmNegProbe = TRUE)

setwd("/Users/upton6/Documents/Nanostring/projects/NS_msWTA/StandR/SKC")
# setwd(dataDir)

# ```
# ## QC
# ### metadata visualization 
# Based on the description of the data, we know that all glomerulus are classified 
# as abnormal and healthy, and tubule are classified as neg and PanCK. 
# We therefore merge the region-related annotations to avoid collinearity, which 
# can affect the process of batch correction.
# ```{r}

# Merge annotations to avoid colinearity (can adversely affect batch correction)

# colData(spe)$regions <- paste0(colData(spe)$region,"_",colData(spe)$SegmentLabel) |> 
#   (\(.) gsub("_Geometric Segment","",.))() |>
#   paste0("_",colData(spe)$pathology) |>
#   (\(.) gsub("_NA","_ns",.))()

colData(spe)$regions <- paste0(colData(spe)$Strain,"_",colData(spe)$Sex,"_",colData(spe)$KO,"_",colData(spe)$Diet)


# todo: get variables to plot from options

# plotSampleInfo(spe, column2plot = c("SlideName","disease_status","regions"))
pdfName <- 'SampleInfo_spe.pdf'
pdf( pdfName, width = 8 , height = 4 )
plotSampleInfo(spe, column2plot = c("SlideName","Strain","Sex","KO","Diet"))
dev.off()


### Gene level QC
spe <- addPerROIQC(spe, rm_genes = TRUE)
pdfName <- 'GeneQC_spe.pdf'
pdf( pdfName, width = 8 , height = 4 )
plotGeneQC(spe, ordannots = "regions", col = regions, point_size = 2)
dev.off()

# 
# Using the `plotGeneQC` function, we can have a look at which were the genes 
# removed and the overall distribution of percentage of non-expressed genes in 
# all ROIs. By default, top 9 genes are plotted here (arranging by mean 
# expression), user can increase the number of plotted genes by changing 
# the `top_n` parameter.
# 
# In this case we don't see any specific biological pattern for the samples 
# expressing this genes (figure above).
# 
# 
# ### ROI level QC
# 
# In the ROI level QC, we first aim to identify (if any) ROI(s) that have 
# relatively low library size and low cell count because they are considered as 
# low quality samples due to insufficient sequencing depth or lack of RNA in the 
# chosen region. 
# 
# In this case, looking at the distribution plots of library size and nuclei 
# count, we don't see any particular spike in the low ends, rather the 
# distributions are relatively smooth. Looking at the dot plot, library sizes are 
# mostly positively correlate with the nuclei count, with some data have 
# relatively low library size while the nuclei count is reasonable. We therefore 
# can try to draw an filtering threshold at the low end of the library size, in 
# this case 50000. By coloring the dot with their slide names, we find that the 
# ROIs below the threshold are all from slide disease1B, suggesting the reason 
# for this might be some technical issues of slide disease1B.
# 


pdfName <- 'ROIQC_spe.pdf'
pdf( pdfName, width = 8 , height = 4 )
plotROIQC(spe, y_threshold = 50000, col = SlideName)
dev.off()

# 
# Since library size of 50000 seems to be a reasonable threshold, here we subset 
# the spatial experiment object based on the library size in `colData`.
# 

spe <- spe[,rownames(colData(spe))[colData(spe)$lib_size > 50000]]

## Inspection of variations on ROI level
### RLE
# Here we can see obvious variation from slides to slides, and small variations 
# are also observed within each slide.
# 

pdfName <- 'RLE_spe.pdf'
pdf( pdfName, width = 8 , height = 4 )
plotRLExpr(spe, ordannots = "SlideName", assay = 2, col = SlideName)
dev.off()

### PCA
# 
# Here we color the PCA with slide information, and shape by regions (tissue). 
# We can see that PC1 is mainly spread out by regions, especially glomerulus and 
# tubule. And grouping based on slide within each tissue are observed. 
# The subtypes in tubule are clearly separated, but different subtypes of 
# glomerulus is still grouping together. Moreover, diseased tissues and control 
# tissues are mixed as well (disease slides and normal slides).
# 
pdfName <- 'PCAA_spe.pdf'
pdf( pdfName, width = 8 , height = 4 )
drawPCA(spe, assay = 2, col = SlideName, shape = regions)
dev.off()


# Data normalization
# 
# As we observed the technical variations in the data in both RLE and PCA plots. 
# It is necessary to perform normalization in the data.
# 
# In the `standR` package, we offer normalization options including TMM, RPKM, 
# TPM, CPM, upperquartile and sizefactor. Among them, RPKM and TPM required gene 
# length information (add `genelength` column to the `rowData` of the object). 
# For TMM, upperquartile and sizefactor, their normalized factor will be stored 
# their `metadata`.
# 
# Here we used TMM to normalize the data.
# colData(spe)$biology <- paste0(colData(spe)$disease_status, "_", colData(spe)$regions)

# x = "Strain_Diet"
# thisVector = strsplit(x, split = "_")[[1]]
# thisVector
# thisVector[1]
# str(thisVector[1])
# 
# colData(spe)$Strain
# colData(spe)$"Strain"
# 
# thisVector[1]
# 
# spe[colData(spe)[,"Strain"] == "SKC"]
# 
# colData(spe)[,thisVector[1]]
# colData(spe)[,thisVector[2]]
# 
# paste0(colData(spe)[,thisVector[1]],"_",colData(spe)[,thisVector[2]])
# 
# 
# paste(thisVector, collapse = "_")
# 
# vector1 <- c("Data Science", "is", "fun")
# vector1
# paste(vector1, collapse = "-")

colData(spe)$biology <- paste0(colData(spe)$Strain, "_", colData(spe)$Diet, "_", colData(spe)$KO, "_", colData(spe)$Sex)

colData(spe)$biology <- paste0(colData(spe)$Strain)
colData(spe)$biology <- paste0(colData(spe)$Diet)
colData(spe)$biology <- paste0(colData(spe)$KO)
colData(spe)$biology <- paste0(colData(spe)$Sex)

spe_tmm <- geomxNorm(spe, method = "TMM")

# 
# # Batch correction
# 
# In the Nanostring's GeoMX DSP protocol, due to the fact that one slide is only 
# big enough for a handful of tissue segments (ROIs), it is common that we see 
# the DSP data being confounded by the batch effect introduced by different 
# slides. In order to establish fair comparison between ROIs later on, it is 
# necessary to remove this batch effect from the data.
# 
# To run RUV4 batch correction, we need to provide a list of "negative control 
# genes (NCGs)".
# 
# The function `findNCGs` allows identifying the NCGs from the data. In this case, 
# since the batch effect is mostly introduced by slide, we therefore want to 
# identify NCGs across all slides, so here we set the `batch_name` to "SlideName", 
# and select the top 500 least variable genes across different slides as NCGs. 
# 

# spe <- findNCGs(spe, batch_name = "SlideName", top_n = 500)
spe <- findNCGs(spe, top_n = 500)

metadata(spe) |> names()

# Here we use k of 5 to perform RUV-4 normalization.

spe_ruv <- geomxBatchCorrection(spe, factors = "biology", 
                   NCGs = metadata(spe)$NCGs, k = 5)

# We can then inspect the PCA of the corrected data with annotations, to inspect 
# the removal of batch effects, and the retaining of the biological factors.
# 

# points <- c(0,1,2,3,4,7,8,9)
#                 # 10,11,15,16,17,18,19,20,21))

pdfName <- 'PCA_spe_ruv.pdf'
# ToDo: Increase number of shapes in plot. change colours to colour by KO or diet
pdf( pdfName, width = 8 , height = 4 )
# plotPairPCA(spe_ruv, assay = 2, color = Strain, shape = regions, title = "RUV4")
plotPairPCA(spe_ruv, assay = 2, color = Diet, shape = regions, title = "RUV4")
plotPairPCA(spe_ruv, assay = 2, color = KO, shape = regions, title = "RUV4")
plotPairPCA(spe_ruv, assay = 2, color = Sex, shape = regions, title = "RUV4")
# # plotPairPCA(spe_ruv, assay = 2, color = KO, shape = regions, title = "RUV4", pch=points)
dev.off()


# 
# Moreover, we can also have a look at the RLE plots of the normalized count.
# 
pdfName <- 'RLE_spe_ruv.pdf'
pdf( pdfName, width = 8 , height = 4 )
plotRLExpr(spe_ruv, assay = 2, color = SlideName) + ggtitle("RUV4")
# plotRLExpr(spe_ruv, assay = 2, color = Strain) + ggtitle("RUV4")
plotRLExpr(spe_ruv, assay = 2, color = Diet) + ggtitle("RUV4")
plotRLExpr(spe_ruv, assay = 2, color = KO) + ggtitle("RUV4")
plotRLExpr(spe_ruv, assay = 2, color = Sex) + ggtitle("RUV4")
dev.off()
# 
# 
# **For more detailed analysis pipeline and usage of the standR package, please see https://github.com/DavisLaboratory/GeoMXAnalysisWorkflow**
# 

# Write RUV corrected results to file
spe_RUV_out <- assay(spe_ruv, i=2)
spe_RUV_out <- 2**spe_RUV_out

write.table(spe_RUV_out, file = "RUV_corrected_values_SKC.csv", sep = ",", row.names = TRUE, col.names = TRUE)



# # SessionInfo
sessionInfo()
































