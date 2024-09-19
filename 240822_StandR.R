################################################################################
###############################      Author : Kyle R. Upton      ###############
################################################################################

### R version 4.3.2 (2023-10-31)

## Attributions: Adapted from StandR Demo scripts

################################################################################
################################### Required R packages ########################

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

# countFile <- "/Users/upton6/Documents/Nanostring/projects/NS_msWTA/DSP_Output/countFile_SKC.tsv"
# sampleAnnoFile <- "/Users/upton6/Documents/Nanostring/projects/NS_msWTA/DSP_Output/sampleAnnoFile_SKC.tsv"
# featureAnnoFile <- "/Users/upton6/Documents/Nanostring/projects/NS_msWTA/DSP_Output/featureAnnoFile.tsv"

countFile <- "/Users/upton6/Documents/Nanostring/Projects/NS_JCU_WTA/DSP_Output/countFile_Epithelial.tsv"
sampleAnnoFile <- "/Users/upton6/Documents/Nanostring/Projects/NS_JCU_WTA/DSP_Output/sampleAnnoFile_Epithelial.tsv"
featureAnnoFile <- "/Users/upton6/Documents/Nanostring/Projects/NS_JCU_WTA/DSP_Output/featureAnnoFile.tsv"

countFile <- "/Users/upton6/Documents/Nanostring/Projects/NS_JCU_WTA/DSP_Output/countFile_Luminal.tsv"
sampleAnnoFile <- "/Users/upton6/Documents/Nanostring/Projects/NS_JCU_WTA/DSP_Output/sampleAnnoFile_Luminal.tsv"
featureAnnoFile <- "/Users/upton6/Documents/Nanostring/Projects/NS_JCU_WTA/DSP_Output/featureAnnoFile.tsv"

################################################################################
#####################################  Run through STANDR ######################


# Set up spatial experiment object

# Negative probe was removed manually in python so rmNegProbe has to be FALSE
# spe <- readGeoMx(countFile, sampleAnnoFile, 
#                  featureAnnoFile = featureAnnoFile, rmNegProbe = TRUE)
spe <- readGeoMx(countFile, sampleAnnoFile, 
                 featureAnnoFile = featureAnnoFile, rmNegProbe = FALSE)

# setwd("/Users/upton6/Documents/Nanostring/projects/NS_msWTA/StandR/SKC")
setwd("/Users/upton6/Documents/Nanostring/projects/NS_JCU_WTA/DSP_Output/StandR/Epithelial")
setwd("/Users/upton6/Documents/Nanostring/projects/NS_JCU_WTA/DSP_Output/StandR/Luminal")
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

# colData(spe)$regions <- paste0(colData(spe)$Strain,"_",colData(spe)$Sex,"_",colData(spe)$KO,"_",colData(spe)$Diet)


# colData(spe)$regions <- paste0(colData(spe)$Core,"_",colData(spe)$Clinical_outcome,"_",colData(spe)$Time_point,"_",colData(spe)$Group)
colData(spe)$regions <- paste0(colData(spe)$Clinical_outcome,"_",colData(spe)$Time_point,"_",colData(spe)$Group)


# todo: get variables to plot from options

# plotSampleInfo(spe, column2plot = c("SlideName","disease_status","regions"))
pdfName <- 'SampleInfo_spe.pdf'
pdf( pdfName, width = 8 , height = 4 )
plotSampleInfo(spe, column2plot = c("ScanLabel","Clinical_outcome","Group","Time_point","Core"))
dev.off()


### Gene level QC
spe <- addPerROIQC(spe, rm_genes = TRUE)
pdfName <- 'GeneQC_spe.pdf'
pdf( pdfName, width = 8 , height = 4 )
plotGeneQC(spe, ordannots = "regions", col = regions, point_size = 2, bin_num = 30, top_n=90)
dev.off()


pdfName <- 'ROIQC_spe.pdf'
pdf( pdfName, width = 8 , height = 4 )
# plotROIQC(spe, y_threshold = 1000000, col = SlideName)
plotROIQC(spe, y_threshold = 100000, col = colData(spe)$ScanLabel)
dev.off()


# Remove samples with a library size lower than a chosen threshold value
# spe <- spe[,rownames(colData(spe))[colData(spe)$lib_size > 50000]]
spe <- spe[,rownames(colData(spe))[colData(spe)$lib_size > 100000]]

## Inspection of variations on ROI level
### RLE
# Here we can see obvious variation from slides to slides, and small variations 
# are also observed within each slide.
# 

pdfName <- 'RLE_spe.pdf'
pdf( pdfName, width = 8 , height = 4 )
plotRLExpr(spe, ordannots = "SegmentLabel", assay = 2, col = colData(spe)$SegmentLabel)
# plotRLExpr(spe, ordannots = "Core", assay = 2, col = colData(spe)$Core)
plotRLExpr(spe, ordannots = "Clinical_outcome", assay = 2, col = colData(spe)$Clinical_outcome)
plotRLExpr(spe, ordannots = "Time_point", assay = 2, col = colData(spe)$Time_point)
plotRLExpr(spe, ordannots = "Group", assay = 2, col = colData(spe)$Group)
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

# Set regions for PCA plots with fewer discrete values for shapes
colData(spe)$smlRegions <- paste0(colData(spe)$Clinical_outcome,"_",colData(spe)$Group)

pdfName <- 'PCAA_spe.pdf'
pdf( pdfName, width = 20 , height = 20 )
# drawPCA(spe, assay = 2, col = Core, shape = smlRegions)
drawPCA(spe, assay = 2, col = Clinical_outcome, shape = smlRegions)
drawPCA(spe, assay = 2, col = Time_point, shape = smlRegions)
drawPCA(spe, assay = 2, col = Group, shape = smlRegions)
dev.off()


# Set regions for actual analysis
colData(spe)$regions <- paste0(colData(spe)$Clinical_outcome,"_",colData(spe)$Time_point,"_",colData(spe)$Group)


set.seed(100)
spe <- scater::runPCA(spe)
pca_results <- reducedDim(spe, "PCA")
pdfName <- 'PCAA_Scree_spe.pdf'
pdf( pdfName, width = 8 , height = 4 )
plotScreePCA(spe, precomputed = pca_results)
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

# colData(spe)$biology <- paste0(colData(spe)$Strain, "_", colData(spe)$Diet, "_", colData(spe)$KO, "_", colData(spe)$Sex)
# colData(spe)$biology <- paste0(colData(spe)$Diet, "_", colData(spe)$KO)

# Changing the order of the factors does not appear to have any effect on normalisation
# colData(spe)$biology <- paste0(colData(spe)$KO, "_", colData(spe)$Diet)
# colData(spe)$biology <- paste0(colData(spe)$KO, "_", colData(spe)$Diet, "_", colData(spe)$Sex)

# colData(spe)$biology <- paste0(colData(spe)$SegmentLabel, "_", colData(spe)$Core, "_", colData(spe)$Clinical_outcome, "_", colData(spe)$Time_point, "_", colData(spe)$Group)
colData(spe)$biology <- paste0(colData(spe)$Clinical_outcome, "_", colData(spe)$Time_point, "_", colData(spe)$Group)


# colData(spe)$biology <- paste0(colData(spe)$Strain)
# colData(spe)$biology <- paste0(colData(spe)$Diet)
# colData(spe)$biology <- paste0(colData(spe)$KO)
# colData(spe)$biology <- paste0(colData(spe)$Sex)

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
spe <- findNCGs(spe, batch_name = "ScanLabel", top_n = 500)

metadata(spe) |> names()

# Here we use k of 5 to perform RUV-4 normalization.

pdfName <- 'PCA_series_Core-Clinical_outcome-Time_point-Group_spe.pdf'
pdf( pdfName, width = 8 , height = 4 )

for(i in seq(20)){
  spe_ruv <- geomxBatchCorrection(spe, factors = "biology", 
                                  NCGs = metadata(spe)$NCGs, k = i)
  
  # print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Core, title = paste0("k = ", i)))
  # print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Clinical_outcome, title = paste0("k = ", i)))
  # print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Time_point, title = paste0("k = ", i)))
  # print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i)))
  plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = SegmentLabel, title = paste0("k = ", i))
  # plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Core, title = paste0("k = ", i))
  plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Clinical_outcome, title = paste0("k = ", i))
  plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Time_point, title = paste0("k = ", i))
  plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Group, title = paste0("k = ", i))
  
}
dev.off()



spe_ruv2 <- geomxBatchCorrection(spe, factors = "biology", 
                                NCGs = metadata(spe)$NCGs, k = 2)
spe_ruv3 <- geomxBatchCorrection(spe, factors = "biology", 
                                 NCGs = metadata(spe)$NCGs, k = 3)
spe_ruv4 <- geomxBatchCorrection(spe, factors = "biology", 
                                 NCGs = metadata(spe)$NCGs, k = 4)
spe_ruv5 <- geomxBatchCorrection(spe, factors = "biology", 
                                NCGs = metadata(spe)$NCGs, k = 5)
spe_ruv6 <- geomxBatchCorrection(spe, factors = "biology", 
                                 NCGs = metadata(spe)$NCGs, k = 6)
spe_ruv7 <- geomxBatchCorrection(spe, factors = "biology", 
                                 NCGs = metadata(spe)$NCGs, k = 7)
spe_ruv8 <- geomxBatchCorrection(spe, factors = "biology", 
                                 NCGs = metadata(spe)$NCGs, k = 8)
spe_ruv9 <- geomxBatchCorrection(spe, factors = "biology", 
                                 NCGs = metadata(spe)$NCGs, k = 9)
spe_ruv10 <- geomxBatchCorrection(spe, factors = "biology", 
                                NCGs = metadata(spe)$NCGs, k = 10)
spe_ruv15 <- geomxBatchCorrection(spe, factors = "biology", 
                                NCGs = metadata(spe)$NCGs, k = 15)
spe_ruv20 <- geomxBatchCorrection(spe, factors = "biology", 
                                NCGs = metadata(spe)$NCGs, k = 20)


pdfName <- 'plotClusterEvalStats_spe_ruv.pdf'
# ToDo: Increase number of shapes in plot. change colours to colour by KO or diet
pdf( pdfName, width = 20 , height = 10 )

spe_list <- list(spe, spe_ruv2,spe_ruv3,spe_ruv4,spe_ruv5,spe_ruv6,spe_ruv7,spe_ruv8,spe_ruv9,spe_ruv10,spe_ruv15,spe_ruv20)
plotClusterEvalStats(spe_list = spe_list,
                     bio_feature_name = "biology",
                     batch_feature_name = "ScanLabel",
                     data_names = c("Raw","RUV4-2","RUV4-3","RUV4-4","RUV4-5","RUV4-6","RUV4-7","RUV4-8","RUV4-9","RUV410","RUV415","RUV420"))


spe_list <- list(spe, spe_ruv2,spe_ruv3,spe_ruv4,spe_ruv5,spe_ruv6,spe_ruv7,spe_ruv8,spe_ruv9,spe_ruv10,spe_ruv15,spe_ruv20)
plotClusterEvalStats(spe_list = spe_list,
                     bio_feature_name = "Clinical_outcome",
                     batch_feature_name = "ScanLabel",
                     data_names = c("Raw","RUV4-2","RUV4-3","RUV4-4","RUV4-5","RUV4-6","RUV4-7","RUV4-8","RUV4-9","RUV410","RUV415","RUV420"))

# spe_list <- list(spe, spe_ruv2,spe_ruv3,spe_ruv4,spe_ruv5,spe_ruv6,spe_ruv7,spe_ruv8,spe_ruv9,spe_ruv10,spe_ruv15)
# plotClusterEvalStats(spe_list = spe_list,
#                      bio_feature_name = "Core",
#                      batch_feature_name = "ScanLabel",
#                      data_names = c("Raw","RUV4-2","RUV4-3","RUV4-4","RUV4-5","RUV4-6","RUV4-7","RUV4-8","RUV4-9","RUV410","RUV415","RUV420"))

spe_list <- list(spe, spe_ruv2,spe_ruv3,spe_ruv4,spe_ruv5,spe_ruv6,spe_ruv7,spe_ruv8,spe_ruv9,spe_ruv10,spe_ruv15,spe_ruv20)
plotClusterEvalStats(spe_list = spe_list,
                     bio_feature_name = "Time_point",
                     batch_feature_name = "ScanLabel",
                     data_names = c("Raw","RUV4-2","RUV4-3","RUV4-4","RUV4-5","RUV4-6","RUV4-7","RUV4-8","RUV4-9","RUV410","RUV415","RUV420"))

spe_list <- list(spe, spe_ruv2,spe_ruv3,spe_ruv4,spe_ruv5,spe_ruv6,spe_ruv7,spe_ruv8,spe_ruv9,spe_ruv10,spe_ruv15,spe_ruv20)
plotClusterEvalStats(spe_list = spe_list,
                     bio_feature_name = "Group",
                     batch_feature_name = "ScanLabel",
                     data_names = c("Raw","RUV4-2","RUV4-3","RUV4-4","RUV4-5","RUV4-6","RUV4-7","RUV4-8","RUV4-9","RUV410","RUV415","RUV420"))

dev.off()



### Cut down list of comparisons by removing worst performers for each
# 
# spe_list <- list(spe, spe_ruv2,spe_ruv3,spe_ruv4,spe_ruv5,spe_ruv6,spe_ruv7,spe_ruv8,spe_ruv9,spe_ruv10,spe_ruv15)
# plotClusterEvalStats(spe_list = spe_list,
#                      bio_feature_name = "Sex",
#                      batch_feature_name = "biology",
#                      data_names = c("Raw","RUV4-2","RUV4-3","RUV4-4","RUV4-5","RUV4-6","RUV4-7","RUV4-8","RUV4-9","RUV410","RUV415"))
# 
# 
# spe_list <- list(spe, spe_ruv2,spe_ruv3,spe_ruv4,spe_ruv5,spe_ruv6,spe_ruv7,spe_ruv8,spe_ruv9,spe_ruv10,spe_ruv15)
# plotClusterEvalStats(spe_list = spe_list,
#                      bio_feature_name = "Diet",
#                      batch_feature_name = "biology",
#                      data_names = c("Raw","RUV4-2","RUV4-3","RUV4-4","RUV4-5","RUV4-6","RUV4-7","RUV4-8","RUV4-9","RUV410","RUV415"))
# 
# spe_list <- list(spe, spe_ruv2,spe_ruv3,spe_ruv4,spe_ruv5,spe_ruv6,spe_ruv7,spe_ruv8,spe_ruv9,spe_ruv10,spe_ruv15)
# plotClusterEvalStats(spe_list = spe_list,
#                      bio_feature_name = "KO",
#                      batch_feature_name = "biology",
#                      data_names = c("Raw","RUV4-2","RUV4-3","RUV4-4","RUV4-5","RUV4-6","RUV4-7","RUV4-8","RUV4-9","RUV410","RUV415"))
# 


spe_ruv <- geomxBatchCorrection(spe, factors = "biology", 
                                NCGs = metadata(spe)$NCGs, k = 6)


# We can then inspect the PCA of the corrected data with annotations, to inspect 
# the removal of batch effects, and the retaining of the biological factors.
# 

# points <- c(0,1,2,3,4,7,8,9)
#                 # 10,11,15,16,17,18,19,20,21))



# Core
# Clinical_outcome
# Time_point
# Group




pdfName <- 'PCA_spe_ruv.pdf'
# ToDo: Increase number of shapes in plot. change colours to colour by KO or diet
pdf( pdfName, width = 20 , height = 20 )
# plotPairPCA(spe_ruv, assay = 2, color = Strain, shape = regions, title = "RUV4")
plotPairPCA(spe_ruv, assay = 2, color = SegmentLabel, shape = smlRegions, title = "RUV4")
# plotPairPCA(spe_ruv, assay = 2, color = Core, shape = smlRegions, title = "RUV4")
plotPairPCA(spe_ruv, assay = 2, color = Clinical_outcome, shape = smlRegions, title = "RUV4")
plotPairPCA(spe_ruv, assay = 2, color = Time_point, shape = smlRegions, title = "RUV4")
plotPairPCA(spe_ruv, assay = 2, color = Group, shape = smlRegions, title = "RUV4")
dev.off()


# 
# Moreover, we can also have a look at the RLE plots of the normalized count.
# 
pdfName <- 'RLE_spe_ruv.pdf'
pdf( pdfName, width = 8 , height = 4 )
plotRLExpr(spe_ruv, assay = 2, color = SegmentLabel) + ggtitle("RUV4")
# plotRLExpr(spe_ruv, assay = 2, color = Core) + ggtitle("RUV4")
plotRLExpr(spe_ruv, assay = 2, color = Clinical_outcome) + ggtitle("RUV4")
plotRLExpr(spe_ruv, assay = 2, color = Time_point) + ggtitle("RUV4")
plotRLExpr(spe_ruv, assay = 2, color = Group) + ggtitle("RUV4")
dev.off()
# 
# 
# **For more detailed analysis pipeline and usage of the standR package, please see https://github.com/DavisLaboratory/GeoMXAnalysisWorkflow**
# 

# Write RUV corrected results to file
spe_RUV_out <- assay(spe_ruv, i=2)
spe_RUV_out <- 2**spe_RUV_out

write.table(spe_RUV_out, file = "RUV_corrected_values.csv", sep = ",", row.names = TRUE, col.names = TRUE)



# # SessionInfo
sessionInfo()
