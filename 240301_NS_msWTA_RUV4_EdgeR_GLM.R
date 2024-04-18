### To Do.

# Read section 3 (Specific Experimental Designs) in User guide
# Run through with basic group-group comparisons (using exact test)

# Build up a descent glm approach



library(edgeR)


setwd("/Users/upton6/Documents/Nanostring/projects/NS_msWTA/StandR/SKC")

# counts <- read.table(file = "RUV_corrected_values.csv", header = TRUE, sep=",")

counts <- read.csv(file = "RUV_corrected_values_SKC.csv", header = TRUE, sep=",")
# counts

min(counts)
# counts = counts - min(counts) ## ToDo: Update this to perform exp transform, then log2+1

head(counts)
dim(counts)


setwd("/Users/upton6/Documents/Nanostring/projects/NS_msWTA/DSP_Output/")
targets <- read.delim("sampleAnnoFile.tsv", sep="\t", header=TRUE)
# targets$SegmentDisplayName

rownames(targets) <- targets$SegmentDisplayName

head(targets)
dim(targets)
targets

### Need to subset targets and counts to PFAC samples only
targets <- targets[targets["Strain"] == "SKC" , ]

head(targets)
dim(targets)


## Extract row names to subset counts
keep <- as.list(rownames(targets))
keep <- gsub("\\|", ".", keep)
keep <- gsub(" ", ".", keep)

keep
# list(keep)
keepCounts <- counts[, keep, drop = FALSE]
head(keepCounts)
dim(keepCounts)


setwd("/Users/upton6/Documents/Nanostring/projects/NS_msWTA/StandR/SKC/EdgeR/")



# group <- factor(paste0(targets$Grade_b, ".", targets$Grade_c))
# group <- factor(paste0("group", ".", targets$Grade_b))


# ToDo: Add Sex as a factor for groups!!
group <- factor(paste0("group", ".", targets$KO, ".", targets$Diet))
group <- factor(paste0("group", ".", targets$KO, ".", targets$Diet, ".", targets$Sex))
group

y <- DGEList(keepCounts, group = group)
y$samples

# y$norm.factors

#### Plot first sample
plotMD(cpm(y, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

#### Plot Multi Dimensional Scaling
#### Points and Colours need work
#set up to use colours for grade and symbols for label
points <- c(0,1,2,15,16,17)
points <- c(0,1,2,3,4,7,8,9,10,11,15,16,17,18,19,20,21)

colors <- c("blue", "blue", "blue", "blue", "blue", 
            "darkgreen", "darkgreen", "darkgreen", "darkgreen", "darkgreen", 
            "red", "red", "red", "red", "red")
points <- rep(c(0,1,2,3,4,5,6),3)

plotMDS(y, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), col=colors[group], pch=points, ncol=3)



#### Is this the best design set-up? Double check the use of 0 vs other options

design <- model.matrix(~0 + group)

# LipidVacoule	Inclusion



colnames(design) <- levels(group)
design

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit)


group

### Set up a list of different contrasts to make here

# # Tumour - TME
# con <- makeContrasts(Tumour - TME, levels=design)
# con <- makeContrasts(Tumour - Full_ROI, levels=design)
# con <- makeContrasts(TME - Full_ROI, levels=design)

# help(make.names)


# Comparisons
con <- makeContrasts(group.KO.TAA - group.WT.TAA, levels=design)
pdf( "MD_plot_KO_TAA_vs_WT_TAA.pdf" , 
     width = 4 , height = 4 ) # in inches
filename <- "MD_plot_KO_TAA_vs_WT_TAA.csv"

con <- makeContrasts(group.KO.TAA - group.KO.Standard, levels=design)
pdf( "MD_plot_KO_TAA_vs_KO_Standard.pdf" , 
     width = 4 , height = 4 ) # in inches
filename <- "MD_plot_KO_TAA_vs_KO_Standard.csv"

con <- makeContrasts(group.KO.Standard - group.WT.Standard, levels=design)
pdf( "MD_plot_KO_Standard_vs_WT_Standard.pdf" , 
     width = 4 , height = 4 ) # in inches
filename <- "MD_plot_KO_Standard_vs_WT_Standard.csv"

con <- makeContrasts(group.WT.TAA - group.WT.Standard, levels=design)
pdf( "MD_plot_WT_TAA_vs_WT_Standard.pdf" , 
     width = 4 , height = 4 ) # in inches
filename <- "MD_plot_WT_TAA_vs_WT_Standard.csv"




# Comparisons
con <- makeContrasts(group.KO.TAA.Male - group.WT.TAA.Male, levels=design)
pdf( "MD_plot_KO_TAA_Male_vs_WT_TAA_Male.pdf" , 
     width = 4 , height = 4 ) # in inches
filename <- "MD_plot_KO_TAA_Male_vs_WT_TAA_Male.csv"

con <- makeContrasts(group.KO.TAA.Female - group.WT.TAA.Female, levels=design)
pdf( "MD_plot_KO_TAA_Female_vs_WT_TAA_Female.pdf" , 
     width = 4 , height = 4 ) # in inches
filename <- "MD_plot_KO_TAA_Female_vs_WT_TAA_Female.csv"


con <- makeContrasts(group.KO.TAA.Male - group.KO.Standard.Male, levels=design)
pdf( "MD_plot_KO_TAA_Male_vs_KO_Standard_Male.pdf" , 
     width = 4 , height = 4 ) # in inches
filename <- "MD_plot_KO_TAA_Male_vs_KO_Standard_Male.csv"

con <- makeContrasts(group.KO.TAA.Female - group.KO.Standard.Female, levels=design)
pdf( "MD_plot_KO_TAA_Female_vs_KO_Standard_Female.pdf" , 
     width = 4 , height = 4 ) # in inches
filename <- "MD_plot_KO_TAA_Female_vs_KO_Standard_Female.csv"


con <- makeContrasts(group.KO.Standard.Male - group.WT.Standard.Male, levels=design)
pdf( "MD_plot_KO_Standard_Male_vs_WT_Standard_Male.pdf" , 
     width = 4 , height = 4 ) # in inches
filename <- "MD_plot_KO_Standard_Male_vs_WT_Standard_Male.csv"

con <- makeContrasts(group.KO.Standard.Female - group.WT.Standard.Female, levels=design)
pdf( "MD_plot_KO_Standard_Female_vs_WT_Standard_Female.pdf" , 
     width = 4 , height = 4 ) # in inches
filename <- "MD_plot_KO_Standard_Female_vs_WT_Standard_Female.csv"


con <- makeContrasts(group.WT.TAA.Male - group.WT.Standard.Male, levels=design)
pdf( "MD_plot_WT_TAA_Male_vs_WT_Standard_Male.pdf" , 
     width = 4 , height = 4 ) # in inches
filename <- "MD_plot_WT_TAA_Male_vs_WT_Standard_Male.csv"

con <- makeContrasts(group.WT.TAA.Female - group.WT.Standard.Female, levels=design)
pdf( "MD_plot_WT_TAA_Female_vs_WT_Standard_Female.pdf" , 
     width = 4 , height = 4 ) # in inches
filename <- "MD_plot_WT_TAA_Female_vs_WT_Standard_Female.csv"








qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf, 15)
summary(decideTests(qlf))

## Save plots and write results to file?
# pdf( "MD_plot_Tumour_1_Tumour_0.pdf" , width = 4 , height = 4 ) # in inches
par( mfrow=c(1 ,1) )
plotMD(qlf)
dev.off() # this tells [R] to close and stop writing to the pdf.

resultsbyP <- topTags(qlf, n = nrow(qlf$table))$table
wh.rows.glm <- match( rownames( resultsbyP ) , rownames( y$counts ) )
results2.tbl <- cbind (resultsbyP, 
                       "Tgw.Disp"=y$tagwise.dispersion[wh.rows.glm], 
                       "UpDown" = decideTestsDGE(qlf)[wh.rows.glm,], 
                       y$counts[wh.rows.glm,] )
head (results2.tbl)
write.table(results2.tbl, file = filename, sep = ",", row.names = TRUE)


