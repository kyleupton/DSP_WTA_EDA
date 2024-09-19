### To Do.

# Read section 3 (Specific Experimental Designs) in User guide
# Run through with basic group-group comparisons (using exact test)

# Build up a descent glm approach



library(edgeR)


segment = 'Luminal'
# segment = 'Epithelial'

if (segment == 'Luminal'){
  setwd("/Users/upton6/Documents/Nanostring/projects/NS_JCU_WTA/DSP_Output/StandR/Luminal")
}else if (segment == 'Epithelial'){
  setwd("/Users/upton6/Documents/Nanostring/projects/NS_JCU_WTA/DSP_Output/StandR/Epithelial")
}

# counts <- read.table(file = "RUV_corrected_values.csv", header = TRUE, sep=",")

counts <- read.csv(file = "RUV_corrected_values.csv", header = TRUE, sep=",")
# counts

min(counts)
# counts = counts - min(counts) ## ToDo: Update this to perform exp transform, then log2+1

head(counts)
dim(counts)


setwd("/Users/upton6/Documents/Nanostring/projects/NS_JCU_WTA/DSP_Output")

if (segment == 'Luminal'){
  targets <- read.delim("sampleAnnoFile_Luminal.tsv", sep="\t", header=TRUE)
}else if (segment == 'Epithelial'){
  targets <- read.delim("sampleAnnoFile_Epithelial.tsv", sep="\t", header=TRUE)
}


rownames(targets) <- targets$SegmentDisplayName

head(targets)
dim(targets)
targets





if (segment == 'Luminal'){
  # setwd("/Users/upton6/Documents/Nanostring/projects/NS_JCU_WTA/DSP_Output/StandR/Luminal")
  setwd("/Users/upton6/Documents/Nanostring/projects/NS_JCU_WTA/DSP_Output/StandR/Coeliac_Gluten_Chalenge/Luminal")
  # setwd("/Users/upton6/Documents/Nanostring/projects/NS_JCU_WTA/DSP_Output/StandR/Gut_Rwesponse_HW_Infection/Luminal")
} else if (segment == 'Epithelial'){
  # setwd("/Users/upton6/Documents/Nanostring/projects/NS_JCU_WTA/DSP_Output/StandR/Epithelial")
  setwd("/Users/upton6/Documents/Nanostring/projects/NS_JCU_WTA/DSP_Output/StandR/Coeliac_Gluten_Chalenge/Epithelial")
  # setwd("/Users/upton6/Documents/Nanostring/projects/NS_JCU_WTA/DSP_Output/StandR/Gut_Rwesponse_HW_Infection/Epithelial")
}


# group <- factor(paste0("group", ".", targets$Clinical_outcome, ".", targets$Time_point, ".", targets$Group))
group <- factor(paste0("group", ".", targets$Time_point, ".", targets$Group))


group

y <- DGEList(counts, group = group)
y$samples

#### Plot first sample
plotMD(cpm(y, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

#### Plot Multi Dimensional Scaling
#### Points and Colours need work
#set up to use colours for grade and symbols for label
# points <- c(0,1,2,15,16,17)
# points <- c(0,1,2,3,4,7,8,9,10,11,15,16,17,18,19,20,21)

colors <- c("blue","darkgreen","pink", "darkgreen", "orange","yellow", 
            "yellow", "yellow", "yellow", "yellow", "darkgreen", 
            "yellow", "yellow", "yellow", "yellow", "red")
# points <- rep(c(0,1,2,3,4,5,6),3)
points <- c(0,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6)

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

# 
# # group.Responder.Post challenge.Hookworm
# con <- makeContrasts(group.Responder.PostChallenge.Hookworm - group.NonResponder.PostChallenge.Hookworm, levels=design)
# pdf( "Responder_vs_NonResponder.pdf" , 
#      width = 4 , height = 4 ) # in inches
# filename <- "MD_plot_Responder_vs_NonResponder.csv"
# 
# 
#


# # group.Responder.PostChallenge.Placebo
# con <- makeContrasts(group.Responder.PostChallenge.Placebo - group.Responder.Baseline.Placebo, levels=design)
# pdf( "Responder_PostChallenge_vs_Responder_Baseline.pdf" ,
#      width = 4 , height = 4 ) # in inches
# filename <- "Responder_PostChallenge_vs_Responder_Baseline.csv"


# 
# # group.Responder.PostHookworm.Hookworm
# con <- makeContrasts(group.Responder.PostHookworm.Hookworm - group.Responder.Baseline.Hookworm, levels=design)
# pdf( "Responder_PostHookworm_vs_Responder_Baseline.pdf" , 
#      width = 4 , height = 4 ) # in inches
# filename <- "Responder_PostHookworm_vs_Responder_Baseline.csv"







# group.PostChallenge.Placebo
con <- makeContrasts(group.PostChallenge.Placebo - group.Baseline.Placebo, levels=design)
pdf( "PostChallenge_vs_Baseline.pdf" ,
     width = 4 , height = 4 ) # in inches
filename <- "PostChallenge_vs_Baseline.csv"



# # group.PostHookworm.Hookworm
# con <- makeContrasts(group.PostHookworm.Hookworm - group.Baseline.Hookworm, levels=design)
# pdf( "PostHookworm_vs_Baseline.pdf" , 
#      width = 4 , height = 4 ) # in inches
# filename <- "PostHookworm_vs_Baseline.csv"



qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf, 15)
summary(decideTests(qlf))

## Save plots and write results to file?
# pdf( "MD_plot_Responder_vs_NonResponder.pdf" , width = 4 , height = 4 ) # in inches
# pdf( "PostChallenge_vs_Baseline.pdf" , width = 4 , height = 4 ) # in inches
# pdf( "PostHookworm_vs_Baseline.pdf" , width = 4 , height = 4 ) # in inches
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

