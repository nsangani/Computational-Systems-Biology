# Install Packages
##install.packages("HelpersMG")
##install.packages("WGCNA")

# The GO.db package will not let me load WGCNA, so let's download it
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GO.db", version = "3.8")

# Load the packages
library(HelpersMG)
library(WGCNA)
library(tidyverse)
library(phyloseq)

#Load data
data = read.csv("C:\\Users\\Sangani\\Desktop\\Course Work\\Spring 2022\\Comp. Methods in System Biology\\Assignments\\Assignment 3\\HW03_expression.csv") # select gene expression data
trait = read.csv("C:\\Users\\Sangani\\Desktop\\Course Work\\Spring 2022\\Comp. Methods in System Biology\\Assignments\\Assignment 3\\HW03_Traits.csv") # select trait data
data = t(data)

colnames(data) <- data[1,]
data <- data[-1, ]

rownames(trait) <- trait[,1]
trait <- trait[,-1]

str(data)
mode(data) <- "numeric"

library(dplyr)
data1 = inner_join(data, trait)
dim(datOTUClr)
dim(datContinuous)

datContinuous = datContinuous[rownames(datOTUClr),]

#Remove null values
datOTUClr = na.omit(data)


head(datContinuous)
dim(datContinuous)

# Pre Analysis ------------------------------------------------------------


# Check for OTUs and samples with too many missing values
good_otus <- goodSamplesGenes(datOTUClr, verbose = 3);
good_otus$allOK

# Re-cluster samples
sampleTree2 <- hclust(dist(datOTUClr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(datContinuous, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datContinuous),
                    main = "Sample dendrogram and trait heatmap")



# WGCNA processes ---------------------------------------------------------

powers <- c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft <- pickSoftThreshold(datOTUClr, powerVector = powers, verbose = 5)


# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))

# Set some parameters
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", 
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")




adjacency <- adjacency(datOTUClr, power = 3)
TOMadj <- TOMsimilarity(adjacency)


dissTOMadj <- 1- TOMadj


# Clustering using TOM
# Call the hierarchical clustering function 
hclustGeneTree <- hclust(as.dist(dissTOMadj), method = "average")

# Plot the resulting clustering tree (dendogram)
sizeGrWindow(12, 9)
plot(hclustGeneTree, xlab = "", sub = "", 
     main = "Gene Clustering on TOM-based disssimilarity", 
     labels = FALSE, hang = 0.04)

minModuleSize <- 12

# Module ID using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = hclustGeneTree, 
                             distM = dissTOMadj,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)


table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

sizeGrWindow(8,6)
plotDendroAndColors(hclustGeneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")


dynamic_MEList <- moduleEigengenes(datOTUClr, colors = dynamicColors)
dynamic_MEs <- dynamic_MEList$eigengenes

# Calculate dissimilarity of module eigengenes
dynamic_MEDiss <- 1-cor(dynamic_MEs)
dynamic_METree <- hclust(as.dist(dynamic_MEDiss))
# Plot the hclust
sizeGrWindow(7,6)
plot(dynamic_METree, main = "Dynamic Clustering of module eigengenes",
     xlab = "", sub = "")

######################## MERGE SIMILAR MODULES
dynamic_MEDissThres <- 0.25

# Plot the cut line
#abline(h = dynamic_MEDissThres, col = "red")

# Call an automatic merging function
merge_dynamic_MEDs <- mergeCloseModules(datOTUClr, dynamicColors, cutHeight = dynamic_MEDissThres, verbose = 3)
# The Merged Colors
dynamic_mergedColors <- merge_dynamic_MEDs$colors

# Eigen genes of the new merged modules
mergedMEs <- merge_dynamic_MEDs$newMEs
mergedMEs
table(dynamic_mergedColors)


sizeGrWindow(12,9)
plotDendroAndColors(hclustGeneTree, cbind(dynamicColors, dynamic_mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Rename Module Colors 
moduleColors <- dynamic_mergedColors

# Construct numerical labels corresponding to the colors 
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs

save(MEs, moduleLabels, moduleColors, hclustGeneTree, file = "data/WGCNA/dynamic-prodOTU-02-networkConstruction.RData")




# Trait- module -----------------------------------------------------------

# define numbers of genes and samples
nOTUs <- ncol(datOTUClr)
nSamples <- nrow(datOTUClr)

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datOTUClr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

names(MEs) <- substring(names(MEs), 3)


moduleTraitCor <- cor(MEs, datContinuous, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# PLOT
sizeGrWindow(10,6)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values within a heatmap
labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = names(datContinuous),
               yLabels = names(MEs), 
               ySymbols = names(MEs), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, 
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait Relationships"))




# OTU in relation with modules and trait(tot_bacprod) ---------------------

bacprod <- as.data.frame(datContinuous)
names(bacprod) <- "bacprod" # rename

# Names (colors of the modules
#modNames <- substring(names(MEs), 3) # Remove the first two letters in string

# Calculate the correlations between modules
geneModuleMembership <- as.data.frame(WGCNA::cor(datOTUClr, MEs, use = "p"))

# What are the p-values for each correlation?
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# What's the correlation for the trait: bacterial production?
geneTraitSignificance <- as.data.frame(cor(datOTUClr, bacprod, use = "p"))

# What are the p-values for each correlation?
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples = nSamples))

names(geneTraitSignificance) <- paste("GS.", names(bacprod), sep = "")
names(GSPvalue) <- paste("p.GS.", names(bacprod), sep = "")
# Visualization -----------------------------------------------------------

dissTOM <- 1-TOMsimilarityFromExpr(datOTUClr, power = 3) 
plotTOM <- dissTOM^7

# Set diagnol to NA for a nicer plot
diag(plotTOM) <- NA

######################## ALL GENES - even though it looks like subsampling 
nSelect = 482
nGenes = ncol(datOTUClr)
nSamples = nrow(datOTUClr)

# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")





# Isolate weight from the clinical traits
bacprod <- as.data.frame(datContinuous)
names(bacprod) = "bacprod"

# Add the weight to existing module eigengenes
MET <- orderMEs(cbind(MEs, bacprod))
#names(MET) <- substring(names(MET), 3)
sizeGrWindow(6,6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendograms = TRUE, xLabeles = 90)
