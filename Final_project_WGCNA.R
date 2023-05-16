# network analysis of liver expression data in female mice

getwd()
setwd('C:/Users/Sangani/Desktop/Course Work/Spring 2022/Comp. Methods in System Biology/Final Project/WGCNA_example/')

library(WGCNA)
options(stringsAsFactors = FALSE) #important step but don't know why???

#Allow multi-threading
enableWGCNAThreads()
 
lnames = load(file = "FemaleLiver-01-dataInput.RData") # How to load both of these files in .RData?
# datExpr (sample x Expr) - 134x3600 & datTraits (sample x Traits) - 134x26


# Step 1*: Network Construction and Module Detection

# Step 2: Choosing soft-thresholding power - network topology
powers = c(c(1:10), seq(from=12, to=20, by=2)) # 1:10; 12:20 by=2

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) #power 6 has highest truncated.R.sq(0.962)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") # assigned cutoff 0.90

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

########################### Method 1: Step-by-step network ######################################

# Step 3: Co-expression similarity and adjacency(nearness)

softPower = 6; # Note identified using beta soft-thresholding power
adjacency = adjacency(datExpr, power = softPower); # file list of gene with connectivity of 6

# Step 4: Topological Overlap Matrix (TOM) # to minimize effects of noise and spurious association

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


# Step 5: Clustering using TOM # hierarchical clustering or dendrogram of genes

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram) # vertical line = gene; interconnected = highly co-expressed genes
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;

# Module identification using dynamic tree cut: #Note: Module = Cluster
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Results: Label 0 = unassigned genes, 1:22 modules are largest to smallest
# Plot these 22 modules using gene dendrogram

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors) #uniq color assigned to each modules

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


# Step 6: Merging of modules based on similar expression profiles

#To quantify co-expression similarity of entire modules, we calculate their eigengenes and cluster them on their correlation
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average"); #23 clusters

# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25 # Corresponds to 0.75 correlation

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors;

# Eigengenes of the new merged modules: #merged dendrogram
mergedMEs = merge$newMEs;

# Replot to view the updated dendrogram
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors 

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")

########################### Method 2: Blockwise network construction ######################################

# Step 2: Block-wise network construction and module detection

# Suppose total genes are 5000 and your computer has limited memory
# A 16GB workstation should handle up to BlockSize of 20000 probes
# BlockSize default 5000
bwnet = blockwiseModules(datExpr, maxBlockSize = 2000,
                         power = 6, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "femaleMouseTOM-blockwise",
                         verbose = 3)


# Load the results of single-block analysis # From Method 1
load(file = "FemaleLiver-02-networkConstruction-auto.RData");

# Relabel blockwise modules
bwLabels = matchLabels(bwnet$colors, moduleLabels);
table(bwLabels) #all the modules and their sizes

# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)

# open a graphics window # Note dendrogram for two blocks
sizeGrWindow(6,6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#  would like to change some of the tree cut; use function recutBlockwiseTrees

# Step 3: Comparing single block and block-wise network analysis

sizeGrWindow(12,9)
plotDendroAndColors(geneTree,
                    cbind(moduleColors, bwModuleColors),
                    c("Single block", "2 blocks"),
                    main = "Single block gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Visual inspection to confirm agreement between the single-block and the block-wise module assignment.
# Next, verify module eigengenes of modules that correspond to one another in the single-block and block-wise approaches are extremely similar
singleBlockMEs = moduleEigengenes(datExpr, moduleColors)$eigengenes;
blockwiseMEs = moduleEigengenes(datExpr, bwModuleColors)$eigengenes;

# calculate the correlations of the corresponding eigengenes
single2blockwise = match(names(singleBlockMEs), names(blockwiseMEs))
signif(diag(cor(blockwiseMEs[, single2blockwise], singleBlockMEs)), 3)

####################### Relating Modules to External Clinical Traits #################################

# Step 1: Quantifying module–trait associations

#correlate eigengenes with external traits and look for most signigicant association

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Step 2: Gene relationship to trait and important modules: Gene Significance (GS) and Module Membership (MM)

# quantify the similarity of all genes on the array to every module by defining Gene Significance GS
# GS = correlation between the gene and the trait
# MM = correlation of the module eigengene and the gene expression profile
# Define variable weight containing the weight column of datTrait

weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

# Step 3: Intramodular analysis: identifying genes with high GS and MM

# brown module has highest association with weight
# Scatterplot of GS vs MM in the brown module

module = "brown" #try other highly associated gene vs trait
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# Step 4: Summary output of network analysis results
# After finding modules with high association with our trait of interest, merge statistical information with gene annotation

names(datExpr) #910
names(datExpr)[moduleColors=="brown"] #409

# connect ProbeIDs to gene_names and Entrez code ...different in our case
annot = read.csv(file = "GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo.csv")









