# Weighted Gene Co-expression Network Analysis (WGCNA)

The following scripts were made following the WGCNA tutorials.

## 1. Filtering low count transcripts, normalisation and batch effect removal

```r
set.seed(0)
library(DESeq2) 
library(vroom)
library(tximport)
library(RColorBrewer)
library(ggplot2)
library(doParallel)
library(tidyverse)
library(ggdendro)
library(viridis)
library(limma)
library(edgeR)
library(ggrepel)
library(devtools)
library(sva)
library(WGCNA)
library(SummarizedExperiment)
library(BioNERO)

################################################################################
################################################################################
#############         1.a . Loading expression data        #####################
################################################################################
################################################################################

PHYLUM = "Annelida"
SPECIES="Arenicola_marina"
SPECIES_CODE = "COCO"

#Paths to give as input
species_path = paste0("/", PHYLUM, "/", SPECIES)
WGCNA_path = paste(species_path,"/WGCNA/", sep="") #Path for saving WGCNA output files
filtering_path= paste(WGCNA_path,"/1.Filtering/", sep="") #Path to save the preprocessing plots (PCA and heatmaps)

QuantsPath=paste0(species_path, "/quants", sep="") #Path to take the Salmon quants files

directories <- list.dirs(QuantsPath, full.names = TRUE, recursive = FALSE)
# Find quant.sf files in each directory
quant_sf_paths <- unlist(
  Filter(Negate(is.null), 
         lapply(directories, function(dir) {
           # Construct the path to quant.sf
           path <- file.path(dir, "quant.sf")
           # Check if the file exists
           if (file.exists(path)) {
             return(path)  # Return the path if file exists
           } else {
             return(NULL)  # Return NULL if file doesn't exist
           }
         }))
)
# Verify quant_sf_paths
print(quant_sf_paths)
# Extract the sample names and remove "_quant"
sampleNames <- sub("_quant$", "", 
                   basename(unlist(Filter(Negate(is.null), 
                                          lapply(directories, function(dir) {
                                            if (file.exists(file.path(dir, "quant.sf"))) dir else NULL
                                          })))))
sampleNames
names(quant_sf_paths) <- sampleNames
# Extract conditions (everything after COCO_head_ or COCO_tail_ and before _quant, without the replicate number)
conditions <- sub("COCO_(head|tail)_([A-Za-z0-9]+?)[0-9]$", "\\2", sampleNames)
conditions

# Extract the "head" or "tail" part
body_part <- sub("COCO_(head|tail)_.*", "\\1", sampleNames)
body_part

# Create metadata data frame
metadata_COCO <- data.frame(
  sampleName = sampleNames,
  condition = factor(conditions),
  body_part = body_part,
  cond_bp = paste0(factor(conditions),"_", body_part)
)
metadata_COCO

#Create one matrix from the RNA-seq quant data 
txi_COCO <- tximport(quant_sf_paths, type = "salmon", txOut=TRUE, countsFromAbundance = "no")
head(txi_COCO$abundance)

################################################################################
################################################################################
##### 1.b. Filter low counts and Variance Stabilizing Transformation   #########
################################################################################
################################################################################

ddsTxi <- DESeqDataSetFromTximport(txi_COCO,
                                   colData = metadata_COCO,
                                   design = ~condition)

#Keep only genes that have at least 10 counts in 90% of the samples
percentage <- dim(txi_COCO$counts)[2]*0.9
keep <- rowSums(counts(ddsTxi) >= 10) >= percentage
ddsTxi_keep <- ddsTxi[keep,]
dim(ddsTxi_keep)

# rlog normalisation
rld <- rlog(ddsTxi_keep,blind=FALSE)
# Plot and save PCA and heatmap of this step 
# Before running the following lines you may need to change 
# plot_PCA  and plot_heatmap funtions for better visualisatiion
salmonSE_rld <- SummarizedExperiment(assays=list(counts=assay(rld)), colData=metadata_COCO)
pdf(paste(filtering_path,SPECIES_CODE, "_rld_DESEQ2.pdf",sep=""),         # File name
    width = 20, height = 12, # Width and height in inches
    bg = "white")#,          # Background color
    # Color model (cmyk is required for most publications)
    #paper = "A4r")          # Paper size
par(mar = c(5, 2, 2, 3))
BioNERO::plot_PCA(salmonSE_rld)
par(mar = c(5, 5, 5, 5))
BioNERO::plot_heatmap(salmonSE_rld)
dev.off()

################################################################################
################################################################################
###########             1.c. limma Remove Batch Effect          ###############
################################################################################
################################################################################

#batch <- metadata_COCO$sampleID
group <- metadata_COCO$condition

df_group <- as.data.frame(group)
rownames(df_group) <- metadata_COCO$sampleName
colnames(df_group) <- "Condition"

df_group

#create the full model matrix 
mod <- model.matrix(~as.factor(Condition), data = df_group)

#Create the null model matrix
mod0 <- model.matrix(~1,data=df_group)
n.sv = num.sv(assay(rld),mod,method="be")
n.sv

sva_res <- svaseq(assay(rld),mod,mod0,n.sv=n.sv)
names(sva_res)

adjusted_counts_rld_limma <- limma::removeBatchEffect(assay(rld), covariates = sva_res$sv, group=metadata_COCO$condition)

# Plot and save PCA and heatmap of this step 
salmonSE_rld_limma <- SummarizedExperiment(assays=list(counts=adjusted_counts_rld_limma), colData=metadata_COCO)
pdf(paste(filtering_path,SPECIES_CODE, "_rld_DESEQ2_limma.pdf",sep=""),         # File name
        width = 20, height = 12, # Width and height in inches
        bg = "white")#,          # Background color
    # Color model (cmyk is required for most publications)
    #paper = "A4r")          # Paper size
par(mar = c(5, 2, 2, 3))
BioNERO::plot_PCA(salmonSE_rld_limma)
par(mar = c(5, 5, 5, 5))
BioNERO::plot_heatmap(salmonSE_rld_limma)
dev.off()

################################################################################
################################################################################
#########      1.d. Checking data for excessive missing values     #############
#########            and identification of outlier samples         #############
################################################################################
################################################################################

# WGCNA works with matrices where rows are the samples and column the genes
adjusted_counts_rld_limma_trans <- t(adjusted_counts_rld_limma)

#We first check for genes and samples with too many missing values:
gsg = goodSamplesGenes(adjusted_counts_rld_limma_trans, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(colnames(adjusted_counts_rld_limma_trans)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(adjusted_counts_rld_limma_trans)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  adjusted_counts_rld_limma_trans_GG = adjusted_counts_rld_limma_trans[gsg$goodSamples, gsg$goodGenes]
} else if (gsg$allOK) 
{
    # If all genes and samples are good, retain the original matrix:
    adjusted_counts_rld_limma_trans_GG <- adjusted_counts_rld_limma_trans
}

# Check if there are genes removed
dim(adjusted_counts_rld_limma_trans_GG)
blocksize <- dim(adjusted_counts_rld_limma_trans_GG)[2]

wgcna_tree_cluster <- function(matrix_trans){
  sampleTree = hclust(dist(matrix_trans), method = "average");
  # Plot the sample tree: 
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
}

salmonSE_rld_limma_GG <- SummarizedExperiment(assays=list(counts=t(adjusted_counts_rld_limma_trans_GG)), colData=metadata_COCO)
pdf(paste(filtering_path,SPECIES_CODE, "_rld_DESEQ2_limma_GG.pdf",sep=""),         # File name
    width = 20, height = 12, # Width and height in inches
    bg = "white")#,          # Background color
# Color model (cmyk is required for most publications)
#paper = "A4r")          # Paper size
par(mar = c(5, 2, 2, 3))
BioNERO::plot_PCA(salmonSE_rld_limma_GG)
par(mar = c(5, 5, 5, 5))
BioNERO::plot_heatmap(salmonSE_rld_limma_GG)
par(mar = c(5, 5, 5, 5))
wgcna_tree_cluster(adjusted_counts_rld_limma_trans_GG)
dev.off()

# The salmonSE_rld_limma_GG is the final matrix that will be used 
# to the network construction.

save(txi_COCO, metadata_COCO, adjusted_counts_rld_limma_trans_GG, salmonSE_rld_limma_GG, blocksize, 
     file = paste0(filtering_path, "1.Filtering.RData"))

```

## 2. Network construcion

```r
set.seed(0)
library(WGCNA)
library(SummarizedExperiment)

#Paths to give as input
PHYLUM = "Annelida"
SPECIES="Arenicola_marina"
SPECIES_CODE = "COCO"

species_path = paste0("/", PHYLUM, "/", SPECIES)
WGCNA_path = paste(species_path,"/WGCNA", sep="") #Path for saving WGCNA output files
filtering_path= paste(WGCNA_path,"/1.Filtering/", sep="") #Path to save the preprocessing plots (PCA and heatmaps)
network_construction_path <-  paste(WGCNA_path,"/2.Network_construction/", sep="")

#load(paste0(filtering_path, "1.Filtering.RData"))
load("../1.Filtering.RData")
temp_cor <- cor 
cor <- WGCNA::cor

load(paste(network_construction_path,SPECIES_CODE,"_sft_and_net.RData", sep=""))

# Choose a set of soft-thresholding powers
powers = seq(4,20,by=1);
# Call the network topology analysis function
sft = pickSoftThreshold(adjusted_counts_rld_limma_trans_GG, powerVector = powers, verbose = 5, networkType = "signed hybrid", RsquaredCut = 0.8, 
                        blockSize = blocksize, corFnc = bicor, corOptions = list(use = 'pairwise.complete.obs', maxPOutliers=0.1, robustY=FALSE))

# Plot the results:
pdf(paste(network_construction_path,SPECIES_CODE, "_SFT_plot.pdf",sep=""),         # File name
    width = 10, height = 8, # Width and height in inches
    bg = "white")          # Background color
    # Color model (cmyk is required for most publications)
    #paper = "A4r")          # Paper size
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence",SPECIES));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity", SPECIES))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

k <- softConnectivity(adjusted_counts_rld_limma_trans_GG,
                      corFnc = "bicor", corOptions = list(use = 'pairwise.complete.obs', maxPOutliers = 0.1, robustY=FALSE), 
                      weights = NULL,
                      type = "signed hybrid",
                      power = sft$powerEstimate,
                      blockSize = blocksize, 
                      minNSamples = NULL, 
                      verbose = 5, indent = 0)

# And plot them for visual inspection
# Plot a histogram of k and a scale free topology plot
pdf(paste(network_construction_path,SPECIES_CODE, "_Soft_connectivity2.pdf",sep=""),         # File name
    width = 10, height = 8, # Width and height in inches
    bg = "white")         # Background color
    # Color model (cmyk is required for most publications)
    #paper = "A4r")    
sizeGrWindow(12,7)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main= paste("Check scale free topology for" , SPECIES, "\n"))
dev.off()

# Build the network
net = blockwiseModules(adjusted_counts_rld_limma_trans_GG,
                       # Data checking options
                       checkMissingData=TRUE, 
                       # Options for splitting data into blocks
                       maxBlockSize = blocksize,
                       # Network construction arguments: correlation options
                       corType = "bicor", corOptions = list(use = 'pairwise.complete.obs', maxPOutliers = 0.1, robustY=FALSE),
                       # Adjacency function options
                       power = sft$powerEstimate, networkType = "signed hybrid", replaceMissingAdjacencies = TRUE, 
                       # Adjacency function options
                       TOMType = "signed", 
                       # Saving or returning TOM
                       saveTOMs = TRUE, saveTOMFileBase = paste(SPECIES_CODE,"_blockwiseTOM",sep=""),
                       # Basic tree cut options
                       deepSplit = 4, minModuleSize = 30,
                       # Gene reassignment, module trimming, and module "significance" criteria
                       reassignThreshold =1e-6, mergeCutHeight = 0.15,
                       # Output options
                       numericLabels = TRUE, pamRespectsDendro = TRUE,
                       nThreads = 30, verbose = 5)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
unmergedColors = labels2colors(net$unmergedColors)

# Plot the dendrogram and the module colors underneath
pdf(paste(network_construction_path,SPECIES_CODE, "_Dendrogram.pdf",sep=""),         # File name
    width = 10, height = 8, # Width and height in inches
    bg = "white")          # Background color
    # Color model (cmyk is required for most publications)
    #paper = "A4r")    
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", main = paste("Cluster Dendrogram of", SPECIES),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(net$dendrograms[[1]], cbind(unmergedColors, mergedColors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)
dev.off()

save(sft, net , file = paste0(network_construction_path, SPECIES_CODE,"_sft_and_net.RData"))
```

## 3. Relate modules to stress experiments

```r
set.seed(0)
library(WGCNA)
library(SummarizedExperiment)
library(dplyr)

#Paths to give as input
PHYLUM = "Annelida"
SPECIES = "Arenicola_marina"
SPECIES_CODE = "COCO"

species_path = paste0("/", PHYLUM, "/", SPECIES)
WGCNA_path = paste0(species_path,"/WGCNA/") #Path for saving WGCNA output files
filtering_path= paste0(WGCNA_path,"/1.Filtering/") #Path to save the preprocessing plots (PCA and heatmaps)
network_construction_path <-  paste0(WGCNA_path,"/2.Network_construction/")
mod2traits_path <- paste0(WGCNA_path,"/3.Relate_modules_to_external_traits/")

# Load the expression and trait data 

load(paste(filtering_path,"1.Filtering.RData", sep=""))
load(paste(network_construction_path, SPECIES_CODE,"_sft_and_net.RData", sep=""))

################################################################################
################################################################################
########           Relating modules to external traits         #################
################################################################################
################################################################################

################################################################################
########      Quantifying module–trait associations      #######################
################################################################################
# In this analysis we would like to identify modules that are significantly 
# associated with each experiment by correlating the module eigengene 
# with the experiments. 

moduleLabels = net$colors #Genes with the number of module they belong 

#Convert the number of each module to a color
moduleColors = labels2colors(net$colors) # module color per gene, same positions as the moduleLabels class: character
MEs = net$MEs; # A df where columns are the modules (as numbers) and rows are the replicates
geneTree = net$dendrograms[[1]];

# Define numbers of genes and samples
nGenes = ncol(adjusted_counts_rld_limma_trans_GG)
nSamples = nrow(adjusted_counts_rld_limma_trans_GG);
datTraits <- binarizeCategoricalVariable(metadata_COCO$condition, includePairwise = FALSE,
                                         includeLevelVsAll = TRUE,
                                         minCount = 1)

datTraits <- cbind(datTraits, metadata_COCO$condition)
datTraits
datTraits <- datTraits[,-13] #Change this accordingly
datTraits
datTraits <- as.data.frame(datTraits)
colnames(datTraits) <- sub("\\.vs.*", "", colnames(datTraits)) #Remove the .vs.all from colnames
datTraits

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(adjusted_counts_rld_limma_trans_GG, moduleColors)$eigengenes # Calculate module eigengenes, same as MEs before, but now we have colors instead of numbers
#Reorder given (eigen-)vectors such that similar ones (as measured by correlation) are next to each other.
MEs = orderMEs(MEs0)
moduleTraitCor = WGCNA::cor(MEs, datTraits, use = "p", method = "pearson"); #Correlation of each module with the Trait 
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
dim(textMatrix) 
dim(moduleTraitCor)

# Display the correlation values within a heatmap plot
pdf(paste(mod2traits_path, SPECIES_CODE, "_Module-traitCorrelations_Heatmap.pdf",sep=""),         # File name
    width = 15, height = 29, # Width and height in inches
    bg = "white" )        # Background color
par(mar = c(6, 10,3,3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships of Arenicola marina"))
dev.off()

################################################################################
########      Gene relationship to trait and important modules:    #############
########        Gene Significance and Module Membership            #############
################################################################################

# We quantify associations of individual genes with our trait of interest (eg. 
# Hypoxia) by defining Gene Significance GS as (the absolute value of) the 
# correlation between the gene and the trait. For each module, we also define a 
# quantitative measure of module membership MM as the correlation of the module 
# eigengene and the gene expression profile. This allows us to quantify the 
# similarity of all genes on the array to every module.

# names (in colors) of the modules without the ME prefix
modNames = substring(names(MEs), 3)
# Correlation of expression with the modules
geneModuleMembership = as.data.frame(cor(adjusted_counts_rld_limma_trans_GG, MEs, use = "p"));
# Calculate Student asymptotic p-value for the above correlation
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
# Assign names 
names(geneModuleMembership) = paste("MM", modNames, sep=""); #Change instead of MEcyan to MMcyan 
names(MMPvalue) = paste("p.MM", modNames, sep=""); #Change instead of MEcyan to p.MMcyan 

# We need to generate the following dataframes for stress experiment

# But, since we are going to have many traits,
# this function returns a list with these 3 dataframes per trait:
# 1. {trait}_df
# 2. geneTraitSignificance_{trait}
# 3. GSPvalue_{trait}

calculate_all_traits_associations <- function(datExpr, datTraits, nSamples) {
  #datExpr = adjusted_counts_rld_limma_trans_GG

  # Initialize lists to store the results
  results <- list()
  
  for (trait_name in colnames(datTraits)) {
    # Extract the current trait
    trait <- as.numeric(datTraits[[trait_name]])
    trait_df <- as.data.frame(trait)
    
    # Set the name of the trait column to the current trait name
    names(trait_df) <- trait_name
  
    
    # Calculate gene trait significance
    geneTraitSignificance <- as.data.frame(cor(datExpr, trait_df, use = "p"))
    GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
    
    # Rename columns
    names(geneTraitSignificance) <- paste("GS.", names(trait_df), sep = "")
    names(GSPvalue) <- paste("p.GS.", names(trait_df), sep = "")
    
    # Dynamically name the trait_df and store the results in the list
    trait_df_name <- paste(trait_name, "_df", sep = "")
    results[[trait_df_name]] <- trait_df
    results[[paste("geneTraitSignificance_", trait_name, sep = "")]] <- geneTraitSignificance
    results[[paste("GSPvalue_", trait_name, sep = "")]] <- GSPvalue
  }
  
  return(results)
}
  
# Dataframe containing the 3 variables mentioned before per trait 
all_traits_associations <- calculate_all_traits_associations(adjusted_counts_rld_limma_trans_GG, datTraits,nSamples)

################################################################################
# Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high 
# significance for each experiment as well as high module membership 
# in interesting modules. 
################################################################################

# Plot a scatterplot per experiment of Gene Significance vs. Module Membership
# for the modules that are significant for each trait. 
# The following function (get_significant_modules_df) will return only the 
# statistically significant modules for each trait (p-value<5e-2).

#############      Get significant modules for each experiment      #############
# This function will return a df where columns are the traits and rows are the significant 
# modules per trait with p-value <5e-2 and also ordered based on p-value 
# (more to the less significant)
find_significant_modules_per_trait_ordered <- function(MEs, all_traits_associations, 
                                                       moduleTraitCor, moduleTraitPvalue, 
                                                       modNames, pval = 0.05) 
{  
  # Initialize lists to store significant modules
  significant_modules_with_corr <- list()
  significant_modules_without_corr <- list()
  
  # Loop through each trait
  for (trait in colnames(moduleTraitCor)) {
    trait_df <- paste0(trait, "_df")
    modOrder <- order(-abs(cor(MEs, all_traits_associations[[trait_df]], use = "p")))
    
    significant_mods_with_corr <- c()
    significant_mods_without_corr <- c()
    
    for (mod in modOrder) {
      sig_mod <- names(MEs)[mod] 
      if (moduleTraitPvalue[sig_mod, trait] < pval && sig_mod != "MEgrey") {
        mod_with_corr <- paste0(modNames[mod], " (", round(moduleTraitCor[sig_mod, trait], 2), ")")
        significant_mods_with_corr <- c(significant_mods_with_corr, mod_with_corr)
        
        mod_without_corr <- modNames[mod]
        significant_mods_without_corr <- c(significant_mods_without_corr, mod_without_corr)
      }
    }
    
    significant_modules_with_corr[[trait]] <- significant_mods_with_corr
    significant_modules_without_corr[[trait]] <- significant_mods_without_corr
  }
  
  # Determine max length for formatting
  max_length <- max(sapply(significant_modules_with_corr, length))
  
  # Create data frames
  combined_df_with_corr <- as.data.frame(matrix(NA, nrow = max_length, ncol = length(significant_modules_with_corr)))
  combined_df_without_corr <- as.data.frame(matrix(NA, nrow = max_length, ncol = length(significant_modules_without_corr)))
  
  colnames(combined_df_with_corr) <- names(significant_modules_with_corr)
  colnames(combined_df_without_corr) <- names(significant_modules_without_corr)
  
  # Fill data frames
  for (trait in names(significant_modules_with_corr)) {
    combined_df_with_corr[1:length(significant_modules_with_corr[[trait]]), trait] <- significant_modules_with_corr[[trait]]
    combined_df_without_corr[1:length(significant_modules_without_corr[[trait]]), trait] <- significant_modules_without_corr[[trait]]
  }
  
  return(list(with_corr = combined_df_with_corr, without_corr = combined_df_without_corr))
}

signif_modules <- find_significant_modules_per_trait_ordered(MEs, all_traits_associations, moduleTraitCor, moduleTraitPvalue, modNames)
signif_modules$with_corr
signif_modules$without_corr
dim(signif_modules$with_corr)
significant_modules <- signif_modules$without_corr
  

# Define the output file name
sigMods2file <-paste0(mod2traits_path, "Significant_modules_per_trait.txt")
# Write the common rows (all columns) to a new output file
write.table(significant_modules, sigMods2file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

################################################################################
######           Summary output of network analysis results               ######        
################################################################################

# We have found modules with high association with our trait of interest, and 
# we have identified their central players by the Module Membership measure. 
# We will now merge this statistical information with gene annotation and write out 
# a file that summarizes the most important results and can be inspected in standard spreadsheet software.

#Path to save gene lists of significant modules
dir.create(paste(mod2traits_path, "/2.signifMods2genes/", sep=""))
sigMods2genes <- paste(mod2traits_path, "/2.signifMods2genes/", sep="")

extract_significant_module_genes <- function(significant_modules, datExpr, 
                                             moduleColors, sigMods2genes, SPECIES_CODE) {
  for (trait in colnames(significant_modules)) {
    for (module in significant_modules[[trait]]) {
      if (!is.na(module)) {
        # Extract the genes associated with the current module (will be used for enrichment analysis)
        significantModule_Genes <- as.data.frame(colnames(datExpr)[moduleColors==module])
        
        # Define the file path
        file_path <- paste(sigMods2genes, SPECIES_CODE, "_", module,"_significant_module_geneList_for_", trait,".txt",  sep = "")
        
        # Write the genes to the file
        write.table(significantModule_Genes, 
                    file = file_path, 
                    quote = FALSE, 
                    row.names = FALSE, 
                    col.names = FALSE)
      }
    }
  }
}

# Call the function to write the geneLists
extract_significant_module_genes(significant_modules, adjusted_counts_rld_limma_trans_GG, moduleColors, sigMods2genes, SPECIES_CODE)

################################################################################
######         WGCNA Tutorials Simulated-06-RelatingToExt         ############## 
###### 6 Relating modules and module eigengenes to external data  ##############
################################################################################

################################################################################
# 6a.Representing modules by eigengenes and relating eigengenes to one another #
################################################################################

# To get a sense of how related the modules are one can summarize each module 
# by its eigengene (first principal component).

datExpr <- adjusted_counts_rld_limma_trans_GG
mergedColors = labels2colors(net$colors)
colorh1 <- mergedColors

datME=moduleEigengenes(datExpr,colorh1)$eigengenes
signif(cor(datME, use="p"), 2)

# We define a dissimilarity measure between the module eigengenes that keeps 
# track of the sign of the correlation between the module eigengenes, 
# and use it to cluster the eigengene:

dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
pdf(paste(mod2traits_path, SPECIES_CODE, "_Eigengenes_Clustering_Tree.pdf",sep=""),         # File name
    width = 13, height = 5, # Width and height in inches
    bg = "white" )        # Background color
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")
dev.off()

################################################################################
###############    Displaying module heatmap and the eigengene    ##############
################################################################################
# Plot with gene expression heatmap and the barplot with the representative expression
# of the module eigenegene
dir.create(paste(mod2traits_path, "3.Module_expresssion_heatmap/", sep = ""))
mod2expr_path = paste(mod2traits_path, "3.Module_expresssion_heatmap/", sep = "")

sigMods_unique
plot_expression_profiles <- function(sigMods_unique, mod2expr_path,datME){
    
  for (which.module in sigMods_unique){

      ME=datME[, paste("ME",which.module, sep="")]
  
      pdf(paste(mod2expr_path, SPECIES_CODE, "_", which.module, "_module_expression_profile.pdf", sep=""),         # File name
           width = 7, height = 9, # Width and height in inches
           bg = "white")        # Background color
      #sizeGrWindow(8,7);
      par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
      plotMat(t(scale(datExpr[,colorh1==which.module ])),
            nrgcols=30,rlabels=F,rcols=which.module,
            main=which.module, cex.main=2)
      par(mar=c(7, 4.2, 0, 0.7))
      bp <- barplot(ME, col=which.module, main="", cex.main=2,axes=TRUE,
            ylab="eigengene expression") #,xlab="Samples")
      # Add rotated labels with a smaller font size and 45-degree angle
      text(x = bp, y = par("usr")[3] - 0.01,  # Position below bars
           labels = rownames(datME), srt = 90,  # Rotate text 45 degrees
          adj = 1, xpd = TRUE,          # Adjust text alignment and plotting outside the plot region
          cex = 0.4)   
      mtext("Samples", side = 1, line = 5)
      dev.off()
  }
}

# Call the function to plot the heatmap and the expression profiles of the 
# significant modules
plot_expression_profiles(sigMods_unique,mod2expr_path,datME)

################################################################################
################################################################################
########################       HUB GENES        ################################
################################################################################
################################################################################

################################################################################
# a. Intramodular connectivity
# We begin by calculating the intramodular connectivity for each gene. 
# (In network literature, connectivity is often referred to as ”degree”.) 
# The function intramodularConnectivity computes:
#             the whole network connectivity kTotal,
#             the within module connectivity kWithin, 
#             kOut=kTotal-kWithin, and 
#             kDiff=kIn-kOut=2*kIN-kTotal

bpower = sft$powerEstimate
ADJ1=abs(cor(datExpr,use="p"))^bpower
Alldegrees1=intramodularConnectivity(ADJ1, colorh1)
head(Alldegrees1)

################################################################################
# b. Relationship between gene significance and intramodular connectivity

dir.create(paste(mod2traits_path, "5.GS_VS_Connectivity/", sep=""))
geneSignificanceVSconnectivity_path = paste(mod2traits_path, "5.GS_VS_Connectivity/", sep="")

plot_GS_vs_Connectivity_scatterplots <- function(significant_modules,all_traits_associations, Alldegrees1,geneSignificanceVSconnectivity_path){
  # Set the graphics window size for a large plot
  sizeGrWindow(16, 12)  # Increase window size as needed
  
  # Set up the layout for a grid of plots (adjust number of rows and columns as needed)
  num_traits = length(colnames(significant_modules))
  par(mfrow = c(4, ceiling(num_traits / 4)))  # Adjust based on the number of traits
  par(mar = c(5, 6, 8, 4))   # Adjust margins
  
  # Iterate through each trait and save the plots in a separate PDF file
  for (trait in colnames(significant_modules)) {
    print(trait)
    
    # Dynamically get the GeneSignificance for the trait
    GeneSignificance1 = abs(eval(parse(text = paste0("all_traits_associations$geneTraitSignificance_", trait)))[, 1])
    
    # Filter out NA values from significant_modules
    non_na_modules = significant_modules[trait][!is.na(significant_modules[trait][, 1]) & significant_modules[trait][, 1] != "", 1]
    
    # Print the number of modules for the trait
    num_modules = length(non_na_modules)
    print(num_modules)
    
    # Open a new PDF device for each trait
    pdf(paste0(geneSignificanceVSconnectivity_path,"scatterplots_", trait, ".pdf"),
        width = 50, height = 15)  # Specify file name, width, and height
    
    
    # Set up the layout for a grid of plots (adjust rows/columns based on the number of non-NA modules)
    par(mfrow = c(4, ceiling(num_modules / 4)))  # Adjust based on the number of non-NA modules
    par(mar = c(5, 6, 8, 4))  # Adjust margins
    
    # Create scatter plots for each module of the trait
    for (i in seq_along(non_na_modules)) {
      whichmodule = non_na_modules[i]
      #print(whichmodule)  # Print the module name being processed
      
      restrict1 = (colorh1 == whichmodule)
      
      # Create scatter plot for each module
      verboseScatterplot(Alldegrees1$kWithin[restrict1],
                         GeneSignificance1[restrict1],
                         col = colorh1[restrict1],
                         main = whichmodule,
                         xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
    }
    
    # Close the PDF device to save the file
    dev.off()
  }
}

plot_GS_vs_Connectivity_scatterplots(significant_modules,all_traits_associations, Alldegrees1, geneSignificanceVSconnectivity_path)

################################################################################
# c. Generalizing intramodular connectivity for all genes on the array
# Extracting the hub genes

datKME=signedKME(datExpr, datME, outputColumnName="MM.")
# Display the first few rows of the data frame
head(datKME)

# d. Finding genes with high gene significance and high intramodular connectivity 
# in interesting modules (the statistically significant ones identified above)
dir.create(paste(mod2traits_path, "/6.HubGenes/", sep=""))
hubGenes_path = paste(mod2traits_path, "/6.HubGenes/", sep="")

HubGenes_extraction <- function(significant_modules, all_traits_associations, datKME, datExpr, output_dir, proteome_fasta) {
  
  # Read the FASTA file headers (proteome) into a vector
  proteome_headers <- readLines(proteome_fasta)
  proteome_headers <- proteome_headers[substr(proteome_headers, 1, 1) == ">"]   # Filter out only the lines that start with '>'
  protein_names <- gsub("^>", "", proteome_headers)  # Remove ">" from the protein names
  
  # Function to match hub genes to protein names
  match_hub_to_protein <- function(hub_genes, protein_names) {
    # Initialize a vector to store the matching protein names
    matching_proteins <- sapply(hub_genes, function(hub_gene) {
      # Search for the hub gene in protein names
      matched <- grep(paste0("^", hub_gene, "$"), gsub("\\.p[0-9]+$", "", protein_names))
      
      # If a match is found, return the full protein name with .p* suffix
      if (length(matched) > 0) {
        return(protein_names[matched])
      } else {
        return("")  # If no match is found, return empty string
      }
    })
    
    return(matching_proteins)
  }
  
  
  hub_counts <- significant_modules  # Copy structure of significant_modules
  for (trait in colnames(significant_modules)) {
    cat("Processing trait:", trait, "\n")
    
    # Filter out NA values from significant_modules
    non_na_modules = significant_modules[trait][!is.na(significant_modules[trait][, 1]) & significant_modules[trait][, 1] != "", 1]
    
    # Get gene significance values
    GS1 <- eval(parse(text = paste0("all_traits_associations$geneTraitSignificance_", trait)))[, 1]
    
    result_list <- list()
    hubs_per_module <- rep(0, length(non_na_modules))
    proteins_per_module <- rep(0, length(non_na_modules))  # Add for proteins
    for (i in seq_along(non_na_modules)) {
      whichmodule <- non_na_modules[i]
      
      # Identify genes that meet filtering criteria
      FilterGenes <- abs(GS1) > 0.2 & abs(datKME[[paste0("MM.", whichmodule)]]) > 0.8
      
      # Get gene names
      hub_genes <- dimnames(data.frame(datExpr))[[2]][FilterGenes]
      
      # Match hub genes to proteins using the match_hub_to_protein function
      protein_names_for_hub_genes <- match_hub_to_protein(hub_genes, protein_names)
      
      # Store results in list
      if (length(hub_genes) > 0) {
        result_list[[whichmodule]] <- data.frame(Gene = hub_genes, Module = whichmodule, Protein = protein_names_for_hub_genes)
        hubs_per_module[i] <- length(hub_genes)  # Count occurrences of each module
        proteins_per_module[i] <- sum(nchar(protein_names_for_hub_genes) > 0)  # Count the number of proteins for the hub genes
      }
    }
    
    # Store hub counts next to modules, including protein counts
    hub_counts[[trait]] <- c(paste0(non_na_modules, " (", hubs_per_module, "/", proteins_per_module, ")"), 
                             rep("", nrow(significant_modules) - length(non_na_modules)))
    
    # Add "(Hubs/Proteins)" to column names
    colnames(hub_counts)[which(colnames(hub_counts) == trait)] <- paste0(trait, " (#Hubs/#Proteins)")
    
    # Combine results for the trait
    final_result <- do.call(rbind, result_list)
    
    if (!is.null(final_result)) {
      output_file <- file.path(output_dir, paste0(trait, "_HubGenes.txt"))
      write.table(final_result, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      cat("Results for", trait, "saved to", output_file, "\n")
    }
  }
  
  # Write the hub_counts to a tab-delimited file
  hub_counts_file <- file.path(output_dir, "HubCounts.txt")
  write.table(hub_counts, file = hub_counts_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("Hub counts saved to", hub_counts_file, "\n")
  
  # Return updated significant_modules with hub counts
  return(hub_counts)
}

Hubs <- HubGenes_extraction(significant_modules, all_traits_associations, datKME, datExpr, hubGenes_path, proteome_fasta)

save(MEs, datTraits, moduleTraitCor, moduleTraitPvalue, 
     all_traits_associations, significant_modules,sigMods_unique,colorh1, ADJ1, Alldegrees1, Hubs, file = paste0(mod2traits_path, "3.Mods2Traits.RData"))

```