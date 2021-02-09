### GO enrichment analysis script
### Part of the publication Schreiber M, Orr J, Barakate A, Waugh R (2021) Barley (Hordeum vulgare) anther and meiocyte RNA sequencing: mapping sequencing reads and downstream data analyses.

##########################################################################
### Required R libraries (require Bioconductor)
library(topGO)
library(GOplot)


##########################################################################
### Read in required files into R

## GO annotation look-up table
geneID2GO <- readMappings("")

## Example table for the barley anther transcriptome:
# geneID2GO <- readMappings("BAnTR_GO_Annotation.txt")

## differential expressed gene output table from 3DRNASeqApp
differentialGenes <- read.csv("DE gene testing statistics.csv")


##########################################################################
### Decide which contrast to look at
t.contrast <- 
## Example contrast for the barley anther transcriptome
#t.contrast <- "Anther.Premeiosis-Anther.Leptotene.Zygotene"

### Decide on p-value
p.value <- 0.01

### Decide on log2FC (positive for upregulated, e.g. 2; negative for downregulated, e.g. -2)
logfold <- 2

### Decide on studying up- or downregulation and uncomment the line below

## upregulation
OneCondition <- differentialGenes[(differentialGenes$contrast==t.contrast & differentialGenes$adj.pval <= p.value & differentialGenes$log2FC >= logfold),]

## downregulation
# OneCondition <- differentialGenes[(differentialGenes$contrast==t.contrast & differentialGenes$adj.pval <= p.value & differentialGenes$log2FC <= logfold),]



##########################################################################
## Gene ontology enrichment pipeline using topGO

## allocate GO terms to the gene subset chosen
geneNames <- names(geneID2GO)
candidates <- factor(as.integer(geneNames %in% OneCondition$target))
names(candidates) <- geneNames


## GO enrichment for the molecular function (MF)
GO_MF <- new("topGOdata", ontology="MF", allGenes=candidates, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=5)
## To test for encrichment this uses fisher statistics and the weight01 algorithm. 
resultFisherWei <- runTest(GO_MF, algorithm = "weight01", statistic="fisher")
# For significance using a p-value below 0.001. If no significant terms can be identified it will still print the two top terms.
mysummary <- summary(attributes(resultFisherWei)$score <= 0.001)
if(length(mysummary) == 2){
	allResMF <- data.frame(colnames(GO_MF))
} else{
	numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001
	if(numsignif == 1){
		allResMF <- GenTable(GO_MF, weight01Fisher=resultFisherWei, topNodes=2)
		allResMF <- allResMF[1,]
	} else{
		allResMF <- GenTable(GO_MF, weight01Fisher=resultFisherWei, topNodes=numsignif)
	}
}
	
category <- as.vector(rep("MF", dim(allResMF)[1]))
output1 <- cbind(allResMF, category)

## GO enrichment for the biological process (BP)
GO_BP <- new("topGOdata", ontology="BP", allGenes=candidates, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=5)
## To test for encrichment this uses fisher statistics and the weight01 algorithm.
resultFisherWei <- runTest(GO_BP, algorithm = "weight01", statistic="fisher")
# For significance using a p-value below 0.001. If no significant terms can be identified it will still print the two top terms.
mysummary <- summary(attributes(resultFisherWei)$score <= 0.001)
if(length(mysummary) == 2){
	allResBP <- data.frame(colnames(GO_BP))
} else{
	numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001
	if(numsignif == 1){
		allResBP <- GenTable(GO_BP, weight01Fisher=resultFisherWei, topNodes=2)
		allResBP <- allResBP[1,]
	} else{
		allResBP <- GenTable(GO_BP,  weight01Fisher=resultFisherWei, topNodes=numsignif)
	}
}

category <- as.vector(rep("BP", dim(allResBP)[1]))
output2 <- cbind(allResBP, category)

## GO enrichment for the cellular component (CC)
GO_CC <- new("topGOdata", ontology="CC", allGenes=candidates, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=5)
## To test for encrichment this uses fisher statistics and the weight01 algorithm.
resultFisherWei <- runTest(GO_CC, algorithm = "weight01", statistic="fisher")
# For significance using a p-value below 0.001. If no significant terms can be identified it will still print the two top terms.
mysummary <- summary(attributes(resultFisherWei)$score <= 0.001)
if(length(mysummary) == 2){
	allResCC <- data.frame(colnames(GO_CC))
} else {
	numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001
	if(numsignif == 1){
		allResCC <- GenTable(GO_CC, weight01Fisher=resultFisherWei, topNodes=2)
		allResCC <- allResCC[1,]
	} else{
		allResCC <- GenTable(GO_CC, weight01Fisher=resultFisherWei, topNodes=numsignif)
	}
}
	
category <- as.vector(rep("CC", dim(allResCC)[1]))
output3 <- cbind(allResCC, category)



##########################################################################
### Preparing final dataframe for writing to file by merging all three gene ontology terms output
GOMerged<- rbind(output1,output2,output3)
	
### Defining z-score as percentage of signifcant genes in comparison to total terms assigned, adjust header and write to a csv output file
GOResults <- cbind(GOMerged, GOMerged$Significant/GOMerged$Annotated)
colnames(GOResults) <- c("ID", "term", "Annotated", "count", "Expected", "adj_pval", "category", "zscore")
GOResults$zscore <- GOResults$zscore*100
write.csv(GOResults, file=paste(t.contrast, "GeneOntologyEnrichment.csv",sep="_"))


##########################################################################
### Plotting gene ontology enrichment results using GOplot
### Need to change p <1e-30 as it can not be plotted and therefore will be set to 1e-30
GOResults$adj_pval <- sub("< 1e-30", "1e-30", GOResults$adj_pval)
GOResultsPlotting <- transform(GOResults, adj_pval=as.numeric(adj_pval))

### Plot and save the results
tiff(paste(t.contrast, "GeneOntologyEnrichment.tiff",sep="_"))
GOBubble(GOResultsPlotting, labels=5, display = 'multiple')
dev.off()


