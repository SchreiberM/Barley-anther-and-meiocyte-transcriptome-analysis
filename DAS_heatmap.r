### GO enrichment analysis script
### Part of the publication Schreiber M, Orr J, Barakate A, Waugh R (2021) Barley (Hordeum vulgare) anther and meiocyte RNA sequencing: mapping sequencing reads and downstream data analyses.

##########################################################################
### Required R libraries (require Bioconductor)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gridExtra)

## differential alternative spliced output table plus the epression for all transcripts and genes in TPM from 3DRNASeqApp
differentialAS <- read.csv("Significant DAS genes list and statistics.csv")
tpmExpression <- read.csv("Transcript TPM.csv")
geneExpression <- read.csv("Gene TPM.csv")


## load in the transcript to gene identification file which was one of the input files for 3D RNA-seq app
transcriptToGene <- read.csv("TranscriptToGene.csv")


TPMplusGeneId <- cbind(transcriptToGene[,2], tpmExpression)

##########################################################################
### Decide which contrast to look at
t.contrast <- 
## Example contrast for the barley anther transcriptome
#t.contrast <- "Anther.Premeiosis-Anther.Leptotene.Zygotene"

### Decide on p-value
p.value <- 0.01

### Decide on the percentage of differences between the isoforms, can be positiv or negative
maxDelta <- 0.5

## maxDelta positive
OneCondition <- differentialAS[(differentialAS$contrast==t.contrast & differentialAS$adj.pval <= p.value  & differentialAS$maxdeltaPS >= maxDelta),]

## maxDelta negative
# OneCondition <- differentialAS[(differentialAS$contrast==t.contrast & differentialAS$adj.pval <= p.value  & differentialAS$maxdeltaPS <= maxDelta),]


###### Choose one gene for plotting a heatmap of differential alternative spliced transcripts

## In this case it takes the top gene on the list
GeneOfInterest <- OneCondition$target[1]

## Merge the transcript expression with the corresponding gene expression
CombinedResults <- rbind(geneExpression[geneExpression[,1]==GeneOfInterest, ] ,TPMplusGeneId[TPMplusGeneId[,1]==GeneOfInterest, 2:dim(TPMplusGeneId)[2]])


#### define short names for the x and y axis labelling - long names for the experimental design are unfortunately often unreadable
## Example labeling for the barley anther transcriptome
labelsStaging <- c()
# labelsStaging <- c("A.PRE.1", "A.PRE.2", "A.PRE.3", "A.LEP-ZYG.1", "A.LEP-ZYG.2", "A.LEP-ZYG.3", "A.PAC-DIP.1", "A.PAC-DIP.2",  "A.PAC-DIP.3", "A.MET.TET.1", "A.MET.TET.2",  "A.MET.TET.3", "M.LEP-ZYG.1", "M.LEP-ZYG.2", "M.LEP-ZYG.3", "M.PAC-DIP.1", "M.PAC-DIP.2",  "M.PAC-DIP.3")


### change dataframe structure for plotting
meltTranscripts <- melt(CombinedResults, id.vars = "X")


############ Plotting
## p1 generates the first figure which is a heatmap
p1 <- ggplot(meltTranscripts, aes(variable, X, fill=value)) + 
  geom_tile() + scale_fill_distiller(palette = "GnBu") + theme_classic() + labs(x=NULL, y=NULL, fill="TPM", title=paste("Heatmap", GeneOfInterest, sep=" "), tag="A") + scale_x_discrete(labels=labelsStaging) + 
  theme(legend.position="bottom", axis.text.x = element_text(size=7, angle=90), axis.text.y=element_text(size=7, angle=60,hjust=0.5), title=element_text(size=7, face="bold"), legend.title=element_text(size=7)) 

## p2 generates the second figure which is a line graph
p2 <- ggplot(meltTranscripts, aes(variable, value, group=X)) + geom_line(aes(color=X), size=0.6) + geom_point(aes(color=X), size=0.6) 
p2 <- p2 + scale_colour_brewer(palette="Dark2")+ theme_classic()  + labs(x=NULL, y="TPM", colour=NULL, title=paste("Line plot", GeneOfInterest, sep=" "), tag="B") + scale_x_discrete(labels=labelsStaging) + 
    theme(legend.position="bottom", axis.text.x = element_text(size=7, angle=90), axis.text.y=element_text(size=7), title=element_text(size=7, face="bold"), legend.text=element_text(size=6)) + 
    guides(colour = guide_legend(ncol = 2))

## arranging the two figures next to each other for the output
p3 <- grid.arrange(p1, p2, nrow = 1)

## saving the figure as .tiff in 300 dpi and 7.5 x 4.5 inches
ggsave(paste(GeneOfInterest, "DifferentialAlternativeSplicingPlot.tiff", sep="_"), plot=p3, width=7.5, height=4.5, dpi=300, compression="lzw")



 ####################### Example on plotting the top 10 differential expressed genes in a heatmap
## differential expressed gene output table from 3DRNASeqApp
differentialExpressedGenes <- read.csv("Significant DE genes list and statistics.csv")
geneExpression <- read.csv("Gene TPM.csv")

## Select the top differential expressed genes. The table is sorted by adjusted p.value in ascending order. So the top genes are the most significant.
TopdiffGenes <- geneExpression[geneExpression[,1]%in%differentialExpressedGenes[1:10,1], ]

## change dataframe structure for plotting
meltGenes <- melt(TopdiffGenes, id.vars = "X")

## log2 transform the tpm values for better visualisation of the results
meltGenes$value <- log2(meltGenes$value+1)

#### define short names for the x and y axis labelling - long names for the experimental design are unfortunately often unreadable
## Example labeling for the barley anther transcriptome
labelsStaging <- c()
# labelsStaging <- c("A.PRE.1", "A.PRE.2", "A.PRE.3", "A.LEP-ZYG.1", "A.LEP-ZYG.2", "A.LEP-ZYG.3", "A.PAC-DIP.1", "A.PAC-DIP.2",  "A.PAC-DIP.3", "A.MET.TET.1", "A.MET.TET.2",  "A.MET.TET.3", "M.LEP-ZYG.1", "M.LEP-ZYG.2", "M.LEP-ZYG.3", "M.PAC-DIP.1", "M.PAC-DIP.2",  "M.PAC-DIP.3")

# Plotting the figure
ggplot(meltGenes, aes(variable, X, fill=value)) + 
  geom_tile() +
  scale_fill_distiller(palette = "GnBu") + theme_classic() + labs(x=NULL, y=NULL, fill="logs2(TPM)", title=paste("Heatmap differential expressed genes", sep=" ")) + scale_x_discrete(labels=labelsStaging) + 
  theme(legend.position="bottom", axis.text.x = element_text(size=7, angle=90), axis.text.y=element_text(size=7), title=element_text(size=7, face="bold"), legend.title=element_text(size=7)) 

## saving the figure as .tiff in 300 dpi and 7.5 x 4.5 inches
ggsave("HeatmapTopDiffGenes.tiff", width=7.5, height=4.5, dpi=300, compression="lzw")
