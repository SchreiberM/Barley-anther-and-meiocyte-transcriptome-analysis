# This script is an example of the preparation of publishable heat and line 
# plots from 3D-RNA-seq output.
#
# Part of the publication Schreiber M, Orr J, Barakate A, Waugh R (2021) Barley
# (Hordeum vulgare) anther and meiocyte RNA sequencing: mapping sequencing reads
# and downstream data analyses.

# Required R libraries (require Bioconductor) ---------------------------------
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gridExtra)

# Required data ---------------------------------------------------------------
# Here we will need the differential alternatively spliced output table
# plus the epression for all transcripts and genes in TPM from 3DRNASeqApp
differentialAS <- read.csv("3D output/result/Significant DAS genes list and statistics.csv")
tpmExpression <- read.csv("3D output/result/Transcript TPM.csv")
geneExpression <- read.csv("3D output/result/Gene TPM.csv")
differentialExpressedGenes <- read.csv("3D output/result/Significant DE genes list and statistics.csv")

# We will also need the transcript to gene identification file which was one of
# the input files for 3D RNA-seq app
transcriptToGene <- read.csv("3DBAnTr/TranscriptToGeneID.csv")

# Reshape the data for plotting -----------------------------------------------
TPMplusGeneId <- cbind(transcriptToGene[,2], tpmExpression)

# Decide which contrast to look at
#
# You can remind yourself of these contrasts by looking at the levels of the 
# DAS conrast column
levels(differentialAS$contrast)

# In the example dataset we will look at anthers at premeiosis compared to 
# anthers at leptotene-zygotene
t.contrast <- "A.Pre-A.Lep.Zyg"

# Decide on p-value
p.value <- 0.01

# Decide the percentage of differences between the isoforms (maxDelta), this can 
# be positive or negative
maxDelta <- 0.5

# we can subset the data to pull out only the transcritps meeting the criteria
# we have set out above with a positive maxdelta value...
OneCondition <- differentialAS[(differentialAS$contrast==t.contrast &
                                differentialAS$adj.pval <= p.value &
                                differentialAS$maxdeltaPS >= maxDelta),]

# # ...or a negative maxDelta negative
# OneCondition <- differentialAS[(differentialAS$contrast==t.contrast & 
#                                 differentialAS$adj.pval <= p.value &
#                                 differentialAS$maxdeltaPS <= maxDelta),]

# From this list we can pull out a gene of interest, here we select the gene
# correspoding to the maximum deltaPS value in the above subset
GeneOfInterest <- as.character(OneCondition[which.max(OneCondition$maxdeltaPS), 1])

# We can then merge the transcript expression with the corresponding gene expression
# to pull out the expression values for all isoforms and the gene as a whole 
# in all samples.
CombinedResults <- rbind(geneExpression[geneExpression[,1]==GeneOfInterest, ],
                         TPMplusGeneId[TPMplusGeneId[,1]==GeneOfInterest,
                                       2:dim(TPMplusGeneId)[2]])

# For plotting our heatmap we need to convert the data from a wide format
# to a long format. We can do that using the melt function from reshape2
meltTranscripts <- melt(CombinedResults,
                        id.vars = "X")

# Plotting --------------------------------------------------------------------

# Define short names for the x and y axis labelling 
#
# long names for the experimental design are unfortunately often unreadable
labelsStaging <- c("A.PRE.1", "A.PRE.2", "A.PRE.3",
                   "A.LEP-ZYG.1", "A.LEP-ZYG.2", "A.LEP-ZYG.3",
                   "A.PAC-DIP.1", "A.PAC-DIP.2",  "A.PAC-DIP.3",
                   "A.MET.TET.1", "A.MET.TET.2",  "A.MET.TET.3",
                   "M.LEP-ZYG.1", "M.LEP-ZYG.2", "M.LEP-ZYG.3",
                   "M.PAC-DIP.1", "M.PAC-DIP.2",  "M.PAC-DIP.3")

# Generate the first figure which is a heatmap
p1 <- ggplot(meltTranscripts,
             aes(variable, X, fill=value)) + 
      geom_tile() +
      scale_fill_distiller(palette = "GnBu") +
      theme_classic() +
      labs(x=NULL,
           y=NULL,
           fill="TPM",
           title=paste("Heatmap", GeneOfInterest, sep=" "), tag="A") +
      scale_x_discrete(labels=labelsStaging) + 
      theme(legend.position="bottom",
            axis.text.x = element_text(size=7, angle=90),
            axis.text.y=element_text(size=7, angle=60,hjust=0.5),
            title=element_text(size=7, face="bold"),
            legend.title=element_text(size=7)) 

# Generate the second figure which is a line graph
p2 <- ggplot(meltTranscripts,
             aes(variable, value, group=X)) +
      geom_line(aes(color=X), size=0.6) +
      geom_point(aes(color=X), size=0.6) +
      scale_colour_brewer(palette="Dark2") +
      theme_classic() +
      labs(x=NULL,
           y="TPM",
           colour=NULL,
           title=paste("Line plot", GeneOfInterest, sep=" "), tag="B") +
      scale_x_discrete(labels=labelsStaging) + 
      theme(legend.position="bottom",
            axis.text.x = element_text(size=7, angle=90),
            axis.text.y=element_text(size=7),
            title=element_text(size=7, face="bold"),
            legend.text=element_text(size=6)) + 
      guides(colour = guide_legend(ncol = 2))

# Arrange the two figures next to each other for the output
p3 <- grid.arrange(p1, p2, nrow = 2)

# saving the figure as .tiff at 600 dpi and 4.5 inches wide
#
# Nb: this sizing is appropriate for a book with a maximum usable width of 
# around 4.5 inches. Journal figures are typically around 180mm in width
# maximum. You can change the units of ggsave to mm by adding unit="mm".
ggsave(paste(GeneOfInterest,
             "DifferentialAlternativeSplicingPlot.tiff",
             sep="_"),
       plot=p3,
       width=4.5,
       dpi=600,
       compression="lzw")

# Example of plotting the top 10 differential expressed genes -----------------

# Select the top differentially expressed genes. The table is sorted by adjusted
# p.value in ascending order. So the top genes are the most significant.
TopdiffGenes <- geneExpression[geneExpression[,1] %in% differentialExpressedGenes[1:10,1], ]

# change dataframe structure for plotting as above
meltGenes <- melt(TopdiffGenes,
                  id.vars = "X")

# log2 transform the tpm values for better visualisation of the results
meltGenes$value <- log2(meltGenes$value+1)

# Plotting the heatmap
ggplot(meltGenes,
       aes(variable, X, fill=value)) + 
       geom_tile() +
       scale_fill_distiller(palette = "GnBu") +
       theme_classic() +
       labs(x=NULL,
            y=NULL,
            fill="log2(TPM)",
            title=paste("Heatmap differential expressed genes", sep=" ")) + 
       scale_x_discrete(labels=labelsStaging) + 
       theme(legend.position="bottom",
             axis.text.x = element_text(size=7, angle=90),
             axis.text.y=element_text(size=7),
             title=element_text(size=7, face="bold"),
             legend.title=element_text(size=7)) 

# saving the figure as .tiff in 600 dpi and 4.5 inches wide
ggsave("HeatmapTopDiffGenes.tiff",
       width=4.5,
       dpi=600)
