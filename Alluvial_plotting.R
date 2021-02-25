# This script can be used to make alluvial plots  
# from 3D RNAseq output using a number of available R packages.
# Part of the publication Schreiber M, Orr J, Barakate A, Waugh R (2021)
# Barley (Hordeum vulgare) anther and meiocyte RNA sequencing: mapping
# sequencing reads and downstream data analyses.

# load required packages ------------------------------------------------------
library(ggalluvial) # You will need this to make the alluvial plot

# Getting the data ready ------------------------------------------------------
# First we must get the 3D RNAseq output statistics file into a format
# that these packages can interpret. 

# load in the 3D RNA seq differential expression statistics data
differentialGenes <- read.csv("3D output/result/DE gene testing statistics.csv")

# Take a look at how the testing statistics are currently organised
colnames(differentialGenes)

# ggaluvial takes data in either wide or long format. Happily, this is 
# simpler to prepare from 3D RNA-seq. output than for venn or upset 
# plotting.
#
# We want to plot the change in expression of each gene in each comparison
# in order of meiotic stages in anthers. 
#
# the only change to the input data format that is needed is the addition of
# an extra categorical column that allows us to group DEGs into
# up-regulated (UP), down-regulated (DOWN), and no significant 
# difference (NSD)
#
# We can create a new column called "change" and set its value for
# each row based on the adj.pval and log2FC columns. 
#
# The code below uses a nested ifelse statement (an ifelse statment with 
# another ifelse statement inside it).
# 
# First, if the adj.pval column meets our DEG cutoff (<0.01) and the log2FC
# column is > 1 we consider this gene significantly up-regulated so the value
# in our new column is set to "UP".
# 
# If these conditions are not satisfied we check the adj.pval and log2FC
# columns again. This time, if the adj.pval cutoff is met and the log2FC
# value is < -1 we consider this gene significantly down-regulated so the
# value of the row in the new column is set to "DOWN".
#
# If this second test is not met either there is no significant
# change in expression and the value of the row in the new column is 
# set to "NSD" meaning "no significant difference"
differentialGenes$change <- ifelse(differentialGenes$adj.pval < 0.01 &
                                differentialGenes$log2FC > 1, "UP",
                                ifelse(differentialGenes$adj.pval < 0.01 &
                                       differentialGenes$log2FC < -1, "DOWN",
                                       "NSD"))

# As only the anther-to-anther and meiocyte-to-meiocyte contrast groups 
# could be considered sequential we will next restrict the data to the
# anther-to-anther contrast groups. 

# first check the contrast names
levels(differentialGenes$contrast)

# The "-" in the contrast names is a headache so we'll replace it with 
# "_"
differentialGenes$contrast <- gsub("-", # substitute "-"
                                   "_", # for "_"
                                   differentialGenes$contrast) # in this column

# the above changes the class of values in this column to "character"
# we will need to change it back to "factor"
differentialGenes$contrast <- as.factor(differentialGenes$contrast)

# Then subset the desired groups
DEG <- subset(differentialGenes, # data frame to subset
              contrast == "A.Pre_A.Lep.Zyg" | # rows for this contrast or...
              contrast == "A.Lep.Zyg_A.Pac.Dip" | # or...
              contrast == "A.Pac.Dip_A.Met.Tet")

# Next, we will reduce the genes to those that are differentially expressed in
# any of these three comparisons only. 
#
# First, we can get a data frame with only those rows corresponding to
# significant differential expression in any of our anther comparisons
sig_DEG <- subset(DEG, 
                  change != 'NSD')

# Then we can use the vector of gene names in the "target" column of the
# above sig_DEG frame to reduce the anther-to-anther comparisons to only those
# that are significant in at least one comparison. 
DEG <- DEG[DEG$target %in% sig_DEG$target, ]

# Making the alluvial plot ----------------------------------------------------
#
# First, choose a colourblind firendly pallette.  
# This is the IBM pallette recovered from: 
# https://davidmathlogic.com/colorblind/#%23648FFF-%23785EF0-%23DC267F-%23FE6100-%23FFB000
IBM <- c("#648FFF",
         "#785EF0",
         "#DC267F",
         "#FE6100",
         "#FFB000")

# ggalluvial is an extension of the ggplot2 library which will plot factors
# in order of their level by default. Unless specified this is alphabetical.
#
# We can change the order in which the "change" value will appear in our plot
# by converting this column to a "factor" and setting its levels in the 
# order we would like.
#
# Check what the class of the values in the "change" column are already.
class(DEG$change)

# Convert this to factor using the as.factor function
DEG$change <- as.factor(DEG$change)

# Look at the default levels of this factor column. 
levels(DEG$change)

# Chage the order from "DOWN", "NSD", "UP" to "Up", "NSD", "DOWN"
DEG$change <- factor(DEG$change,
                     levels(DEG$change)[c(3,2,1)])

# Similarly, we want to plot the contrasts in order of meiotic progression
# from A.Pre to A.Met.Tet
#
# We can check and change the order as above.
#
# Note that list of levels for this factor will include all the original values
# of this column, including the ones we have excluded already.
levels(DEG$contrast)
DEG$contrast <- factor(DEG$contrast,
                       levels(DEG$contrast)[c(6,2,4)])

levels(DEG$contrast)

# All that remains is to make the plot. 
#
# Note that this will plot the flow of each gene (alluvia) between contrast
# groups (axis) depending on the DEG status (strata). 
#
# Because each gene is plotted, this can take some time to comlplete.
ggplot(DEG, # use the DEG data frame
      aes(x = contrast, # plot contast group on the x axis
          stratum = change, # group the strata by DE status
          alluvium = target, # chart the flow of genes
          fill = change, # colour the strata by DE status
          label = change, # label by the DE status
          color=change)) + # colour the alluvia by DE status
scale_fill_manual(values = IBM[c(1,3,5)]) + # use the 1st, 3rd, and 5th IBM colours
scale_color_manual(values = IBM[c(1,3,5)]) + # for both fill and colour
geom_flow(stat = "alluvium", # plot flow between strata
         lode.guidance = "frontback", # the position of each gene in a strata (lode) is able to move back and forth
         alpha=0.1) + # set the alluvia between strata to 10% opacity
geom_stratum(alpha=0.5) + # set the strata to 50% opacity
theme(legend.position = "bottom", # place the legend at the bottom
      text = element_text(size = 20)) + 
xlab("Comparison") + # label the x axis
ylab("Differentially Expressed Genes") # label the y axis

