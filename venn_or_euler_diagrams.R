# This script can be used to make venn or euler diagrams 
# form 3D RNAseq output using a number of available R packages.
# Part of the publication Schreiber M, Orr J, Barakate A, Waugh R (2021)
# Barley (Hordeum vulgare) anther and meiocyte RNA sequencing: mapping
# sequencing reads and downstream data analyses.

# load required packages ------------------------------------------------------
library(tidyverse) # You will need this package to reshape the data

# You will need one of these packages to draw a venn diagram
library(VennDiagram)
library(limma)
library(ggvenn)
library(ggVennDiagram)

# you will need one of these packages to draw a Euler diagram
library(venneuler)
library(eulerr)

# Getting the data ready ------------------------------------------------------
# First we must get the 3D RNAseq output statistics file into a format
# that these packages can interpret. 

# load in the 3D RNA seq differential expression statistics data
differentialGenes <- read.csv("3D output/result/DE transcripts testing statistics.csv")

# We need to convert the data frame from long format (where both genes
# in the target column and the contrast groups in the contrast column) 
# to a wider format (where each gene is represented by a single row, 
# and contrast group scoring statistics are represented as columns)

# Take a look at how the testing statistics are currently organised
colnames(differentialGenes)

# We will use the pivot_wider function (part of the tidyverse collection
# of R packages) to reshape the data.
DEG_wide <- pivot_wider(differentialGenes, # input data frame
                        names_from = contrast, # new column name suffix
                        values_from = 3:4) # values the new columns will contain

# Take a look at the new, wider structure of the data
colnames(DEG_wide)

# Most intersection plotting functions in R require one of two data structures:
#   1. A data frame where each group (circle in your venn diagram) is a column 
#      and the row values are a binary representation of whether a member of the 
#      group (genes in our case) are present or absent. This can be represented 
#      numerically (1 = present, 0 = absent) or with a logical (TRUE = present,
#      FALSE = absent).
#    2. An R list object which contains the names of your groups and its members.
#
# In this example we will plot the intersectons of differentially expressed genes
# between contrast groups. 
# 
# First we will convert the testing statistics (adjusted P value and LFC) to 
# an indication of whether or not we consider each gene to be differentially 
# expressed in each contrast group.
#
# The loop below iterates through each adjusted P value column (2 to 8)
# If the value in the column is <0.01 and the corresponding LFC value 
# (7 columns along; x+7) is > 1 
# OR (|) 
# If the value in the column is <0.01 and the corresponding LFC value
# is < -1 then the value is replaced with the logical binary value "TRUE".
# Otherwise, if none of these conditions is met the value is replaced with "FALSE".
for (x in 2:8) {
  DEG_wide[,x] <- ifelse(DEG_wide[,x]<0.01 & DEG_wide[,x+7] > 1 |
                         DEG_wide[,x]<0.01 & DEG_wide[,x+7] < -1, TRUE, FALSE)
}

# Take a quick look at what this has done
DEG_wide[1:10,] # print the first ten rows

# We no longer need the LFC rows so we'll get rid of those
# For ease of plotting we'll keep only the five direct anther to anther
# and anther to meioctye comparisons (columns 1:6)
DEG_wide <- DEG_wide[,1:6]

# Each column now represents a contrast group we can remove the measurement prefix
# from the column names as it's no longer relevant. We will use the gsub function for this.
names(DEG_wide) <- gsub(x = names(DEG_wide), # in the column names in DEG_wide
                        pattern = "adj.pval_", # replace this string in quotes
                        replacement = "") # with nothing ("")

# For our practical reasons the "-" in the contrast column is a nuisance so we'll
# replace this with "_"
names(DEG_wide) <- gsub(x = names(DEG_wide), # in the column names in DEG_wide
                        pattern = "-", # replace this string in quotes
                        replacement = "_") # with nothing ("")

# This is now ready for input into most of our plotting functions.
# However, one of our functions requires the R list object format.
# We can use the format we've prepared above to enerate this.
# First, we have to make a vector for each contrast group we want to
# plot that contains the names of all it's members. Because most venn diagram
# plotting functions preclude groups of more than four we'll just create a list
# object for each of the anther to anther comparisons. 
# The code below returns the gene names (DEG_wide$target) if the gene is differentially
# expressed (indicated as TRUE in DEG_wide as prepared above) in the relevant contrast group
# column (DEG_wide$`A.LEP.ZYG-A.PRE`).
A.Pre_A.Lep.Zyg <- DEG_wide$target[DEG_wide$A.Pre_A.Lep.Zyg] 
A.LepZyg_A.PacDip <- DEG_wide$target[DEG_wide$A.Lep.Zyg_A.Pac.Dip]
A.Pac.Dip_A.Met.Tet <- DEG_wide$target[DEG_wide$A.Pac.Dip_A.Met.Tet]

# These vectors can then be combined as an R list object
x <- list(A.Pre_A.Lep.Zyg,
          A.LepZyg_A.PacDip,
          A.Pac.Dip_A.Met.Tet)

# Plotting Venn diagrams ------------------------------------------------------
# Now that we have our data in the required format we can try out different 
# plotting functions. 
#
# First, store the name of each contrast group (column headers) as a vector
# This will save you typing out your group names later. 
setnames <- colnames(DEG_wide[,2:6])

# Next, choose a colourblind firendly pallette for your figures.
# This is optional as all plotting functions below have an inbuilt
# defualt colour pallette. Although most are not colourblind friendly.  
# This is the IBM pallette recovered from: 
# https://davidmathlogic.com/colorblind/#%23648FFF-%23785EF0-%23DC267F-%23FE6100-%23FFB000
IBM <- c("#648FFF",
         "#785EF0",
         "#DC267F",
         "#FE6100",
         "#FFB000")

# Plot a venn diagram with limma.
# 
# Limma's vennDiagram function lets us plot five contrast groups.
vennDiagram(DEG_wide[2:6], # plot columns 2:6
            circle.col = IBM) # use the IBM colour palette
# This shows a lot of information but is very difficult to read

# Plot a venn diagram with ggvenn.
#
# ggvenn allows ven diagram plotting with the popular ggplot2 graphics
# format. It allows intersects bertween up to four groups to be plotted
ggvenn(DEG_wide,
  setnames[1:3], # plot only the first three contrasts (anther comparisions)
  fill_color = IBM[c(1,3,5)],
  stroke_size = 0.5, # the thickness of the circle outlines
  set_name_size = 6, # the size of the text of group labels
  text_size = 6) # the size of numerical labels

# Plot with ggVennDiagram.
#
# This is similar to ggvenn, using ggplot2 formatting. However, instead
# of colouring by group this function colours by the size of the intersection.
# 
# This is the one function in this script that uses the R list format data only.
# It is also limited to a maximum of four groups.
ggVennDiagram(x, # use the R list data as input
        category.names = setnames[1:3])

# Plot Euler diagrams ---------------------------------------------------------
#
# Euler diagrams have two main advantages over venn diagrams for larger groups.
# The size of the circles representing each group is proportional to their total
# size, so they convey this information more intuitively. This also makes 
# plotting a large number of groups somewhat less cluttered. 

# Plot a euler diagram using venneuler
# 
# The venneuler  function first creates a list object from the input data
vd <- venneuler(DEG_wide[,2:6])
# Which can then be plotted
plot(vd, # Using the above list object as input
     col = IBM,
     cex = 2)
# Venneuler has somewhat limited flexibility in the plot format.

# Plot a euler diagram using eulerr 
#
# Eulerr is based on venneuler with added functionality.
# As above, first an R list object is created by the euler function
vd <- euler(DEG_wide[,2:6])
# Which can then be plotted
plot(vd,
legend = list(cex = 1.5), # include a figure legend at 1.5x default size
quantities = list(cex = 1.5), # print the numbers in each group and intersect at 1.5x default size
fills = IBM,
alpha = 0.6)
