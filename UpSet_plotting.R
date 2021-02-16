# This script can be used to make upset plots  
# from 3D RNAseq output using a number of available R packages.
# Part of the publication Schreiber M, Orr J, Barakate A, Waugh R (2021)
# Barley (Hordeum vulgare) anther and meiocyte RNA sequencing: mapping
# sequencing reads and downstream data analyses.

# load required packages ------------------------------------------------------
library(tidyverse) # You will need this package to reshape the data
library(UpSetR) # You will need this to make the upset plot

# Getting the data ready ------------------------------------------------------
# First we must get the 3D RNAseq output statistics file into a format
# that these packages can interpret. 

# load in the 3D RNA seq differential expression statistics data
differentialGenes <- read.csv("3D_output/DE_gene_testing_statistics.csv")

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

# UpsetR requires a data frame where each group is a column and the row values
# are a numerical binary representation of whether a member of the group (genes
# in our case) are present or absent (1 = present, 0 = absent).
#
# In this example we will plot the intersectons of differentially expressed genes
# between contrast groups. 
# 
# First we will convert the testing statistics (adjusted P value and LFC) to 
# an indication of whether or not we consider each gene to be differentially 
# expressed in each contrast group.
#
# The loop below iterates through each adjusted P value column (2 to 7)
# If the value in the column is <0.01 and the corresponding LFC value 
# (6 columns along; x+6) is > 1 
# OR (|) 
# If the value in the column is <0.01 and the corresponding LFC value
# is < -1 then the value is replaced with 1. Otherwise, if none of these
# conditions is met the value is replaced with 0.
for (x in 2:7) {
  DEG_wide[,x] <- ifelse(DEG_wide[,x]<0.01 & DEG_wide[,x+6] > 1 |
                           DEG_wide[,x]<0.01 & DEG_wide[,x+6] < -1, 1, 0)
}

# Take a quick look at what this has done
DEG_wide[1:10,] # print the first ten rows

# We no longer need the LFC rows so we'll get rid of those
DEG_wide <- DEG_wide[,1:6] 

# Each column now represents a contrast group. We can remove the measurement
# prefix from the column names as it's no longer relevant. We will use the
# gsub function for this.
names(DEG_wide) <- gsub(x = names(DEG_wide), # in the column names in DEG_wide
                        pattern = "adj.pval_", # replace this string in quotes
                        replacement = "") # with nothing ("")

# UpstR cannot cope with rows (here, genes) returning 0 (here, not 
# differentially expressed) in any group/set.
# 
# To account for this we can add a new column which contains the sum of all rows
DEG_wide$total_sig_contrasts <- rowSums(DEG_wide[,2:6])

# This column can be useful for subsetting particular genes later on. 
# 
# For now we can use the subset function to remove genes which are not
# differentially expressed in any contrast group
sig_DEG <- subset(DEG_wide, # input data frame
                  total_sig_contrasts > 0) # return only the rows whose value in our new column is > 0

# UpsetR also strugles to handle the tibble R data format so we will
# need to convert it to a data frame.
sig_DEG <- as.data.frame(sig_DEG)

# Plotting intersections with UpsetR ------------------------------------------
# Now that we have our data in the required format we can create an upset plot.
#
# First, store the name of each contrast group (column headers) as a vector
# This will save you typing out your group names later. 
setnames <- colnames(sig_DEG[,2:6])

# Then we can generate our UpSet Plot
upset(sig_DEG, # data frame
      order.by = "freq", # sort by largest to smallest
      nintersects = 25, # only plot the top 25 intersects
      sets = setnames, # sets to look for intersects
      keep.order = TRUE) # Keep the order of sets in setnames

# With upset you can add or highlight information in the plot by 
# annotating or colouring different aspects.
# 
# As an example we can add the coding status of our genes to the plot
# by setting the colour of the bars in the histogram above the intersect
# matrix.

# load in the results of our gene coding status analysis
codingstat <- read.csv("GeneCodingStatus.csv")
colnames(codingstat)

# We can combine the coding status and DEG data by the gene names but first
# the name of the column containing the names needs to be identical in 
# both data frames.
#
# We can do this using the rename function from dplyr (which is part of the 
# tidyverse collection of R packages)

# currently the column containg gene names in the DEG data frame is called
# "target"
colnames(sig_DEG)

# Lets change it to BAnTr
sig_DEG <- rename(sig_DEG, # data frame in which to rename
                  "BAnTr" = "target") # new column name = old column name

# check that this worked 
colnames(sig_DEG)

# now that the column names match we can merge them
sig_DEG_coding <- merge(sig_DEG, # data frame 1
                        codingstat, # data frame 2
                        by="BAnTr") # shared column name

# This has created a new data frame which merges both data sets
colnames(sig_DEG_coding)

# The coding.noncoding column contains the data we want to add to the upset
# intersect plot. It has three possible values: "coding", "noncoding", and 
# "undefined".
levels(sig_DEG_coding$coding.noncoding)

# Next, choose a colourblind firendly pallette.  
# This is the IBM pallette recovered from: 
# https://davidmathlogic.com/colorblind/#%23648FFF-%23785EF0-%23DC267F-%23FE6100-%23FFB000
IBM <- c("#648FFF",
         "#785EF0",
         "#DC267F",
         "#FE6100",
         "#FFB000")

# now we can add this information to the upset plot
# 
# This is slightly complicated. 
#
# In the code below, we generate the upset plot as before.
# Then we are colouring each bar in the intersection size histogram
# by the value in the coding.noncoding column.
#
# We have to layer this colouring operation to add more than one colour
# to each bar.
# 
# First, we colour the entire bar (all three values).
#
# Second, we colour only undefined and noncoding values leaving the remainder
# of the bar in the first colour to represent coding genes.
#
# Third, we colour only noncoding values leaving the section in between the
# first and last colour to represent "undefined" values.
upset(sig_DEG_coding,
      order.by = "freq",
      sets=setnames,
      nintersects = 25,
      keep.order = TRUE,
      queries = list(list(query = elements,
                          params = list("coding.noncoding",
                                        c("coding", "undefined", "noncoding")),
                          color = IBM[1],
                          active = T,
                          query.name = "coding"),
                     list(query = elements,
                          params = list("coding.noncoding",
                                        c("undefined", "noncoding")),
                          color = IBM[3],
                          active = T,
                          query.name = "undefined"),
                     list(query = elements,
                          params = list("coding.noncoding",
                                        "noncoding"),
                          color = IBM[5],
                          active = T,
                          query.name = "non-coding")),
      query.legend = "bottom")


