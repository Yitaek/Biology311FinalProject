## Biology 311 Final Project
# Exploring VNG0258H Gene

# Read in the Gene Expression Data
ccdata = read.delim("GEdata.txt", header=TRUE, sep="\t", row.names=1)

# Find the VNG0258H subset
RosR <- subset(ccdata, grepl("^VNG0258H", row.names(ccdata)))

# We want to replicate Figure 1 in the Sharma et al paper

# halo.cor <- cor(ccdata, use="pairwise.complete.obs")
# halo.dist <-as.dist(1-halo.cor)
# halo.tree <-hclust(halo.dist, method="complete")

# plot(halo.tree)