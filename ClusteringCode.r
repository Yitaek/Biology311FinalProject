## Biology 311 Final Project
# Exploring RosR - Clustering Analysis

## Loading all the required packages
library(dendextend)
library(gplots)

# Read in the Gene Expression Data
gedata = read.csv("GEdata.csv", header=T, row.names=1)

# Subset of the Parent Strain
parent = t(gedata[1:8])

# Subset of the RosR Deletion
RosRDel = t(gedata[9:16])

# Heatmap for the Parent Strain	
parent.T = t(parent)
parent.max.times <- apply(parent.T, 1, which.max)
parent.gene.order <- rev(order(parent.max.times))
reordered.parent <- parent.T[parent.gene.order, 1:8]

png("heatmap_reordered_parent.png", width=600, height=600)
heatmap(reordered.parent, Colv=NA, Rowv=NA, col=greenred(64), labRow=NA, cexCol=0.5)
dev.off()

# Clustering for the Parent Strain

my.dist <- function(x) {as.dist(1 - cor(x, use="pairwise.complete.obs"))}
parent.dist <- my.dist(parent)
parent.tree <- hclust(parent.dist, method="complete")
parent.dend <- as.dendrogram(parent.tree)

png("heatmap_with_dendrogram.png", width=600, height=600)
heatmap(parent.T, Colv=NA, Rowv=parent.dend, col=greenred(64), labRow=NA, cexCol=0.5)
dev.off()


colored.parent <- set(parent.dend, "branches_k_color", h=1)

parent.clusters <- cutree(parent.dend, k=7)
table(parent.clusters)

# VNG0258H gene is in cluster 6
parent.clust6 <- names(parent.clusters[parent.clusters==6])

