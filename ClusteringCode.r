ccdata = read.delim("ChIP-chip.txt", header=TRUE, sep="\t", row.names=1)

halo.cor <- cor(ccdata, use="pairwise.complete.obs")
halo.dist <-as.dist(1-halo.cor)
halo.tree <-hclust(halo.dist, method="complete")

plot(halo.tree)