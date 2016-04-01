## Biology 311 Final Project
# Exploring RosR - Network Analysis

# Load the network functions
source("NetworkFxn.r")

# Read in the Gene Expression Data
ccdata = read.csv("ChIP-chip.csv", header=T)

# Get the TF Names
tf.names <- ccdata$Gene

# Get the experiment conditions
conditions <- names(ccdata)[2:5]

# Filter out all edges below p = 0.001 for each condition
tft <- get.all.targets(ccdata[,2:5], tf.names)

# Get counts for each condition
cts <- lapply(tft, length)
tf.cts <- unlist(cts)

# Get common genes across all of them
Reduce(intersect, tft)

## Creating incidence matrix

ccdata.mtx <- as.matrix(ccdata[,2:5])
rownames(ccdata.mtx) <-as.character(ccdata$Gene)
colnames(ccdata.mtx) <-names(ccdata)[2:5]

ccdata.mtx[ccdata.mtx < 0.001] <- -1
ccdata.mtx[ccdata.mtx > 0] <- 0
ccdata.mtx[is.na(ccdata.mtx)] <- 0
ccdata.mtx <- -1 * ccdata.mtx

# Creating Graph
library(igraph)

g <- graph_from_incidence_matrix(t(ccdata.mtx), direct = T, mode="out")

g <- delete.vertices(g, degree(g) == 0)

graph.tfs <- V(g)[V(g)$type == F]

V(g)$shape <- c("circle", "square")[V(g)$type+1]
V(g)$color <- c("steelblue", "lightcoral")[V(g)$type+1]

big.cc <- giant.component(g)

png("network_attempt.png")
plot(big.cc, vertex.label.cex=0.65, vertex.label.degree=180, vertex.label.dist=.1)
dev.off()

big.cc.tfs <- V(big.cc)[V(big.cc)$type == F]
tf.shared <- number.in.common(big.cc, big.cc.tfs$name)

g.shared <- graph_from_adjacency_matrix(tf.shared, weighted=T, mode="undirected")

V(g.shared)$size <- sqrt(degree(g.shared))
E(g.shared)$width <- sqrt(E(g.shared)$weight)

l <- layout_with_fr(g.shared, weights=E(g.shared)$weight)
plot(g.shared, vertex.label=NA, layout=l)