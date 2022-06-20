# Loading required packages
library(NetCoMi)
library(phyloseq)
library(WGCNA)

# Setting the working directory
setwd("F:/Research/Methodology")

# Loading the wet zone phyloseq object
Wet_data <- readRDS("Wet Zone.rds")

# Constructing the network
net_single <- netConstruct(Wet_data,
                           measure = "sparcc",
                           measurePar = list(iter=20),
                           sparsMethod = "t-test",
                           alpha = 0.05,
                           adjust = "adaptBH",
                           verbose = 3,
                           seed = 123456)

saveRDS(net_single, "Wet Zone/Wet network.rds")


# Analyzing the network
props_single <- netAnalyze(net_single, 
                           centrLCC = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree","eigenvector"),
                           weightDeg = FALSE,
                           normDeg = FALSE,
                           verbose = 3)

saveRDS(props_single, "Wet Zone/Wet analysis.rds")

# Summary of the analysis
net.summary <- summary(props_single)


# Preparing the network for plotting
p <- plot(props_single,
     shortenLabels = "none",
     labelScale = FALSE,
     rmSingles = "all",
     nodeSize = "degree",
     nodeColor = "cluster",
     hubBorderCol = "blue",
     cexNodes = 1,
     cexLabels = 1,
     edgeWidth = 1,
     highlightHubs = TRUE,
     cexHubs = 1.5,
     title1 = "Wet Zone Network on genus level with SparCC Method.", 
     showTitle = TRUE,
     cexTitle = 1.5)

saveRDS(p, "Wet Zone/Wet plot.rds")

# Extracting required information for exporting the network
from <- p$labels$labels1[p$q1$Edgelist$from]
to <- p$labels$labels1[p$q1$Edgelist$to]
direction <- vector(mode = "integer")

for (x in 1:length(from)){
        if (net_single$assoMat1[from[x],to[x]] < 0){
                direction[x] <- -1 
        }
        else if (net_single$assoMat1[from[x],to[x]] > 0){
                direction[x] <- 1 
        }
} 

edges <- data.frame(row.names = c(1:length(p$q1$Edgelist$from)))
edges$from <- p$q1$Edgelist$from
edges$to <- p$q1$Edgelist$to
edges$weight <- p$q1$Edgelist$weight
edges$association <- direction
        
write.csv(edges, file = "Wet Zone/Wet edge data.csv", row.names = FALSE)

hubs <- props_single$hubs$hubs1
write(hubs, "Wet Zone/Wet Zone Hubs.txt")

node.lables <- p$labels$labels1
clust <- props_single$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1


nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs
degrees <- props_single$centralities$degree1[nodes$lable]
nodes$degree <- degrees

write.csv(nodes, file = "Wet Zone/Wet node data.csv", row.names = FALSE)

# exporting data of the keystone species
hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Wet Zone/Wet hub data.csv", row.names = FALSE)


hubinteractions <- edges[(edges$from %in% hubdata$index) | (edges$to %in% hubdata$index),]
from <- p$labels$labels1[hubinteractions$from]
to <- p$labels$labels1[hubinteractions$to]
assoc <- vector(mode = "integer")

for (x in 1:length(from)){
                assoc[x] <- net_single$assoMat1[from[x],to[x]] 

}
hubinteractions$association <- assoc

positives <- subset(hubinteractions, association>0.75, select = -weight)
positives$from_cl <- nodes$cluster[positives$from]
positives$to_cl <- nodes$cluster[positives$to]
positives$from <- node.lables[positives$from]
positives$to <- node.lables[positives$to]

negatives <- subset(hubinteractions, association < (-0.7), select = -weight)
negatives$from_cl <- nodes$cluster[negatives$from]
negatives$to_cl <- nodes$cluster[negatives$to]
negatives$from <- node.lables[negatives$from]
negatives$to <- node.lables[negatives$to]

func <- read.csv("Wet Zone/func result.csv",header = FALSE)
func$cluster <- nodes$cluster[func$V1]
for (x in 1:length(func$V1)){
        func$cluster[x] <- nodes$cluster[which(nodes$lable == func$V1[x])] 
        
}
