library(NetCoMi)
library(phyloseq)
library(WGCNA)

setwd("F:/Research/Methodology")

Int_data <- readRDS("Intermediate Zone.rds")

top_Int <- prune_taxa(names(sort(taxa_sums(Int_data),TRUE)[1:50]), Int_data)
plot_heatmap(top_Int)

# ?netConstruct
net_single <- netConstruct(Int_data,
                           measure = "sparcc",
                           measurePar = list(iter=20),
                           sparsMethod = "t-test",
                           alpha = 0.05,
                           adjust = "adaptBH",
                           trueNullMethod = "convest",
                           verbose = 3,
                           seed = 123456)

saveRDS(net_single, "Intermediate Zone/Intermediate network.rds")


# ?netAnalyze
props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "degree",
                           weightDeg = FALSE,
                           normDeg = TRUE,
                           verbose = 3)

saveRDS(props_single, "Intermediate Zone/Intermediate analysis.rds")

#?summary.microNetProps
net.summary <- summary(props_single)


# ?plot.microNetProps

p <- plot(props_single,
          shortenLabels = "none",
          # labelLength = 16,
          # charToRm = "g__",
          labelScale = FALSE,
          rmSingles = "all",
          nodeSize = "degree",
          nodeColor = "cluster",
          hubBorderCol = "blue",
          cexNodes = 1,
          cexLabels = 1,
          edgeWidth = 1,
          # negDiffCol = TRUE,
          # posCol = '#009900',
          # negCol = 'red',
          highlightHubs = TRUE,
          cexHubs = 1.5,
          # edgeFilter = "threshold",
          # edgeFilterPar = 0.3,
          # cexHubLabels = 2,
          title1 = "Intermediate Zone Network on genus level with SparCC Method.", 
          showTitle = TRUE,
          cexTitle = 1.5)

saveRDS(p, "Intermediate Zone/Intermediate plot.rds")

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

write.csv(edges, file = "Intermediate Zone/Intermediate edge data.csv", row.names = FALSE)

hubs <- props_single$hubs$hubs1
write(hubs, "Intermediate Zone/Intermediate Zone Hubs.txt")

node.lables <- p$labels$labels1
clust <- props_single$clustering$clust1[node.lables]
node.hubs <- integer(length((node.lables)))
node.hubs[match(hubs,node.lables)] <- 1
degrees <- props_single$centralities$degree1[nodes$lable]
nodes$degree <- degrees

nodes <- data.frame(row.names = 1:length(node.lables))
nodes$index <- 1:length(node.lables)
nodes$lable <- node.lables
nodes$cluster <- clust
nodes$hubs <- node.hubs

write.csv(nodes, file = "Intermediate Zone/Intermediate node data.csv", row.names = FALSE)

hubdata <- subset(nodes, hubs==1, select = -hubs)
write.csv(hubdata, file = "Intermediate Zone/Intermediate hub data.csv", row.names = FALSE)


hubinteractions <- edges[(edges$from %in% hubdata$index) & (edges$to %in% hubdata$index),]

positives <- subset(hubinteractions, association==1, select = -association)
positives$from_cl <- nodes$cluster[positives$from]
positives$to_cl <- nodes$cluster[positives$to]
# positives$from <- node.lables[positives$from]
# positives$to <- node.lables[positives$to]

negatives <- subset(hubinteractions, association== (-1), select = -association)
negatives$from_cl <- nodes$cluster[negatives$from]
negatives$to_cl <- nodes$cluster[negatives$to]
# negatives$from <- node.lables[negatives$from]
# negatives$to <- node.lables[negatives$to]