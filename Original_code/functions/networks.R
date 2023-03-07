#######

library(igraph)

anc <- simulation$anc
species_number <- max(anc[,3], na.rm = T)

pres <- simulation$pres
pres_real <- pres[1:nsteps,1:species_number]
extant_sp <- pres[nsteps, 1:species_number]
to_keep <- which(extant_sp == 1)

#djacency_mats <- test$L_list
#adj_mat <- adjacency_mats[[349]]
#actual_adj_mat <- adj_mat[1:species_number, 1:species_number]
L <- simulation$L

actual_adj_mat_ex <- L[c(to_keep), c(to_keep)]
colnames(actual_adj_mat_ex) = to_keep
rownames(actual_adj_mat_ex) = to_keep

net2 <- graph_from_adjacency_matrix(actual_adj_mat_ex, mode = "undirected", weighted = TRUE)
plot.igraph(net2)

deg <-degree(net2, mode="all")
V(net2)$size <- deg*0.5
V(net2)$label <- NA

plot.igraph(net2, vertex.color="lightgreen")

l <-layout_in_circle(net2)
plot(net2, layout=l, vertex.color="lightgreen")


dist_mat_tree <- get_dist_mat(anc, nsteps, species_number)
test_phylo <- dist_mat_tree[c(to_keep), c(to_keep)]
colnames(test_phylo) = to_keep
rownames(test_phylo) = to_keep

phylo <- hclust(as.dist(test_phylo))
groups_sp <- cutree(phylo, k = 2)

V(net2)$color <- colrs[groups_sp]

layout <- layout_in_circle((net2), order=order(groups_sp))

radian.rescale <- function(x, start=0, direction=1) {
   c.rotate <- function(x) (x + start) %% (2 * pi) * direction
   c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
 }

lab.locs <- radian.rescale(x=1:species_number, direction=-1, start=0)

plot(net2, layout=layout, vertex.label.dist=1, vertex.label.degree=lab.locs)
V(net2)$label <- NA
plot(net2, layout=layout)


colrs <-c("gray50", "tomato", "gold", "lightgreen", "lightblue", "steel blue", "aquamarine")#"yellowgreen", "black", "purple")

plot(hclust(as.dist(test_phylo), method = "average"), hang = - 1, color = colrs[groups_sp])

colrs <-c("gray50", "tomato", "gold", "lightgreen", "lightblue", "steel blue") # "yellowgreen"
V(net2)$color <- colrs[groups_sp] # groups = 6
plot(net2, layout=layout_in_circle)
plot(net2)

inc.edges <-incident_edges(net2,V(net2)[groups_sp==1], mode="all")
ecol <-rep("lightgrey",ecount(net2))
ecol[unlist(inc.edges)] <- "gold"
plot(net2, layout=l, edge.color=ecol)
plot(net2, layout=layout_in_circle, edge.color=ecol)
plot(net2, edge.color=ecol)


##########
# Degree #
##########

# For one simulation
pres <- simulation$pres
degree <- simulation$degree

sp_pres <- pres[68,] == 1
hist(degree[68, sp_pres], breaks = 100, xlim = c(0,30))

sp_pres <- pres[1,] == 1
hist(degree[1, sp_pres], breaks = 100, xlim = c(0,30), ylim = c(0,25), main = paste0("time step =", i))
for(i in 2:nrow(degree)){
	Sys.sleep(0.05)
	sp_pres <- pres[i,] == 1
	hist(degree[i, sp_pres], breaks = 100, xlim = c(0,30), ylim = c(0,25), main = paste0("time step =", i))
}

# Relative degree distribution (relative to the number of species into the community)
sp_pres <- pres[1,] == 1
hist(degree[1, sp_pres], breaks = 10, xlim = c(0,1), ylim = c(0,25), main = paste0("time step =", i))
for(i in 2:nrow(degree)){
	Sys.sleep(0.05)
	sp_pres <- pres[i,] == 1
	hist(degree[i, sp_pres], breaks = 10, xlim = c(0,1), ylim = c(0,25), main = paste0("time step =", i))
}


## Degree distribution over time

library(vioplot)

degree_time_sp <- matrix(NA, nr = nrow(degree), nc = ncol(degree))

for(i in 1:nsteps){
	sp_pres <- pres[i,] == 1
	degree_time_sp[i, sp_pres] <- degree[i, sp_pres]
}

cdat <- as.list(as.data.frame(t(as.matrix(degree_time_sp))))
vioplot(cdat, col="lightgreen")
boxplot(cdat, col="lightgreen")



## Modularity

library(igraph)
pres <- simulation$pres
sp_pres <- pres[nsteps,] == 1

adj_mat <- simulation$L[sp_pres, sp_pres]
network <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE)

test_mod <- modularity(network)

par(mfrow = c(1,1))
ebc <- edge.betweenness.community(network)
membership(ceb)
ceb <- cluster_edge_betweenness(network) #Community detection based on edge betweenness (Newman-Girvan)
dendPlot(ceb, mode="hclust")
plot(ceb, network)
modularity(ceb)

cfg <- cluster_fast_greedy(as.undirected(network)) #fg <- fastgreedy.community(network, merges=TRUE, modularity=TRUE)
plot(cfg, as.undirected(network))
modularity(cfg)

spg <- spinglass.community(network)
plot(spg, network)
modularity(spg)

wtc <- walktrap.community(network)
plot(wtc, network)
modularity(wtc)

coords = layout_with_fr(net2)

par(mfrow = c(1,2))
memb <- membership(cfg)
plot(net2, vertex.color=memb, vertex.label = NA, vertex.size = 4, layout = coords)

memb2 = groups_sp
plot(net2, vertex.color=memb2, vertex.label = NA, vertex.size = 4, layout = coords)

plot(net2, vertex.color=colrs[memb], vertex.label = NA, vertex.size = 4, layout = layout_in_circle((net2), order = order(memb)))
plot(net2, vertex.color=colrs[memb2], vertex.label = NA, vertex.size = 4, layout = layout_in_circle((net2), order = order(memb)))

# Range of interactions
range_interaction <- matrix(traits[,2], nr = nsteps, nc = ncol(pres), byrow=TRUE)
range_pres <- pres*range_interaction
to_keep <- which(range_pres[1,] != 0)
hist(2*range_pres[1, to_keep], breaks = 10, xlim = c(0,1), ylim = c(0,25), main = paste0("time step =", i))
for(i in 2:nrow(range_pres)){
	Sys.sleep(0.05)
	to_keep <- which(range_pres[i,] != 0)
	hist(2*range_pres[i, to_keep], breaks = 10, xlim = c(0,1), ylim = c(0,25), main = paste0("time step =", i))
}
