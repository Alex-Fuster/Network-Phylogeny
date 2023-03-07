
anc <- readRDS("anc.rds") #AF

library(ape)
library(apTreeshape)

source("tree_functions.R")
source("phylo_functions.R")

#anc <- read.csv("ancestor.txt", header = TRUE, sep = '\t') # les colonnes ne sont pas dans le bon ordre

species_number <- max(anc[,3], na.rm = T)

dist_mat_tree <- get_dist_mat(anc, nsteps, species_number)
#plot(hclust(as.dist(dist_mat_tree), method = "average"), hang = -1)
#plot(hclust(as.dist(dist_mat_tree), method = "average"), label = FALSE)

tree_gendist <- as.phylo(hclust(as.dist(dist_mat_tree), method = "average"))
# Quelque chose à faire ici pour
ttgen <- timetree(tree_gendist)
gamma_stat <- gammaStat((ttgen))

ttgen_s <- as.treeshape(ttgen)
sackin_index <- sackin(ttgen_s)
sackin_norm_index <- sackin(ttgen_s, norm = "yule")

#tree.gendist=as.phylo(hclust(as.dist(gendist), method="average"))
#ttgen <-timetree(tree.gendist)
#gamma.values[i,1]=gammaStat((ttgen))
#ttgen.s <- as.treeshape(ttgen)
#sackin.values[i,1]=sackin(ttgen.s)
#sackin.values[i,2]=sackin(ttgen.s,norm="yule")


hclust = hclust(as.dist(dist_mat_tree))
# Check dendrogram looks sensible
plot(hclust)
class(hclust) # check that class is hclust
# Save to Newick file
my_tree <- as.phylo(hclust)
#write.tree(phy=my_tree, file="ExampleTree.newick") # Writes a Newick file


#### Until here  I think is what I want - the tree ###

# I could ignore the code below









# Produce tree
plot(unroot(my_tree),type="unrooted",cex=1.5,label = FALSE,
     use.edge.length=TRUE,lab4ut="axial",
     edge.width=2,
     no.margin=TRUE)

plot(unroot(my_tree),type="unrooted",cex=0.6,
    use.edge.length=FALSE,lab4ut="axial",
    no.margin=TRUE)

plotTree(my_tree,type="fan",fsize=0.7,lwd=1,
    ftype="i")

plot(my_tree,type="cladogram")

# Plot phylogenetic tree of existing species at the end of the simulation
pres <- test$pres
extant_sp <- pres[nsteps, 1:species_number]
to_keep <- which(extant_sp == 1)
test_phylo <- dist_mat_tree[c(to_keep), c(to_keep)]
colnames(test_phylo) = to_keep
rownames(test_phylo) = to_keep
plot(hclust(as.dist(test_phylo), method = "average"), hang = - 1)


plot(unroot(as.phylo(hclust(as.dist(test_phylo)))),type="cladogram",cex=0.6,
    use.edge.length=FALSE,lab4ut="axial",
    no.margin=TRUE)

plot(unroot(as.phylo(hclust(as.dist(test_phylo)))),type="phylogram",cex=0.6,
    use.edge.length=FALSE,lab4ut="axial",
    no.margin=TRUE)

plot(unroot(as.phylo(hclust(as.dist(test_phylo)))),type="radial",cex=0.6,
    use.edge.length=FALSE,lab4ut="axial",
    no.margin=TRUE)

plot(unroot(as.phylo(hclust(as.dist(test_phylo)))),type="fan",cex=0.6,
    use.edge.length=FALSE,lab4ut="axial",
    no.margin=TRUE)


tree_gendist <- as.phylo(hclust(as.dist(test_phylo), method = "average"))
# Quelque chose à faire ici pour
ttgen <- timetree(tree_gendist)
gamma_stat <- gammaStat((ttgen))

ttgen_s <- as.treeshape(ttgen)
sackin_index <- sackin(ttgen_s)
sackin_norm_index <- sackin(ttgen_s, norm = "yule")

# dendrogram with colored species by group
library(dendextend)

phylo <- hclust(as.dist(test_phylo))

dend <- as.dendrogram(phylo)
dend <- dend %>%
   color_branches(k=5) %>%
   color_labels
plot(dend)

colrs <-c("gray50", "tomato", "gold", "lightgreen", "lightblue", "steel blue", "aquamarine", "yellowgreen", "black", "purple")
groups_sp <- cutree(phylo, k = 5)

dend <- as.dendrogram(phylo)
labels_colors(dend) <- colrs[groups_sp][order.dendrogram(dend)]
plot(dend)
