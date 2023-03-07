library("treespace")
library("adegenet")
library("adegraphics")
library("rgl")

library(viridis)
library(wesanderson)

class(tree_space_an) <- "multiPhylo"
class(tree_space_an)
tree_space_an

write.tree(tree_space_an, file = "Data/Positive_interactions/pos_trees.Rdata")
pos_trees <- read.tree("Data/Positive_interactions/pos_trees.Rdata")

write.tree(tree_space_an, file = "Data/Negative_interactions/neg_trees.Rdata")
neg_trees <- read.tree("Data/Negative_interactions/neg_trees.Rdata")

test_treespace <- treespace(tree_space_an, nf = 3, lambda = 1)
test_treespace <- treespace(tree_space_an, nf = 3, lambda = 0)
test_treespace <- treespace(tree_space_an, method = "nNodes", nf = 3, lambda = 0)

names(tree_space_an) <- c(1:nsim)
plotGrovesD3(test_treespace$pco, treeNames = names(tree_space_an))#, color.palette=viridis)
table.image(test_treespace$D, nclass=30)

table.value(test_treespace$D, nclass=5, method="color", symbol="circle", col=redpal(6))
plotGroves(test_treespace$pco, lab.show=TRUE, lab.cex=1.5)

#faire des groupes
treespace_group <- findGroves(test_treespace, nclust=6)
plotGrovesD3(treespace_group)
# alternative with improved legend and tooltip text, giving the tree numbers:
plotGrovesD3(treespace_group, tooltip_text=paste0("Tree ",1:length(tree_space_an)), legend_width=50, col_lab="Cluster")

# prepare a colour palette:
colours <- fac2col(treespace_group$groups, col.pal=viridis)
colours <- fac2col(treespace_group$groups, col.pal=funky)


#plot en 3D
plot3d(treespace_group$treespace$pco$li[,1],
       treespace_group$treespace$pco$li[,2],
       treespace_group$treespace$pco$li[,3],
       col=colours, type="s", size=1.5,
       xlab="", ylab="", zlab="")

pos_trees <- read.tree("Data/Positive_interactions/pos_trees.Rdata")
neg_trees <- read.tree("Data/Negative_interactions/neg_trees.Rdata")
test_tree_pn <- c(pos_trees, neg_trees)

test_total_treespace <- treespace(test_tree_pn, nf = 3, lambda = 1)
test_total_treespace <- treespace(test_tree_pn, method = "nNodes", nf = 3, lambda = 0)
plotGroves(test_total_treespace$pco, lab.show=TRUE, lab.cex=1.5)
plotGrovesD3(test_total_treespace$pco, colors = c(rgb(77,146,33, max = 255), rgb(197,27,125, max = 255)))

treespace_group <- findGroves(test_total_treespace, nclust=2)
treespace_group$groups <- c(rep("positive",100), rep("negative",100))

plotGrovesD3(treespace_group, colors = c(rgb(77,146,33, max = 255), rgb(197,27,125, max = 255)), col_lab="Interactions")

#plot en 3D
plot3d(treespace_group$treespace$pco$li[,1],
       treespace_group$treespace$pco$li[,2],
       treespace_group$treespace$pco$li[,3],
       col=c(rgb(77,146,33, max = 255), rgb(197,27,125, max = 255)), type="s", size=1.5,
       xlab="x axis", ylab="y axis", zlab="z axis")


test_pcoa1 <- dudi.pco(d=test_total_treespace$D,scannf=is.null(3),nf=3)
test_pcoa2 <- pcoa(test_total_treespace$D)

biplot(test_pcoa1)
x11()
biplot(test_pcoa2)
