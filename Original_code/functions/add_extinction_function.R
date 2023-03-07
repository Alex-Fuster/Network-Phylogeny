## Function add information about species extinction into tree data changing the
##length of le last branch of each species depending on its time to extinction

add_extinction <- function(distance_matrix, extinction, anc){
  # Which species are extinct?
  extinct_sp <- which(!is.na(extinction[,2]))

  tree <- ape::as.phylo(hclust(as.dist(distance_matrix)))

  #How have to give back the real branches length that have been divided by two
  tree$edge.length = tree$edge.length*2

  for(i in 1:length(extinct_sp)){
    # Check if the species that extinct gave birth to (an) other species
    ancestor <- which(anc[,2] == extinct_sp[i])

    # Gives you the last line in anc where the species appear, it would be the species last node, which means that we'll modify the branch related to this pĥylogenetic node and the species
    if(length(ancestor) < 1){
      # The species never gave a new species
      last_speciation <- max(which(anc[,3] == extinct_sp[i]))
    }else if(length(ancestor) > 0){
      last_speciation <- max(ancestor)
    }

    # Gives you the phylogenetic node related to the species
    pos_1 <- which(tree$edge[,2] == extinct_sp[i])

    # Modification of the last branch related to the species
    tree$edge.length[pos_1] = extinction[extinct_sp[i],1] - anc[last_speciation,1]
  }
  return(tree)
}


#### Test du package TreeSpace ###
#
##install.packages("treespace")
#library("treespace")
## Library nécessaires pour les graphiques
#library("adegenet")
#library("adegraphics")
#library("rgl")
#
## generate list of trees
#suppressWarnings(RNGversion("3.5.0"))
#set.seed(1)
#x <- rmtree(10, 20)
#names(x) <- paste("tree", 1:10, sep = "")
#
## use treespace
#res <- treespace(x, nf=3)
#names(res)
#res
#
## table.image
#table.image(res$D, nclass=30)
#
## table.value with some customization
#table.value(res$D, nclass=5, method="color",
#            symbol="circle", col=redpal(5))
#
#plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5)
#
#plotGrovesD3(res$pco, treeNames=1:10)
#
#data(woodmiceTrees)
#wm.res <- treespace(woodmiceTrees,nf=3)
#
## PCs are stored in:
#head(wm.res$pco$li)
## plot results
#plotGrovesD3(wm.res$pco)
#
