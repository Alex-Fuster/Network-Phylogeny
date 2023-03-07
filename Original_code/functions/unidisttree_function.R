unidisttree <- function(distance_matirx){
  tree <- ape::as.phylo(hclust(as.dist(distance_matrix)))
  uni_edge_length <- rep(1, length(tree$edge.length))
  tree$edge.length = 1
}

## FONCTIONNE PAS ! CAR LES TIPS NE CORRESPONDENT PAS AUX ESPECES, DONC LA MATRICE DE DISTANCES OBTENUE N'A RIEN à VOIR DES LA MATRICE DE DISTANCE DES TRAITS
