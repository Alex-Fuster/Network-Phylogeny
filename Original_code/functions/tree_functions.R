######### Function to build a phylogenetic tree from the anc table #############

# Function to append a list into an existing list
list_append <- function(lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

# Function to find the cell that are going to be filled
cell_position <- function(n_row, x, y) n_row * (x - 1) + y

# Function to assign "time-distance" between species
fill_distance <- function(new_group, dist_mat, anc_rev, nsteps, loop_position){
  # Create all possible pairwise combinations between the species of the group
  pairs <- expand.grid(new_group, new_group)
  sp1 <- pairs[,1]
  sp2 <- pairs[,2]
  to_fill <- cell_position(nrow(dist_mat), sp1, sp2)
  to_test <- dist_mat[to_fill] == 0
  dist_mat[to_fill[to_test]] = nsteps - anc_rev[loop_position,1]
  return(dist_mat)
}

# Function to fill the time-distance matrix between all species of the system
group_distance <- function(anc_rev, groups, nsteps, dist_mat){
  for(i in 1:nrow(anc_rev)){
    ancestor <- anc_rev[i,2]
    new_species <- anc_rev[i,3]

    if(i==1){
      # Create a new group of realated species and add it to the existing list of group
      new_group <- c(ancestor, new_species)
      groups <- list_append(groups, new_group)
      # Record the time-distance of the related species
      dist_mat <- fill_distance(new_group, dist_mat, anc_rev, nsteps, i)
    }

    if(i>1){
      # See if any of the two species (ancestor or new species) are already involved in a group
      test_anc <- which(sapply(groups, function(e) is.element(ancestor, e)))
      test_new_spe <- which(anc_rev[1:i-1, 2] == new_species)

      if(length(test_anc)<1 && length(test_new_spe)<1){
        # The two species are note yet related to other species
        # Create a new group of realated species and add it to the existing list of group
        new_group <- c(ancestor, new_species)
        groups <- list_append(groups, new_group)
        # Record the time-distance of the related species
        dist_mat <- fill_distance(new_group, dist_mat, anc_rev, nsteps, i)
      }
      else if(length(test_anc)==1 && length(test_new_spe)<1){
        # The ancestor is already related to one group
        # Create a new group of realated species and add it to the existing list of group
        new_group <- c(new_species, unlist(groups[test_anc]))
        groups <- list_append(groups, new_group)
        # Record the time-distance of the related species
        dist_mat <- fill_distance(new_group, dist_mat, anc_rev, nsteps, i)
      }
      else if(length(test_anc)>1 && length(test_new_spe)<1){
        # The ancestor is already related to several groups
        # Create a new group of realated species and add it to the existing list of group
        new_group <- unique(c(new_species, unlist(groups[test_anc])))
        groups <- list_append(groups, new_group)
        # Record the time-distance of the related species
        dist_mat <- fill_distance(new_group, dist_mat, anc_rev, nsteps, i)
      }
      else if(length(test_anc)<1 && length(test_new_spe)==1){
        # The new species is already related to an other group
        # Create a new group of realated species and add it to the existing list of group
        new_group <- c(ancestor, unlist(groups[test_new_spe]))
        groups <- list_append(groups, new_group)
        # Record the time-distance of the related species
        dist_mat <- fill_distance(new_group, dist_mat, anc_rev, nsteps, i)
      }
      else if(length(test_anc)==1 && length(test_new_spe)==1){
        # The ancestor AND the new species are already related to one other group (one each)
        # Create a new group of realated species and add it to the existing list of group
        new_group <- c(unlist(groups[test_new_spe]), unlist(groups[test_anc]))
        groups <- list_append(groups, new_group)
        # Record the time-distance of the related species
        dist_mat <- fill_distance(new_group, dist_mat, anc_rev, nsteps, i)
      }
      else if(length(test_anc)>1 && length(test_new_spe)==1){
        # The ancestor is realted to several other groups, the new species is related to one other group
        # Create a new group of realated species and add it to the existing list of group
        new_group <- unique(c(unlist(groups[test_new_spe]), unlist(groups[test_anc])))
        groups <- list_append(groups, new_group)
        # Record the time-distance of the related species
        dist_mat <- fill_distance(new_group, dist_mat, anc_rev, nsteps, i)
      }
      else if(length(test_anc)<1 && length(test_new_spe)>1){
        # The new species is already related to several other groups
        # Create a new group of realated species and add it to the existing list of group
        new_group <- unique(c(ancestor, unlist(groups[test_new_spe])))
        groups <- list_append(groups, new_group)
        # Record the time-distance of the related species
        dist_mat <- fill_distance(new_group, dist_mat, anc_rev, nsteps, i)
      }
      else if(length(test_anc)==1 && length(test_new_spe)>1){
        # The ancestor is related to one other group and the new species to several other groups
        # Create a new group of realated species and add it to the existing list of group
        new_group <- unique(c(unlist(groups[test_new_spe]), unlist(groups[test_anc])))
        groups <- list_append(groups, new_group)
        # Record the time-distance of the related species
        dist_mat <- fill_distance(new_group, dist_mat, anc_rev, nsteps, i)
      }
      else if(length(test_anc)>1 && length(test_new_spe)>1){
        # The ancestor and the new species are both related to several groups
        # Create a new group of realated species and add it to the existing list of group
        new_group <- unique(c(unlist(groups[test_new_spe]), unlist(groups[test_anc])))
        groups <- list_append(groups, new_group)
        # Record the time-distance of the related species
        dist_mat <- fill_distance(new_group, dist_mat, anc_rev, nsteps, i)
      }
    }
  }
  return(dist_mat)
}

# Function to get the final time-distance matrix between species and build an approximative phylogenetic tree (branch length are not correct)
get_dist_mat <- function(anc, nsteps, species_number){
  anc_df <- anc[1:species_number,]
  # Put the matrix upside-down!
  anc_rev <- apply(anc_df, 2, rev)

  dist_mat <- matrix(0, species_number, species_number)
  groups <- list()

  end_sim <- anc[species_number,1]

  filled_dist_mat <- group_distance(anc_rev, groups, end_sim, dist_mat)

  #plot(hclust(as.dist(filled_dist_mat)), method = "average")

  return(filled_dist_mat)
}










#anc <- read.csv("ancestor.txt", header = TRUE, sep = '\t')
#nsteps <- 250
#species_number <- max(anc[,1])
#
#test_tree <- draw_the_tree(anc, nsteps, species_number)

#library(profvis)

#profvis({
#  dist_mat_test2 <- get_dist_mat_vec(anc, nsteps, species_number)
#})
