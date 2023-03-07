# Function to append a list into an existing list
list_append <- function(lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

## Function to find the cell that are going to be filled
#cell_position <- function(n_row, x, y) n_row * (x - 1) + y

# Function to regroup species into groups of the most recent common ancestor (groups of 2 or 3)
groups_creation <- function(anc, groups){
  first_groups <- list()
  #Starts at i = 2 because the first row is always NA
  for(i in 2:nrow(anc)){
    new_group <- c(anc[i,2], anc[i,3])
    first_groups <- list_append(first_groups, new_group)
  }

  #Starts at 2 because we use des i-1
  for(i in 2:length(first_groups)){
    if(first_groups[[i]][1] == first_groups[[i-1]][1]){
      new_group <- unique(c(unlist(first_groups[[i]]), unlist(first_groups[[i-1]])))
      groups <- list_append(groups, new_group)
    }

    if(first_groups[[i]][1] != first_groups[[i-1]][1]){
      if(i == length(first_groups)){
        new_group <- first_groups[[i]]
        groups <- list_append(groups, new_group)
      } else if(first_groups[[i]][1] != first_groups[[i+1]][1]){
        if(first_groups[[i-1]][1] == 1){
          new_group <- first_groups[[i-1]]
          groups <- list_append(groups, new_group)
          new_group <- first_groups[[i]]
          groups <- list_append(groups, new_group)
        } else {
          new_group <- first_groups[[i]]
          groups <- list_append(groups, new_group)
        }
      } else if(first_groups[[i]][1] == first_groups[[i+1]][1]){
        if(first_groups[[i-1]][1] == 1){
          new_group <- first_groups[[i-1]]
          groups <- list_append(groups, new_group)
        }
      }
    }
  }
  return(groups)
}

# Function to count the number of nodes separating two species.
## The tree si divided into two main sides, species 1 is the root and species 2 enhance a first side where species 3 enhance the other one.
## If the two species are from the same side, the count end when they are lineary directly related
## If the two species are not from the same side, we count every group untill we reach the group (1,2,3).
## The number of nodes is calculated as number of node from species 1 + number of nodes from species 2 - 1
## Data should be a two colomns data frame with all possible species combinations
node_count <- function(data, groups){
  sp1 <- data[1]
  sp2 <- data[2]
  # Ordered species from the youngest to the oldest
  youngest <- max(sp1, sp2)
  oldest <- min(sp1,sp2)

  # Assomption the species are not from the same side of the tree
  same_side <- FALSE
  node_distance <- 999

  # Potential common ancestor
  pot_anc1 <- NULL

  # Start counting with the youngest species
  sp_cible <- youngest

  if(youngest == oldest & youngest == 1){
    node_distance <- 0
  } else {
    while(sp_cible[1]!= 1){
      if(youngest == oldest){
        node_distance <- 0
        break
      }

      groups_sp1 <- grep(TRUE, lapply(groups, '%in%', sp_cible))
      group_cible <- groups[[groups_sp1[1]]]

      # Stop if species are from the same side of the tree
      if(oldest %in% group_cible & group_cible[1] != 1){
        pot_anc1 <- cbind(pot_anc1, group_cible[1])
        same_side <- TRUE
        break
      }

      # Continue if we have to go through species 1 to reach the other species
      pot_anc1 <- cbind(pot_anc1, group_cible[1])
      sp_cible <- group_cible[1]
    }
  }


  if(same_side == TRUE){
    node_distance <- ncol(pot_anc1)
  } else if(youngest == oldest & youngest == 1){
    node_distance <- 0
  } else {
    sp_cible <- oldest

    # Potential common ancestor
    pot_anc2 <- NULL

    if(sp_cible[1] == 1){
      node_distance <- ncol(pot_anc1)
    }

    while(sp_cible[1]!= 1){
      groups_sp1 <- grep(TRUE, lapply(groups, '%in%', sp_cible))
      group_cible <- groups[[groups_sp1[1]]]

      if(group_cible[1] %in% pot_anc1 & group_cible[1] != 1){
        pot_anc1 <- pot_anc1[1:which(pot_anc1 == group_cible[1])]
        node_distance <- length(pot_anc1) + ncol(pot_anc2)
        break
      } else if(group_cible[1] %in% pot_anc1 & group_cible[1] == 1){
        pot_anc1 <- pot_anc1[1:which(pot_anc1 == group_cible[1])]
        pot_anc2 <- cbind(pot_anc2, group_cible[1])
        node_distance <- length(pot_anc1) + ncol(pot_anc2) - 1
      break
      } else {
        pot_anc2 <- cbind(pot_anc2, group_cible[1])
        sp_cible <- group_cible[1]
      }
    }
  }
  return(node_distance)
}

make_dist_mat <- function(time, pres, groups){
  # Sélectionner les espèces présentent à un temps t, ici, pour le test, on va dire t = 100
  #pres <- list_res[[1]]$presence_matrix
  esp_to_test <- grep(1, pres[time,])
  #esp_to_test

  if(length(esp_to_test) == 1){
    distance_matrix <- 0
  } else {
    # Préparer les paires d'especes pour les distances
    n_sp <- length(esp_to_test)

    #library(gtools)
    sp_pairs <- combinations(n_sp, 2, 1:n_sp)

    for(i in 1:n_sp){
      sp_pairs[sp_pairs == i] <- esp_to_test[i]
    }

    # calculer la distance entre les paires d'especes
    node_distance <- NULL

    for(i in 1:nrow(sp_pairs)){
      esp_distance <- node_count(c(sp_pairs[i, 1], sp_pairs[i, 2]), groups)
      node_distance <- rbind(node_distance, esp_distance[1])
    }

    # Ranger les distances sous forme de matrix
    distance_matrix <- matrix(0, n_sp, n_sp)

    position_into_matrix <- combinations(n_sp, 2, 1:n_sp)

    for(i in 1:nrow(position_into_matrix)){
      distance_matrix[position_into_matrix[i,2],position_into_matrix[i,1]] = node_distance[i]
    }

    distance_matrix <- as.dist(distance_matrix)

  }

  return(distance_matrix)
}

##### TEST ####
#library(gtools)
#
#load(file = "Data/Positive_interactions/list_res_pos.Rdata")
#anc <- list_res[[1]]$parentage_matrix
#pres <- list_res[[1]]$presence_matrix
#
##distance_matrix <- matrix(0, nrow = nrow(anc), ncol = nrow(anc))
##diag(distance_matrix) = 1
#
#groups <- list()
#groups <- groups_creation(anc, groups)
#
#time <- 30
#
#make_dist_mat(time, pres, groups)
