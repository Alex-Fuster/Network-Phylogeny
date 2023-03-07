




# Compute phylogenetic distances for each timestep from the anc matrix


compute_phylo_dist_all_steps <- function(df_anc, timesteps, initial_timestep) {
  
  
  
  phylo_dist_list <- list()
  
  
  for (i in initial_timestep:timesteps) {
    
    
    n_spp <- max(df_anc[1:i,3], na.rm = T)
    
    dist_mat_tree <- get_dist_mat(anc = df_anc, nsteps = i, species_number = n_spp)
    
    phylo_dist_list[[i]] <- dist_mat_tree
    
    timestep_anc.matrix <- df_anc[i,1]
    
    names(phylo_dist_list)[[i]] <- paste("mat_timestep", timestep_anc.matrix)
    
    
  }
  
  
  return(phylo_dist_list)
  
}










# Compute phylogenetic trees for each timestep from the anc matrix


compute_phylog_all_timesteps <- function(df_anc, timesteps, initial_timestep) {
  
  
  
  Phylo_list <- list()
  
  
  for (i in initial_timestep:timesteps) {
    
    
    n_spp <- max(df_anc[1:i,3], na.rm = T)
    
    dist_mat_tree <- get_dist_mat(anc = df_anc, nsteps = i, species_number = n_spp)
    
    Phylo_list[[i]] <- as.phylo(hclust(as.dist(dist_mat_tree), method = "average"))
    

    
  }
  
  
  return(Phylo_list)
  
}


# Same but without hclust

compute_phylog_all_timesteps_nohclust <- function(df_anc, timesteps, initial_timestep) {
  
  
  
  Phylo_list <- list()
  
  
  for (i in initial_timestep:timesteps) {
    
    
    n_spp <- max(df_anc[1:i,3], na.rm = T)
    
    dist_mat_tree <- get_dist_mat(anc = df_anc, nsteps = i, species_number = n_spp)
    
    Phylo_list[[i]] <- as.dist(dist_mat_tree)
    
    Phylo_list[[i]]
    
    
  }
  
  
  return(Phylo_list)
  
}




# Plot trees of the trees list obtained in the simulation

plot_trees <- function(trees_list, plot_type) { # cladogram, phylogram, radial, fan
  
  list_plots <- list()
  
  for (i in 3:length(trees_list)) {
    
    list_plots[i] <- plot(unroot(trees_list[[i]]),type= plot_type,cex=0.6,
                         use.edge.length=FALSE,lab4ut="axial",
                         no.margin=TRUE)
    
    
  }
  

  
  return(list_plots)
  
}





# set column and row names as "spp_n", n being their position

set_sppNames_icolrows <-  function(matrix) {
  
  mat <- matrix
  rownames(mat) <- rownames(mat, do.NULL = FALSE, prefix = "spp_")
  colnames(mat) <- colnames(mat, do.NULL = FALSE, prefix = "spp_")
  
  return(mat)
  
}


# compute dissimilarity matrix computting jaccard index between interacting species

compute_jaccard <- function(matrix) {
  
  mat <- as.matrix(vegdist(matrix[-which(rowSums(matrix) == 0),], method = "jac", binary = T))
  
  return(mat)
  
}

