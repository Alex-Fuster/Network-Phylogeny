




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

set_sppNames_icolrows_spp <-  function(matrix) {
  
  mat <- matrix
  rownames(mat) <- rownames(mat, do.NULL = FALSE, prefix = "spp_")
  colnames(mat) <- colnames(mat, do.NULL = FALSE, prefix = "spp_")
  
  return(mat)
  
}






set_sppNames_numbers <- function(mat) {
  
  rownames(mat) <- seq(1:1000)
  colnames(mat) <- seq(1:1000)
  
  return(mat)  
}





set_sppNames_icolrows <-  function(matrix) {
  
  mat <- matrix
  rownames(mat) <- rownames(mat, do.NULL = FALSE, prefix = "")
  colnames(mat) <- colnames(mat, do.NULL = FALSE, prefix = "")
  
  return(mat)
  
}


convert_sppnames_toletters <-  function(mat) {
  
  
  rownames(mat) <- chartr("0123456789", "ABCDEFGHIJ", rownames(mat))
  colnames(mat) <- chartr("0123456789", "ABCDEFGHIJ", rownames(mat))
  
  return(mat)
  
}




change_sppnames_letters_ancdist.table <- function(mat) {
  
  mat[,"ancestor"] <- chartr("0123456789", "ABCDEFGHIJ", mat[,"ancestor"])
  
  mat[,"spp"] <- chartr("0123456789", "ABCDEFGHIJ", mat[,"spp"])
  
  return(mat)
  
}




# compute dissimilarity matrix computting jaccard index between interacting species

compute_jaccard <- function(matrix, distance) {
  
  mat <- as.matrix(vegdist(matrix[-which(rowSums(matrix) == 0),], method = "jac", binary = T))
  
  return(mat)
  
}


compute_chisq <- function(matrix, distance) {
  
  mat <- as.matrix(vegdist(matrix[-which(rowSums(matrix) == 0),], method = "chisq", binary = T))
  
  return(mat)
  
}



# Obtain tree from anc_table 




sister.group<-function(sis.names,sis.dist){
  n<-length(sis.dist)
  distances<-paste0(rep(":",n),sis.dist)
  sis<-paste0(sis.names,distances)
  res<-paste0(sis,collapse=",")
  res<-paste0("(",res,")")
  res
}

ToPhylo<-function(data){
  data.2<-data
  data.2$repr<-data$spp
  sisters<-levels(as.factor(data$spp))
  mothers<-levels(as.factor(data$ancestor))
  tips<-setdiff(sisters,mothers)
  root<-setdiff(mothers,sisters)
  foc.nodes<-unique(data[which(data$spp%in%tips),"ancestor"])
  n<-length(foc.nodes)
  data.2$repr[data.2$spp%in%tips]<-data.2$repr[data.2$spp%in%tips]
  while(n>1){
    foc.nodes2<-unique(data.2[which(data.2$spp%in%foc.nodes),"ancestor"])
    for(i in 1:n){
      daughters<-data.2[which(data.2$ancestor==foc.nodes[i]),"repr"]
      #print(daughters)
      daughters.dist<-data.2[which(data.2$ancestor==foc.nodes[i]),"distance"]
      data.2$repr[data.2$spp==foc.nodes[i]]<-paste0(sister.group(daughters,daughters.dist),foc.nodes[i])
    }
    tips<-foc.nodes
    foc.nodes<-foc.nodes2
    n<-length(foc.nodes)
  }
  daughters<-data.2[which(data.2$ancestor==foc.nodes[1]),"repr"]
  #print(daughters)
  daughters.dist<-data.2[which(data.2$ancestor==foc.nodes[1]),"distance"]
  paste0(sister.group(daughters,daughters.dist),root)
}



## Order alhpabetically rows and columns of dataframe

order_col.row_names <- function(df) {
  
  df_ordered <- df[order(names(df)) , order(names(df))]
  
  return(df_ordered)
  
}




















