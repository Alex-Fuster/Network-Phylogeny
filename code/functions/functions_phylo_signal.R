




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




# set colnames and rownames as numbers


set_sppNames_numbers <- function(mat) {
  
  rownames(mat) <- seq(1:nrow(mat))
  colnames(mat) <- seq(1:ncol(mat))
  
  return(mat)  
}





#set_sppNames_icolrows <-  function(matrix) {
  
 # mat <- matrix
 # rownames(mat) <- rownames(mat, do.NULL = FALSE, prefix = "")
 # colnames(mat) <- colnames(mat, do.NULL = FALSE, prefix = "")
  
#  return(mat)
  
#}



# convert colnames and rownames from numbers to letters

convert_sppnames_toletters <-  function(mat) {
  
  
  rownames(mat) <- chartr("0123456789", "ABCDEFGHIJ", rownames(mat))
  colnames(mat) <- chartr("0123456789", "ABCDEFGHIJ", colnames(mat))
  
  return(mat)
  
}



# convert species names of ancestry table from numbers to letters

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


compute_chisq_pred <- function(matrix, distance) {
  
  #Eliminate those species with no interactions (none as predators & none as preys)
  
  if(length(which(rowSums(matrix)==0 & colSums(matrix)==0)) > 0){
    
    mat <- matrix[-which(rowSums(matrix)==0 & colSums(matrix) ==0),
    -which(rowSums(matrix)==0 & colSums(matrix)==0)]

#compute dista ce between spp(columns) as predators

mat_dist <- as.matrix(vegdist(mat, method = distance, binary = T))
    
  }else{
    mat_dist <- as.matrix(vegdist(matrix, method = distance, binary = T))
  }
  
  return(mat_dist)
}


compute_chisq_prey <- function(matrix, distance) {
  
  #Eliminate those species with no interactions (none as predators & none as preys)
  
  if(length(which(rowSums(matrix)==0 & colSums(matrix)==0)) > 0){
    
    mat <- matrix[-which(rowSums(matrix)==0 & colSums(matrix) ==0),
                  -which(rowSums(matrix)==0 & colSums(matrix)==0)]
    
    #compute dista ce between spp(columns) as predators
    
    mat_dist <- as.matrix(vegdist(t(mat), method = distance, binary = T))
    
  }else{
    mat_dist <- as.matrix(vegdist(t(matrix), method = distance, binary = T))
  }
  
  return(mat_dist)
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




#####################
#Compute NMI
####################


## NMI is for similarity. This function makes it a distance metric (1 - MNI)

nmi_to_distance <- function(value) {
  
  result <- 1 - value
  
  return(result)
  
}


compute_nmi_igraph_pred <- function(matrix) {
  
  
  #  Vector source for column combinations
  n <- seq_len( ncol(matrix) )
  
  #  Make combinations
  id <- expand.grid( n , n )
  
  #  Get result
  matrix_nmi <- matrix(0, ncol = ncol(matrix), nrow = nrow(matrix))
  
  for (i in 1:nrow(id)) {
    
    matrix_nmi[id[i,1], id[i,2]] <- compare(
      matrix[,id[i,1]],
      matrix[,id[i,2]],
      method = "nmi")
    
    matrix_nmi[id[i,2], id[i,1]] <- compare(
      matrix[,id[i,1]],
      matrix[,id[i,2]],
      method = "nmi")
    
  }
  
  
  colnames(matrix_nmi) <- colnames(matrix)
  rownames(matrix_nmi) <- rownames(matrix)
  
  matrix_nmi_dist <- nmi_to_distance(matrix_nmi)
  
  return(matrix_nmi_dist)
  
  
}



compute_nmi_igraph_prey <- function(matrix) {
  
  
  #  Vector source for column combinations
  n <- seq_len( ncol(matrix) )
  
  #  Make combinations
  id <- expand.grid( n , n )
  
  #  Get result
  matrix_nmi <- matrix(0, ncol = ncol(matrix), nrow = nrow(matrix))
  
  for (i in 1:nrow(id)) {
    
    matrix_nmi[id[i,1], id[i,2]] <- compare(
      matrix[id[i,1],],
      matrix[id[i,2],],
      method = "nmi")
    
    matrix_nmi[id[i,2], id[i,1]] <- compare(
      matrix[id[i,1],],
      matrix[id[i,2],],
      method = "nmi")
    
  }
  
  
  colnames(matrix_nmi) <- colnames(matrix)
  rownames(matrix_nmi) <- rownames(matrix)
  
  matrix_nmi_dist <- nmi_to_distance(matrix_nmi)
  
  return(matrix_nmi_dist)
  
  
}






compute_nmi_aricode_pred <- function(matrix) {
  
  
  #  Vector source for column combinations
  n <- seq_len( ncol(matrix) )
  
  #  Make combinations
  id <- expand.grid( n , n )
  
  #  Get result
  matrix_nmi <- matrix(0, ncol = ncol(matrix), nrow = nrow(matrix))
  
  for (i in 1:nrow(id)) {
    
    matrix_nmi[id[i,1], id[i,2]] <- NMI(
      matrix[,id[i,1]],
      matrix[,id[i,2]])
    
    matrix_nmi[id[i,2], id[i,1]] <- NMI(
      matrix[,id[i,1]],
      matrix[,id[i,2]])
    
  }
  
  
  colnames(matrix_nmi) <- colnames(matrix)
  rownames(matrix_nmi) <- rownames(matrix)
  
  matrix_nmi_dist <- nmi_to_distance(matrix_nmi)
  
  return(matrix_nmi_dist)
  
  
}




compute_nmi_aricode_prey <- function(matrix) {
  
  
  #  Vector source for column combinations
  n <- seq_len( ncol(matrix) )
  
  #  Make combinations
  id <- expand.grid( n , n )
  
  #  Get result
  matrix_nmi <- matrix(0, ncol = ncol(matrix), nrow = nrow(matrix))
  
  for (i in 1:nrow(id)) {
    
    matrix_nmi[id[i,1], id[i,2]] <- NMI(
      matrix[id[i,1],],
      matrix[id[i,2],])
    
    matrix_nmi[id[i,2], id[i,1]] <- NMI(
      matrix[id[i,1],],
      matrix[id[i,2],])
    
  }
  
  
  colnames(matrix_nmi) <- colnames(matrix)
  rownames(matrix_nmi) <- rownames(matrix)
  
  matrix_nmi_dist <- nmi_to_distance(matrix_nmi)
  
  return(matrix_nmi_dist)
  
  
}





# Compute mean matrix from a pair of matrix 

compute_mean_two_mat_from_list<- function(list){
  y<-apply(array(unlist(list), c(dim(list[[1]]), dim(list[[2]]), length(list))), 
           c(1,2), mean)
  colnames(y)<-colnames(list[[1]])
  rownames(y)<-rownames(list[[1]])
  return(y)
}



## Check that colnames (number and order) coincide in two lists


check_colnames_lists <- function(list1, list2) {
  
  vec_check <- c()
  
  for (i in 1:length(list1)) {
    
    if(length(which((colnames(list1[[i]]) == colnames(list2[[i]])) == FALSE)) == 0) {
      
      vec_check[i] <- "good"
      
    } else {
      
      vec_check[i] <- "PROBLEM"
      
      
    }
    
  }
  
  return(vec_check)
  
}


# Set diagonal elements of matrix to 1



diag_to0 <- function(matrix) {
  
  diag(matrix)<-0
  
  return(matrix)
  
}



# convert NaN to 1


convert_nan_to_0 <- function(vector) {
  
  vector[is.nan(vector)] <- 0
  
  return(vector)
  
}



# convert nan to 1 for all columns or rows


convet_nan_to_0_matrix <- function(matrix, marg) {
  
  m <- apply(matrix,
        MARGIN = marg,
        FUN = convert_nan_to_0)
  
  return(m)
  
}


check_all.finite <- function(matrix) {
  
  result <- all( is.finite( matrix ) )
  
  return(result)
  
}




################################################################


#Run all the functions together to compute a dataframe with the correlation values 
#for the predator, prey, and mean(pred-prey) interaction distances from the object 
#containing the results from 1 simulation

compute_df_signal_time <- function(path) {
  
  list_res <- readRDS(path)
  
  # count number of timesteps where there were spp
  
  n_steps <- length(list_res$network_list)
  
  
  
  # Crop presence matric to the number of timesteps of the simulation
  
  #(The presence_matrix (rows = timesteps, cols = species, entries = 1,0) has extra timesteps in case the simulation is longer (max = 350))
  
  presence_matrix <- list_res$presence_matrix
  
  presence_matrix <- presence_matrix[1:n_steps,]
  
  # identify timesteps with less than 3 spp (It requires at leas 3 species to compute the phylogenetic distances.
  
  non.valid_timesteps_phylo_distance <- c(which(rowSums(presence_matrix) < 3))
  
  
  
  #homogenize elements to start from valid timesteps
  
  list_anc_dist <- list_res$list_anc_dist[-non.valid_timesteps_phylo_distance]
  
  network_list <- list_res$network_list[-non.valid_timesteps_phylo_distance]
  
  
  
  # Convert spp names from numbers to letters
  
  
  ## ancestry-distances table
  
  list_anc_dist_letters <- lapply(list_anc_dist, change_sppnames_letters_ancdist.table)
  
  
  ## Network list 
  
  list_networks_sppnames_numbers <- lapply(network_list, set_sppNames_numbers)
  
  ## convert numbers to letters
  
  list_networks_sppnames_letters <- lapply(list_networks_sppnames_numbers, convert_sppnames_toletters)
  
  
  print("1/4 obtaining phylogenetic distances")
  
  # Loop for obtaining phylogenetic distances:
  
  list_newick <- list()
  list_trees <- list()
  list_newick_tails <- list()
  list_dist.phylo <- list()
  
  
  
  
  for (i in 1:length(list_anc_dist_letters)) {
    
    
    print(sprintf("step %s", i))
    
    list_newick[[i]] <- ToPhylo(list_anc_dist_letters[[i]])
    
    list_newick_tails[[i]] <- paste(list_newick[[i]], "root")
    
    list_trees[[i]] <- read.tree(text = sub("A root",";",list_newick_tails[[i]]))
    
    list_dist.phylo[[i]] <- cophenetic.phylo(list_trees[[i]])
    
    
  }
  
  
  
  
  
  
  # Retain only present species in network matrices
  
  ## Set the same spp names for the presence_matrix than for the interacion matrices
  
  
  colnames(presence_matrix) <- seq(1:1000)
  
  colnames(presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(presence_matrix))
  
  
  ## Discard same timesteps (rows) than the discarted phylogenetic distance matrices
  
  presence_matrix <- presence_matrix[-non.valid_timesteps_phylo_distance,]
  
  ## crop the interaction matrix with present spp
  
  list_net_present_spp.letters <- list()
  
  for (i in 1:length(list_networks_sppnames_letters)) {
    
    list_net_present_spp.letters[[i]] <- list_networks_sppnames_letters[[i]][names(which(presence_matrix[i,] == 1)), names(which(presence_matrix[i,] == 1))]
    
  }
  
  
  print("2/4 obtaining interaction distances")
  
  # Compute interaction distances (NMI)
  
  
  
  ## DISTANCES AS PREDATORS (columns)
  
  
  list_interact_distances_pred <- lapply(list_net_present_spp.letters, FUN = compute_nmi_aricode_pred)
  
  
  ## DISTANCES AS PREYS (rows)
  
  list_interact_distances_prey <- lapply(list_net_present_spp.letters, FUN = compute_nmi_aricode_prey)
  
  
  ## convert NAN in diagonal to 0
  
  list_interact_distances_pred <- lapply(list_interact_distances_pred, FUN = diag_to0)
  
  list_interact_distances_prey <- lapply(list_interact_distances_prey, FUN = diag_to0)
  
  
  ## convert the NAN to 0
  
  list_interact_distances_pred_corrected <- list()
  
  for (i in 1:length(list_interact_distances_pred)) {
    
    list_interact_distances_pred_corrected[[i]] <- convet_nan_to_0_matrix(matrix = list_interact_distances_pred[[i]],
                                                                          marg = 2) # 1 (rows), 2 (col), or c(1,2)
    
  }
  
  
  list_interact_distances_prey_corrected <- list()
  
  for (i in 1:length(list_interact_distances_prey)) {
    
    list_interact_distances_prey_corrected[[i]] <- convet_nan_to_0_matrix(matrix = list_interact_distances_prey[[i]],
                                                                          marg = 1) # 1 (rows), 2 (col), or c(1,2)
    
  }
  
  
  
  # For phylogenetic distance matrices, retain only those species present
  
  list_dist.phylo_pres <- list()
  
  for (i in 1:length(list_dist.phylo)) {
    
    spp_present <- colnames(list_interact_distances_pred[[i]])
    list_dist.phylo_pres[[i]] <- list_dist.phylo[[i]][spp_present ,spp_present] 
    
  }
  
  
  
  # Compute mean interaction distances between distances as preys and predators
  
  list_interact_distances_mean_corrected <- list()
  
  
  for (i in 1:length(list_interact_distances_pred_corrected)) {
    
    pair_mat <- list(list_interact_distances_pred_corrected[[i]], list_interact_distances_prey_corrected[[i]])
    
    list_interact_distances_mean_corrected[[i]] <- compute_mean_two_mat_from_list(list = pair_mat)
    
  }
  
  
  
  print("3/4 obtaining MDS")
  
  # Compute MDS
  
  
  ## First, discard those timesteps where all distances are 0
  
  matrices_all0 <- which(sapply(list_interact_distances_prey_corrected, function(x) all(x == 0))==TRUE)
  
  if(length(matrices_all0) > 0) {
    
    list_interact_distances_pred_corrected <- list_interact_distances_pred_corrected[-matrices_all0]
    
    list_interact_distances_prey_corrected <- list_interact_distances_prey_corrected[-matrices_all0]
    
    list_interact_distances_mean_corrected <- list_interact_distances_mean_corrected[-matrices_all0]
    
    list_dist.phylo_pres <- list_dist.phylo_pres[-matrices_all0]
    
  }
  
  
  
  
  # compute mds
  
  list_mds.phy <- list()
  list_mds.int_pred <- list()
  list_mds.int_prey <- list()
  list_mds.int_mean <- list()
  
  
  
  for (i in 1:length(list_dist.phylo_pres)) {
    
    list_mds.phy[[i]] <- monoMDS(list_dist.phylo_pres[[i]], y = cmdscale(list_dist.phylo_pres[[i]])) 
    
    list_mds.int_pred[[i]] <- monoMDS(list_interact_distances_pred_corrected[[i]], y = cmdscale(list_interact_distances_pred_corrected[[i]]))
    
    list_mds.int_prey[[i]] <- monoMDS(list_interact_distances_prey_corrected[[i]], y = cmdscale(list_interact_distances_prey_corrected[[i]]))
    
    list_mds.int_mean[[i]] <- monoMDS(list_interact_distances_mean_corrected[[i]], y = cmdscale(list_interact_distances_mean_corrected[[i]]))
  }
  
  
  
  
  print("4/4 obtaining correlation")
  
  
  ############# # Compute correlation - Procrustes test ##############
  
  
  
  protest_pred <- list()
  procrustes_pred <- list()
  
  protest_prey <- list()
  procrustes_prey <- list()
  
  protest_mean <- list()
  procrustes_mean <- list()
  
  
  for (i in 1:length(list_mds.phy)) {
    
    protest_pred[[i]] <- protest(list_mds.int_pred[[i]], list_mds.phy[[i]])
    
    procrustes_pred[[i]] <- procrustes(list_mds.int_pred[[i]],list_mds.phy[[i]])
    
    
  }
  
  for (i in 1:length(list_mds.phy)) {
    
    protest_prey[[i]] <- protest(list_mds.int_prey[[i]], list_mds.phy[[i]])
    
    procrustes_prey[[i]] <- procrustes(list_mds.int_prey[[i]],list_mds.phy[[i]])
    
    
  }
  
  for (i in 1:length(list_mds.phy)) {
    
    protest_mean[[i]] <- protest(list_mds.int_mean[[i]], list_mds.phy[[i]])
    
    procrustes_mean[[i]] <- procrustes(list_mds.int_mean[[i]],list_mds.phy[[i]])
    
    
  }
  
  
  
  
  
  
  
  #### store correlation and significance in dataframe
  
  
  protest_pval_pred <- c()
  protest_corr_pred <- c()
  protest_t_pred <- c()
  
  protest_pval_prey <- c()
  protest_corr_prey <- c()
  protest_t_prey <- c()
  
  protest_pval_mean <- c()
  protest_corr_mean <- c()
  protest_t_mean <- c()
  
  for (i in 1:length(protest_mean)) {
    
    protest_pval_pred[i] <- protest_pred[[i]]$signif
    
    protest_corr_pred[i] <-protest_pred[[i]]$t0
    
    protest_t_pred[i] <-mean(protest_pred[[i]]$t)
    
    
    protest_pval_prey[i] <- protest_prey[[i]]$signif
    
    protest_corr_prey[i] <-protest_prey[[i]]$t0
    
    protest_t_prey[i] <-mean(protest_prey[[i]]$t)
    
    
    protest_pval_mean[i] <- protest_mean[[i]]$signif
    
    protest_corr_mean[i] <-protest_mean[[i]]$t0
    
    protest_t_mean[i] <-mean(protest_mean[[i]]$t)
    
    
  }
  
  discarted_timesteps <- length(list_res$network_list) - length(protest_mean)
  
  timesteps = (discarted_timesteps+1):length(list_res$network_list)
  
  
  df_signal_time <- data.frame(timesteps,
                               protest_pval_pred,protest_corr_pred,protest_t_pred,
                               protest_pval_prey,protest_corr_prey,protest_t_prey,
                               protest_pval_mean,protest_corr_mean,protest_t_mean)
  
  df_signal_time$sign_pred <- with(df_signal_time, ifelse(protest_pval_pred < 0.051, 'sign', 'non.sign'))
  df_signal_time$sign_prey <- with(df_signal_time, ifelse(protest_pval_prey < 0.051, 'sign', 'non.sign'))
  df_signal_time$sign_mean <- with(df_signal_time, ifelse(protest_pval_mean < 0.051, 'sign', 'non.sign'))
  
  
  list_results <- list(df_signal_time,
         list_interact_distances_pred_corrected,list_dist.phylo_pres)
  
  return(list_results)
  
  
}




# Take the dataframe with the correlations and plot the results

plot_phylo.sign_time <- function(df_corr) {
  
  p_pred <- ggplot(df_corr, aes(x=timesteps, y=protest_corr_pred)) +
    geom_line(color="black", linetype="twodash")+
    #geom_point(aes(fill = sign_mean, color=sign_mean), size=2, alpha = 0.7) +
    theme_bw()+
    my.theme+
    scale_colour_manual(values = c("black", "red"))+
    ggtitle("as predators")+
    xlab("timesteps")+
    ylab("Phylogenetic signal")
  
  
  p_prey <- ggplot(df_corr, aes(x=timesteps, y=protest_corr_prey)) +
    geom_line(color="black", linetype="twodash")+
    #geom_point(aes(fill = sign_mean, color=sign_mean), size=2, alpha = 0.7) +
    theme_bw()+
    my.theme+
    scale_colour_manual(values = c("black", "red"))+
    ggtitle("as preys")+
    xlab("timesteps")+
    ylab("Phylogenetic signal")
  
  p_mean <- ggplot(df_corr, aes(x=timesteps, y=protest_corr_mean)) +
    geom_line(color="black", linetype="twodash")+
    #geom_point(aes(fill = sign_mean, color=sign_mean), size=2, alpha = 0.7) +
    theme_bw()+
    my.theme+
    scale_colour_manual(values = c("black", "red"))+
    ggtitle("mean prey-pred")+
    xlab("timesteps")+
    ylab("Phylogenetic signal")
  
  
  
  result_list <- list(p_pred,p_prey,p_mean)
  
  return(result_list)
  
}




compute_df_signal_time_v2 <- function(path) {
  
  list_res <- readRDS(path)
  
  # count number of timesteps where there were spp
  
  n_steps <- length(list_res$network_list)
  
  
  
  # Crop presence matric to the number of timesteps of the simulation
  
  #(The presence_matrix (rows = timesteps, cols = species, entries = 1,0) has extra timesteps in case the simulation is longer (max = 350))
  
  presence_matrix <- list_res$presence_matrix
  
  presence_matrix <- presence_matrix[1:n_steps,]
  
  # identify timesteps with less than 3 spp (It requires at leas 3 species to compute the phylogenetic distances.
  
  non.valid_timesteps_phylo_distance <- c(which(rowSums(presence_matrix) < 3))
  
  
  
  #homogenize elements to start from valid timesteps
  
  list_anc_dist <- list_res$list_anc_dist[-non.valid_timesteps_phylo_distance]
  
  network_list <- list_res$network_list[-non.valid_timesteps_phylo_distance]
  
  
  
  # Convert spp names from numbers to letters
  
  
  ## ancestry-distances table
  
  list_anc_dist_letters <- lapply(list_anc_dist, change_sppnames_letters_ancdist.table)
  
  
  ## Network list 
  
  list_networks_sppnames_numbers <- lapply(network_list, set_sppNames_numbers)
  
  ## convert numbers to letters
  
  list_networks_sppnames_letters <- lapply(list_networks_sppnames_numbers, convert_sppnames_toletters)
  
  
  print("1/4 obtaining phylogenetic distances")
  
  # Loop for obtaining phylogenetic distances:
  
  list_newick <- list()
  list_trees <- list()
  list_newick_tails <- list()
  list_dist.phylo <- list()
  
  
  
  
  for (i in 1:length(list_anc_dist_letters)) {
    
    
    print(sprintf("step %s", i))
    
    list_newick[[i]] <- ToPhylo(list_anc_dist_letters[[i]])
    
    list_newick_tails[[i]] <- paste(list_newick[[i]], "root")
    
    list_trees[[i]] <- read.tree(text = sub("A root",";",list_newick_tails[[i]]))
    
    list_dist.phylo[[i]] <- cophenetic.phylo(list_trees[[i]])
    
    
  }
  
  
  
  
  
  
  # Retain only present species in network matrices
  
  ## Set the same spp names for the presence_matrix than for the interacion matrices
  
  
  colnames(presence_matrix) <- seq(1:1000)
  
  colnames(presence_matrix) <- chartr("0123456789", "ABCDEFGHIJ", colnames(presence_matrix))
  
  
  ## Discard same timesteps (rows) than the discarted phylogenetic distance matrices
  
  presence_matrix <- presence_matrix[-non.valid_timesteps_phylo_distance,]
  
  ## crop the interaction matrix with present spp
  
  list_net_present_spp.letters <- list()
  
  for (i in 1:length(list_networks_sppnames_letters)) {
    
    list_net_present_spp.letters[[i]] <- list_networks_sppnames_letters[[i]][names(which(presence_matrix[i,] == 1)), names(which(presence_matrix[i,] == 1))]
    
  }
  
  
  print("2/4 obtaining interaction distances")
  
  # Compute interaction distances (NMI)
  
  
  
  ## DISTANCES AS PREDATORS (columns)
  
  
  list_interact_distances_pred <- lapply(list_net_present_spp.letters, FUN = compute_nmi_aricode_pred)
  
  
  ## DISTANCES AS PREYS (rows)
  
  list_interact_distances_prey <- lapply(list_net_present_spp.letters, FUN = compute_nmi_aricode_prey)
  
  
  ## convert NAN in diagonal to 0
  
  list_interact_distances_pred <- lapply(list_interact_distances_pred, FUN = diag_to0)
  
  list_interact_distances_prey <- lapply(list_interact_distances_prey, FUN = diag_to0)
  
  
  ## convert the NAN to 0
  
  list_interact_distances_pred_corrected <- list()
  
  for (i in 1:length(list_interact_distances_pred)) {
    
    list_interact_distances_pred_corrected[[i]] <- convet_nan_to_0_matrix(matrix = list_interact_distances_pred[[i]],
                                                                          marg = 2) # 1 (rows), 2 (col), or c(1,2)
    
  }
  
  
  list_interact_distances_prey_corrected <- list()
  
  for (i in 1:length(list_interact_distances_prey)) {
    
    list_interact_distances_prey_corrected[[i]] <- convet_nan_to_0_matrix(matrix = list_interact_distances_prey[[i]],
                                                                          marg = 1) # 1 (rows), 2 (col), or c(1,2)
    
  }
  
  
  
  # For phylogenetic distance matrices, retain only those species present
  
  list_dist.phylo_pres <- list()
  
  for (i in 1:length(list_dist.phylo)) {
    
    spp_present <- colnames(list_interact_distances_pred[[i]])
    list_dist.phylo_pres[[i]] <- list_dist.phylo[[i]][spp_present ,spp_present] 
    
  }
  
  
  
  # Compute mean interaction distances between distances as preys and predators
  
  list_interact_distances_mean_corrected <- list()
  
  
  for (i in 1:length(list_interact_distances_pred_corrected)) {
    
    pair_mat <- list(list_interact_distances_pred_corrected[[i]], list_interact_distances_prey_corrected[[i]])
    
    list_interact_distances_mean_corrected[[i]] <- compute_mean_two_mat_from_list(list = pair_mat)
    
  }
  
  
  
  print("3/4 obtaining MDS")
  
  # Compute MDS
  
  
  ## First, discard those timesteps where all distances are 0
  
  matrices_all0 <- which(sapply(list_interact_distances_prey_corrected, function(x) all(x == 0))==TRUE)
  
  if(length(matrices_all0) > 0) {
    
    list_interact_distances_pred_corrected <- list_interact_distances_pred_corrected[-matrices_all0]
    
    list_interact_distances_prey_corrected <- list_interact_distances_prey_corrected[-matrices_all0]
    
    list_interact_distances_mean_corrected <- list_interact_distances_mean_corrected[-matrices_all0]
    
    list_dist.phylo_pres <- list_dist.phylo_pres[-matrices_all0]
    
  }
  
  
  
  
  
  
  print("4/4 obtaining correlation")
  
  
  ############# # Compute correlation - Procrustes test ##############
  
  
  
  protest_pred <- list()
  procrustes_pred <- list()
  
  protest_prey <- list()
  procrustes_prey <- list()
  
  protest_mean <- list()
  procrustes_mean <- list()
  
  
  for (i in 1:length(list_dist.phylo_pres)) {
    
    protest_pred[[i]] <- protest(list_interact_distances_pred_corrected[[i]], list_dist.phylo_pres[[i]])
    
    procrustes_pred[[i]] <- procrustes(list_interact_distances_pred_corrected[[i]],list_dist.phylo_pres[[i]])
    
    
  }
  
  for (i in 1:length(list_dist.phylo_pres)) {
    
    protest_prey[[i]] <- protest(list_interact_distances_prey_corrected[[i]], list_dist.phylo_pres[[i]])
    
    procrustes_prey[[i]] <- procrustes(list_interact_distances_prey_corrected[[i]],list_dist.phylo_pres[[i]])
    
    
  }
  
  for (i in 1:length(list_dist.phylo_pres)) {
    
    protest_mean[[i]] <- protest(list_interact_distances_mean_corrected[[i]], list_dist.phylo_pres[[i]])
    
    procrustes_mean[[i]] <- procrustes(list_interact_distances_mean_corrected[[i]],list_dist.phylo_pres[[i]])
    
    
  }
  
  
  
  
  
  
  
  #### store correlation and significance in dataframe
  
  
  protest_pval_pred <- c()
  protest_corr_pred <- c()
  protest_t_pred <- c()
  
  protest_pval_prey <- c()
  protest_corr_prey <- c()
  protest_t_prey <- c()
  
  protest_pval_mean <- c()
  protest_corr_mean <- c()
  protest_t_mean <- c()
  
  for (i in 1:length(protest_mean)) {
    
    protest_pval_pred[i] <- protest_pred[[i]]$signif
    
    protest_corr_pred[i] <-protest_pred[[i]]$t0
    
    protest_t_pred[i] <-mean(protest_pred[[i]]$t)
    
    
    protest_pval_prey[i] <- protest_prey[[i]]$signif
    
    protest_corr_prey[i] <-protest_prey[[i]]$t0
    
    protest_t_prey[i] <-mean(protest_prey[[i]]$t)
    
    
    protest_pval_mean[i] <- protest_mean[[i]]$signif
    
    protest_corr_mean[i] <-protest_mean[[i]]$t0
    
    protest_t_mean[i] <-mean(protest_mean[[i]]$t)
    
    
  }
  
  discarted_timesteps <- length(list_res$network_list) - length(protest_mean)
  
  timesteps = (discarted_timesteps+1):length(list_res$network_list)
  
  
  df_signal_time <- data.frame(timesteps,
                               protest_pval_pred,protest_corr_pred,protest_t_pred,
                               protest_pval_prey,protest_corr_prey,protest_t_prey,
                               protest_pval_mean,protest_corr_mean,protest_t_mean)
  
  df_signal_time$sign_pred <- with(df_signal_time, ifelse(protest_pval_pred < 0.051, 'sign', 'non.sign'))
  df_signal_time$sign_prey <- with(df_signal_time, ifelse(protest_pval_prey < 0.051, 'sign', 'non.sign'))
  df_signal_time$sign_mean <- with(df_signal_time, ifelse(protest_pval_mean < 0.051, 'sign', 'non.sign'))
  
  
  
  return(df_signal_time)
  
  
}



eliminate_basals <- function(matrix, nbasals) {
  
   mat <- matrix[(nbasals+1):nrow(matrix),]
  
   return(mat)
}




































