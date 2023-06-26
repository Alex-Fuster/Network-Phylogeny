

# Create pool by combining traits dataframes of simulations


create_pool <- function(list_file_paths) {
  
  pool <- data.frame(matrix(ncol=4,nrow=0, 
                            dimnames=list(NULL, c("community", "n", 
                                                  "r","o")))) 
  
  for (i in 1:length(file_paths_fac)) {
    
    sim <- readRDS(file_paths_fac[[i]])
    
    present_spp_f <- which(sim$presence_matrix[length(sim$network_list),] == 1)
    
    df_traits <- sim$traits_df[present_spp_f,]
    
    res.add <- data.frame("community" = rep(i, times = nrow(df_traits)),
                          "n" = df_traits[,"n"],
                          "r" = df_traits[,"r"],
                          "o" = df_traits[,"o"])
    
    pool<-rbind(pool,res.add)
    
  }
  
  return(pool)
  
}