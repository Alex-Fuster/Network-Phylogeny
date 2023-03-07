rm(list = ls())

pars = list()

# General parameters
#pars$int = 0 # 0 for competition, 1 for facilitation, 2 for predation
pars$Smax = 10000 # Maximal number of species in the system

# Related to interaction range trait
pars$av_r = 0.1 # Half range of the niche of the first species
pars$sd = 0.1*pars$av_r + 0.0001 # Standard deviation of the normal distribution used to calculate the niche optimum trait of a new species

pars$u_max = 0.15 # Speciation probability

# Extinction probability, negative interactions communities
pars$a_eneg = 0.025 # Shape of the exponential decay of the negative extinction - interaction relationship
pars$e_0neg = 0.5 # Asymptotic extinction probability with infinite negative interactions
pars$e_1neg = 1  # Extinction probability with absence of interactions

# Extinction probability, positive interactions communities
pars$a_epos = 1.2  # Shape of the exponential decay of the positive extinction - interaction relationship
pars$e_0pos = 0.075 # Asymptotic extinction probability with infinite positive interactions
pars$e_1pos = 5.19 # 1 - pars$e_0pos

# Establishment probability, negative interactions communities
pars$a_uneg = 0.075 # Shape of the exponential decay of the colonization - interaction relationship
pars$u_0neg = 0.075 # Asymptotic establishment probability with infinite competition interactions
pars$u_1neg = 2 # Establishment probability with absence of competition interaction

# Establishment probability, positive interactions communities
pars$a_u = 0.45 # Shape of the exponential decay of the colonization - interaction relationship
pars$u_0 = 1 # Asymptotic establishment probability with infinite facilitation interactions
pars$u_1 = -1 # Establishment probability with absence of facilitation interactions
pars$d = 0.5 # Decrease speed of the establishment probability

# Extinction and establishment probabilitie for both, positive and negative interactions communities
pars$Bspe = 4 # Constant minimal number of interaction per species (extablishment prob)
pars$Bext = 4 # Constant minimal number of interaction per species (extinction prob)
pars$I_max = 40 # Maximal number of interactioning species

source("code/functions/functions_M-C_bifurcation.R")
source("code/functions/list_append_function.R")
source("code/figures/traits_overtime_plot_function.R")


nsteps = 350 # Set the maximum timestep per simulation
nsim = 1 # Set the number of simulations

# List of used seed per simulation batch
reccorded_seed_pos_list <- list()
reccorded_seed_neg_list <- list()

seed_list <- NULL # Set an object to record seed that will be used to launch simulations

total_nb_sim <- 0
list_res <- list() # Set the list to store the results

#par(mfrow = c(4,5))

#profvis::profvis({
start_time <- Sys.time()

for(int in 0:1){
	pars$int = int # Define the interaction type: 0 = competition ; 1 = facilitation
	successful_sim <- 1 # Set the count of simulations

	if(pars$int == 1){
		print("Positive interctions")
	} else {
		print("Negative interctions")
	}

	seed_record <- c(rep(NA, nsim)) # Set the vector to record the seed that have been used

	while(successful_sim < nsim+1){

		seed <- sample.int(80000, 1) # Pick a random number, this number will be the seed for one simulation
		print(seed)

		# If a seed have already been tested, pick another one
		if(length(which(seed_list == seed)) != 0){
			while(seed %in% seed_list){
  			seed <- sample.int(80000, 1)
			}
		}

		seed_list <- c(seed_list, seed) # Add the seed

		# SIMULATION ; from the file code/functions/functions_M-C_bifurcation.R
		simulation <- sim_model_bif(seed = seed, pars, nsteps = nsteps)

		total_nb_sim <- total_nb_sim +1 # Count the total aount of simulations

		# Test if we have enough species at the timestep 150 in the simulation
		success <- sum(simulation$pres[150,])
		#print(paste0("Success = ", success))

		if(success >= 20){
			#print("the simulation as more than 19 species at the time step 150")
			traits_overtime_plot(simulation)

  		seed_record[successful_sim] <- seed # Record seed which alows us to have "good" simulation
  		res_sim <- list(simulation_number = paste0("simulation", successful_sim), seed = seed, parameters = pars, presence_matrix = simulation$pres, traits_df = simulation$traits, parentage_matrix = simulation$anc, extinxction_matrix = simulation$extinct) # Record results from the simulation
  		list_res <- list_append(list_res, res_sim)

  		print(paste("simulation ", successful_sim, "sur", nsim))
  		successful_sim <- successful_sim + 1 # Count the amount of "good" simulations
		}
	}
}

end_time <- Sys.time()
end_time - start_time
#})

x11()

par(mfrow = c(1,1))
com_dynamic <- matrix(nr = nsteps-1, nc = length(list_res))

for(sim in 1:length(list_res)){
	parameters <- as.list(unlist(list_res[[sim]][3]))
  names(parameters) <- gsub("parameters.","",names(parameters))
  pres <- matrix(unlist(list_res[[sim]][4]), nr = nsteps, nc = parameters$Smax)

	species_richness <- rowSums(pres)
	species_richness <- replace(species_richness, species_richness == 0, NA)
	com_dynamic[, sim] <- species_richness[-350]
}

plot(com_dynamic[,1], type = "l", col = "red", ylim = c(0,50))
for(i in 2:(length(list_res)/2)){
	#Sys.sleep(0.01)
	lines(com_dynamic[,i], col = "red")
}
for(i in ((length(list_res)/2)+1):length(list_res)){
	#Sys.sleep(0.01)
	lines(com_dynamic[,i])
}


MeanPlotNegPos <- function(measure, thing_we_measure){
	mean_measure_neg <- apply(measure[,1:nsim/2],1, mean, na.rm = TRUE)
	quantiles_neg <- apply(measure[,1:nsim/2], 1, quantile, probs = c(0.05, 0.95), na.rm = T)

	y_min_neg <- min(quantiles_neg, na.rm = T)
	y_max_neg <- max(quantiles_neg, na.rm = T)

	if(y_min_neg < 0 | y_max_neg > 1){
		y_range_neg = c(y_min_neg, y_max_neg)
	}else{
		y_range_neg = c(0, 1)
	}

	mean_measure_pos <- apply(measure[,((nsim/2)+1):nsim],1, mean, na.rm = TRUE)
	quantiles_pos <- apply(measure[,((nsim/2)+1):nsim], 1, quantile, probs = c(0.05, 0.95), na.rm = T)

	y_min_pos <- min(quantiles_pos, na.rm = T)
	y_max_pos <- max(quantiles_pos, na.rm = T)

	if(y_min_pos < 0 | y_max_pos > 1){
		y_range_pos = c(y_min_pos, y_max_pos)
	}else{
		y_range_pos = c(0, 1)
	}

	## Plot for Positive interactions
	par(mar = c(3,3,0.5,0.5), mgp = c(1.5, 0.3, 0), tck = -.008, las = 1)
	plot(mean_measure_pos[1:200], pch=19, xlab = "Time", xlim = c(0,200), ylab = paste("Mean", thing_we_measure), ylim = y_range_pos, col = c(rgb(0, 158, 115, max = 255)))
	polygon(c(1:200,rev(1:200)), c(quantiles_pos[1,1:200], rev(quantiles_pos[2,1:200])), col = adjustcolor(rgb(161,215,106, max = 255), alpha = 0.5), border = NA)
	points(mean_measure_neg[1:200], pch=19, col = c(rgb(230, 159, 0, max = 255)))
	polygon(c(1:200,rev(1:200)), c(quantiles_neg[1,1:200], rev(quantiles_neg[2,1:200])), col = adjustcolor(rgb(230, 159, 0, max = 255), alpha = 0.3), border = NA)

}


source("code/functions/functions_M-C_bifurcation.R")
source("code/functions/tree_functions.R")
source("code/functions/phylo_functions.R")
source("code/functions/add_extinction_function.R")
source("code/functions/node_dist_function.R")

nsim <- length(list_res)
nsteps <- 350
selected_time <- c(25, 50, 75, 100, 125, 150)


#############################################
## Matrices vides pour les figures de base ##
#############################################

# Connectance matrix
connectance <- matrix(nr = nsteps, nc = nsim)
# Mean degree
mean_degree <- matrix(NA, nr = nsteps, nc = nsim)
# Modularity matrix
L_n <- matrix(0, nr = sp_max, nc = sp_max)
# Gamma statistics matrix
alpha_stats <- matrix(nr = nsteps, nc = nsim)
# Tree space
tree_space_an <- list()
# Community dynamic matrix
com_dynamic <- matrix(nr = nsteps-1, nc = nsim)
#time_to_equi <- c(rep(NA, nsim))
## Interaction range
int_range <- matrix(nr = nsteps, nc = nsim)
# Correlation between phylogenetic and trait distances
correlation_pdtd <- matrix(0, nrow = nsim, ncol = length(selected_time))

Species <- matrix(NA, nr = nsteps, nc = nsim)
Links <- matrix(NA, nr = nsteps, nc = nsim)
Links_in <- matrix(NA, nr = nsteps, nc = nsim)
Links_out <- matrix(NA, nr = nsteps, nc = nsim)

degree_in_tot <- NULL
degree_out_tot <- NULL
age_tot <- NULL

sim_to_remove <- NULL

#############################################################
## Remplissage des matrices pour faire les figures de base ##
#############################################################

for(sim in 1:length(list_res)){
  print(paste("simulation", sim, "on", nsim))
  parameters <- as.list(unlist(list_res[[sim]][3]))
  names(parameters) <- gsub("parameters.","",names(parameters))
  pres <- matrix(unlist(list_res[[sim]][4]), nr = nsteps, nc = parameters$Smax)
  traits <- as.data.frame(list_res[[sim]][5])
  names(traits) <- gsub("traits_df.","",names(traits))
  anc <- matrix(unlist(list_res[[sim]][6]), nr = parameters$Smax, nc = 3) #time, ancestor, new species
  extinction <- matrix(unlist(list_res[[sim]][7]), nr = parameters$Smax, nc = 2) #time, species extinct

  # Phylogenetic statistics
	species_total <- max(anc[,3], na.rm = T)

	sum(pres[150,])

	if(sum(pres[150,]) >= 15){
		distance_mat <- get_dist_mat(anc, nsteps, species_total)

		for(p in 1:nsteps){
			species_number <- sum(pres[p,])
			if(species_number > 2){
				to_keep <- which(pres[p,] == 1)
				last_timestep_dist <- distance_mat[c(to_keep), c(to_keep)]
				colnames(last_timestep_dist) <- to_keep
				rownames(last_timestep_dist) <- to_keep

	    	tree_dist <- ape::as.phylo(hclust(as.dist(last_timestep_dist), method = "average"))
				ttdist <- timetree(tree_dist)
				alpha_stats[p, sim] <- alpha_stat(ttdist) #function in phylo_functions.R file
			}
		}
  	extinction_tree <- add_extinction(distance_mat, extinction, anc)
  	tree_space_an <- list_append(tree_space_an, extinction_tree)
	} else {
		print(paste("Not enough species for the", list_res[[sim]]$simulation_number, "of interaction type =", list_res[[1]]$parameters$int," to calculate phylogenetic mesures such as alpha stat (1 = positive, 0 = negative)"))
		sim_to_remove  <- c(sim_to_remove, i)
	}

  # Network structure indices
  # Connectance, age and degree

	age_sp <- matrix(NA, nr = nrow(pres), nc = ncol(pres))
  degree_in <- matrix(NA, nr = nrow(pres), nc = ncol(pres))
  degree_out <- matrix(NA, nr = nrow(pres), nc = ncol(pres))

	for(m in 1:nsteps){

		# Age and Degree
		if(m == 1){
      age_sp[m,] <- pres[m,]
    } else {
      age_sp[m,] <- colSums(pres[1:m,])
    }

    sp_pres <- pres[m,] == 1
    not_pres <- which(sp_pres == FALSE)
    age_sp[m, not_pres] = 0
		#

		pres_vec <- pres[m,]
		cooc <- matrix(pres_vec, nr=parameters$Smax, nc=parameters$Smax, byrow=TRUE)*
		matrix(pres_vec, nr=parameters$Smax, nc=parameters$Smax, byrow=FALSE)
		L_n <- get_L_mat(parameters, traits)
		L_n <- L_n*cooc

		#sp_pres <- pres[m,] == 1
		adj_mat <- L_n[sp_pres, sp_pres]

		Links_nb <- sum(rowSums(L_n))
		Species_nb <- sum(pres_vec)

		connectance[m, sim] <- Links_nb / Species_nb^2

		# Age and Degree
		Species[m, sim] <- sum(pres_vec)
    Links[m, sim] <- sum(L_n)

    degree_in[m,] <- rowSums(L_n)
    degree_out[m,] <- colSums(L_n)

   # Degree
     if(length(adj_mat)>1){
       mean_degree[m, sim] <- mean(rowSums(adj_mat), na.rm = T)
     }

   # Interaction ranges
     r <- matrix(traits[,2], nr = nsteps, nc = ncol(pres), byrow=TRUE)
     range_pres <- pres[1:nsteps,]*r
     for(rg in 1:nsteps){
       sp_pres <- pres[rg,] == 1
       int_range[rg, sim] <- mean(range_pres[rg, sp_pres]*2)
     }
	}

	# Age and Degree
  degree_in_tot <- c(degree_in_tot, degree_in)
  degree_out_tot <- c(degree_out_tot, degree_out)
  age_tot <- c(age_tot, age_sp)

  # Community Dynamic
	species_richness <- rowSums(pres)
	species_richness <- replace(species_richness, species_richness == 0, NA)
  com_dynamic[, sim] <- species_richness[-350]

	##print(paste("simulation successfull", successful_sim, "sur", nsim))
	#print(paste("données de la simulation", i, "enregistrées ! :)"))
}

## Record Data
#save(com_dynamic, file = "Data/species_richness_dynamic.csv")
#save(alpha_stats, file = "Data/alpha_statistics.csv")
#save(connectance, file = "Data/connectance.csv")
#save(int_range, file = "Data/interaction_range.csv")
#save(degree, file = "Data/mean_degree.csv")


# Remove unvalid simulations
com_dynamic[,sim_to_remove] = NA
alpha_stats[,sim_to_remove] = NA
connectance[,sim_to_remove] = NA
int_range[,sim_to_remove] = NA
mean_degree[,sim_to_remove] = NA


################################################################################
################### FIGURE 1 - Model's Vital Check #############################
################################################################################

x11()

#path <- paste0("Figures/basic_figures_", pourcentage, ".png")

#png(path)
#par(mfrow = c(3,2))
layout(matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), 2, byrow = TRUE))

MeanPlotNegPos(com_dynamic, "species richness")
MeanPlotNegPos(alpha_stats, "alpha value")
MeanPlotNegPos(connectance, "connectance")
MeanPlotNegPos(int_range, "niche range")
MeanPlotNegPos(mean_degree, "of mean degree")

#dev.off()


################################################################################
############# Dynamic Figure : Species turn-over ratio #########################
################################################################################

turnover_rate <- matrix(NA, nrow = nsteps, ncol = nsim)

for(sim in 1:length(list_res)){
  print(paste("simulation", sim, "on", nsim))
  parameters <- as.list(unlist(list_res[[sim]][3]))
  names(parameters) <- gsub("parameters.","",names(parameters))
  pres <- matrix(unlist(list_res[[sim]][4]), nr = nsteps, nc = parameters$Smax)

	sp_turnover <- NULL

	for(sp in 1:ncol(pres)){
		life_time <- expand.grid(sp, which(pres[, sp] == 1))
		sp_turnover <- rbind(sp_turnover, life_time)
	}

	abundances <- as.integer(rep(1, nrow(sp_turnover)))
	sp_turnover <- cbind(sp_turnover, abundances)

	relative_total_turnover <- codyn::turnover(df= sp_turnover,
               	time.var = "Var2",
               	species.var = "Var1",
						 	 	abundance.var = "abundances")

	turnover_rate[2:(nrow(relative_total_turnover)+1), sim] <- relative_total_turnover[,1]
}

mean_turnover_neg <- rowMeans(turnover_rate[1:150,1:10], na.rm = TRUE)
mean_turnover_pos <- rowMeans(turnover_rate[1:150,11:20], na.rm = TRUE)

x11()
#path <- paste0("Figures/community_turnover_", pourcentage, ".png")
#png(path)

par(mfrow=c(2,1))

MeanPlotNegPos(com_dynamic, "species richness")
mtext("A.", 2, adj=3, las=1, padj=-10)


plot(mean_turnover_pos, type = "l", lwd = 3, ylim = c(0,0.5), col = rgb(0, 158, 115, max = 255), xlab = "Time", ylab = "Mean Species Turnover")
lines(mean_turnover_neg, lwd = 3, col = rgb(230, 159, 0, max = 255))
mtext("B.", 2, adj=3, las=1, padj=-10)

#dev.off()

################################################################################
########### Dynamic Figure : Origination and Extinction rates ##################
################################################################################

origination_matrix <- matrix(NA, nrow = nsteps, ncol = nsim)
extinction_matrix <- matrix(NA, nrow = nsteps, ncol = nsim)
sp_pres <- matrix(NA, nrow = nsteps, ncol = nsim)

for(sim in 1:length(list_res)){
  print(paste("simulation", sim, "on", nsim))
  parameters <- as.list(unlist(list_res[[sim]][3]))
  names(parameters) <- gsub("parameters.","",names(parameters))
  pres <- matrix(unlist(list_res[[sim]][4]), nr = nsteps, nc = parameters$Smax)
	anc <- matrix(unlist(list_res[[sim]][6]), nr = parameters$Smax, nc = 3) #time, ancestor, new species
  ext <- matrix(unlist(list_res[[sim]][7]), nr = parameters$Smax, nc = 2) #time, species extinct

	# Number of species at each timesteps
	sp_pres[, sim] <- rowSums(pres)

	time <- as.data.frame(c(1:nsteps))
	names(time) <- c("time")

	# Number of new species establisment at each timestep
	origination_freq <- as.data.frame(table(anc[, 1]))
	names(origination_freq) <- c("time", "freq")
	origination_freq$time <- as.numeric(as.character(origination_freq$time))

	origination <- dplyr::left_join(time, origination_freq, by = "time")

	last_time_origi <- max(anc[,1], na.rm = TRUE)
	no_origi <- which(is.na(origination$freq[1:last_time_origi]))
	origination$freq[no_origi] <- 0

	origination_matrix[, sim] <- origination[, 2]

	# Number of species extinction at each timestep
	extinction_freq <- as.data.frame(table(ext[, 1]))
	names(extinction_freq) <- c("time", "freq")
	extinction_freq$time <- as.numeric(as.character(extinction_freq$time))

	extinction <- dplyr::left_join(time, extinction_freq, by = "time")

	last_time_extinct <- max(ext[,1], na.rm = TRUE)
	no_extinct <- which(is.na(extinction$freq[1:last_time_extinct]))
	extinction$freq[no_extinct] <- 0

	extinction_matrix[, sim] <- extinction[, 2]
}

sp_tmoins1 <- sp_pres[1:(nsteps - 1),]

origination_t <- origination_matrix[2:nsteps,]

origination_neg <- rowMeans(origination_t[,1:10] / sp_tmoins1[,1:10], na.rm = TRUE)
origination_pos <- rowMeans(origination_t[,11:20] / sp_tmoins1[,11:20], na.rm = TRUE)

extinction_t <- extinction_matrix[2:nsteps,]

extinction_neg <- rowMeans(extinction_t[,1:10] / sp_tmoins1[,1:10], na.rm = TRUE)
extinction_pos <- rowMeans(extinction_t[,11:20] / sp_tmoins1[,11:20], na.rm = TRUE)

x11()
#path <- paste0("Figures/orignation_extinction_rates", pourcentage, ".png")
#png(path)

par(mfrow=c(2,1))

y_max <- max(c(origination_pos, origination_neg, extinction_pos, extinction_neg), na.rm = TRUE)

plot(origination_pos, type = "l", lwd = 3, col = rgb(0, 158, 115, max = 255), ylim = c(0,y_max), xlim = c(0,150), xlab = "Time", ylab = "Oigination rate")
lines(origination_neg, lwd = 3, col = rgb(230, 159, 0, max = 255))
mtext("A.", 2, adj=3, las=1, padj=-10)

plot(extinction_pos, type = "l", lwd = 3, col = rgb(0, 158, 115, max = 255), ylim = c(0,y_max), xlim = c(0,150), xlab = "Time", ylab = "Extinction rate")
lines(extinction_neg, lwd = 3, col = rgb(230, 159, 0, max = 255))
mtext("B.", 2, adj=3, las=1, padj=-10)

#dev.off()


################################################################################
################### FIGURE 2 - Degrees and the rest ############################
################################################################################

# Différenciate point between positive and negative interactions communities values
midd <- length(degree_in_tot)/2

# Positive interactions communities
# Create a table to store values of degree in (first row), degree out (second row) and related species age (third row)
degree_age_neg <- rbind(degree_in_tot[1:midd], degree_out_tot[1:midd], age_tot[1:midd])
## Remove the ones with zero interactions
to_keep_neg <- which(degree_age_neg[3,] == 1 | degree_age_neg[3,] > 1)
degree_age_neg <- degree_age_neg[,to_keep_neg]
#test <- degree_age_pos[, !duplicated(t(degree_age_pos))]

# Negative interactions communities
# Create a table to store values of degree in (first row), degree out (second row) and related species age (third row)
degree_age_pos <- rbind(degree_in_tot[(midd+1):length(degree_in_tot)], degree_out_tot[(midd+1):length(degree_in_tot)], age_tot[(midd+1):length(degree_in_tot)])
## Remove the ones with zero interactions
to_keep_pos <- which(degree_age_pos[3,] == 1 | degree_age_pos[3,] > 1)
degree_age_pos <- degree_age_pos[,to_keep_pos]

# Plot age vs. degree
#path <- paste0("Figures/age-degree_figures_", pourcentage, ".png")
#png(path)

x11()
#par(mfrow= c(1,1), mar = c(3,3,0.5,0.5), mgp = c(1.5, 0.3, 0), tck = -.008, las = 1)

plot(degree_age_pos[1,], degree_age_pos[3,], col = rgb(0, 158, 115, max = 255, alpha = 75), pch = 16, xlim = c(0,50), ylim = c(0,80), xlab = "Interactions number", ylab = "Brench length")
points(degree_age_neg[1,], degree_age_neg[3,], col = rgb(230, 159, 0, max = 255, alpha = 75), pch = 16, xlim = c(0,50), ylim = c(0,80))

#dev.off()

# Figure with separated positive and negative simulations
#path <- paste0("Figures/age-degree_figures_pos-neg_", pourcentage, ".png")
#png(path)

x11()
par(mfrow= c(2,1), mar = c(3,3,0.5,0.5), mgp = c(1.5, 0.3, 0), tck = -.008, las = 1)

# Only positive interactions simulations
plot(degree_age_pos[1,], degree_age_pos[3,], col = rgb(0, 158, 115, max = 255, alpha = 50), pch = 16, xlim = c(0,50), ylim = c(0,80), cex = 1.5, xlab = "Interactions number", ylab = "Species age", las = 1)
# Only negative interactions simulations
plot(degree_age_neg[1,], degree_age_neg[3,], col = rgb(230, 159, 0, max = 255, alpha = 50), pch = 16, xlim = c(0,50), ylim = c(0,80), cex = 1.5, xlab = "Interactions number", ylab = "Species age", las = 1)

#dev.off()

################################################################################
################# Nombre d'interactions total à l'équilibre ####################
################################################################################

## In degree
eq_time <- 175

n_sp <- c(rep(NA, nsim))
n_interact <- c(rep(NA, nsim))

degree_eq <- NULL
age_eq <- NULL

sp_max <- list_res[[1]][3]$parameters$Smax

for(sim in 1:length(list_res)){
	pres <- list_res[[sim]]$presence_matrix
	trait_mat <- list_res[[sim]]$traits_df

	n_sp[sim] <- sum(pres[eq_time,])

	pres_vec = pres[eq_time,]
	cooc = matrix(pres_vec, nr=sp_max, nc=sp_max, byrow=TRUE) * matrix(pres_vec, nr=sp_max, nc=sp_max, byrow=FALSE)
	L = get_L_mat(list_res[[sim]]$parameters, trait_mat)
	L = L*cooc
	in_degree <- rowSums(L)
	n_interact[sim] <- mean(in_degree[which(in_degree != 0)])

	sp_pres_eq <- pres[eq_time,]
	sp_pres <- which(sp_pres_eq != 0)
	age_sp_pres <- colSums(pres[1:eq_time, sp_pres])
	age_eq <- append(age_eq, age_sp_pres)
	degree_eq <- append(degree_eq, in_degree[which(sp_pres_eq != 0)])
}

#path <- paste0("Figures/interaction_nb_at_eq_", pourcentage, ".png")
#png(path)

x11()

par(mfrow= c(1,1), mar = c(3,3,0.5,0.5), mgp = c(1.5, 0.3, 0), tck = -.008, las = 1)

x_lim <- c(min(n_interact, na.rm = T), max(n_interact, na.rm = T))
y_lim <- c(min(n_sp, na.rm = T), max(n_sp, na.rm = T))

plot(n_interact[((length(list_res)/2)+1):length(list_res)], n_sp[((length(list_res)/2)+1):length(list_res)], col = rgb(0, 158, 115, max = 255), pch = 16, xlim = x_lim,  ylim = y_lim, ylab = "Species richness", xlab = "Mean interactions number")
points(n_interact[1:length(list_res)/2], n_sp[1:length(list_res)/2], col = rgb(230, 159, 0, max = 255), pch = 16)

#dev.off()


# Species age at the equilibrium

neg <- sum(n_sp[1:(length(list_res)/2)])

x11()
plot(degree_eq[neg+1:length(degree_eq)], age_eq[neg+1:length(degree_eq)], col = rgb(0, 158, 115, max = 255), pch = 16)
points(degree_eq[1:neg], age_eq[1:neg], col = rgb(230, 159, 0, max = 255), pch = 16)


interaction_type <- c(rep("Negative", neg), rep("Positive", (length(age_eq) - neg)))
data_age_eq <- data.frame(age_eq, interaction_type)


# Violin plot of species age at the equilibrium
#install.packages("vioplot")

#path <- paste0("Figures/sp_age_at_eq_", pourcentage, ".png")
#png(path)

x11()

with(data_age_eq , vioplot::vioplot(
  age_eq[interaction_type=="Positive"] , age_eq[interaction_type=="Negative"],
  col = c(adjustcolor(rgb(0, 158, 115, max = 255), alpha = 0.5), adjustcolor(rgb(230, 159, 0, max = 255), alpha = 0.5)),
  border = c(rgb(0, 158, 115, max = 255), rgb(230, 159, 0, max = 255)),
  lineCol = c(rgb(0, 158, 115, max = 255), rgb(230, 159, 0, max = 255)),
  rectCol = c(rgb(0, 158, 115, max = 255), rgb(230, 159, 0, max = 255)),
  names=c("Positive interactions", "Negative interactions"),
	ylab="Species age"
)
)

#dev.off()

#path <- paste0("Figures/sp_degree_age_at_eq_", pourcentage, ".png")
#png(path)

################################################################################
##################### Ages des espèces dans le temps ###########################
################################################################################

# Test avec une simulation
x11()
par(mfrow = c(4,5))

for(sim in 1:length(list_res)){
	parameters <- as.list(unlist(list_res[[sim]][3]))
	names(parameters) <- gsub("parameters.","",names(parameters))
	pres <- matrix(unlist(list_res[[sim]][4]), nr = nsteps, nc = parameters$Smax)

	sp_age <- colSums(pres)
	plot(sp_age, type = "l", ylim = c(0,50))
}

################################################################################
################### Occupation de niche dans le temps ##########################
################################################################################

source("code/functions/list_append_function.R")

x11()
par(mfrow = c(4,5))

for(sim in 1:length(list_res)){
	density_time <- list()
	pres <- list_res[[sim]]$presence_matrix
	trait_mat <- list_res[[sim]]$traits_df

	for(time in selected_time){
		pres_sp <- which(pres[time,] == 1)
		pres_traits <- trait_mat[pres_sp, 1]

		if(sum(pres[time,]) > 1){
			d <- density(pres_traits)
			density_time <- list_append(density_time, d)

		}
	}

	plot(density_time[[1]], xlim = c(0,1))
	for(plot in 2:length(density_time)){
		lines(density_time[[plot]], col = plot)
	}
}


################################################################################
############ Probabiliés en fonction de la richesse specifique #################
################################################################################

origination_t
extinction_t

sp_richness <-  matrix(NA, nrow = nsteps, ncol = nsim)

for(sim in 1:length(list_res)){
	pres <- list_res[[sim]]$presence_matrix
	sp_richness[,sim] <- rowSums(pres)

}

x11()
par(mfrow = c(1,2))
plot(sp_richness[-350,1:10], origination_t[,1:10], col = rgb(0, 158, 115, max = 255))
plot(sp_richness[-350,11:20], origination_t[,11:20], col =  rgb(230, 159, 0, max = 255))
title("origination")

x11()
par(mfrow = c(1,2))
plot(sp_richness[-350,1:10], extinction_t[,1:10], col = rgb(0, 158, 115, max = 255))
plot(sp_richness[-350,11:20], extinction_t[,11:20], col =  rgb(230, 159, 0, max = 255))
title("extinction")


# Plot at specific timesteps, from selected_time object
x11()
par(mfrow = c(2,2))

# Origination
plot(mean(sp_richness[1,1:10]), mean(origination_t[1,1:10]), pch = 16, col = rgb(230, 159, 0, max = 255), xlim = c(0, max(sp_richness)), ylim = c(0,max(origination_t, na.rm = T)))
for(time in selected_time){
	points(mean(sp_richness[time,1:10]), mean(origination_t[time,1:10]), pch = 16, col = time)
}

plot(mean(sp_richness[1,11:20]), mean(origination_t[1,11:20]), pch = 16, col = rgb(0, 158, 115, max = 255), xlim = c(0, max(sp_richness)), ylim = c(0,max(origination_t, na.rm = T)))
for(time in selected_time){
	points(mean(sp_richness[time,11:20]), mean(origination_t[time,11:20]), pch = 16, col = time)
}

# Extinction
plot(mean(sp_richness[1,1:10]), mean(extinction_t[1,1:10]), pch = 16, col = rgb(230, 159, 0, max = 255), xlim = c(0, max(sp_richness)), ylim = c(0,max(origination_t, na.rm = T)))
for(time in selected_time){
	points(mean(sp_richness[time,1:10]), mean(extinction_t[time,1:10]), pch = 16, col = time)
}

plot(mean(sp_richness[1,11:20]), mean(extinction_t[1,11:20]), pch = 16, col = rgb(0, 158, 115, max = 255), xlim = c(0, max(sp_richness)), ylim = c(0,max(origination_t, na.rm = T)))
for(time in selected_time){
	points(mean(sp_richness[time,11:20]), mean(extinction_t[time,11:20]), pch = 16, col = time)
}
