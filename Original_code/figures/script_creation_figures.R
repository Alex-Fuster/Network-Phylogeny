### Script to creat automatically figures ###
### Figures that are build by this script:
# mean species richness over time
# mean alpha value over time
# mean connectance over time
# mean species degree over time
# mean interaction range over time
# degree vs. branch length (age of species)
# centrality vs. branch length (age of species)
# MSD (PCoA) of phylogenetic trees topology
## All these figures are made for positive and negative interactions

################################################################################
############################# Data Extraction ##################################
################################################################################

#install.packages("ape")

source("code/functions/functions_M-C_bifurcation.R")
source("code/functions/tree_functions.R")
source("code/functions/phylo_functions.R")
source("code/functions/add_extinction_function.R")
source("code/functions/node_dist_function.R")

# Load output from the simulations
load(file = "Data/Positive_interactions/list_res_pos_100.Rdata")
list_res_pos = list_res
load(file = "Data/Negative_interactions/list_res_neg_100.Rdata")
list_res_neg = list_res

list_res <- c(list_res_pos, list_res_neg)
total_number_sim <- length(list_res)

# Pourcentage of simulation you want to take to make the figures
pourcentage <- 100

if(pourcentage != 100){
	# Random sampling of simulations
	sim_pos <- sample(1:(total_number_sim/2), pourcentage)
	sim_neg <- sample(1:(total_number_sim/2), pourcentage)
	list_res <- c(list_res_pos[sim_pos], list_res_neg[sim_neg])
}

nsim <- length(list_res)
nsteps <- 350

sp_max <- list_res[[1]][3]$parameters$Smax
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
		print(paste("Not enough species for the", list_res[[sim]]$simulation_number, "of interaction type =", list_res[[1]]$parameters$int,"(1 = positive, 0 = negative)"))
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

source("code/figures/mean_plot_function.R")

path <- paste0("Figures/basic_figures_", pourcentage, ".png")

png(path)
#par(mfrow = c(3,2))
layout(matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), 2, byrow = TRUE))

MeanPlotPosNeg(com_dynamic, "species richness")
MeanPlotPosNeg(alpha_stats, "alpha value")
MeanPlotPosNeg(connectance, "connectance")
MeanPlotPosNeg(int_range, "niche range")
MeanPlotPosNeg(mean_degree, "of mean degree")

dev.off()


################################################################################
################### FIGURE 2 - Degrees and the rest ############################
################################################################################

# Différenciate point between positive and negative interactions communities values
midd <- length(degree_in_tot)/2

# Positive interactions communities
# Create a table to store values of degree in (first row), degree out (second row) and related species age (third row)
degree_age_pos <- rbind(degree_in_tot[1:midd], degree_out_tot[1:midd], age_tot[1:midd])
## Remove the ones with zero interactions
to_keep_pos <- which(degree_age_pos[3,] == 1 | degree_age_pos[3,] > 1)
degree_age_pos <- degree_age_pos[,to_keep_pos]
#test <- degree_age_pos[, !duplicated(t(degree_age_pos))]

# Negative interactions communities
# Create a table to store values of degree in (first row), degree out (second row) and related species age (third row)
degree_age_neg <- rbind(degree_in_tot[(midd+1):length(degree_in_tot)], degree_out_tot[(midd+1):length(degree_in_tot)], age_tot[(midd+1):length(degree_in_tot)])
## Remove the ones with zero interactions
to_keep_neg <- which(degree_age_neg[3,] == 1 | degree_age_neg[3,] > 1)
degree_age_neg <- degree_age_neg[,to_keep_neg]

# Plot age vs. degree
path <- paste0("Figures/age-degree_figures_", pourcentage, ".png")
png(path)

par(mfrow= c(1,1), mar = c(3,3,0.5,0.5), mgp = c(1.5, 0.3, 0), tck = -.008, las = 1)

plot(degree_age_pos[1,], degree_age_pos[3,], col = rgb(0, 158, 115, max = 255, alpha = 75), pch = 16, xlim = c(0,40), ylim = c(0,80), xlab = "Interactions number", ylab = "Brench length")
points(degree_age_neg[1,], degree_age_neg[3,], col = rgb(230, 159, 0, max = 255, alpha = 75), pch = 16, xlim = c(0,50), ylim = c(0,80))

dev.off()

# Figure with separated positive and negative simulations
path <- paste0("Figures/age-degree_figures_pos-neg_", pourcentage, ".png")
png(path)

par(mfrow= c(2,1), mar = c(3,3,0.5,0.5), mgp = c(1.5, 0.3, 0), tck = -.008, las = 1)

# Only positive interactions simulations
plot(degree_age_pos[1,], degree_age_pos[3,], col = rgb(0, 158, 115, max = 255, alpha = 50), pch = 16, xlim = c(0,35), ylim = c(0,80), cex = 1.5, xlab = "Interactions number", ylab = "Species age", las = 1)
# Only negative interactions simulations
plot(degree_age_neg[1,], degree_age_neg[3,], col = rgb(230, 159, 0, max = 255, alpha = 50), pch = 16, xlim = c(0,35), ylim = c(0,80), cex = 1.5, xlab = "Interactions number", ylab = "Species age", las = 1)

dev.off()


# test of what ?
test_pres <- list()

for(i in 1:nsim){
	v <- (nsteps*sp_max) * i
	pres <- matrix(age_tot[(v-((nsteps*sp_max)-1)):v], nrow = nsteps, ncol = sp_max, byrow = FALSE)
	pres[pres != 0] <- 1
	pres[pres == 0] <- NA
	test_pres <- list_append(test_pres, pres)
}

test <- list()

for(i in 1:nsim){
	v <- (nsteps*sp_max) * i
	deg <- matrix(degree_in_tot[(v-((nsteps*sp_max)-1)):v], nrow = nsteps, ncol = sp_max, byrow = FALSE)
	deg_pres <- deg * test_pres[[i]]
	test <- list_append(test, deg_pres)
}

for(i in 1:length(test)){
	Sys.sleep(0.5)
	mean_degree_time <- rowMeans(test[[i]], na.rm = TRUE)
	plot(mean_degree_time, type = "l")
}


test <- test[-sim_to_remove]

# Do we have some groups into the way degree are evolving over time?

medoid_groups <- pam(test, 2, metric = "euclidean", stand = FALSE)


################################################################################
################# Nombre d'interactions total à l'équilibre ####################
################################################################################

## In degree
eq_time <- 175

n_sp_pos <- c(rep(NA, nsim))
n_interact_pos <- c(rep(NA, nsim))

n_sp_neg <- c(rep(NA, nsim))
n_interact_neg <- c(rep(NA, nsim))

degree_eq_pos <- NULL
age_eq_pos <- NULL

degree_eq_neg <- NULL
age_eq_neg <- NULL


for(sim in 1:length(list_res_pos)){
	if(!(sim %in% sim_to_remove)){
		pres <- list_res_pos[[sim]]$presence_matrix
		trait_mat <- list_res_pos[[sim]]$traits_df

		n_sp_pos[sim] <- sum(pres[eq_time,])

		pres_vec = pres[eq_time,]
		cooc = matrix(pres_vec, nr=sp_max, nc=sp_max, byrow=TRUE) * matrix(pres_vec, nr=sp_max, nc=sp_max, byrow=FALSE)
		L = get_L_mat(list_res_pos[[i]]$parameters, trait_mat)
		L = L*cooc
		in_degree <- rowSums(L)
		n_interact_pos[i] <- mean(in_degree[which(in_degree != 0)])

		sp_pres_eq <- pres[eq_time,]
		sp_pres <- which(sp_pres_eq != 0)
		age_sp_pres <- colSums(pres[1:eq_time, sp_pres])
		age_eq_pos <- append(age_eq_pos, age_sp_pres)
		degree_eq_pos <- append(degree_eq_pos, in_degree[which(sp_pres_eq != 0)])
	}

  pres <- list_res_neg[[sim]]$presence_matrix
  trait_mat <- list_res_neg[[sim]]$traits_df

  n_sp_neg[sim] <- sum(pres[eq_time,])

	pres_vec <- pres[eq_time,]
	cooc <- matrix(pres_vec, nr=sp_max, nc=sp_max, byrow=TRUE) * matrix(pres_vec, nr=sp_max, nc=sp_max, byrow=FALSE)
	L <- get_L_mat(list_res_pos[[i]]$parameters, trait_mat)
	L <- L*cooc
	in_degree <- rowSums(L)
  n_interact_neg[sim] <- mean(in_degree[which(in_degree != 0)])

	sp_pres_eq <- pres[eq_time,]
	sp_pres <- which(sp_pres_eq != 0)
	age_sp_pres <- colSums(pres[1:eq_time, sp_pres])
	age_eq_neg <- append(age_eq_neg, age_sp_pres)
	degree_eq_neg <- append(degree_eq_neg, in_degree[which(sp_pres_eq != 0)])
}

path <- paste0("Figures/interaction_nb_at_eq_", pourcentage, ".png")
png(path)

par(mfrow= c(1,1), mar = c(3,3,0.5,0.5), mgp = c(1.5, 0.3, 0), tck = -.008, las = 1)

plot(n_interact_pos, n_sp_pos, col = rgb(0, 158, 115, max = 255), pch = 16, xlim = c(0,20),  ylim = c(20,50), ylab = "Species richness", xlab = "Mean interactions number")
points(n_interact_neg, n_sp_neg, col = rgb(230, 159, 0, max = 255), pch = 16, ylim = c(20,50), ylab = "Species richness", xlab = "Mean interactions number")

dev.off()


# Species age at the equilibrium
x11()
plot(degree_eq_pos, age_eq_pos, col = rgb(0, 158, 115, max = 255), pch = 16)
points(degree_eq_neg, age_eq_neg, col = rgb(230, 159, 0, max = 255), pch = 16)

age_at_eq <- c(age_eq_pos, age_eq_neg)
interaction_type <- c(rep("Positive", length(age_eq_pos)), rep("Negative", length(age_eq_neg)))
data <- data.frame(age_at_eq, interaction_type)


# Violin plot of species age at the equilibrium
library(vioplot)

path <- paste0("Figures/sp_age_at_eq_", pourcentage, ".png")
png(path)

with(data , vioplot(
  age_at_eq[interaction_type=="Positive"] , age_at_eq[interaction_type=="Negative"],
  col = c(adjustcolor(rgb(0, 158, 115, max = 255), alpha = 0.5), adjustcolor(rgb(230, 159, 0, max = 255), alpha = 0.5)),
  border = c(rgb(0, 158, 115, max = 255), rgb(230, 159, 0, max = 255)),
  lineCol = c(rgb(0, 158, 115, max = 255), rgb(230, 159, 0, max = 255)),
  rectCol = c(rgb(0, 158, 115, max = 255), rgb(230, 159, 0, max = 255)),
  names=c("Positive interactions", "Negative interactions"),
	ylab="Species age"
)
)

dev.off()

path <- paste0("Figures/sp_degree_age_at_eq_", pourcentage, ".png")
png(path)

layout(matrix(c(1,0,1,3,2,3,2,0), 4, byrow = TRUE))
# figure age_degree pour pos communities
# figure age_degree pour neg communities
# violin plot de l'age des espèces a l'equilibre

dev.off()


# Violin plot age of species at different time steps

time_steps <- seq(10, 160, by = 10)

mean_age_time <- matrix(NA, nrow = length(time_steps), ncol = length(list_res_pos))

for(sim in 1:length(list_res_pos)){
	row_count <- 1
	for(time in time_steps){
		if(!(sim %in% sim_to_remove)){
			# récupération des données de présences des espèces
			pres <- list_res_pos[[sim]]$presence_matrix
			pres_time <- pres[time, ]
			# Quelles espèce est présente et quel est son àge à la tine step i
			sp_pres <- which(pres_time == 1)
			if(length(sp_pres) > 1){
				sp_age <- colSums(pres[1:time, sp_pres])
				# Moyenne de cet âge
				mean_age_time[row_count, sim] <- mean(sp_age)
			} else {
				sp_age <- sum(pres[1:time, sp_pres])
				mean_age_time[row_count, sim] <- sp_age
			}
			row_count <- row_count + 1
		}
	}
}

age_sp_pos <- rowMeans(mean_age_time, na.rm = TRUE)
plot(age_sp_pos, type = "l", xlim = c(0,16), ylim = c(0,20))

library(tidyr)

data <- data.frame(time_steps, mean_age_time)
data_long <- gather(data, simulation, values, X1:X100)

vioplot(
  data_long[time_steps=="10",],
  col = c(adjustcolor(rgb(0, 158, 115, max = 255), alpha = 0.5)),
  border = c(rgb(0, 158, 115, max = 255)),
  lineCol = c(rgb(0, 158, 115, max = 255)),
  rectCol = c(rgb(0, 158, 115, max = 255)),
  names=c("Pegative interactions"),
	ylab="Species age"
)


mean_age_time <- matrix(NA, nrow = length(time_steps), ncol = length(list_res_neg))

for(k in 1:length(list_res_neg)){
	row_count <- 1
	for(i in time_steps){
		# récupération des données de présences des espèces
		pres <- list_res_neg[[k]]$presence_matrix
		pres_time <- pres[i, ]
		# Quelles espèce est présente et quel est son àge à la tine step i
		sp_pres <- which(pres_time == 1)
		if(length(sp_pres) > 1){
			sp_age <- colSums(pres[1:i, sp_pres])
			# Moyenne de cet âge
			mean_age_time[row_count, k] <- mean(sp_age)
		} else {
			sp_age <- sum(pres[1:i, sp_pres])
			mean_age_time[row_count, k] <- sp_age
		}
		row_count <- row_count + 1
	}
}

plot(rep(1, 100), mean_age_time[1,])
age_sp_neg <- rowMeans(mean_age_time, na.rm = TRUE)
x11()
plot(age_sp_neg, type = "l", xlim = c(0,16), ylim = c(0,20))


###

dens_pos <- density(age_eq_pos, na.rm = TRUE)
dens_neg <- density(age_eq_neg, na.rm = TRUE)

plot(dens_pos, col = rgb(0, 158, 115, max = 255), xlim = c(0, 60), ylim = c(0,0.15), xlab = "Species age", main = "")
lines(dens_neg, col = rgb(230, 159, 0, max = 255))

##
eq_time <- 150

degree_eq <- NA
age_eq <- NA

for(i in 1:length(list_res_pos)){
	if(!(i %in% sim_to_remove)){
		pres <- list_res_pos[[i]]$presence_matrix
		sp_pres_eq <- pres[eq_time,]
		sp_pres <- which(sp_pres != 0)
		age_sp_pres <- colSums(pres[1:eq_time, sp_pres])
	}
}

cell_position <- function(n_row, x, y) n_row * (x - 1) + y
position_cell <- function(n_row, x){
	col_nb <- ceiling(x/n_row)
	row_nb <- x - (n_row * (col_nb-1))
	return(c(row_nb,col_nb))
}

test_position <- which(pres == 1)
test_coord_cell <- mapply(position_cell, 350, test_position) #[1,] => time, [2,] => species
plot(test_coord_cell[1,], rev(test_coord_cell[2,]))

test_coord_cell[,1:10]

#species age

age_sp <- NULL
sp_arrival <- NULL

for(i in 1:length(list_res_pos)){
	if(!(i %in% sim_to_remove)){
		pres <- list_res_pos[[i]]$presence_matrix
		age_sp <- append(age_sp, colSums(pres))

		anc <- list_res_pos[[i]]$parentage_matrix
		sp_arrival <- append(sp_arrival, anc[,1])
	}
}

#CHECk for NA
#arrival_age[which(is.na(sp_arrival)),]

#for(i in 1:length(list_res_pos)){
#	pres <- list_res_pos[[i]]$presence_matrix
#	test_first_sp <- colSums(pres)
#	print(test_first_sp[1])
#}


sp_arrival[which(is.na(sp_arrival))] <- 1
arrival_age <- as.data.frame(cbind(sp_arrival, age_sp))
test_aa <- arrival_age[order(sp_arrival, age_sp),]

par(mfrow = c(1,2))
plot(test_aa[,1],test_aa[,2], col = rgb(0, 158, 115, max = 255, alpha = 75), pch = 16)


age_sp <- NULL
sp_arrival <- NULL

for(sim in 1:length(list_res_neg)){
	pres <- list_res_neg[[sim]]$presence_matrix
	age_sp <- append(age_sp, colSums(pres))

	anc <- list_res_neg[[sim]]$parentage_matrix
	sp_arrival <- append(sp_arrival, anc[,1])
}

sp_arrival[which(is.na(sp_arrival))] <- 1
arrival_age <- as.data.frame(cbind(sp_arrival, age_sp))
test_aa <- arrival_age[order(sp_arrival, age_sp),]
plot(test_aa[,1],test_aa[,2], col = rgb(230, 159, 0, max = 255, alpha = 75), pch = 16)


################################################################################
######### FIGURE 3 : Distance between species in the niche space ###############
################################################################################

unselected_sim <- matrix(NA, nrow = 200, ncol = length(selected_time))
col_unselected_sim <- 1

path <- paste0("Figures/Dist_btw_sp", pourcentage, ".png")

png(path)

par(mfrow=c(2,round(length(selected_time)/2))) #Partition de la figure en 2 rangs et n colonnes (nombre de colonne relative aux nombre de pas de temps observé)

for(time in selected_time){
	traits_eq_list <- list()

	for(sim in 1:length(list_res)){ #for each simulation
	  #print(paste("simulation", i, "on", nsim))
		# Harvest selected parameters (species presence and traits)
	  parameters <- as.list(unlist(list_res[[sim]][3]))
	  names(parameters) <- gsub("parameters.","",names(parameters))
	  pres <- matrix(unlist(list_res[[sim]][4]), nr = nsteps, nc = parameters$Smax)
		pres_eq <- pres[time, ]
	  traits <- as.data.frame(list_res[[sim]][5])
	  names(traits) <- gsub("traits_df.","",names(traits))

		# Distance between species niche optimum trait value
		traits_to_keep <- sort(traits[which(pres_eq == 1), 1])
		delta_traits <- diff(traits_to_keep)

		# Frenquency of distance between species niche optimum trait value
		if(length(delta_traits)>1){
			freq_dist_traits <- density(delta_traits)
		} else {
			unselected_sim[i, col_unselected_sim] <- 1
		}
		traits_eq_list <- list_append(traits_eq_list, freq_dist_traits)
	}

	# Count of removed simulations - for the plot
	line_to_remove_pos <- sum(unselected_sim[1:100,col_unselected_sim], na.rm = TRUE)
	line_to_remove_neg <- sum(unselected_sim[101:200,col_unselected_sim], na.rm = TRUE)

	col_unselected_sim <- col_unselected_sim + 1

	mean_distrib_pos <- mean(traits_eq_list)

	# Plot
	plot(traits_eq_list[[1]], type = "l", xlim = c(0,0.5), ylim = c(0,43), col = rgb(77,146,33, max = 255), main = paste("t =", t), xlab = "Distance between traits", ylab = "Density")
	for(i in 2:(100-line_to_remove_pos)){
		lines(traits_eq_list[[i]], col = rgb(77,146,33, max = 255))
	}
	for(i in (101-line_to_remove_pos):length(traits_eq_list)){
		lines(traits_eq_list[[i]], col = rgb(250, 159, 0, max = 255))
	}
}

dev.off()

### Répartition des valeurs des traits des communautés à l'équilibre

par(mfrow=c(2,round(length(selected_time)/2))) #Partition de la figure en 2 rangs et n colonnes (nombre de colonne relative aux nombre de pas de temps observé)

diff_to_one_pos <- matrix(NA, nrow = 1, ncol = length(selected_time))
diff_to_one_neg <- matrix(NA, nrow = 1, ncol = length(selected_time))
iterateur <- 0

# Une figure par pas de temps sélectionné (voir en haut du script)
for(t in selected_time){
	traits_eq <- matrix(NA, nr = length(list_res), nc = 50) #Objet pour stocker les traits de façon temporaire
	iterateur <- iterateur+1

	# Récupération des valeurs de traits pour le pas de temps t
	for(i in 1:length(list_res)){
	  #print(paste("simulation", i, "on", nsim))
		## Renommer les paramètres
	  parameters <- as.list(unlist(list_res[[i]][3]))
	  names(parameters) <- gsub("parameters.","",names(parameters))
	  pres <- matrix(unlist(list_res[[i]][4]), nr = nsteps, nc = parameters$Smax) #Récupération des espèces présentes dans la simulation
		pres_eq <- pres[t, ] # Sélection du vecteur d'espèce présentes au pas de temps désiré
	  traits <- as.data.frame(list_res[[i]][5]) #Récupération du data frame contenant les traits des espèces de la simulation
	  names(traits) <- gsub("traits_df.","",names(traits)) #Changement des noms de colonne

		to_keep <- sort(traits[which(pres_eq == 1), 1]) #Sélection des valeurs de trait n des espèces présentes au temps t
		traits_eq[i, 1:length(to_keep)] <- to_keep #insertion des valeurs de trait dans la matrice temporaire

	}

	# Mean trait value for positive interactions communities
	mean_traits_pos <- colMeans(traits_eq[1:100,], na.rm = TRUE)

	sp_tot <- sum(!is.na(mean_traits_pos))
	y <- seq(0, 1, length.out=sp_tot)

	diff_to_one_pos[iterateur] <- sum(mean_traits_pos[1:sp_tot] - y)
	print(paste("diff pos", sum(mean_traits_pos[1:sp_tot] - y)))

	plot(mean_traits_pos, type = "l", ylim = c(0,1),col = "red")
	lines(y, lty = 2, col = "red")

	mean_traits_neg <- colMeans(traits_eq[101:200,], na.rm = TRUE)
	lines(mean_traits_neg, type = "l")
	sp_tot <- sum(!is.na(mean_traits_neg))
	y <- seq(0, 1, length.out=sp_tot)
	lines(y, lty = 2)

	diff_to_one_neg[iterateur] <- sum(mean_traits_neg[1:sp_tot] - y)
	print(paste("diff neg", sum(mean_traits_neg[1:sp_tot] - y)))
}

path <- paste0("Figures/Niche_occupancy", pourcentage, ".png")

png(path)

par(mfrow=c(1,1))

plot(diff_to_one_pos[1,], type = "l", col = rgb(0, 158, 115, max = 255), xaxt = "n", xlab = "Time", ylim = c(0,10), ylab = "`Distance` à la courbe 1:1")
axis(1, at = 1:length(selected_time), labels = selected_time)
lines(diff_to_one_neg[1,], col = rgb(230, 159, 0, max = 255))

dev.off()


# Spatial analysis of niche trait values

library(spatstat)

#plot(1, type="n", xlab="species sorted by nearest neighbor", ylab="mean nearest neighbor", xlim=c(0, 50), ylim=c(0, 1))
#line_type = 1

pseudo_pvalue <- matrix(NA, nrow = nsim, ncol = length(selected_time))
time_pos <- 1

for(time in selected_time){
	print(time)
	# Récupération des valeurs de traits pour le pas de temps t

	nearest_neighbor <- matrix(NA, nrow = 50, ncol = length(list_res))

	for(sim in 1:length(list_res)){
		## Renommer les paramètres
	  parameters <- as.list(unlist(list_res[[sim]][3]))
	  names(parameters) <- gsub("parameters.","",names(parameters))
	  pres <- matrix(unlist(list_res[[sim]][4]), nr = nsteps, nc = parameters$Smax) #Récupération des espèces présentes dans la simulation
		pres_select <- pres[time, ] # Sélection du vecteur d'espèce présentes au pas de temps désiré
	  traits <- as.data.frame(list_res[[sim]][5]) #Récupération du data frame contenant les traits des espèces de la simulation
	  names(traits) <- gsub("traits_df.","",names(traits)) #Changement des noms de colonne

		to_keep <- sort(traits[which(pres_select == 1), 1]) #Sélection des valeurs de trait n des espèces présentes au temps t
		traits_pres <- to_keep #insertion des valeurs de trait dans la matrice temporaire

		y <- rep(0, length(traits_pres))
		nearest_neighbor[1:length(traits_pres), sim] <- sort(nndist(traits_pres, y, k = 1, by = NULL))

		# Test for clustering/dispersion
		ann.p <- mean(nndist(traits_pres, y, k=1))
		ann.p

		n <- 999               # Number of simulations
		ann.r <- vector(length = n) # Create an empty object to be used to store simulated ANN values
		for (ite in 1:n){
  		rand.p   <- rpoint(n=length(traits_pres))  # Generate random point locations ; n = number of points
			rand.p$y <- rep(0, length(rand.p$y)) # Make the point in one dimension
  		ann.r[ite] <- mean(nndist(rand.p, k=1))  # Tally the ANN values
		}

		# Computi a pseudo p-value from the simulation
		N.greater <- sum(ann.r > ann.p) #find the number of simulated ANN values greater than our observed ANN value
		pseudo_pvalue[sim, time_pos] <- min(N.greater + 1, n + 1 - N.greater) / (n +1)
	}

	time_pos <- time_pos + 1
}

head(pseudo_pvalue)
plot(pseudo_pvalue[,6])
plot(rep(1, 100), pseudo_pvalue[1:100,1])
plot(rep(1, 100), pseudo_pvalue[101:200,1])

plot(c(1:100), sort(pseudo_pvalue[1:100,1]))
plot(c(1:100), sort(pseudo_pvalue[101:200,1]))

plot(c(1:100), sort(pseudo_pvalue[1:100,2]))
plot(c(1:100), sort(pseudo_pvalue[101:200,2]))

plot(c(1:100), sort(pseudo_pvalue[1:100,3]))
plot(c(1:100), sort(pseudo_pvalue[101:200,3]))

plot(c(1:100), sort(pseudo_pvalue[1:100,4]))
plot(c(1:100), sort(pseudo_pvalue[101:200,4]))

plot(c(1:100), sort(pseudo_pvalue[1:100,5]))
plot(c(1:100), sort(pseudo_pvalue[101:200,5]))


################################################################################
################# Distribution des espèces sans déscendance ####################
################################################################################

no_desc_sp <- rep(NA,length(list_res))
id_no_desc_sp <- matrix(NA, nrow = length(list_res), ncol = 500)
desc_sp <- rep(NA,length(list_res))

for(sim in 1:length(list_res)){
	print(sim)
  anc <- matrix(unlist(list_res[[sim]][6]), nr = parameters$Smax, nc = 3) #time, ancestor, new species
  extinction <- matrix(unlist(list_res[[sim]][7]), nr = parameters$Smax, nc = 2) #time, species extinct

	no_desc_sp[sim] <- length(setdiff(extinction[,2], anc[,2]))
	id_no_desc_sp[sim, 1:no_desc_sp[sim]] <- setdiff(extinction[,2], anc[,2])
	desc_sp[sim] <- length(unique(anc[,2]))
}

plot(no_desc_sp[1:(length(list_res)/2)], xlim = c(1,length(list_res)/2), ylim = c(1,max(no_desc_sp)), col = rgb(0, 158, 115, max = 255))
points(no_desc_sp[(length(list_res)/2 + 1):length(list_res)], col = rgb(230, 159, 0, max = 255))

plot(desc_sp[1:(length(list_res)/2)], xlim = c(1,length(list_res)/2), ylim = c(1,max(desc_sp)), col = rgb(0, 158, 115, max = 255))
points(desc_sp[(length(list_res)/2 + 1):length(list_res)], col = rgb(230, 159, 0, max = 255))

################################################################################
########################### Simulations Stats ##################################
################################################################################

load(file = "Data/Positive_interactions/needed_sim_pos.Rdata")
load(file = "Data/Negative_interactions/needed_sim_neg.Rdata")

load(file = "Data/Positive_interactions/reccorded_seed_pos_list.Rdata")
load(file = "Data/Negative_interactions/reccorded_seed_neg_list.Rdata")

nsim_test <- 10

pourcentage_of_success <- function(success_nb, total){
	return(success_nb * 100 / total)
}

success_pourc_pos <- mapply(pourcentage_of_success, nsim_test, needed_sim_pos)
mean_success_pos <- mean(success_pourc_pos)
sd_success_pos <- sd(success_pourc_pos)

success_pourc_neg <- mapply(pourcentage_of_success, nsim_test, needed_sim_neg)
mean_success_neg <- mean(success_pourc_neg)
sd_success_neg <- sd(success_pourc_neg)

################################################################################
###################### Plot MDS on tree topography #############################
################################################################################

#install.packages("treespace")

class(tree_space_an) <- "multiPhylo"
test_treespace <- treespace::treespace(tree_space_an, nf = 3, lambda = 0) # based on topography (Lambda = 0)

names_trees <- c(1:nsim)
names(tree_space_an) <- names_trees[-sim_to_remove]
#plotGrovesD3(test_treespace$pco, treeNames = names(tree_space_an))#, color.palette=viridis)

treespace_group <- treespace::findGroves(test_treespace, nclust=2)
treespace_group$groups <- c(rep("positive",98), rep("negative",100))

# Plot de la MDS sur la topographie des arbres, couleurs en fonction du type d'interaction
plotGrovesD3(treespace_group, colors = c(rgb(0, 158, 115, max = 255), rgb(250, 159, 0, max = 255)), col_lab="Interactions")

# Finding cluster into the cloud of trees
test_treespace_dist <- test_treespace$D
library(cluster)
test_pam <- pam(as.dist(test_treespace_dist), 2, metric = "euclidian", stand = FALSE)
test_pam$interact_groups <- c(rep(1,98), rep(2,100))

# Plot de la MDS sur la topographie des arbres, couleurs en fonction du type d'interaction, synbole en fonction du groupe
plotGrovesD3(treespace_group, colors = c(rgb(0, 158, 115, max = 255), rgb(250, 159, 0, max = 255)), col_lab="Interactions",
	symbol_var = test_pam$cluster, symbols = c("wye", "circle"))
#point types : "circle", "cross", "diamond", "square", "star", "triangle", and "wye"

################################################################################

#if(pars$int == 1){
#  color_int = c(rgb(250, 159, 0, max = 255))0
#} else {
#  color_int = c(rgb(0, 158, 115, max = 255))
#}

################################################################################
############## Correlation trait distance - phylogenetic distance ##############
################################################################################

corr_pos <- read.csv("Data/Positive_interactions/corr_pd-td_pos.csv")
corr_pos <- corr_pos[,-1]
corr_neg <- read.csv("Data/Negative_interactions/corr_pd-td_neg.csv")
corr_neg <- corr_neg[,-1]

mean_distance <- colMeans(corr_pos, na.rm = TRUE)
plot(selected_time, mean_distance, type = "l", ylim = c(0,1), lty = 2, col = rgb(77,146,33, max = 255), xlab = "Time")
#Add quantils
quantiles <- apply(corr_pos[1:nsim,], 2, quantile, probs = c(0.05, 0.95), na.rm = T)
polygon(c(selected_time,rev(selected_time)), c(rep(0,length(selected_time)), rev(quantiles[2,])), col = adjustcolor(rgb(161,215,106, max = 255), alpha = 0.5), border = NA)
lines(selected_time, quantiles[2,], col = adjustcolor(rgb(161,215,106, max = 255), alpha = 0.5))
title("Correlation between phylogenetic and trait distances for positive and negative interactions")

mean_distance <- colMeans(corr_neg, na.rm = TRUE)
#plot(time, mean_distance, type = "l", ylim = c(0,1), lty = 2, col = rgb(197,27,125, max = 255), xlab = "Time")
lines(selected_time, mean_distance, lty = 2, col = rgb(250, 159, 0, max = 255))
#Add quantils
quantiles <- apply(corr_neg[1:nsim,], 2, quantile, probs = c(0.05, 0.95), na.rm = T)
polygon(c(selected_time,rev(selected_time)), c(rep(0,length(selected_time)), rev(quantiles[2,])), col = adjustcolor(rgb(250, 159, 0, max = 255), alpha = 0.2), border = NA)
lines(selected_time, quantiles[2,], col = adjustcolor(rgb(250, 159, 0, max = 255), alpha = 0.4))
#title("Correlation between phylogenetic and trait distances for negative interactions")

## Do some groups with similar simulation
## For positive interactions
# Create distance matrix between simulations
dist_corr_pos <- dist(corr_pos, method = "euclidean")
corr_pos_hclust <- hclust(dist_corr_pos, method = "complete")
plot(corr_pos_hclust)
simsim_groups <- cutree(corr_pos_hclust, k = 4)

library(vegan)
corrpos_dist <- vegdist(corr_pos, method = "euclidean", na.rm = TRUE)
corrpos_dist_complete <- hclust(corrpos_dist, method = "complete")
plot(corrpos_dist_complete)
sim_simpos_groups <- cutree(corrpos_dist_complete, k = 6)

corrneg_dist <- vegdist(corr_neg, method = "euclidean", na.rm = TRUE)
corrneg_dist_complete <- hclust(corrneg_dist, method = "complete")
plot(corrneg_dist_complete)
sim_simneg_groups <- cutree(corrneg_dist_complete, k = 4)
sim_simneg_groups <- sim_simneg_groups + 6

similar_sim_groups <- c(sim_simpos_groups, sim_simneg_groups)
similar_sim_groups <- similar_sim_groups[-sim_to_remove]
treespace_group$groups <- as.factor(similar_sim_groups)
interaction_type <- c(rep("positive",98), rep("negative",100))

# Plot de la MDS sur la topographie des arbres, couleurs en fonction du type d'interaction, synbole en fonction du groupe
plotGrovesD3(treespace_group,
	#colors = c(rgb(84,48,5, max = 255), rgb(140,81,10, max = 255), rgb(191,129,45, max = 255), rgb(223,194,125, max = 255), rgb(246,232,195, max = 255),
	#rgb(199,234,229, max = 255), rgb(128,205,193, max = 255), rgb(53,151,143, max = 255), rgb(1,102,94, max = 255), rgb(0,60,48, max = 255)),
	colors = c(rgb(165,0,38, max = 255), rgb(215,48,39, max = 255),rgb(244,109,67, max = 255),rgb(253,174,97, max = 255),rgb(254,224,139, max = 255),
	rgb(217,239,139, max = 255),rgb(166,217,106, max = 255),rgb(102,189,99, max = 255),rgb(26,152,80, max = 255),rgb(0,104,55, max = 255)),
	col_lab="groups",
	symbol_var = interaction_type, symbols = c("wye", "triangle"))

# Plot de la MDS sur la topographie des arbres, couleurs en fonction du groupe de correlation entre trait distance and phylogeny distance
# Positive interactions

similar_sim_groups <- c(sim_simpos_groups, rep(7, 100))
similar_sim_groups <- similar_sim_groups[-sim_to_remove]
treespace_group$groups <- as.factor(similar_sim_groups)
interaction_type <- c(rep("positive",98), rep("negative",100))

plotGrovesD3(treespace_group,
	#colors = c(rgb(84,48,5, max = 255), rgb(140,81,10, max = 255), rgb(191,129,45, max = 255), rgb(223,194,125, max = 255), rgb(246,232,195, max = 255),
	#rgb(199,234,229, max = 255), rgb(128,205,193, max = 255), rgb(53,151,143, max = 255), rgb(1,102,94, max = 255), rgb(0,60,48, max = 255)),
	colors = c(rgb(199,233,180, max = 255),rgb(127,205,187, max = 255),rgb(65,182,196, max = 255),rgb(29,145,192, max = 255),rgb(34,94,168, max = 255),rgb(37,52,148, max = 255),"lightgrey"),
	col_lab="groups")

# Plot de la MDS sur la topographie des arbres, couleurs en fonction du groupe de correlation entre trait distance and phylogeny distance
# Negative interactions

similar_sim_groups <- c(rep(5, 100), sim_simneg_groups)
similar_sim_groups <- similar_sim_groups[-sim_to_remove]
treespace_group$groups <- as.factor(similar_sim_groups)
interaction_type <- c(rep("positive",98), rep("negative",100))

plotGrovesD3(treespace_group,
	#colors = c(rgb(84,48,5, max = 255), rgb(140,81,10, max = 255), rgb(191,129,45, max = 255), rgb(223,194,125, max = 255), rgb(246,232,195, max = 255),
	#rgb(199,234,229, max = 255), rgb(128,205,193, max = 255), rgb(53,151,143, max = 255), rgb(1,102,94, max = 255), rgb(0,60,48, max = 255)),
	colors = c("lightgrey", rgb(199,233,180, max = 255),rgb(127,205,187, max = 255),rgb(65,182,196, max = 255),rgb(29,145,192, max = 255)),
	col_lab="groups")

###
pca_corrpos <- prcomp(na.omit(corr_pos))
biplot(pca_corrpos, col=simsim_groups)

raw <- pca_corrpos$x[,1:2]
plot(raw[,1], raw[,2], col=simsim_groups, pch=20)

#K-medoids method to find groups
fviz_nbclust(corr_pos, pam, method = "silhouette",k.max=20)
medoid_groups <- pam(corr_pos, 2, metric = "euclidean", stand = FALSE)


# PAM (Partitioning Around Medoid)
library(cluster)
test_pam <- pam(corr_pos, 2, metric = "euclidean", stand = FALSE)

clusplot(test_pam, main = "Cluster plot, k = 2", color = TRUE)

plot(raw[,1], raw[,2], col=medoid_groups$clustering, pch=20)


fviz_pca_ind(pca_corrpos$x[,1:2],
             col.ind = medoid_groups$clustering, # colorer par groupes
             palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Ellipse de concentration
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
             )


################################################################################
###################### Pas de temps pour arriver au Smax #######################
################################################################################

time_to_Smax_pos <- c(rep(NA, nsim))
time_to_Smax_neg <- c(rep(NA, nsim))

for(sim in 1:length(list_res_pos)){
  anc <- list_res_pos[[sim]]$parentage_matrix
  time_to_Smax_pos[sim] <- anc[pars$Smax,1]
  anc <- list_res_neg[[sim]]$parentage_matrix
  time_to_Smax_neg[sim] <- anc[pars$Smax,1]
}

plot(rep(1, length(list_res_pos)),time_to_Smax_pos, col = rgb(77,146,33, max = 255), ylim = c(100, 300), xlim = c(0,3))
points(rep(2, length(list_res_pos)),time_to_Smax_neg, col = rgb(197,27,125, max = 255))

library(vioplot)

time_to_Smax <- c(time_to_Smax_pos, time_to_Smax_neg)
interaction_type <- c(rep("Positive", nsim), rep("Negative", nsim))
data <- data.frame(time_to_Smax, interaction_type)

with(data , vioplot(
  time_to_Smax[interaction_type=="Positive"] , time_to_Smax[interaction_type=="Negative"],
  col = c(adjustcolor(rgb(161,215,106, max = 255), alpha = 0.5), adjustcolor(rgb(233,163,201, max = 255), alpha = 0.5)),
  border = c(rgb(161,215,106, max = 255), rgb(233,163,201, max = 255)),
  lineCol = c(rgb(161,215,106, max = 255), rgb(233,163,201, max = 255)),
  rectCol = c(rgb(161,215,106, max = 255), rgb(233,163,201, max = 255)),
  names=c("Positive interactions", "Negative interactions")
))

dens_pos <- density(time_to_Smax_pos, na.rm = TRUE)
dens_neg <- density(time_to_Smax_neg, na.rm = TRUE)

plot(dens_pos, col = rgb(77,146,33, max = 255), xlim = c(125, 300), ylim = c(0,0.035), xlab = "Time", main = "")
lines(dens_neg, col = rgb(197,27,125, max = 255))


################################################################################
########################## Temps pour arriver à Smax ###########################
################################################################################

time_to_Smax_pos <- c(rep(NA, nsim/2))
time_to_Smax_neg <- c(rep(NA, nsim/2))

Smax <- list_res[[1]][3]$parameters$Smax

for(sim in 1:length(list_res_pos)){
  anc <- list_res_pos[[sim]]$parentage_matrix
  time_to_Smax_pos[sim] <- anc[Smax,1]
  anc <- list_res_neg[[sim]]$parentage_matrix
  time_to_Smax_neg[sim] <- anc[Smax,1]
}

library(vioplot)

time_to_Smax <- c(time_to_Smax_pos, time_to_Smax_neg)
interaction_type <- c(rep("Positive", nsim/2), rep("Negative", nsim/2))
data <- data.frame(time_to_Smax, interaction_type)

with(data , vioplot(
  time_to_Smax[interaction_type=="Positive"] , time_to_Smax[interaction_type=="Negative"],
  col = c(adjustcolor(rgb(0, 158, 115, max = 255), alpha = 0.5), adjustcolor(rgb(230, 159, 0, max = 255), alpha = 0.5)),
  border = c(rgb(0, 158, 115, max = 255), rgb(230, 159, 0, max = 255)),
  lineCol = c(rgb(0, 158, 115, max = 255), rgb(230, 159, 0, max = 255)),
  rectCol = c(rgb(0, 158, 115, max = 255), rgb(230, 159, 0, max = 255)),
  names=c("Positive interactions", "Negative interactions")
))

dens_pos <- density(time_to_Smax_pos, na.rm = TRUE)
dens_neg <- density(time_to_Smax_neg, na.rm = TRUE)

plot(dens_pos, col = rgb(77,146,33, max = 255), xlim = c(125, 300), ylim = c(0,0.035), xlab = "Time", main = "")
lines(dens_neg, col = rgb(197,27,125, max = 255))


################################################################################
######################### PCoA phylogenetic trees ##############################
################################################################################

library(igraph)

source("code/functions/list_append_function.R")

phylo_distances <- list()
interactions_matrices <- list()
pcoa_phylo_list <- list()
pcoa_node_list <- list()


for(sim in 1:length(list_res)){

	print(paste("simulation", sim, "on", nsim))
  parameters <- as.list(unlist(list_res[[sim]][3]))
  names(parameters) <- gsub("parameters.","",names(parameters))
  pres <- matrix(unlist(list_res[[sim]][4]), nr = nsteps, nc = parameters$Smax)
  traits <- as.data.frame(list_res[[sim]][5])
  names(traits) <- gsub("traits_df.","",names(traits))
  anc <- matrix(unlist(list_res[[sim]][6]), nr = parameters$Smax, nc = 3) #time, ancestor, new species
  extinction <- matrix(unlist(list_res[[sim]][7]), nr = parameters$Smax, nc = 2) #time, species extinct

	## Global distance and interaction matrices
  # Phylogenetic statistics
	species_total <- max(anc[,3], na.rm = T)

	# Calculate phylogenetic distance between species
	phylo_distance_mat <- get_dist_mat(anc, nsteps, species_total)
	max_phylo_dist <- max(phylo_distance_mat)
	distance_mat <- phylo_distance_mat / max_phylo_dist

	# PCoA on phylogenetic distances
	pcoa_phylo <- cmdscale(distance_mat,k = 2, eig=TRUE)

	# Record results on lists
	phylo_distances <- list_append(phylo_distances, distance_mat)
	pcoa_phylo_list <- list_append(pcoa_phylo_list, pcoa_phylo)

	# Interactions matrix
	pres_vec <- rep(1, species_total)
	cooc <- matrix(pres_vec, nr=parameters$Smax, nc=parameters$Smax, byrow=TRUE)*
	matrix(pres_vec, nr=parameters$Smax, nc=parameters$Smax, byrow=FALSE)
	adj_mat <- get_L_mat(parameters, traits)

	g <- graph_from_adjacency_matrix(adj_mat)
	node_dist <- distances(g)
	node_dist[node_dist == Inf] = 999
	max_node_dist <- max(node_dist)

	node_distance <- node_dist / max_node_dist

	pcoa_node <- cmdscale(node_distance, k = 2, eig=TRUE)

	# Record results on lists
	interactions_matrices <- list_append(interactions_matrices, adj_mat)
	pcoa_node_list <- list_append(pcoa_node_list, pcoa_node)

}

for(sim in 1:length(list_res)){

  parameters <- as.list(unlist(list_res[[sim]][3]))
  names(parameters) <- gsub("parameters.","",names(parameters))
  pres <- matrix(unlist(list_res[[sim]][4]), nr = nsteps, nc = parameters$Smax)

	phylo_dist <- phylo_distances[[sim]]
	node_dist <- interactions_matrices[[sim]]

	for (time in selected_time){

		pres_time <- pres[time, ]
		# Quelles espèce est présente et quel est son àge à la tine step i
		sp_pres <- which(pres_time == 1)

		# Phylogeny
		distance_mat <- phylo_dist[sp_pres,sp_pres]
		pcoa_phylo_time <- cmdscale(distance_mat,k = 2, eig=TRUE)

		# interactions
		node_distance_time <- node_dist[sp_pres,sp_pres]
		pcoa_node_time <- cmdscale(node_distance_time,k = 2, eig=TRUE)

	}

}


## PCoA plot

x <- pcoa_phylo_time$points[,1]
y <- pcoa_phylo_time$points[,2]
plot(x,y)

x2<- pcoa_node$points[,1]
y2 <- pcoa_node$points[,2]
plot(x2,y2)


g <- graph_from_adjacency_matrix(adj_mat)
g_small <- graph_from_adjacency_matrix(adj_mat[1:30,1:30])

test_dist <- distances(g_small)
test_path <- shortest.paths(g_small)

################################################################################
######################## Community dynamic figure ##############################
################################################################################

## Community dynamic matrix
com_dynamic <- matrix(nr = nsteps-1, nc = nsim)

# Data extraction
for(sim in 1:length(list_res)){
  print(paste("simulation", sim, "on", nsim))
  parameters <- as.list(unlist(list_res[[sim]][3]))
  names(parameters) <- gsub("parameters.","",names(parameters))
  pres <- matrix(unlist(list_res[[sim]][4]), nr = nsteps, nc = parameters$Smax)

	species_richness <- rowSums(pres)
	species_richness <- replace(species_richness, species_richness == 0, NA)
	com_dynamic[, i] <- species_richness[-350]
}

# Growth rate matrix
growth_rate <- matrix(nr = nsteps, nc = nsim)

# Growth rate calculation (r = (Nt+1 - Nt) / Nt)
for(j in 1:ncol(com_dynamic)){
	for(i in 2:(nrow(com_dynamic)-1)){
		growth_rate[i, j] <- (com_dynamic[i, j] - com_dynamic[i-1, j]) / com_dynamic[i-1, j]
	}
}

# Mean growth rates
mean_growth_rate_pos <- rowMeans(growth_rate[,1:100], na.rm = T)
mean_growth_rate_neg <- rowMeans(growth_rate[,101:200], na.rm = T)

# Plot
plot(mean_growth_rate_pos, type = "l", xlim = c(0,350), ylim = c(0,1))
lines(mean_growth_rate_neg,col = "red")


## Species richness and growth rate relationship
growth_rate[, 1]
com_dynamic[, 1]

mean_com_dynamic_pos <- rowMeans(com_dynamic[,1:100], na.rm = T)
mean_com_dynamic_neg <- rowMeans(com_dynamic[,101:200], na.rm = T)

par(mfrow = c(1,2))
plot(mean_growth_rate_pos[1:150], mean_com_dynamic_pos[1:150], xlim = c(-0.2,0.3), ylim = c(0,50))
plot(mean_growth_rate_neg[1:150], mean_com_dynamic_neg[1:150], xlim = c(-0.2,0.3), ylim = c(0,50))

plot(mean_com_dynamic_pos[1:150], mean_growth_rate_pos[1:150], xlim = c(0,50), ylim = c(-0.2,0.3))
plot(mean_com_dynamic_neg[1:150], mean_growth_rate_neg[1:150], xlim = c(0,50), ylim = c(-0.2,0.3))

par(mfrow = c(1,2))
plot(growth_rate[1:150, 1:100], com_dynamic[1:150, 1:100])
plot(growth_rate[1:150, 101:200], com_dynamic[1:150, 101:200])


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

mean_turnover_pos <- rowMeans(turnover_rate[1:150,1:100], na.rm = TRUE)
mean_turnover_neg <- rowMeans(turnover_rate[1:150,101:200], na.rm = TRUE)

source("code/figures/mean_plot_function.R")

path <- paste0("Figures/community_turnover_", pourcentage, ".png")
png(path)

par(mfrow=c(2,1))

MeanPlotPosNeg(com_dynamic, "species richness")
mtext("A.", 2, adj=3, las=1, padj=-10)


plot(mean_turnover_pos, type = "l", lwd = 3, ylim = c(0,0.5), col = rgb(0, 158, 115, max = 255), xlab = "Time", ylab = "Mean Species Turnover")
lines(mean_turnover_neg, lwd = 3, col = rgb(230, 159, 0, max = 255))
mtext("B.", 2, adj=3, las=1, padj=-10)

dev.off()


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

origination_pos <- rowMeans(origination_t[,1:100] / sp_tmoins1[,1:100], na.rm = TRUE)
origination_neg <- rowMeans(origination_t[,101:200] / sp_tmoins1[,101:200], na.rm = TRUE)

extinction_t <- extinction_matrix[2:nsteps,]

extinction_pos <- rowMeans(extinction_t[,1:100] / sp_tmoins1[,1:100], na.rm = TRUE)
extinction_neg <- rowMeans(extinction_t[,101:200] / sp_tmoins1[,101:200], na.rm = TRUE)

path <- paste0("Figures/orignation_extinction_rates", pourcentage, ".png")
png(path)

par(mfrow=c(2,1))

y_max <- max(c(origination_pos, origination_neg, extinction_pos, extinction_neg), na.rm = TRUE)

plot(origination_pos, type = "l", lwd = 3, col = rgb(0, 158, 115, max = 255), ylim = c(0,y_max), xlim = c(0,150), xlab = "Time")
lines(origination_neg, lwd = 3, col = rgb(230, 159, 0, max = 255))
mtext("A.", 2, adj=3, las=1, padj=-10)

plot(extinction_pos, type = "l", lwd = 3, col = rgb(0, 158, 115, max = 255), ylim = c(0,y_max), xlim = c(0,150), xlab = "Time")
lines(extinction_neg, lwd = 3, col = rgb(230, 159, 0, max = 255))
mtext("B.", 2, adj=3, las=1, padj=-10)

dev.off()
