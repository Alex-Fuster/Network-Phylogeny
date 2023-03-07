
# PARAMETERS
# u_max
# u_0
# u_1
# e_max
# e_0neg
# e_1neg
# e_0pos
# e_1pos
# a_u
# a_eneg
# a_epos
# av_r: average range
# beta_r$ : shape parameter for the range
# gam_0: intercept of the n-o relationship
# gam_1: slope of the n-o relationship
# var_o: variance in the niche optimum around the expected niche optimum
# var_n: inheritance parameter (variance of the niche around the ancestor trait)

#setwd("~/Bureau/Documents/Manuscripts/Inprep/ms_gravel_diversification_networks/code")
rm(list = ls())

#pars = list()
#pars$u_max = 0.075 # mutation probability
#pars$u_0 = 1  # 1 for predation & facilitation, 0 for competition
#pars$u_1 = -1 # -1 for predation & facilitation, 1 for competition
#pars$e_0neg = 0.15 # asymptotic extinction probability with infinite negative interactions
#pars$e_1neg = -pars$e_0neg  # Extinction probability with absence of interactions
#pars$e_0pos = 0.01 # asymptotic extinction probability with infinite positive interactions
#pars$e_1pos = 1 - pars$e_0pos # 1 - e0
#pars$a_u = 0.5 # shape of the exponential decay of the colonization - interaction relationship
#pars$a_eneg = 0.025 # shape of the exponential decay of the negative extinction - interaction relationship
#pars$a_epos = 0.5  # shape of the exponential decay of the positive extinction - interaction relationship
#pars$av_r = 0.5 # range of the niche
#pars$beta_n = 1 # parameter of the beta distribution
#pars$int = 1 # 0 for competition, 1 for facilitation, 2 for predation
##pars$Sbasal = 0 # number of basal species
#pars$Smax = 1000 # Maximal number of species allowed
#pars$Bspe = 15 # Basal species impact on the speciation probality
#pars$Bext = 12 # Basal species impact on extinction probability


#source("functions.R")
source("functions_M-C.R")

##############################
# Run the model
# attach(pars)
#seed = 1
seed = c(6877, 2643, 4728, 8219, 4808, 344, 9521, 6433, 8339, 205, 1740, 6562, 7044, 7911, 2525, 9360, 5058, 6128, 7016, 6832)
seed2 = c(1718, 3006, 3259, 1918, 3672, 7510,  473, 9867, 1438, 7594)
nsteps = 250
test = sim_model_bif(seed = 8301, pars, nsteps = nsteps) #seed2[5] = 3672

plot(test$dynamic[,3], type = "l", col = "black", ylim = c(0,50), xlab = "Time", ylab = "Species number")
lines(test025$dynamic[,3], col = "red")
legend("topleft",bty = "n", legend = c("P(ext) = 0.05", "P(ext) = 0.025"), lty = 1, col = c("black","red"))


##############################
# Plot mean trait value
pres = test$pres
traits = test$traits
t=matrix(traits[,1],nr = nsteps, nc = ncol(pres),byrow=TRUE)
mean_trait = apply(pres*t,1,sum,na.rm=T)/apply(pres,1,sum)
plot(c(1:nsteps),mean_trait, type = "l")

##############################
# Plot diversification dynamics
pres = test$pres
traits = test$traits
Stot = ncol(pres)
S = apply(pres,1,sum)

dev.new(width = 5, height = 4)
par(mar = c(5,6,2,1))
plot(c(1:nsteps),S,type = "l", xlab = "Time", ylab = "Species richness",cex.axis = 1.25, cex.lab = 1.5,lwd = 1.5)

##############################
# Plot the trait distribution over time

traits_overtime_plot <- function(res){
	pres = res$pres
	Stot = ncol(pres)
	traits = res$traits

	Time = c(1:nrow(pres))

	plot(Time[pres[,1]==1],rep(traits[1,1],sum(pres[,1])),xlim = c(1,nsteps),ylim = c(0,1), pch = 19, cex = 0.3, xlab = "Time", ylab = "Niche position")
	for(i in 1:Stot) points(Time[pres[,i]==1],rep(traits[i,1],sum(pres[,i])),cex = 0.3,pch = 19)
}

##############################
# Number of speciation and extinction events
t0 = pres[1:(nsteps-1),]
t1 = pres[2:nsteps,]
spec_mat = pres[1:(nsteps-1),]*0
ext_mat = pres[1:(nsteps-1),]*0
spec_mat[t1-t0==1] = 1
ext_mat[t1-t0==-1] = 1
spec = apply(spec_mat ,1,sum)
ext = apply(ext_mat, 1, sum)

S = apply(pres,1,sum)[2:nsteps]

plot(Time[2:nsteps], spec/S, type = "l", ylim = range(spec/S,ext/S),xlab = "Time", ylab = "Rate", col = "darkblue")
lines(Time[2:nsteps], ext/S, col = "darkred")
legend("topleft",bty = "n", legend = c("Speciation", "Extinction"), lty = 1, col = c("darkblue","darkred"))


##############################
# Diversification-richness dependence
plot(S[2:nsteps],spec/S, xlab = "Species richness", ylab = "Rate per species", pch = 21, bg = "darkblue",cex = 0.8, cex.axis = 1.25, cex.lab = 1.5)

points(S[2:nsteps],ext/S, pch = 21, bg = "darkred",cex = 0.8)
legend("topright",bty = "n", legend = c("Speciation", "Extinction"), pch = 21, pt.bg = c("darkblue","darkred"))


##############################
# Reconstruct networks over time

n = numeric(nsteps)
r = numeric(nsteps)
o = numeric(nsteps)
nL = numeric(nsteps)

# Number of links
L = test$L_list[[250]]

for(t in 1:nsteps) {

	# Subset only the species that are present
	n = traits[pres[t,]==1,1]
	r = traits[pres[t,]==1,2]
	o = traits[pres[t,]==1,3]

	# Compute the average
	n[t] = mean(n)
	r[t] = mean(r)
	o[t] = mean(o)

	# Subset interactions between present Species
	subL = L[pres[t,]==1,pres[t,]==1]

	nL[t] = sum(subL)
}
S = apply(pres,1,sum)

dev.new(width = 5, height = 4)
par(mar = c(5,5,2,1))
plot(c(1:nsteps),nL/S/S, type = "l", xlab = "Time", ylab = "Connectance",cex.axis = 1.25, cex.lab = 1.5,lwd = 1.5)


##########################

#pars = list()
#pars$u_max = 0.075 # mutation probability
#pars$u_0 = 1  # 1 for predation & mutualism, 0 for competition
#pars$u_1 = -1 # -1 for predation & mutualism, 1 for competition
#pars$e_0neg = 0.15 # asymptotic extinction probability with infinite negative interactions
#pars$e_1neg = -pars$e_0neg  # Extinction probability with absence of interactions
#pars$e_0pos = 0.05 # asymptotic extinction probability with infinite positive interactions
#pars$e_1pos = 1 - pars$e_0pos # 1 - e0
#pars$a_u = 0.5 # shape of the exponential decay of the colonization - interaction relationship
#pars$a_eneg = 0.025 # shape of the exponential decay of the negative extinction - interaction relationship
#pars$a_epos = 0.5  # shape of the exponential decay of the positive extinction - interaction relationship
#pars$av_r = 0.2 # range of the niche
#pars$beta_n = 1 # parameter of the beta distribution
#pars$int = 1 # 0 for competition, 1 for mutualism, 2 for predation
#pars$Smax = 1000 # Maximal number of species allowed
#
#source("functions_M-C.R")
#
##out_put <- matrix(0, 7*8, 3)
##out_put <- as.data.frame(out_put)
##names(out_put) = c("Bext","Bspe","S")
#
#colors_p <- c(rgb(0, 0, 0, max = 255), #noir
#							rgb(86, 180, 233, max = 255), #bleu ciel
#							rgb(213, 94, 0, max = 255), #rouge
#							rgb(0, 158, 115, max = 255), #vert sapin
#							rgb(0, 114, 178, max = 255), #bleu
#							rgb(204, 121, 167, max = 255), #rose
#							rgb(230, 159, 0, max = 255), #orange
#							rgb(240, 228, 66, max = 255), #jaune
#							"#a6cee3",
#							"#1f78b4",
#							"#b2df8a",
#							"#33a02c",
#							"#fb9a99",
#							"#e31a1c",
#							"#fdbf6f",
#							"#ff7f00",
#							"#cab2d6",
#							"#6a3d9a",
#							"#ffff99",
#							"#b15928"
#						)
#
#dynamic_list = list()
#
#par(mfrow=c(1,1))
#
#k = 1
#
#for(i in 10:25){
#	for(j in 5:5){
#		print(paste("Bext = ", i, ", Bspe =", j, "et ligne =", k))
#
#		pars$Bext = i # Basal species impact on extinction probability
#		pars$Bspe = j # Basal species impact on the speciation probality
#
#		seed = 1
#		nsteps = 250
#		test = sim_model(seed = seed, pars, nsteps = nsteps)
#
#		dynamic_list[[k]] = test$dynamic
#
#		#pres <- test$pres
#		#L_list <- test$L_list
#		#n_sp <- length(L_list)
#		#S = sum(pres[n_sp,])
#
#		#out_put[k, 1] = i
#		#out_put[k, 2] = j
#		#out_put[k, 3] = S
#
#		k <- k+1
#	}
#}
#
## Plot les dynamics
#plot(dynamic_list[[1]][,1], dynamic_list[[1]][,3], type ="l", xlim = c(0,350), ylim = c(0, 100), xlab = "Time", ylab = "Species number", main = paste("Bext = ", i))
#for(p in 2:15){
#	lines(dynamic_list[[p]][,1], dynamic_list[[p]][,3], col = colors_p[p])
#}
#
#
#n_b = 9
#
#ymax = max(out_put[,3])
#ymin = min(out_put[,3])
#
##pdf("test_Bspe-Bext_Smax1000.pdf")
#plot(c(12:20), out_put[1:n_b,3], type = "l", ylim = c(ymin, ymax), xlab = "Bspe", ylab = "S at the end", col = rgb(0, 0, 0, max = 255))
#lines(c(12:20), out_put[(n_b+1):(2*n_b),3], col = rgb(86, 180, 233, max = 255))
#lines(c(12:20), out_put[(2*n_b+1):(3*n_b),3], col = rgb(213, 94, 0, max = 255))
#lines(c(12:20), out_put[(3*n_b+1):(4*n_b),3], col = rgb(0, 158, 115, max = 255))
#lines(c(12:20), out_put[(4*n_b+1):(5*n_b),3], col = rgb(0, 114, 178, max = 255))
##lines(c(12:20), out_put[(5*n_b+1):(6*n_b),3], col = rgb(204, 121, 167, max = 255))
#lines(c(12:20), out_put[(6*n_b+1):(7*n_b),3], col = rgb(230, 159, 0, max = 255))
##lines(c(12:20), out_put[(7*n_b+1):(8*n_b),3], col = rgb(240, 228, 66, max = 255))
##lines(c(12:20), out_put[(8*n_b+1):(9*n_b),3], col = rgb(86, 180,233, max = 255))
##lines(c(12:20), out_put[(9*n_b+1):(10*n_b),3], col = rgb(0, 158,115, max = 255))
##lines(c(12:20), out_put[(10*n_b+1):(11*n_b),3], col = rgb(0, 114,178, max = 255))
##lines(c(12:20), out_put[(11*n_b+1):(12*n_b),3], col = rgb(0, 114,178, max = 255))
##lines(c(12:20), out_put[(12*n_b+1):(13*n_b),3], col = rgb(213, 94,0, max = 255))
##lines(c(12:20), out_put[(13*n_b+1):(14*n_b),3], col = rgb(86, 180,233, max = 255))
##lines(c(12:20), out_put[(14*n_b+1):(15*n_b),3], col = rgb(0, 114,178, max = 255))
#legend("bottomright",
#			legend = c("Bext = 6", "Bext = 7", "Bext = 8", "Bext = 9", "Bext = 10", "Bext = 12", "Bext > 13"),
#			col=c(rgb(230, 159,0, max = 255), rgb(86, 180,233, max = 255), rgb(213, 94,0, max = 255), rgb(0, 158,115, max = 255), rgb(0, 114,178, max = 255), rgb(0,0,0, max = 255), rgb(230, 159,0, max = 255)),
#			bty = "n",
#			lty = 1
#		)
#title(paste0("Smax = ", pars$Smax))
##dev.off()


#############

#start_time <- Sys.time()
#for(i in 1:20){
#	test = sim_model(seed = 1, pars, nsteps = nsteps)
#}
#end_time <- Sys.time()
#end_time - start_time

##########################
##########################
## TEST DES SIMULATIONS ##
##########################
##########################

pars = list()
pars$u_max = 0.15 # mutation probability
pars$av_r = 0.11 # range of the niche
pars$sd = 2*pars$av_r + 0.0001 #standard deviation of the normal distribution used to calculate the trait n of a new species

pars$a_eneg = 0.025 # shape of the exponential decay of the negative extinction - interaction relationship
pars$e_0neg = 0.15 # asymptotic extinction probability with infinite negative interactions
pars$e_1neg = -pars$e_0neg  # Extinction probability with absence of interactions
pars$a_epos = 1.2  # shape of the exponential decay of the positive extinction - interaction relationship
pars$e_0pos = 0.05
pars$e_1pos = 5.19
pars$Bext = 4 # constant minimal number of interaction per species (extinction prob)

pars$a_u = 0.45 # shape of the exponential decay of the colonization - interaction relationship
pars$u_0 = 1
pars$u_1 = -1
pars$d = 0.5 # decrease speed of the establishment probability
pars$Bspe = 4 # constant minimal number of interaction per species (extablishment prob)

pars$I_max = 12 # maximal number of interaction per species
pars$Smax = 2000 # maximal number of species in the system
pars$int = 1 # 0 for competition, 1 for mutualism, 2 for predation


source("functions_M-C.R")
#bext = c(10,12,15)
#bspe = c(10,12,15)
#eopos = 0.05
#u0 = c(0.6,0.8)
#u1 = c(-0.8,-0.9)
#tested_parameters = expand.grid(bext,bspe,eopos,u0,u1)

dyna = matrix(0, 249, 10) #nrow = nstep - 1 ; ncol = nombre total de simulations
#dynamic_list = list()

traits_sim = matrix(0, 2000, 10) #nrow = S_max ; ncol = nombre total de simulations
pres_sim = list()

list_append <- function (lst, ...){
lst <- c(lst, list(...))
  return(lst)
}

# Nombre de simulations total, en comptant les réplicats
sim_number = nrow(tested_parameters) * 10

# Tirage aléatoire des seeds pour les simulations
seed2 = sample.int(10000, 10)

#for(k in 1:nrow(tested_parameters)){
	#cat("Simulation ", k, "sur ",nrow(tested_parameters),'\n')

	#pars$Bext = tested_parameters[k,1] # Basal species impact on extinction probability
	#pars$Bspe = tested_parameters[k,2] # Basal species impact on the speciation probality
	#pars$e_0pos = tested_parameters[k,3]
	#pars$e_1pos = 1 - pars$e_0pos
	#pars$u_0 = tested_parameters[k,4]
	#pars$u_1 = tested_parameters[k,3]

	nsteps = 250

	for(s in 1:10){
		cat("sn = ", s, "sur 10, et le seed = ", seed2[s], '\n')
		test = sim_model(seed = seed2[s], pars, nsteps = nsteps)
		#print(max(test$dynamic[,3]))
		dyna[,s] = test$dynamic[,3]
		#traits_sim[,s] = test$traits[,1]
		#pres_sim = list_append(pres_sim, test$pres)
		# progress
		#cat("sim", s, "of", 20, '\n')
	}
	#dynamic_list[[k]] = dyna
#}

## Plot the
#lot(test$dynamic[,2], xlab = "Time", ylab = "Cummulative number of species", type = "l")
#plot(test$dynamic[,3], xlab = "Time", ylab = "Actual number of species", type = "l")

plot(dyna[,1], type = "l", ylim=c(0,20), xlim = c(0,20))
for(i in 2:ncol(dyna)){
	Sys.sleep(1)
	lines(dyna[,i], col = i)
}

dead_sim <- length(which(is.na(dyna[15,])))

# Species longevity
pres05 = test05$pres

sp_long05 <- colSums(pres)
hist(sp_long, xlab = "species longevity", breaks = 200, main = "P(ext) = 0.05", xlim = c(0,250), ylim = c(0,1200))

pres025 = test025$pres

sp_long025 <- colSums(pres025)
hist(sp_long025, xlab = "species longevity", breaks = 200, main = "P(ext) = 0.025", xlim = c(0,250), ylim = c(0,1200))

####
# Variation coefficient

presence = test025$pres
traits = test025$traits

trait_pres_mat <- matrix(NA, nsteps, Stot)
for(i in 1:ncol(presence)){
	trait <- traits[i, 1]
	for(j in 1:nsteps)
		if(pres[j,i] == 1){
			trait_pres_mat[j,i] = trait
		}
}

cv <- matrix(0, 1, nsteps)

for(i in 1:nsteps){
	if(sum(presence[i,])>1){
		mean <- mean(trait_pres_mat[i,], na.rm=TRUE)
		var <- sd(trait_pres_mat[i,], na.rm=TRUE)
		cv[i] <- var/mean
	}else{
		mean <- mean(trait_pres_mat[i,], na.rm=TRUE)
		cv[i] <- 0/mean
	}
}

plot(cv[1,], xlab = "Time", ylab = "Coefficient of Variation", main = "P(ext) = 0.025", type = "l")


#output_meandynamic = matrix(0,nsteps-1,length(dynamic_list))
#output_sddynamic = matrix(0,nsteps-1,length(dynamic_list))
#
#for(i in 1:length(dynamic_list)){
#	output_meandynamic[,i] = rowMeans(dynamic_list[[i]])
#}
#
#plot(output_meandynamic[,1], ylim = c(0,max(output_meandynamic)), type = "l")
#for(i in 2:ncol(output_meandynamic)){
#	lines(output_meandynamic[,i])
#}

## plot des dynamiques moyennes (avec variance) pour chaque set de paramètres
##par(mfrow = c(3, 3))
#for(i in seq(1, ncol(output_meandynamic), by=9)){
#	x11()
#	#par(mfrow = c(3, 3))
#	for(j in 0:8){
#		plot(NULL, xlim = c(0, nsteps), ylim = c(0, 20), xlab = "Time", ylab = "Species number", main = paste("eopos = ", tested_parameters[i,3], "u0 = ", tested_parameters[i,4],"u1 = ", tested_parameters[i,5]))
#		lines(output_meandynamic[,i], col = rgb(86, 180, 233, max = 255), lty = 1) #bleu ciel = Bspe = 10, bext = 10
#		lines(output_meandynamic[,i+1+j], col = rgb(86, 180, 233, max = 255), lty = 3) #Bext = 12
#		lines(output_meandynamic[,i+2+j], col = rgb(86, 180, 233, max = 255), lty = 5) #Bext = 15
#		lines(output_meandynamic[,i+3+j], col = rgb(213, 94, 0, max = 255), lty = 1)#rouge = Bspe = 12
#		lines(output_meandynamic[,i+4+j], col = rgb(213, 94, 0, max = 255), lty = 3)
#		lines(output_meandynamic[,i+5+j], col = rgb(213, 94, 0, max = 255), lty = 5)
#		lines(output_meandynamic[,i+6+j], col = rgb(0, 158, 115, max = 255), lty = 1)#vert sapin = Bspe = 15
#		lines(output_meandynamic[,i+7+j], col = rgb(0, 158, 115, max = 255), lty = 3)
#		lines(output_meandynamic[,i+8+j], col = rgb(0, 158, 115, max = 255), lty = 5)
#		#legend("topleft",
#		#			legend = c("Bspe = 10", "Bspe = 12", "Bspe = 15"),
#		#			col=c(rgb(86, 180,233, max = 255), rgb(213, 94,0, max = 255), rgb(0, 158,115, max = 255)),
#		#			bty = "n",
#		#			lty = 1
#		#		)
#	}
#}



spe = matrix(0,1,60)
for(i in 1:60){
	spe[1,i] <- (pars$u_0 + pars$u_1 * exp(-pars$a_u * i))#/(1 + exp(pars$d * (i - pars$I_max)))
}

ext = matrix(0,1,60)
for(i in 1:60){
	ext[1,i] <- pars$e_0pos + pars$e_1pos*exp(-pars$a_epos* i)
}

spe_neg = matrix(0,1,60)
for(i in 1:60){
	spe_neg[1,i] <- pars$u_0neg + pars$u_1neg*exp(-pars$a_uneg * i)
}

extneg = matrix(0,1,60)
for(i in 1:60){
	extneg[1,i] <- pars$e_0neg * (1 - exp(-pars$a_eneg * i))
}

par(mfrow = c(2,1))
par(mfrow = c(1,1))

# Probiblities for facilitation interactions
plot(spe[1,], type = "l", xlab = "Number of interaction", ylim = c(0,1), ylab = "Probability", frame.plot = FALSE, main = "Facilitation", xlim = c(0,50))
lines(ext[1,], col = rgb(0, 158, 115, max = 255))
legend("topright",
			legend = c("establishment", "extinction"),
			col = c("black", rgb(0, 158, 115, max = 255)),
			bty = "n",
			lty = 1
		)
*mutation[1,]

# Probiblities for competition interactions
plot(spe_neg[1,], type = "l", xlab = "Number of interaction", ylim = c(0,1), ylab = "Probability", frame.plot = FALSE, main = "Competition", xlim = c(0,50))
lines(extneg[1,], col = rgb(0, 158, 115, max = 255))
legend("topright",
			legend = c("establishment", "extinction"),
			col = c("black", rgb(0, 158, 115, max = 255)),
			bty = "n",
			lty = 1
		)


#png("estab_prb.png")
#plot(spe[1,], type = "l", xlab = "Number of interaction", ylim = c(0,1), ylab = "Probability of establishment", frame.plot = FALSE)
#dev.off()
#
#png("ext_prb_facil.png")
#plot(ext[1,], type = "l", xlab = "Number of interaction", ylim = c(0,1), ylab = "Probability of extinction", frame.plot = FALSE)
#dev.off()
#
#png("ext_prb_compet.png")
#plot(extneg[1,], type = "l", xlab = "Number of interaction", ylim = c(0,1), ylab = "Probability of extinction", frame.plot = FALSE)
#dev.off()

#png("beta_dist.png")
#curve(dbeta(x, 0.5, 1), 0,1)
#dev.off()

#n = runif(1,0,1)
n = 0.5
b = 1/n
a = b*n/(1-n)
#paste0("n = ", n, ", alpha = ", a, ", beta = ",pars$beta_n)
curve(dbeta(x, a, b), ylab = "PDF")
legend("topleft",
			legend = c(paste0("n = ", round(n, digit = 4)), paste0("α = ", round(a, digit = 4)), paste0("β = ",pars$beta_n)),
			bty = "n"
		)

curve(dnorm(x, n, pars$sd, log = FALSE), xlab = "trait n of the new species", ylab = "PDF")
legend("topleft",
			legend = c(paste0("n = ", round(n, digit = 4)), paste0("sd = ", pars$sd), paste0("niche range = ",pars$av_r*2)),
			bty = "n"
		)
abline(v=n, col=rgb(0, 158, 115, max = 255))

functional_seeds <- c(521, 948, 907, 1260, 2786, 9933, 4115, 836, 2627, 4336, 4606, 5250, 7880, 7165, 7108, 5677, 6122, 9525, 8584, 6427, 7396, 144, 3168, 9984, 7954, 331, 6134)

## Tests mutation function

mutation = matrix(0,1,100)
for(i in 1:100){
	mutation[1,i] <- pars$u_max/(1 + exp(pars$d * (i - 40)))
}
plot(mutation[1,], type = "l",xlim = c(0,60), ylab = "Speciation probability", xlab = "Species richness")

mutation_proba <- (pars$u_max/(1 + exp(pars$d * (sum_I - pars$I_max)))
