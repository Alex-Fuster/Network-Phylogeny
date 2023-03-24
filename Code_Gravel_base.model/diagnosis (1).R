setwd("/home/local/USHERBROOKE/grad3002/Bureau/Docs/Manuscripts/Inprep/ms_gravel_diversification_networks/code")
rm(list = ls())
source("functions.R")

pars = list()
pars$u_max = 0.1 # mutation probability
pars$u_0 = 1  # 1 for predation & mutualism, 0 for competition
pars$u_1 = -1 # -1 for predation & mutualism, 1 for competition
pars$e_0neg = 0.0 # asymptotic extinction probability with infinite negative interactions
pars$e_1neg = -pars$e_0neg  # Extinction probability with absence of interactions
pars$e_0pos = 0.01 # asymptotic extinction probability with infinite positive interactions
pars$e_1pos = 1 - pars$e_0pos # 1 - e0
pars$a_u = 0.5 # shape of the exponential decay of the colonization - interaction relationship	
pars$a_eneg = 0.025 # shape of the exponential decay of the negative extinction - interaction relationship 
pars$a_epos = 0.5  # shape of the exponential decay of the positive extinction - interaction relationship 
pars$av_r = 0.2 # range of the niche
pars$beta_n = 1 # parameter of the beta distribution 
pars$int =  2 # 0 for competition, 1 for mutualism, 2 for predation
pars$Sbasal = 25 # number of basal species
pars$Smax = 15 # Maximal number of species allowed

plot_curves = function(pars) {
	with(as.list(pars),{
		I = seq(0,100,0.1)
		eneg = e_0neg + e_1neg*exp(-a_eneg*I)
		epos = e_0pos + e_1pos*exp(-a_epos*I)
		s = u_0 + u_1*exp(-a_u*I)
		dev.new(width = 15, height = 10)
		plot(I, eneg, col = "darkred", ylim = c(0,1),type = "l",lwd = 2)
		lines(I, epos, col = "darkblue", lwd = 2)
		lines(I, eneg+epos, lwd = 2)
		lines(I, s, col = "darkgreen",lwd = 2)
	})
}

#plot_curves(pars)

##############################
# Run the model
# attach(pars)
seed = 1
nsteps = 25
test = sim_model(seed = seed, pars, nsteps = nsteps)

#dev.new(width = 14, height = 20)
#par(mar = c(5,6,2,1),mfrow=c(3,1))

##############################
# Plot diversification dynamics
pres = test$pres
traits = test$traits
Stot = ncol(pres)
S = apply(pres,1,sum)
#plot(c(1:nsteps),S,type = "l", xlab = "Time", ylab = "Species richness",
#cex.axis = 2.25, cex.lab = 2.5, lwd = 1.5 )

##############################
# Plot the trait distribution over time
pres = test$pres
traits = test$traits

Time = c(1:nrow(pres))

#plot(Time[pres[,1]==1],rep(traits[1,1],sum(pres[,1])),xlim = c(1,nsteps),
#ylim = c(0,1), pch = 19, cex = 0.3, xlab = "Time", ylab = "Niche position",
#cex.axis = 2.25, cex.lab = 2.5)
#for(i in 1:pars$Smax) points(Time[pres[,i]==1],rep(traits[i,1],sum(pres[,i])),cex = 0.3,pch = 19)

##############################
# Frequency distribution of traits at the end of simulation

step = 80
#hist(traits[pres[step,]==1,1],cex.axis = 2.25, cex.lab = 2.5,main="",xlab = "Niche position")





