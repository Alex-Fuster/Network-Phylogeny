
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

setwd("~/Bureau/Documents/Manuscripts/Inprep/ms_gravel_diversification_networks/code")
rm(list = ls())
source("functions.R")

pars = list()
pars$u_max = 0.075 # mutation probability
pars$u_0 = 1  # 1 for predation & mutualism, 0 for competition
pars$u_1 = -1 # -1 for predation & mutualism, 1 for competition
pars$e_0neg = 0.15 # asymptotic extinction probability with infinite negative interactions
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
pars$Smax = 1000 # Maximal number of species allowed

##############################
# Run the model
# attach(pars)
seed = 1
nsteps =  250
test = sim_model(seed = runif(1), pars, nsteps = nsteps)

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
pres = test$pres
traits = test$traits

Time = c(1:nrow(pres))

plot(Time[pres[,1]==1],rep(traits[1,1],sum(pres[,1])),xlim = c(1,nsteps),ylim = c(0,1), pch = 19, cex = 0.3, xlab = "Time", ylab = "Niche position")
for(i in 1:Stot) points(Time[pres[,i]==1],rep(traits[i,1],sum(pres[,i])),cex = 0.3,pch = 19)

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
	subL = L[c(rep(1,25),pres[t,])==1,pres[t,]==1]

	nL[t] = sum(subL)
}
S = apply(pres,1,sum)+25

dev.new(width = 5, height = 4)
par(mar = c(5,5,2,1))
plot(c(1:nsteps),nL/S/S, type = "l", xlab = "Time", ylab = "Connectance",cex.axis = 1.25, cex.lab = 1.5,lwd = 1.5)










