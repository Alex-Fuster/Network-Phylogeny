# Plot traits of existing species at n = 350
pres <- simulation$pres
pres_real <- pres[1:nsteps,1:species_number]
extant_sp <- pres[350, 1:species_number]
to_keep <- which(extant_sp == 1)

traits <- simulation$traits
sp_traits <- traits[to_keep, 1]
y <- c(rep(0,length(sp_traits)))

#colrs <-c("gray50", "tomato", "gold", "lightgreen", "lightblue", "steel blue") #, "yellowgreen"

plot(x = sp_traits, y = y, pch = 19, col = colrs[groups_sp])


#par(mfrow = c(3,1))

# Plot the species life span
pres <- simulation$pres
sp_lifespan <- apply(pres, 2, sum)
hist(sp_lifespan)


# Plot the trait distribution over time
pres = simulation$pres
traits = simulation$traits
Stot = ncol(pres)

Time = c(1:nrow(pres))

plot(Time[pres[,1]==1],rep(traits[1,1],sum(pres[,1])),xlim = c(1,nsteps),ylim = c(0,1), pch = 19, cex = 0.3, xlab = "Time", ylab = "Niche position")
for(i in 1:Stot) points(Time[pres[,i]==1],rep(traits[i,1],sum(pres[,i])),cex = 0.3,pch = 19)
