traits_overtime_plot <- function(res){
	pres = res$pres
	Stot = ncol(pres)
	traits = res$traits

	Time = c(1:nrow(pres))

	plot(Time[pres[,1]==1],rep(traits[1,1],sum(pres[,1])),xlim = c(1,nsteps),ylim = c(0,1), pch = 19, cex = 0.3, xlab = "Time", ylab = "Niche position")
	for(i in 1:Stot) points(Time[pres[,i]==1],rep(traits[i,1],sum(pres[,i])),cex = 0.3,pch = 19)
}
