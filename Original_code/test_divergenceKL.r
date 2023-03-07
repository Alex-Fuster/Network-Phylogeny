x <- seq(-5, 5, length=300)
y <- cbind(u=dunif(x,0,2), n1=dnorm(x), n2=dnorm(x,1), t=dt(x, df=10))
matplot(x,y[,"n1"],type="l")

KL <- function(p, q, seuil=1e-10 ) {
    step = 1/sum(p)
    use <- (p>seuil) & (q>seuil) #pour éviter des pbms numériques
    k <- sum( p[use] * log(p[use]/q[use]) )
    return( k #* step )
    }

DivKL <- function(p, q){
    kl <-
    for(i in length(kl)){
        kl[i] <- p[i] * ln(p[i]/q[i])
    }
    KLD <- sum(kl)
    return(KLD)
}

KL(y[,"n1"] , y[,"n1"] )
KL(y[,"n1"] , y[,"n2"] )
KL(y[,"n1"] , y[,"t"] )

#Notons que KL(P,Q) != KL(Q,P)
KL(y[,"t"] , y[,"n1"])

#Attention, les cas suivants ne sont pas corrects si on inverse P et Q (pourquoi ?)
KL(y[,"u"] , y[,"n1"] )
KL(y[,"u"] , y[,"n2"] )

x <- seq(0, 1, length=300)
y <- cbind(u=ppois(x, 0), n1=dnorm(x))
matplot(x,y[,"u"],type="l")

KL(y[,"n1"] , y[,"u"] )

0.005*log(0.005/0.01)+0.995*log(0.995/0.99)
0.01*log(0.01/0.005)+0.99*log(0.99/0.995)

library(FNN)
KL.divergence(y[,"n1"] , y[,"u"])
library(tawny)
divergence.kl(y[,"n1"] , y[,"u"])
library(LaplacesDemon)
KLD(y[,"n1"] , y[,"u"])
