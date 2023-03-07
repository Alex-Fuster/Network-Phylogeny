timetree <- function(tree){
  try(ape::chronoMPL(ape::multi2di(tree)), silent=TRUE)
}

#### Gamma-statistics

gamma_stat<-function(Nspp,alpha0){
  gamma_values=data.frame(0,0)
  colnames(gamma_values)=c("gamma0","gammap0")

  tk=rep(0,Nspp)
  tkp=rep(0,Nspp)
  tk = 0.0
  tkp = 0.0
  tk[2] = 2.0^(1.0-alpha0)
  tkp[2] = -log(2.0)*tk[2]
  for(i in 3:Nspp){
    ai = as.double(i)
    aux = ai^(1.0-alpha0)
    tk[i] = tk[i-1] + aux
    tkp[i] = tkp[i-1] - log(ai)*aux
  }

  deno = tk[Nspp]/sqrt(as.double(Nspp-2)*12.0)

  aux1 = 0.0
  aux2 = 0.0
  for(i in 2:(Nspp-1)){
    aux1 = aux1 + tk[i]
    aux2 = aux2 + tkp[i]
  }

  gamma_values$gamma0 = (aux1/as.double(Nspp-2) - 0.5*tk[Nspp])/deno
  gamma_values$gammap0 = (aux2 - aux1*tkp[Nspp]/tk[Nspp])/(deno*as.double(Nspp-2))

  return(gamma_values)
}


#### Alpha-value
alpha_value <- function(Nspp,gm){
  error = 0.001
  alpha0 = 1.0
  deltaalpha = 1.0

  while(abs(deltaalpha) > error){
    gm_calc = gamma_stat(Nspp,alpha0)
    deltaalpha = (gm - gm_calc$gamma0)/gm_calc$gammap0
    alpha0 = alpha0 + deltaalpha
  }
  alpha = alpha0
  return(alpha)
}

alpha_stat <- function(my_tree){
  Ntips=length(my_tree$tip.label) #number of spp/tips
  gamma_ape=ape::gammaStat(my_tree) #calculating gamma with ape package
  alpha_calc=alpha_value(Ntips, gamma_ape)
  return(alpha_calc)
}
