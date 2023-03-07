########################################
# Function to generate a new set of traits for ancestors
# Each species is characterized by a set of 3 traits: n (niche), o (optimum) and r (niche's range)
rand_traits_anc = function(pars) {
	with(as.list(pars),{
		n = runif(1, 0, 1)
		r = av_r
		o = runif(1,0,n)
		traits = c(n = n, r = r, o = o)
		traits
	})
}

########################################
# Function to generate a new set of traits from mutants
rand_traits_mut = function(traits_anc, pars) {

	with(as.list(c(traits_anc, pars)),{

		if(int == 0) {
			a = beta_n*n/(1-n)
			n_m = rbeta(1, shape1 = a, shape2 = beta_n)
			traits_mut = c(n = n_m, r = r, o = n)
		}

		if(int == 1) {
			a = beta_n*n/(1-n)
			n_m = rbeta(1, shape1 = a, shape2 = beta_n)
			traits_mut = c(n = n_m, r = r, o = n)
		}

		if(int == 2) {
			a = beta_n*n/(1-n)
			n_m = rbeta(1, shape1 = a, shape2 = beta_n)
#			o_m = runif(1,0,n_m)
			o_m = n_m/2
			traits_mut = c(n = n_m, r = r, o = o_m)
		}
		traits_mut
	})
}

########################################
# Function to compute the interaction network from a set of traits
get_L_mat = function(basal, pars, traits_mat) {
	with(as.list(pars),{
		L = matrix(0, nr = pars$Smax+pars$Sbasal, nc = pars$Smax)

		# Lower boundary
		low = traits_mat$o - traits_mat$r
		low_mat = matrix(low, nr = pars$Smax+pars$Sbasal, nc = pars$Smax, byrow = TRUE)

		# Upper boundary
		high = traits_mat$o + traits_mat$r
		high_mat = matrix(high, nr = pars$Smax+pars$Sbasal, nc = pars$Smax, byrow = TRUE)
		S = nrow(traits_mat)

		# Matrix of niche positions
		n_mat = matrix(traits_mat$n, nr = pars$Smax, nc = pars$Smax, byrow = FALSE)

		# Add the basal species
		n_basal = matrix(basal, nr = pars$Sbasal, nc = pars$Smax, byrow = FALSE)
		n_mat = rbind(n_basal, n_mat)

	  	# Test interactions
	  	L[n_mat > low_mat & n_mat < high_mat] = 1
	  	if(pars$Smax > 1) diag(L[(pars$Sbasal+1):(pars$Sbasal+pars$Smax),]) = 0
	  	L
  })
}

########################################
# Function to compute the interactions of a given species
get_L_vec = function(basal, pars, traits_mat, traits_mut) {
	with(as.list(pars),{
		L_vec = numeric(pars$Smax+pars$Sbasal)

		# Lower boundary
		low = traits_mut["o"] - traits_mut["r"]

		# Upper boundary
		high = traits_mut["o"] + traits_mut["r"]

		# Vector of niche positions
		n_vec = c(basal, traits_mat$n)

	  	# Test interactions
	  	L_vec[n_vec > as.numeric(low) & n_vec < as.numeric(high)] = 1
	  	L_vec
  })
}

########################################
sim_model = function(seed, pars, nsteps) {

	with(pars,{

	set.seed(seed)

	# Draw the traits of the producers
#	basal = runif(pars$Sbasal, 0, 1)
	basal = runif(pars$Sbasal, 0, 0.2)

	# Draw the first species traits
	traits_mat = matrix(nr = pars$Smax, nc = 3)
	traits_mat[1,] = rand_traits_anc(pars)
	traits_mat = as.data.frame(traits_mat)
	names(traits_mat) = c("n","r","o")

	# Set the presence/absence matrix
	pres = matrix(0, nr = nsteps, nc = pars$Smax)
	pres[1,1] = 1

 	# Set the ancestry object
	anc = matrix(0,nr = pars$Smax, nc = 3)

	# Set the extinctions object
	extinct = matrix(0,nr = pars$Smax, nc = 2)

	# Record the matrices
	L_list = list()

	# Species count
	S = 1

	##########################
	# MAIN LOOP
	for(step in 2:nsteps) {
		ActualS = sum(pres[step-1,])
		cat("Step = ", step-1," Total S = ", S," Actual S = ", ActualS,'\n')

		# Test for successful speciation probability
		for(i in 1:Smax) {

			if(S >= Smax) break

			# Speciation occurs if the species is present
			if(pres[step-1,i] == 1) {

				# Species is maintained
				pres[step, i] = 1

				# Test if there is mutation
				if(runif(1, 0, 1) < u_max) {

					# Pick new parameters
					traits_mut = rand_traits_mut(traits_mat[i,], pars)

					# Recompute the interactions
					I = get_L_vec(basal, pars, traits_mat, traits_mut)

					# Compute the number of interactions among present species
					sum_I = sum(I*c(rep(1,pars$Sbasal),pres[step,]))

					# Compute the probability of successful speciation
					spec_prob = pars$u_0 + pars$u_1*exp(-pars$a_u*sum_I)

					# Test if there is speciation
					if(runif(1) < spec_prob) {
						S = S + 1
						traits_mat[S,] = traits_mut
						pres[step,S] = 1
						anc[S,] = c(step, i, S) #Time step, ancestor, new_species
					}
				}
			}
		}

		if(S >= Smax) break
		# Test for extinction

		# Compute the interaction matrix among present species
		pres_vec = pres[step,]
		cooc = matrix(pres_vec, nr=pars$Smax, nc=pars$Smax, byrow=TRUE)*
		matrix(pres_vec, nr=pars$Smax, nc=pars$Smax, byrow=FALSE)
		L = get_L_mat(basal, pars, traits_mat)
		L[c((Sbasal+1):(Sbasal+Smax)),]= L[c((Sbasal+1):(Sbasal+Smax)),]*cooc
		L_list[[step]] = L

		if(int == 0) {
			in_I = colSums(L)
			ext_prob = e_0neg + e_1neg*(1 - exp(-a_eneg*in_I))
		}

	 	if(int == 1) {
			in_I = colSums(L)
			ext_prob = e_0pos + e_1pos*exp(-a_epos*in_I)
		}

		if(int == 2) {
			in_I = colSums(L)
			out_I = rowSums(L)[(Sbasal+1):(Sbasal+Smax)]
			ext_prob = e_0neg + e_1neg*exp(-a_eneg*out_I) + e_0pos + e_1pos*exp(-a_epos*in_I)
		}

		# Perform extinctions
		pres[step, pres[step-1,] & runif(Smax,0,1) < ext_prob] = 0

		for(i in 1:Smax) {
			if(pres[step,i] != pres[step-1,i] & pres[step,i] == 0){
				extinct[i,] = c(step, i)
			}
		}
	} # End of main loop

	list(pres = pres, traits = traits_mat, anc = anc, L_list = L_list, basal = basal, extinct = extinct)
	})
}

#### Tester le temps de calcul des fonctions
start_time <- Sys.time()
rand_traits_anc(pars)
end_time <- Sys.time()
end_time - start_time
