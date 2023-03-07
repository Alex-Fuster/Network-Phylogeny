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
			#b = 1/n
			#a = b*n/(1-n)
			#n_m = rbeta(1, shape1 = a, shape2 = b)
			n_m = rnorm(1, mean = n, sd = pars$sd)
			r_m = -1
			while(r_m < 0 || r_m > 1){r_m = r + rnorm(1, mean = 0, sd = 0.02)}
			traits_mut = c(n = n_m, r = r_m, o = n)
		}

		if(int == 1) {
			#b = 1/n
			#a = b*n/(1-n)
			#n_m = rbeta(1, shape1 = a, shape2 = b)
			#P(n)={sqrt[1/(2*Pi*sd)]}*{exp[(n-n_p)^2/2*sd^2]}
			n_m = rnorm(1, mean = n, sd = pars$sd)
			r_m = -1
			while(r_m < 0 || r_m > 1){r_m = r + rnorm(1, mean = 0, sd = 0.02)}
			traits_mut = c(n = n_m, r = r_m, o = n)
		}

		#if(int == 2) {
		#	a = beta_n*n/(1-n)
		#	n_m = rbeta(1, shape1 = a, shape2 = beta_n)
#	#		o_m = runif(1,0,n_m)
		#	o_m = n_m/2
		#	traits_mut = c(n = n_m, r = r, o = o_m)
		#}
		traits_mut
	})
}

########################################
# Function to compute the interaction network from a set of traits
get_L_mat = function(pars, traits_mat) {
	with(as.list(pars),{
		L = matrix(0, nr = pars$Smax, nc = pars$Smax)

		# Lower boundary
		low = traits_mat$n - traits_mat$r
		low_mat = matrix(low, nr = pars$Smax, nc = pars$Smax, byrow = TRUE)

		# Upper boundary
		high = traits_mat$n + traits_mat$r
		high_mat = matrix(high, nr = pars$Smax, nc = pars$Smax, byrow = TRUE)
		S = nrow(traits_mat)

		# Matrix of niche positions
		n_mat = matrix(traits_mat$n, nr = pars$Smax, nc = pars$Smax, byrow = FALSE)

		## Add the basal species
		#n_basal = matrix(basal, nr = pars$Sbasal, nc = pars$Smax, byrow = FALSE)
		#n_mat = rbind(n_basal, n_mat)

	  	# Test interactions
	  	L[n_mat > low_mat & n_mat < high_mat] = 1
	  	if(pars$Smax > 1) diag(L) = 0
	  	L
  })
}

########################################
# Function to compute the interactions of a given species
get_L_vec = function(pars, traits_mat, traits_mut) {
	with(as.list(pars),{
		L_vec = numeric(pars$Smax)

		# Lower boundary
		low = traits_mut["n"] - traits_mut["r"]

		# Upper boundary
		high = traits_mut["n"] + traits_mut["r"]

		# Vector of niche positions
		n_vec = c(traits_mat$n)

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
  #basal = runif(pars$Sbasal, 0, 1)
	#basal = runif(pars$Sbasal, 0, 0.2)

	# Draw the first species traits
	traits_mat = matrix(nr = pars$Smax, nc = 3)
	traits_mat[1,] = rand_traits_anc(pars)
	traits_mat = as.data.frame(traits_mat)
	names(traits_mat) = c("n","r","o")

	# Set the presence/absence matrix
	pres = matrix(0, nr = nsteps, nc = pars$Smax)
	pres[1,1] = 1

 	# Set the ancestry object
	anc = matrix(NA,nr = pars$Smax, nc = 3)

	# Set the extinctions object
	extinct = matrix(NA,nr = pars$Smax, nc = 2)

	# Record the matrices
	#L_list = list()

	# Species count
	S = 1

	Stotact <- matrix(NA, nsteps, 3)
	Stotact <- as.data.frame(Stotact)
	names(Stotact) = c("Step","TotalS","ActualS")

	##########################
	# MAIN LOOP
	for(step in 2:nsteps) {
		ActualS = sum(pres[step-1,])
		#cat("Step = ", step-1," Total S = ", S," Actual S = ", ActualS,'\n')
		Stotact[step-1,] <- c(step-1, S, ActualS)

		# Test for successful speciation probability
		for(i in 1:pars$Smax) {

			if(S >= pars$Smax) break

			# Speciation occurs if the species is present
			if(pres[step-1,i] == 1) {

				# Species is maintained
				pres[step, i] = 1

				# Test if there is mutation
				test_number = runif(1, 0, 1)
				#print(paste("random number =", test_number))
				mutation_prob = pars$u_max/(1 + exp(pars$d * (ActualS - pars$I_max)))
				#print(mutation_prob)
				if(test_number < mutation_prob) {

					# Pick new parameters
					traits_mut <- c(1.5, 0,0)
					while(traits_mut[1] < 0 || traits_mut[1] > 1){
						traits_mut = rand_traits_mut(traits_mat[i,], pars)
						#print(paste("new n trait = ", traits_mut[1]))
					}
					#print(paste("mutant traits =", traits_mut))

					# Recompute the interactions
					I = get_L_vec(pars, traits_mat, traits_mut)

					# Compute the number of interactions among present species
					sum_I = sum(I*pres[step,]) + pars$Bspe
					#print(paste("interaction sum =", sum_I))

					# Compute the probability of successful speciation
					if(pars$int == 0){
						spec_prob = pars$u_0neg + pars$u_1neg*exp(-pars$a_uneg * sum_I)
					}

					if(pars$int == 1){
						spec_prob = (pars$u_0 + pars$u_1*exp(-pars$a_u*sum_I))#/(1 + exp(pars$d * (sum_I - pars$I_max)))
						#print(paste("estab prob =",spec_prob, " sum of interactions =", sum_I))
					}
					#print(spec_prob)

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
		if(ActualS == 0) break

		# Compute the interaction matrix among present species
		pres_vec = pres[step,]
		cooc = matrix(pres_vec, nr=pars$Smax, nc=pars$Smax, byrow=TRUE)*
		matrix(pres_vec, nr=pars$Smax, nc=pars$Smax, byrow=FALSE)
		L = get_L_mat(pars, traits_mat)
		L = L*cooc
		#L_list[[step]] = L

		# Test for extinction
		if(int == 0) {
			in_I = colSums(L) + Bext*pres_vec
			ext_prob = e_0neg * (1 - exp(-a_eneg*in_I))
		}

	 	if(pars$int == 1) {
			in_I = colSums(L) + pars$Bext*pres_vec
			ext_prob = pars$e_0pos + pars$e_1pos*exp(-pars$a_epos*in_I)
			#print(paste("extinction probability =", ext_prob[1]))
		}

		#if(int == 2) {
		#	in_I = colSums(L)
		#	out_I = rowSums(L)[(Sbasal+1):(Sbasal+Smax)]
		#	ext_prob = 0 #e_0neg + e_1neg*exp(-a_eneg*out_I) + e_0pos + e_1pos*exp(-a_epos*in_I)
		#}

		# Perform extinctions
		#test_extprob = runif(S,0,1)
		#test_extprob = c(test_extprob, rep(0,pars$Smax - S))
		#print(paste("random number for the extinction prob = ", test_extprob[1]))
		present_spe <- grep(1,pres[step,])
		test_extprob <- rep(0,pars$Smax)
		random_number <- runif(length(present_spe),0,1)
		test_extprob[present_spe] <- random_number
		#cat("S actuel = ", ActualS," longueur random_number = ", length(random_number),'\n')
		pres[step, pres[step-1,] & test_extprob < ext_prob] = 0

		for(i in 1:pars$Smax) {
			if(pres[step,i] != pres[step-1,i] & pres[step,i] == 0){
				extinct[i,] = c(step, i)
			}
		}
	} # End of main loop

	list(pres = pres, traits = traits_mat, anc = anc, extinct = extinct, dynamic = Stotact, L = L) #pres = pres, traits = traits_mat, anc = anc, L_list = L_list, extinct = extinct, dynamic = Stotact
	})
}
