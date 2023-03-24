pres = test$pres
pres_vec = pres[nsteps,]
S = pars$Smax
anc = test$anc
age = nsteps - anc[,1]

# Function to compute the distance between any pair of species

PD = matrix(0,nr = S, nc = S)
ANC = matrix(0,nr = S, nc = S) 
diag(ANC) = c(1:S)

for(i in 1:S) {
#	if(pres_vec[i]==1) {
		for(j in 1:S) {
#			if(pres_vec[j]==1) {
				if(i<j) {
					anc_i = anc[i,2]
					anc_j = anc[j,2] 

					if(anc_i == j) {
						PD[i,j] = nsteps - anc[anc_i,1]
						ANC[i,j] = j
					}

					else if(anc_j == i) {
						PD[i,j] = nsteps - anc[anc_j,1]
						ANC[i,j] = i
					}

					else if(anc_i == anc_j)  {
						if(anc_i == 0) {
							PD[i,j] = nsteps
							ANC[i,j] = anc_i
						}						
						else {
							PD[i,j] = nsteps - anc[anc_i,1]
							ANC[i,j] = anc_i
						}
					}

					else {
						for(k in 1:nsteps) {
							if(anc_i > 0) anc_i = anc[anc_i,2]
							if(anc_j > 0) anc_j = anc[anc_j,2]

							if(anc_i == anc_j) {
								if(anc_i == 0) {
									PD[i,j] = nsteps
									ANC[i,j] = anc_i
								}
								else {
									PD[i,j] = nsteps - anc[anc_i,1]
									ANC[i,j] = anc_i
								}
								break
							}
							cat(anc_i, " ", anc_j,'\n')
						}
					}
				}			
#			}
			cat(i, " ", j, '\n')			
		}
#	}
}
PD = PD + t(PD)
ANC = ANC + t(ANC)
