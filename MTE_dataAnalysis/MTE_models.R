##########################################################################
# Function to pick model

pickModel <- function(m){
	if(m == 1) model <- model1 else if(m == 2) model <- model2 else if (m==3) model <- model3 else if(m==4) model <- model4 else if (m==5) model <- model5 else if (m==6) model <- model6 else if (m==7) model <- model7 else if (m==8) model <- model8 else if (m==9) model <- model9 else if (m==10) model <- model10   
	return(model)
}

##########################################################################
# Model 1: I I I
##########################################################################
model1 <- function(){
	
	# Priors on hyperparameters
	for(i in 1:nT){
		for(j in 1:3){ #for each parameter mu0, mu1, rho, sigma
			p.mean.log[i,j] ~ dnorm(0, 0.01)
			p.mean[i,j] <- exp(p.mean.log[i,j])
		}}
	
	# Sigma - constant across temperature, but no metabolic theory underpinning
	sigma.log ~ dnorm(0, 0.01)
	for(i in 1:nT){p.mean[i,4] <- exp(sigma.log)}
	
	# Model parameters (u0, u1, rho, sigma) for each temp and replicate
	for(i in 1:nT){
		for(k in 1:nk){
			for(j in 1:4){
				p.rep[i,j,k]  ~ dlnorm(log(p.mean[i,j]), sigp^(-2))}}}
	
	# Starting number
	for(i in 1:nT){
		for(k in 1:nk){
			N0[i,k] ~ dpois(50)}}
	
	#------------------------------------------------------------------
	# Process model
	#------------------------------------------------------------------
	for(i in 1:nT){ # for each temperature
		for(k in 1:nk){ # for each replicate
			# Pre-infectives
			for(j in 1:nt){
				N[i,1,j,k]<-N0[i,k]*exp(-p.rep[i,1,k]*t[j])*(1-phi((log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)/p.rep[i,4,k]))
			}
			# Infectives that develop at time "jj"
			for(j in 1:nt){
				N1_start[i,j,k]<-N0[i,k]/(t[j]*sqrt(2*3.141593*p.rep[i,4,k]^2))*exp(-1/(2*p.rep[i,4,k]^2)*(log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)^2-p.rep[i,1,k]*t[j])
			}
			# ...and survive to time j
			for(j in 1:nt){
				for(jj in 1:j){
					N1_surv[i,j,k,jj]<-N1_start[i,jj,k]*exp(-p.rep[i,2,k]*(t[j]-t[jj]))
				}
				# Summed from 1:j	
				N[i,2,j,k]<-sum(N1_surv[i,j,k,1:j]*dt)
			} #end j
		}}
	
	#----------------------------
	# Data model
	#----------------------------
	for(i in 1:nT){ #for each temperature
		for(k in 1:nk){ # for each replicate
			for(l in 1:3){ # for each of three counts
				n_obs[i,1,k,l,1] ~ dpois(N0[i,k])
				for(s in 1:2){ #for each stage
					for(j in 2:nt_obs[i,k]){
						n_obs[i,j,k,l,s] ~ dpois(max(10^-10, N[i,s,ind[i,j,k],k]))
					}}}}}
} #end model

##########################################################################
# Model 2: Arrhenius for all three
##########################################################################
model2 <- function(){
	for(j in 1:3){
		a[j] ~ dlnorm(log(0.5), pow(1, -2)) 
		E[j] ~ dlnorm(log(0.65), pow(0.5, -2))
	}
	
	Eh <- 0
	Th <- 0
	
	for(i in 1:nT){
		for(j in 1:3){# u0, u1, rho
			p.mean[i,j] <- a[j]*exp(-E[j]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))
		}
	}
	# Sigma - constant across temperature, but no metabolic theory underpinning
	sigma.log ~ dnorm(0, 0.01)
	for(i in 1:nT){p.mean[i,4] <- exp(sigma.log)}
	
	# Model parameters (u0, u1, rho, sigma) for each temp and replicate
	for(i in 1:nT){
		for(k in 1:nk){
			for(j in 1:4){
				p.rep[i,j,k]  ~ dlnorm(log(p.mean[i,j]), sigp^(-2))}}}
	
	# Starting number
	for(i in 1:nT){
		for(k in 1:nk){
			N0[i,k] ~ dpois(50)}}
	
	#------------------------------------------------------------------
	# Process model
	#------------------------------------------------------------------
	for(i in 1:nT){ # for each temperature
		for(k in 1:nk){ # for each replicate
			# Pre-infectives
			for(j in 1:nt){
				N[i,1,j,k]<-N0[i,k]*exp(-p.rep[i,1,k]*t[j])*(1-phi((log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)/p.rep[i,4,k]))
			}
			# Infectives that develop at time "jj"
			for(j in 1:nt){
				N1_start[i,j,k]<-N0[i,k]/(t[j]*sqrt(2*3.141593*p.rep[i,4,k]^2))*exp(-1/(2*p.rep[i,4,k]^2)*(log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)^2-p.rep[i,1,k]*t[j])
			}
			# ...and survive to time j
			for(j in 1:nt){
				for(jj in 1:j){
					N1_surv[i,j,k,jj]<-N1_start[i,jj,k]*exp(-p.rep[i,2,k]*(t[j]-t[jj]))
				}
				# Summed from 1:j	
				N[i,2,j,k]<-sum(N1_surv[i,j,k,1:j]*dt)
			} #end j
		}}
	
	#----------------------------
	# Data model
	#----------------------------
	for(i in 1:nT){ #for each temperature
		for(k in 1:nk){ # for each replicate
			for(l in 1:3){ # for each of three counts
				n_obs[i,1,k,l,1] ~ dpois(N0[i,k])
				for(s in 1:2){ #for each stage
					for(j in 2:nt_obs[i,k]){
						n_obs[i,j,k,l,s] ~ dpois(max(10^-10, N[i,s,ind[i,j,k],k]))
					}}}}}
} #end model

##########################################################################
# Model 3: A A SSU
##########################################################################
model3 <- function(){
	
	# For each parameter (u0, u1, rho), draw base value (a) and activiation 
	# energy (E)
	for(j in 1:3){ #  
		a[j] ~ dlnorm(log(0.5), pow(1, -2)) 
		E[j] ~ dlnorm(log(0.65), pow(0.5, -2))
		}
	
	# Set upper threshold to zero for the two Arrhenius relationships
	for(j in 1:2){
		Eh[j]<-0
		Th[j]<-0
		}
	
	# Draw upper thresholds for the SS relationship of u1
	Eh[3] ~ dlnorm(log(3.25), pow(1, -2)) 
	Th[3] ~ dnorm(23, pow(3, -2))
	
	# Calculate temperature-specific mean parameter value
	for(i in 1:nT){ # for each temperature treatment
		for(j in 1:2){# for the two Arrhenius relationships (u0, u1)
			p.mean[i,j] <- a[j]*exp(-E[j]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))
		}
		
		# For the SS relationships (rho)
		p.mean[i,3] <- a[3]*exp(-E[3]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))*(1+exp(Eh[3]/(8.62*10^-5)*(-1/(T.obs[i]+273.15)+1/(Th[3]+273.15))))^(-1)
	}
	
	# Sigma - constant across temperature; no metabolic theory underpinning
	sigma.log ~ dnorm(0, 0.01)
	for(i in 1:nT){p.mean[i,4] <- exp(sigma.log)}
	
	# Draw random model parameters (u0, u1, rho, sigma) for each temp and replicate
	for(i in 1:nT){ # for each temperature
		for(k in 1:nk){ # for each replicate
			for(j in 1:4){ # for each parameter
				p.rep[i,j,k]  ~ dlnorm(log(p.mean[i,j]), sigp^(-2))}}}
	
	# Draw starting number of eggs in each temp and replicate
	for(i in 1:nT){ # for each temperature
		for(k in 1:nk){ # for each replicate
			N0[i,k] ~ dpois(50)}}
	
	#------------------------------------------------------------------
	# Process model
	#------------------------------------------------------------------
	for(i in 1:nT){ # for each temperature
		for(k in 1:nk){ # for each replicate
			
			# Pre-infectives
			for(j in 1:nt){
				N[i,1,j,k]<-N0[i,k]*exp(-p.rep[i,1,k]*t[j])*(1-phi((log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)/p.rep[i,4,k]))
			}
			
			# Infectives that develop at time "jj"
			for(j in 1:nt){
				N1_start[i,j,k]<-N0[i,k]/(t[j]*sqrt(2*3.141593*p.rep[i,4,k]^2))*exp(-1/(2*p.rep[i,4,k]^2)*(log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)^2-p.rep[i,1,k]*t[j])
			}
			# ...and survive to time j
			for(j in 1:nt){
				for(jj in 1:j){
					N1_surv[i,j,k,jj]<-N1_start[i,jj,k]*exp(-p.rep[i,2,k]*(t[j]-t[jj]))
				}
				# Summed from 1:j	
				N[i,2,j,k]<-sum(N1_surv[i,j,k,1:j]*dt)
			} #end j
		}}
	
	#----------------------------
	# Data model
	#----------------------------
	# n_obs[i,j,k,l,s] is the observed number of  larvae alive in 
	# temperature i, time j, replicate k, count l, and stage s
	
	for(i in 1:nT){ #for each temperature
		for(k in 1:nk){ # for each replicate
			for(l in 1:3){ # for each of three counts
				
				# Likelihood of initial number of eggs
				n_obs[i,1,k,l,1] ~ dpois(N0[i,k])
				
				# Likelihood of preinfective and infective stages
				for(s in 1:2){ # for each stage
					for(j in 2:nt_obs[i,k]){ # for each timestep
						n_obs[i,j,k,l,s] ~ dpois(max(10^-10, N[i,s,ind[i,j,k],k]))
					}}}}}
} #end model

##########################################################################
# Model 4: A C A
##########################################################################
model4 <- function(){
	for(j in 1:3){
		a[j] ~ dlnorm(log(0.5), pow(1, -2)) 
	}
	
	for(j in c(1, 3)){
		E[j] ~ dlnorm(log(0.65), pow(0.5, -2))
	}
	
	E[2] <- 0
	Eh <- 0
	Th <- 0
	
	for(i in 1:nT){
		for(j in c(1,3)){# u0, u1, rho
			p.mean[i,j] <- a[j]*exp(-E[j]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))
		}
		p.mean[i,2] <- a[2]
		}
	# Sigma - constant across temperature; no metabolic theory underpinning
	sigma.log ~ dnorm(0, 0.01)
	for(i in 1:nT){p.mean[i,4] <- exp(sigma.log)}
	
	# Model parameters (u0, u1, rho, sigma) for each temp and replicate
	for(i in 1:nT){
		for(k in 1:nk){
			for(j in 1:4){
				p.rep[i,j,k]  ~ dlnorm(log(p.mean[i,j]), sigp^(-2))}}}
	# Starting number
	for(i in 1:nT){
		for(k in 1:nk){
			N0[i,k] ~ dpois(50)}}
	
	#------------------------------------------------------------------
	# Process model
	#------------------------------------------------------------------
	for(i in 1:nT){ # for each temperature
		for(k in 1:nk){ # for each replicate
			# Pre-infectives
			for(j in 1:nt){
				N[i,1,j,k]<-N0[i,k]*exp(-p.rep[i,1,k]*t[j])*(1-phi((log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)/p.rep[i,4,k]))
			}
			# Infectives that develop at time "jj"
			for(j in 1:nt){
				N1_start[i,j,k]<-N0[i,k]/(t[j]*sqrt(2*3.141593*p.rep[i,4,k]^2))*exp(-1/(2*p.rep[i,4,k]^2)*(log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)^2-p.rep[i,1,k]*t[j])
			}
			# ...and survive to time j
			for(j in 1:nt){
				for(jj in 1:j){
					N1_surv[i,j,k,jj]<-N1_start[i,jj,k]*exp(-p.rep[i,2,k]*(t[j]-t[jj]))
				}
				# Summed from 1:j	
				N[i,2,j,k]<-sum(N1_surv[i,j,k,1:j]*dt)
			} #end j
		}}
	
	#----------------------------
	# Data model
	#----------------------------
	for(i in 1:nT){ #for each temperature
		for(k in 1:nk){ # for each replicate
			for(l in 1:3){ # for each of three counts
				n_obs[i,1,k,l,1] ~ dpois(N0[i,k])
				for(s in 1:2){ #for each stage
					for(j in 2:nt_obs[i,k]){
						n_obs[i,j,k,l,s] ~ dpois(max(10^-10, N[i,s,ind[i,j,k],k]))
					}}}}}
} #end model

##########################################################################
# Model 5: A C SSU
##########################################################################
model5 <- function(){
	for(j in 1:3){
		a[j] ~ dlnorm(log(0.5), pow(1, -2)) 
	}
	
	E[1] ~ dlnorm(log(0.65), pow(0.5, -2))
	Eh[1]<-0
	Th[1]<-0
	
	E[2] <-0
	Eh[2]<-0
	Th[2]<-0
	
	E[3] ~ dlnorm(log(0.65), pow(0.5, -2))
	Eh[3] ~ dlnorm(log(3.25), pow(1, -2)) 
	Th[3] ~ dnorm(23, pow(3, -2))
	
	for(i in 1:nT){
		p.mean[i,1] <- a[1]*exp(-E[1]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))
		p.mean[i,2] <- a[2]
		p.mean[i,3] <- a[3]*exp(-E[3]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))*(1+exp(Eh[3]/(8.62*10^-5)*(-1/(T.obs[i]+273.15)+1/(Th[3]+273.15))))^(-1)
	}
	
	# Sigma - constant across temperature, but no metabolic theory underpinning
	sigma.log ~ dnorm(0, 0.01)
	for(i in 1:nT){p.mean[i,4] <- exp(sigma.log)}
	
	# Model parameters (u0, u1, rho, sigma) for each temp and replicate
	for(i in 1:nT){
		for(k in 1:nk){
			for(j in 1:4){
				p.rep[i,j,k]  ~ dlnorm(log(p.mean[i,j]), sigp^(-2))}}}
	
	# Starting number
	for(i in 1:nT){
		for(k in 1:nk){
			N0[i,k] ~ dpois(50)}}
	
	#------------------------------------------------------------------
	# Process model
	#------------------------------------------------------------------
	for(i in 1:nT){ # for each temperature
		for(k in 1:nk){ # for each replicate
			# Pre-infectives
			for(j in 1:nt){
				N[i,1,j,k]<-N0[i,k]*exp(-p.rep[i,1,k]*t[j])*(1-phi((log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)/p.rep[i,4,k]))
			}
			# Infectives that develop at time "jj"
			for(j in 1:nt){
				N1_start[i,j,k]<-N0[i,k]/(t[j]*sqrt(2*3.141593*p.rep[i,4,k]^2))*exp(-1/(2*p.rep[i,4,k]^2)*(log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)^2-p.rep[i,1,k]*t[j])
			}
			# ...and survive to time j
			for(j in 1:nt){
				for(jj in 1:j){
					N1_surv[i,j,k,jj]<-N1_start[i,jj,k]*exp(-p.rep[i,2,k]*(t[j]-t[jj]))
				}
				# Summed from 1:j	
				N[i,2,j,k]<-sum(N1_surv[i,j,k,1:j]*dt)
			} #end j
		}}
	
	#----------------------------
	# Data model
	#----------------------------
	for(i in 1:nT){ #for each temperature
		for(k in 1:nk){ # for each replicate
			for(l in 1:3){ # for each of three counts
				n_obs[i,1,k,l,1] ~ dpois(N0[i,k])
				for(s in 1:2){ #for each stage
					for(j in 2:nt_obs[i,k]){
						n_obs[i,j,k,l,s] ~ dpois(max(10^-10, N[i,s,ind[i,j,k],k]))
					}}}}}
} #end model

##########################################################################
# Model 6: SSU A A
##########################################################################
model6 <- function(){
	for(j in 1:3){
		a[j] ~ dlnorm(log(0.5), pow(1, -2))
		E[j] ~ dlnorm(log(0.65), pow(0.5, -2))
	}
	
	Eh[1] ~ dlnorm(log(3.25), pow(1, -2)) 
	Th[1] ~ dnorm(23, pow(3, -2))
	
	Eh[2]<-0
	Th[2]<-0
	
	Eh[3]<-0
	Th[3]<-0
	
	
	z<-c(0,1,-1)
	for(i in 1:nT){
		p.mean[i,1] <- a[1]*exp(-E[1]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))*(1+exp(Eh[1]/(8.62*10^-5)*(-1/(T.obs[i]+273.15)+1/(Th[1]+273.15))))^(1)
		for(j in 2:3){
			p.mean[i,j] <- a[j]*exp(-E[j]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))
		}
	}
	
	# Sigma - constant across temperature, but no metabolic theory underpinning
	sigma.log ~ dnorm(0, 0.01)
	for(i in 1:nT){p.mean[i,4] <- exp(sigma.log)}
	
	# Model parameters (u0, u1, rho, sigma) for each temp and replicate
	for(i in 1:nT){
		for(k in 1:nk){
			for(j in 1:4){
				p.rep[i,j,k]  ~ dlnorm(log(p.mean[i,j]), sigp^(-2))}}}
	
	# Starting number
	for(i in 1:nT){
		for(k in 1:nk){
			N0[i,k] ~ dpois(50)}}
	
	#------------------------------------------------------------------
	# Process model
	#------------------------------------------------------------------
	for(i in 1:nT){ # for each temperature
		for(k in 1:nk){ # for each replicate
			# Pre-infectives
			for(j in 1:nt){
				N[i,1,j,k]<-N0[i,k]*exp(-p.rep[i,1,k]*t[j])*(1-phi((log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)/p.rep[i,4,k]))
			}
			# Infectives that develop at time "jj"
			for(j in 1:nt){
				N1_start[i,j,k]<-N0[i,k]/(t[j]*sqrt(2*3.141593*p.rep[i,4,k]^2))*exp(-1/(2*p.rep[i,4,k]^2)*(log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)^2-p.rep[i,1,k]*t[j])
			}
			# ...and survive to time j
			for(j in 1:nt){
				for(jj in 1:j){
					N1_surv[i,j,k,jj]<-N1_start[i,jj,k]*exp(-p.rep[i,2,k]*(t[j]-t[jj]))
				}
				# Summed from 1:j	
				N[i,2,j,k]<-sum(N1_surv[i,j,k,1:j]*dt)
			} #end j
		}}
	
	#----------------------------
	# Data model
	#----------------------------
	for(i in 1:nT){ #for each temperature
		for(k in 1:nk){ # for each replicate
			for(l in 1:3){ # for each of three counts
				n_obs[i,1,k,l,1] ~ dpois(N0[i,k])
				for(s in 1:2){ #for each stage
					for(j in 2:nt_obs[i,k]){
						n_obs[i,j,k,l,s] ~ dpois(max(10^-10, N[i,s,ind[i,j,k],k]))
					}}}}}
} #end model



##########################################################################
# Model 7: SSU A SSU
##########################################################################
model7 <- function(){
	for(j in 1:3){
		a[j] ~ dlnorm(log(0.5), pow(1, -2)) 
		E[j] ~ dlnorm(log(0.65), pow(0.5, -2))
	}
	
	for(j in c(1,3)){
		Eh[j] ~ dlnorm(log(3.25), pow(1, -2)) 
		Th[j] ~ dnorm(23, pow(3, -2))
	}
	
	Eh[2] <- 0
	Th[2] <- 0
	
	z <- c(1, 1, -1) #Exponent in relationship
	for(i in 1:nT){
		for(j in c(1,3)){ # u0, u1, rho
				
				p.mean[i,j] <- a[j]*exp(-E[j]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))*(1+exp(Eh[j]/(8.62*10^-5)*(-1/(T.obs[i]+273.15)+1/(Th[j]+273.15))))^z[j]
				
		}
		
	p.mean[i,2] <- a[2]*exp(-E[2]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))
			
	} # end nT
	
	
	# Sigma - constant across temperature, but no metabolic theory underpinning
	sigma.log ~ dnorm(0, 0.01)
	for(i in 1:nT){p.mean[i,4] <- exp(sigma.log)}
	
	# Model parameters (u0, u1, rho, sigma) for each temp and replicate
	for(i in 1:nT){
		for(k in 1:nk){
			for(j in 1:4){
				p.rep[i,j,k]  ~ dlnorm(log(p.mean[i,j]), sigp^(-2))
			}
		}
	}
	
	# Starting number
	for(i in 1:nT){
		for(k in 1:nk){
			N0[i,k] ~ dpois(50)
		}
	}
	
	#------------------------------------------------------------------
	# Process model
	#------------------------------------------------------------------
	for(i in 1:nT){ # for each temperature
		for(k in 1:nk){ # for each replicate
			# Pre-infectives
			for(j in 1:nt){
				N[i,1,j,k]<-N0[i,k]*exp(-p.rep[i,1,k]*t[j])*(1-phi((log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)/p.rep[i,4,k]))
			}
			# Infectives that develop at time "jj"
			for(j in 1:nt){
				N1_start[i,j,k]<-N0[i,k]/(t[j]*sqrt(2*3.141593*p.rep[i,4,k]^2))*exp(-1/(2*p.rep[i,4,k]^2)*(log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)^2-p.rep[i,1,k]*t[j])
			}
			# ...and survive to time j
			for(j in 1:nt){
				for(jj in 1:j){
					N1_surv[i,j,k,jj]<-N1_start[i,jj,k]*exp(-p.rep[i,2,k]*(t[j]-t[jj]))
				}
				# Summed from 1:j	
				N[i,2,j,k]<-sum(N1_surv[i,j,k,1:j]*dt)
			} #end j
		}}
	
	#----------------------------
	# Data model
	#----------------------------
	for(i in 1:nT){ #for each temperature
		for(k in 1:nk){ # for each replicate
			for(l in 1:3){ # for each of three counts
				n_obs[i,1,k,l,1] ~ dpois(N0[i,k])
				for(s in 1:2){ #for each stage
					for(j in 2:nt_obs[i,k]){
						n_obs[i,j,k,l,s] ~ dpois(max(10^-10, N[i,s,ind[i,j,k],k]))
					}}}}}
} #end model

##########################################################################
# Model 8: SSU C A
##########################################################################
model8 <- function(){
	for(j in 1:3){
		a[j] ~ dlnorm(log(0.5), pow(1, -2))
	}
	
	for(j in c(1,3)){
		E[j] ~ dlnorm(log(0.65), pow(0.5, -2))
	}
	
	Eh[1] ~ dlnorm(log(3.25), pow(1, -2)) 
	Th[1] ~ dnorm(23, pow(3, -2))
	
	E[2] <- 0
	Eh[2] <- 0
	Th[2] <- 0
	
	Eh[3]<-0
	Th[3]<-0
	
	
	for(i in 1:nT){
		p.mean[i,1] <- a[1]*exp(-E[1]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))*(1+exp(Eh[1]/(8.62*10^-5)*(-1/(T.obs[i]+273.15)+1/(Th[1]+273.15))))^(1)
		p.mean[i,2] <- a[2]
		p.mean[i,3] <- a[3]*exp(-E[3]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))
		}
	
	# Sigma - constant across temperature, but no metabolic theory underpinning
	sigma.log ~ dnorm(0, 0.01)
	for(i in 1:nT){p.mean[i,4] <- exp(sigma.log)}
	
	# Model parameters (u0, u1, rho, sigma) for each temp and replicate
	for(i in 1:nT){
		for(k in 1:nk){
			for(j in 1:4){
				p.rep[i,j,k]  ~ dlnorm(log(p.mean[i,j]), sigp^(-2))}}}
	
	# Starting number
	for(i in 1:nT){
		for(k in 1:nk){
			N0[i,k] ~ dpois(50)}}
	
	#------------------------------------------------------------------
	# Process model
	#------------------------------------------------------------------
	for(i in 1:nT){ # for each temperature
		for(k in 1:nk){ # for each replicate
			# Pre-infectives
			for(j in 1:nt){
				N[i,1,j,k]<-N0[i,k]*exp(-p.rep[i,1,k]*t[j])*(1-phi((log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)/p.rep[i,4,k]))
			}
			# Infectives that develop at time "jj"
			for(j in 1:nt){
				N1_start[i,j,k]<-N0[i,k]/(t[j]*sqrt(2*3.141593*p.rep[i,4,k]^2))*exp(-1/(2*p.rep[i,4,k]^2)*(log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)^2-p.rep[i,1,k]*t[j])
			}
			# ...and survive to time j
			for(j in 1:nt){
				for(jj in 1:j){
					N1_surv[i,j,k,jj]<-N1_start[i,jj,k]*exp(-p.rep[i,2,k]*(t[j]-t[jj]))
				}
				# Summed from 1:j	
				N[i,2,j,k]<-sum(N1_surv[i,j,k,1:j]*dt)
			} #end j
		}}
	
	#----------------------------
	# Data model
	#----------------------------
	for(i in 1:nT){ #for each temperature
		for(k in 1:nk){ # for each replicate
			for(l in 1:3){ # for each of three counts
				n_obs[i,1,k,l,1] ~ dpois(N0[i,k])
				for(s in 1:2){ #for each stage
					for(j in 2:nt_obs[i,k]){
						n_obs[i,j,k,l,s] ~ dpois(max(10^-10, N[i,s,ind[i,j,k],k]))
					}}}}}
} #end model

##########################################################################
# Model 9: SSU C SSU
##########################################################################
model9 <- function(){
	for(j in 1:3){
		a[j] ~ dlnorm(log(0.5), pow(1, -2))
	}
	
	for(j in c(1,3)){
		E[j] ~ dlnorm(log(0.65), pow(0.5, -2))
		Eh[j] ~ dlnorm(log(3.25), pow(1, -2)) 
		Th[j] ~ dnorm(23, pow(3, -2))
	}

	E[2] <- 0
	Eh[2] <- 0
	Th[2] <- 0
	
	z<-c(1, 1, -1)
	for(i in 1:nT){
		for(j in c(1,3)){
			p.mean[i,j] <- a[j]*exp(-E[j]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))*(1+exp(Eh[j]/(8.62*10^-5)*(-1/(T.obs[i]+273.15)+1/(Th[j]+273.15))))^(z[j])
			}
		p.mean[i,2] <- a[2]
	}
	
# Sigma - constant across temperature, but no metabolic theory underpinning
sigma.log ~ dnorm(0, 0.01)
for(i in 1:nT){p.mean[i,4] <- exp(sigma.log)}

# Model parameters (u0, u1, rho, sigma) for each temp and replicate
for(i in 1:nT){
	for(k in 1:nk){
		for(j in 1:4){
			p.rep[i,j,k]  ~ dlnorm(log(p.mean[i,j]), sigp^(-2))}}}

# Starting number
for(i in 1:nT){
	for(k in 1:nk){
		N0[i,k] ~ dpois(50)}}

#------------------------------------------------------------------
# Process model
#------------------------------------------------------------------
for(i in 1:nT){ # for each temperature
	for(k in 1:nk){ # for each replicate
		# Pre-infectives
		for(j in 1:nt){
			N[i,1,j,k]<-N0[i,k]*exp(-p.rep[i,1,k]*t[j])*(1-phi((log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)/p.rep[i,4,k]))
		}
		# Infectives that develop at time "jj"
		for(j in 1:nt){
			N1_start[i,j,k]<-N0[i,k]/(t[j]*sqrt(2*3.141593*p.rep[i,4,k]^2))*exp(-1/(2*p.rep[i,4,k]^2)*(log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)^2-p.rep[i,1,k]*t[j])
		}
		# ...and survive to time j
		for(j in 1:nt){
			for(jj in 1:j){
				N1_surv[i,j,k,jj]<-N1_start[i,jj,k]*exp(-p.rep[i,2,k]*(t[j]-t[jj]))
			}
			# Summed from 1:j	
			N[i,2,j,k]<-sum(N1_surv[i,j,k,1:j]*dt)
		} #end j
	}}

#----------------------------
# Data model
#----------------------------
for(i in 1:nT){ #for each temperature
	for(k in 1:nk){ # for each replicate
		for(l in 1:3){ # for each of three counts
			n_obs[i,1,k,l,1] ~ dpois(N0[i,k])
			for(s in 1:2){ #for each stage
				for(j in 2:nt_obs[i,k]){
					n_obs[i,j,k,l,s] ~ dpois(max(10^-10, N[i,s,ind[i,j,k],k]))
				}}}}}
} #end model

##########################################################################
# Model 10: A SSU SSU
##########################################################################
model10 <- function(){
	for(j in 1:3){
		a[j] ~ dlnorm(log(0.5), pow(1, -2)) 
	}
	
	E[1] ~ dlnorm(log(0.65), pow(0.5, -2))
	Eh[1]<-0
	Th[1]<-0
	
	E[2] ~ dlnorm(log(0.65), pow(0.5, -2))
	Eh[2] ~ dlnorm(log(3.25), pow(1, -2))
	Th[2] ~ dnorm(23, pow(3, -2))
	
	E[3] ~ dlnorm(log(0.65), pow(0.5, -2))
	Eh[3] ~ dlnorm(log(3.25), pow(1, -2)) 
	Th[3] ~ dnorm(23, pow(3, -2))
	
	for(i in 1:nT){
		p.mean[i,1] <- a[1]*exp(-E[1]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))
		p.mean[i,2] <- a[2]*exp(-E[2]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))*(1+exp(Eh[2]/(8.62*10^-5)*(-1/(T.obs[i]+273.15)+1/(Th[2]+273.15))))
		p.mean[i,3] <- a[3]*exp(-E[3]/(8.62*10^-5)*(1/(T.obs[i]+273.15)-1/(15+273.15)))*(1+exp(Eh[3]/(8.62*10^-5)*(-1/(T.obs[i]+273.15)+1/(Th[3]+273.15))))^(-1)
	}
	
	# Sigma - constant across temperature, but no metabolic theory underpinning
	sigma.log ~ dnorm(0, 0.01)
	for(i in 1:nT){p.mean[i,4] <- exp(sigma.log)}
	
	# Model parameters (u0, u1, rho, sigma) for each temp and replicate
	for(i in 1:nT){
		for(k in 1:nk){
			for(j in 1:4){
				p.rep[i,j,k]  ~ dlnorm(log(p.mean[i,j]), sigp^(-2))}}}
	
	# Starting number
	for(i in 1:nT){
		for(k in 1:nk){
			N0[i,k] ~ dpois(50)}}
	
	#------------------------------------------------------------------
	# Process model
	#------------------------------------------------------------------
	for(i in 1:nT){ # for each temperature
		for(k in 1:nk){ # for each replicate
			# Pre-infectives
			for(j in 1:nt){
				N[i,1,j,k]<-N0[i,k]*exp(-p.rep[i,1,k]*t[j])*(1-phi((log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)/p.rep[i,4,k]))
			}
			# Infectives that develop at time "jj"
			for(j in 1:nt){
				N1_start[i,j,k]<-N0[i,k]/(t[j]*sqrt(2*3.141593*p.rep[i,4,k]^2))*exp(-1/(2*p.rep[i,4,k]^2)*(log(t[j]*p.rep[i,3,k])+p.rep[i,4,k]^2/2)^2-p.rep[i,1,k]*t[j])
			}
			# ...and survive to time j
			for(j in 1:nt){
				for(jj in 1:j){
					N1_surv[i,j,k,jj]<-N1_start[i,jj,k]*exp(-p.rep[i,2,k]*(t[j]-t[jj]))
				}
				# Summed from 1:j	
				N[i,2,j,k]<-sum(N1_surv[i,j,k,1:j]*dt)
			} #end j
		}}
	
	#----------------------------
	# Data model
	#----------------------------
	for(i in 1:nT){ #for each temperature
		for(k in 1:nk){ # for each replicate
			for(l in 1:3){ # for each of three counts
				n_obs[i,1,k,l,1] ~ dpois(N0[i,k])
				for(s in 1:2){ #for each stage
					for(j in 2:nt_obs[i,k]){
						n_obs[i,j,k,l,s] ~ dpois(max(10^-10, N[i,s,ind[i,j,k],k]))
					}}}}}
} #end model
