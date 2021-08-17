###############################################################################
#
# Code to run simulations of host-parasite PDEs and summarize results
#
# Note that these simulations are not run in parallel because they use a large 
# amount of memory. Output is saved as .rds files that are called in PaperFigs.R
#
# Author: Stephanie Peacock <stephanie.j.peacock@gmail.com
# Date: August 17, 2021
#
###############################################################################

source("simulations/bouSetup.R")
source("simulations/popFunctions.R")

# For three different migration paths
# (Could choose any of the 10 different paths that were simulated and saved.)
for(N in c(6, 1, 3)){
	
	###############################################################################
	# Set up grid and initial conditions
	###############################################################################
	
	# Read in temperature data including CC scenarios
	# Can choose different temperature profiles depending on migration route
	tempDat <- readRDS(paste0("tempData/allTemps_migPath", N, ".rds"))
	
	# Starting abundance of parasites in adults
	# From Bathurst surveys, appprox. 
	initP <- 914
	
	y <- 100
	dt <- 1
	startstop <- "agg"
	
	source("simulations/calcGrid.R")
	
	Bou0 <- rbind(
		calf_mov = rep(0, length(x)),
		yearling_mov = rep(0, length(x)),
		adult_mov = rep(0, length(x)),
		L4_mov = rep(0, length(x)),
		L4A_mov = rep(0, length(x)),
		P_mov = rep(0, length(x)),
		calf_stat = initDist(totPop = bouDist1985['calf'], x),
		yearling_stat = initDist(totPop = bouDist1985['yearling'], x),
		adult_stat = initDist(totPop = c(bouDist1985['cow'] + bouDist1985['bull']), x),
		L4_stat = rep(0, length(x)),
		L4A_stat = rep(0, length(x)),
		P_stat = initP * initDist(totPop = 10^5 * c(bouDist1985['cow'] + bouDist1985['bull'])/sum(bouDist1985), x),
		L0 = rep(0, length(x)),
		L3 = rep(0, length(x)),
		L_uptake = rep(0, length(x))
	)
	
	# Simulate a 100-year "burnin" for each level of transmission (3) and parasite inhibition (4)
	V.init <- list(); length(V.init) <- 12; dim(V.init) <- c(3, 4)
	initBou <- list(); length(initBou) <- 12; dim(initBou) <- c(3, 4)
	
	# Want to store final variables but also mean parasite burdens in October 10
	# (breedDOY - 240) which will influence year 1 pregnancy rates in subsequent simulations
	n <- which(timeDat$DOY == 283 & timeDat$year == 100)
	
	# Only keep the last 20 years, otherwise too much memory!
	for(i in 1:3){ # for 3 levels of transmission
		for(p in 1:4){ # for 3 levels of parasite inhibition
			
			V.dum <- simBou(
				initBou = Bou0, 
				temp = tempDat$current,
				ppnInhibit = pArrestCalc(p),
				transmission = c("low", "base", "high")[i])
			
			V.init[[i, p]] <- V.dum[, , which(timeDat$year > 80)]
			
			initBou[[i, p]] <- list(
				var = V.dum[, , dim(V.dum)[3]],
				OctP = sum(V.dum[c('P_stat', 'P_mov'), , n]) / sum(V.dum[c('adult_stat', 'adult_mov'), , n])
			)
			
		}}
	
	# # Check
	# i <- 2; p <- 1
	# par(mfrow = c(3,1))
	# for(j in 1:3){
	# 	J <- grep(c("adult_", "P_", "L3")[j], dimnames(initBou[[i,p]]$var)[[1]])
	# 	plot(x, initBou[[i,p]]$var[J[1], ], "l", ylim = range(initBou[[i,p]]$var[J, ]), main = c("adult_", "P_", "L3")[j])
	# 	if(length(J) > 1) lines(x, initBou[[i,p]]$var[J[2], ], col = 4)
	# }
	
	saveRDS(initBou, file = paste0("simulations/output/initBou_migPath", N, ".rds"))
	# initBou0 <- readRDS(paste0("simulations/output/initBou_migPath", N, ".rds"))
	
	
	###############################################################################
	# Look at annual dynamics in year 30 of each simulation
	###############################################################################
	# Note that the broad patterns don't change as we run the simulations for longer,
	# just the scale of the 100% inhibition parasite burdens is much higher and
	# thus comparisons are more difficult.
	# Similarly the annual patterns are the same in the peak and trough of the longer
	# term cycles, it's just the scale that shifts.
	y <- 30
	dt <- 1
	startstop <- "agg"
	
	source("simulations/calcGrid.R")
	
	V.annual <- list(); length(V.annual) <- 24; dim(V.annual) <- c(2, 3, 4)
	for(m in 1:2){ # for migratory and non-migratory simulations
		for(i in 1:3){ # for 3 levels of transmission
			for(p in 1:4){ # for 4 levels of parasite inhibition
				V.dum <- simBou(
					initBou = Bou0, 
					temp = tempDat$current,
					ppnInhibit = pArrestCalc(p),
					migSpeed = c(14, 0)[m],
					transmission = c("low", "base", "high")[i])
				
				V.annual[[m, i, p]] <- V.dum[, , which(timeDat$year == 30)]
			}}
	}
	
	#------------------------------------------------------------------------------
	# Summarize annual dynamics
	#------------------------------------------------------------------------------
	
	oneYear <- array(NA, dim = c(2, 3, 4, 5, 365))
	# dimensions = beta (3) * inhibition (4) * variable * day
	
	for(m in 1:2){ # for migratory and resident
		for(i in 1:3){
			for(p in 1:4){ # for each level of ppnInhibit (1, 0.5, 0)
				for(j in 1:365){
					oneYear[m, i, p, 1, j] <- sum(V.annual[[m, i, p]][c('P_mov', 'P_stat'), , j]) / sum(V.annual[[m, i, p]][c('adult_mov', 'adult_stat'), , j])
					oneYear[m, i, p, 2, j] <- sum(V.annual[[m, i, p]][c('L4A_mov', 'L4A_stat'), , j]) / sum(V.annual[[m, i, p]][c('adult_mov', 'adult_stat'), , j])
					oneYear[m, i, p, 3, j] <- sum(V.annual[[m, i, p]][c('L4_mov', 'L4_stat'), , j]) / sum(V.annual[[m, i, p]][c('adult_mov', 'adult_stat'), , j])
					oneYear[m, i, p, 4, j] <- sum(V.annual[[m, i, p]]['L0', , j])
					oneYear[m, i, p, 5, j] <- sum(V.annual[[m, i, p]]['L3', , j])
				} # end day j
				
			}}
	} # end m
	
	saveRDS(oneYear, file = paste0("simulations/output/oneYear30_migPath", N, ".rds"))
	
	###############################################################################
	# Run simulations for subsequent 100 years
	###############################################################################
	
	y <- 100
	dt <- 1
	startstop <- "agg"
	
	source("simulations/calcGrid.R")
	
	# Store AnnualP and Annual N
	# Store full dynamics for the last 20 years, so peak and trough can be further 
	# explored
	
	# Simulate a 100-year simulation for each 
	# r: climate scenario (NULL, RCP2.6, RCP8.5), 
	# i: level of transmission (low, base, high), and 
	# p: parasite inhibition (1, 0.5, 0, variable).
	
	V <- list(); length(V) <- 12; dim(V) <- c(4, 3)
	annualSumm <- list(); length(annualSumm) <- 36; dim(annualSumm) <- c(3, 4, 3)
	
	# Only keep the last 20 years, otherwise too much memory!
	
	for(r in 1:3){ # for three climate scenarios: past, RCP 2.6, and RCP 8.5
		for(i in 1:3){ # for 3 levels of transmission
			for(p in 1:4){ # for 3 levels of parasite inhibition
				V.dum <- simBou(
					initBou = initBou[[i, p]]$var, 
					temp = tempDat[[r]],
					ppnInhibit = pArrestCalc(p),
					transmission = c("low", "base", "high")[i], 
					OctP = initBou[[i, p]]$OctP)
				
				# Store last 20 years of each simulation
				if(i == 2) V[[p, r]] <- V.dum[, , which(timeDat$year == 100)]
				
				# Calculate AnnualP and AnnualN
				annualSumm[[i, p, r]] <- array(NA, dim = c(y, 3), dimnames = list(c(1:y), c("annualP", "annualN", "breedP")))
				
				for(n in 1:y){
					V.n <- V.dum[, , which(timeDat$year == n)]
					Pavg <- numeric(365)
					Ntot <- numeric(365)
					for(j in 1:365){
						# mov <- which(V.n['adult_mov', , j] > 10e-10)
						# stat <- which(V.n['adult_stat', , j] > 10e-10)
						# 
						# if(length(mov) == 0){
						# 	Pavg[j] <- sum(V.n['P_stat', stat, j])/ sum(V.n['adult_stat', stat, j]) 
						# } else if(length(stat) == 0){
						# 	Pavg[j] <- sum(V.n['P_mov', mov, j])/ sum(V.n['adult_mov', mov, j]) 
						# } else {
						# 	Pavg[j] <- sum(V.n['P_mov', mov, j])/ sum(V.n['adult_mov', mov, j]) + sum(V.n['P_stat', stat, j])/ sum(V.n['adult_stat', stat, j])
						# }
						Pavg[j] <- sum(V.n[c('P_mov', "P_stat"), , j])/ sum(V.n[c('adult_mov', "adult_stat"), , j])
						Ntot[j] <- sum(V.n[c('adult_mov', 'adult_stat'), , j])
					} # end day j
					
					annualSumm[[i, p, r]][n, 1] <- sum(Pavg, na.rm = TRUE)
					annualSumm[[i, p, r]][n, 2] <- sum(V.n[c('adult_mov', 'yearling_mov', "calf_mov", "adult_stat", "yearling_stat", "calf_stat"), , 159]) * dx
					annualSumm[[i, p, r]][n, 3] <- sum(Pavg, na.rm = TRUE)
					
				} # end n
				
			} # end r climate scenario
		} # end P
	} # end i
	
	
	saveRDS(annualSumm, file = paste0("simulations/output/annualSumm_migPath", N, ".rds"))
	saveRDS(V, file = paste0("simulations/output/V_baseBeta_migPath", N, ".rds")) # Take a looooong time!!
	
} # End three different migration paths
