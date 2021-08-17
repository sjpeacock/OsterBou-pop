###############################################################################
#
# Code to plot figures that visualize output of simulations from PaperSims.R
#
# Author: Stephanie Peacock <stephanie.j.peacock@gmail.com
# Date: August 17, 2021
#
###############################################################################

source("simulations/bouSetup.R")
source("simulations/popFunctionsGround.R")

library(gplots)
library(ggsci)
library(pals)

# Set colour plots for the figures
colPal <- pal_locuszoom()
cols <- colPal(7)
colPal2 <- pal_locuszoom(alpha = 0.5)
cols2 <- colPal2(7) 

###############################################################################
# Which migration path to plot?
###############################################################################

N <- 6

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

###############################################################################
# Look at annual dynamics in year 30 of each simulation (pre-climate change)
###############################################################################

oneYear <- readRDS(paste0("simulations/output/oneYear30_migPath", N, ".rds"))

#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------

# Plot all levels of transmission or just the base level?
allBeta <- FALSE

# Plot resident populations in dashed
plotRes <- FALSE

if(allBeta == TRUE){
	I <- c(1:3)
	quartz(width = 6.3, height = 6, pointsize = 10)
	layout(mat = cbind(
		matrix(c(rep(1:3, each = 3), 16, rep(4:5, each = 3)), ncol = 1),
		5 + matrix(c(rep(1:3, each = 3), 11, rep(4:5, each = 3)), ncol = 1),
		10 + matrix(c(rep(1:3, each = 3), 6, rep(4:5, each = 3)), ncol = 1)))
	
} else{
	I <- 2
	quartz(width = 4, height = 6, pointsize = 10)
	layout(mat =	matrix(c(rep(1:3, each = 3), 6, rep(4:5, each = 3)), ncol = 1))
}

xDate <- as.Date(paste(30, c(1:365), sep = "-"), format = "%y-%j")
yPerc <- 0.2

par(oma = c(4, 5, 4, 1))
par(mar = c(0,4,0,1))
# [m, i, p, 5, j]

for(i in I){
	for(s in 1:5){ # For each stage
		if(allBeta == FALSE){
			ylim.s <- extendrange(oneYear[1, i, , s, ], f = c(0, yPerc))
		} else {
			ylim.s <- extendrange(oneYear[1, , , s, ], f = c(0, yPerc))
		}
		
		if(s == 3|s == 5){
			plot(xDate, oneYear[1, i, 1, 1, ], "n", ylim = ylim.s, yaxt = "n", ylab = "", xaxs = "i")
		}else{
			plot(xDate, oneYear[1, i, 1, 1, ], "n", ylim = ylim.s, yaxt = "n", ylab = "", xaxt = "n", xaxs = "i")
		}
		
		# Plot migration periods as shaded
		u <- par('usr')
		polygon(x = as.Date(paste(2100, c(110, 153, 153, 110), sep = "-"), format = "%Y-%j"), y = u[c(3,3,4,4)], col = grey(0.9), border = NA) # Spring migration
		polygon(x = as.Date(paste(2100, c(169, 180, 180, 169), sep = "-"), format = "%Y-%j"), y = u[c(3,3,4,4)], col = grey(0.9), border = NA) # Post-calving migration
		polygon(x = as.Date(paste(2100, c(251, 355, 355, 251), sep = "-"), format = "%Y-%j"), y = u[c(3,3,4,4)], col = grey(0.9), border = NA) # Fall migration
		
		
		
		for(p in 1:4){
			lines(xDate, oneYear[1, i, p, s, ], col = cols[p], lwd = 1.5)
			if(plotRes == TRUE) lines(xDate, oneYear[2, i, p, s, ], col = cols[p], lty = 2)
		}
		
		if(i == 1 | allBeta == FALSE) mtext(side = 2, line = 1, c("Adult\nparasites", "Arrested\nlarvae", "Developing\nlarvae", "Pre-infective", "Infective")[s])
		if((i == 1 | allBeta == FALSE) & s == 1){
			# legend("topright", col = cols[c(1:4)], title = "Arresting", lwd = 1.5, legend = c("100%", "50%", "0%", "variable"), bg = "white")
		}
		mtext(side = 3, line = -1.5, paste("", LETTERS[s]), adj = 0)
		
		if(s == 1) legend("topright", col = cols[c(1:4)], title = "Arresting", lwd = 1.5, legend = c("100%", "50%", "0%", "variable"), bg = "white", xpd = NA)
	} #  end s
} #end i

u <- par('usr')
segments(x0 = u[1] - 0.25*(u[2] - u[1]), x1 = u[1] - 0.25*(u[2] - u[1]), y0 = 5.5*u[4], y1 = 2.3*u[4], xpd = NA)
text(u[1] - 0.3*(u[2] - u[1]), mean(c(5.5*u[4], 2.3*u[4])), "Within-host parasite stages", xpd = NA, srt = 90, cex = 1.5)
segments(x0 = u[1] - 0.25*(u[2] - u[1]), x1 = u[1] - 0.25*(u[2] - u[1]), y0 = 2*u[4], y1 = 0, xpd = NA)
text(u[1] - 0.3*(u[2] - u[1]), u[4], "Free-living parasite stages", xpd = NA, srt = 90, cex = 1.5)

segments(x0 = as.Date(paste(2100, L4startDOY, sep = "-"), format = "%Y-%j"), x1 = as.Date(paste(2100, L4startDOY, sep = "-"), format = "%Y-%j"), y0 = 0, y1 = 1.6e+14, lty = 3, xpd = NA)
text(as.Date(paste(2100, L4startDOY, sep = "-"), format = "%Y-%j"), 1.6e+14, pos = 3, "L4 resume", xpd = NA)

# segments(x0 = as.Date(paste(2030, breedDOY, sep = "-"), format = "%Y-%j"), x1 = as.Date(paste(2030, breedDOY, sep = "-"), format = "%Y-%j"), y0 = 0, y1 = 1.6e+14, lty = 3, xpd = NA)
# text(as.Date(paste(2030, breedDOY, sep = "-"), format = "%Y-%j"), 1.6e+14, pos = 3, "Calving", xpd = NA)

# legend(u[2] - 70, 1.65e+14, col = cols[c(1:4)], title = "Arresting", lwd = 1.5, legend = c("100%", "50%", "0%", "variable"), bg = "white", xpd = NA)


###############################################################################
# Run simulations for subsequent 100 years
###############################################################################

y <- 100
dt <- 1
startstop <- "agg"

source("simulations/calcGrid.R")

annualSumm <- readRDS(paste0("simulations/output/annualSumm_migPath", N, ".rds"))
V <- readRDS(paste0("simulations/output/V_baseBeta_migPath", N, ".rds")) # Take a looooong time!!


###############################################################################
# Summarize cycles in average parasite burden and host abundance over 100 years
###############################################################################

# all climate scenarios have the same y-axis scale for comparison
ylims <- list(
	P = list(
		lowTrans <- range(annualSumm[[1, 3, 1]][, 1]),
		baseTrans <- range(annualSumm[[2, 3, 1]][, 1]),
		highTrans <- range(annualSumm[[3, 3, 1]][, 1])
	),
	N = list(
		lowTrans <- range(annualSumm[[1, 3, 1]][, 2]),
		baseTrans <- range(annualSumm[[2, 3, 1]][, 2]),
		highTrans <- range(annualSumm[[3, 3, 1]][, 2])
	))

for(i in 1:3){
	for(r in 1:3){
		for(p in 1:4){
			ylims[[1]][[i]] <- range(c(ylims[[1]][[i]], annualSumm[[i, p, r]][, 1]))
			ylims[[2]][[i]] <- range(c(ylims[[2]][[i]], annualSumm[[i, p, r]][, 2]))
		}
	}
}



# pdf(file = "figures/supplement/annualP_allScenarios_zoomed_GT.pdf", width = 8, height = 8, pointsize = 10)
# # quartz(width = 6.3, height = 5, pointsize = 10)
# par(mfrow = c(3,3), mar = rep(0, 4), oma = c(5, 6,3,8))
# for(i in 1:3){ # for each of three levels of transmission
# 	for(r in 1:3){ # for each climate scenario
# 		
# 		# Plot all three inhibition on the same figure using cols[1-3] (red = 100% inhibition)
# 		plot(1:y, annualSumm[[i, 1, r]][, 1], "n", ylim = c(0, 1e8), yaxt = "n", xaxt = "n")#, ylim = ylims[[1]][[i]]
# 		for(p in 1:4)	lines(1:y, annualSumm[[i, p, r]][, 1], col = cols[p], lwd = 1.2)
# 		if(r == 1) axis(side = 2, las = 1) #else axis(side = 2, labels = FALSE)
# 		if(i == 3) axis(side = 1) #else axis(side = 1, labels = FALSE)
# 		
# 		if(i == 1 & r == 1) legend("topleft", col = cols[1], lty = c(1,3), lwd = 1.2, bty = "n", legend = c("Parasite", "Host"))
# 		
# 		# # Plot host population over top in dashed line
# 		par(new = TRUE)
# 		plot(1:y, annualSumm[[i, 1, r]][, 2], "n", ylim = c(0, 8e7), yaxt = "n", xaxt = "n")
# 		for(p in 1:4)	lines(1:y, annualSumm[[i, p, r]][, 2], col = cols[p], lty = 3, lwd = 1.2)
# 
# 		if(i == 1) mtext(side = 3, line = 1, c("Current", "RCP 2.6", "RCP 8.5")[r])
# 		
# 		if(r == 3){
# 			axis(side = 4, las = 1)
# 			mtext(side = 4, c(expression(paste("Low ", beta)), expression(paste("Base ", beta)), expression(paste("High ", beta)))[i], line = 7)
# 		}
# 		
# 		# if(r == 3 & i == 1) legend("right", col = cols[1:3], title = "Inhibition", legend = c("100%", "50%", "0%"), bty = "n", lwd = 1.2)
# 	}
# 	}
# mtext(side = 1, "Year in simulation", outer = TRUE, line = 4)
# mtext(side = 2, "Average annual parasite pressure", outer = TRUE, line = 5)
# mtext(side = 4, "Total host population size", outer = TRUE, line = 5)
# 
# 
# dev.off()

#------------------------------------------------------------------------------
# Plot of base transmission under current conditions for main text
#-----------------------------------------------------------------------------
# 
# i <- 2 # Transmission rate scenario (2 = base)
# r <- 3 # Climate scenario (1 = current)
# 
# #--------
# # Inset plot of all arrested scenarios
# par(mfrow = c(2,1), mar = c(1,6,1,6), oma = rep(0, 4))
# plot(1:y, annualSumm[[i, 1, r]][, 1], "n", ylim = c(0, max(c(annualSumm[[2, 1, 1]][, 1], annualSumm[[2, 2, 1]][, 1], annualSumm[[2, 3, 1]][, 1], annualSumm[[2, 4, 1]][, 1]))), yaxt = "n", ylab = "", xlab = "", xaxs = "i", xlim = c(40, 100))
# axis(side = 2, at = c(0, 0.5e9, 1e9, 1.5e9, 2e9), labels = c(0, 5000, "10,000", "15,000", "20,000"), las = 1)#c(0, expression(0.5%*%10^9), expression(1%*%10^9), expression(1.5%*%10^9), expression(2%*%10^9)), las = 1)
# for(p in 1:4)	lines(1:y, annualSumm[[i, p, r]][, 1], "o", pch = 19, cex = 0.8, col = cols[p])
# 
# par(new = TRUE)
# plot(1:y, annualSumm[[i, 1, r]][, 2], "n", ylim = c(0, max(c(annualSumm[[2, 1, 1]][, 2], annualSumm[[2, 2, 1]][, 2], annualSumm[[2, 3, 1]][, 2], annualSumm[[2, 4, 1]][, 2]))), yaxt = "n", xaxt = "n", xlab = "", ylab = "", xaxs = "i", xlim = c(40, 100))
# for(p in 1:4)	lines(1:y, annualSumm[[i, p, r]][, 2], "o", pch = 21, bg = "white", cex = 0.8, col = cols[p], lty = 2)
# axis(side = 4, las = 1, at = c(0, 4, 8, 12)*10^8, labels = c(0, 400, 800, 1200))
# 
# # Inset plot of last 3 arrested scenarios
# # par(mfrow = c(1,1), mar = c(4,6,2,6), oma = rep(0, 4))
# plot(1:y, annualSumm[[i, 1, r]][, 1], "n", ylim = c(0, max(c(annualSumm[[2, 2, 1]][, 1], annualSumm[[2, 3, 1]][, 1], annualSumm[[2, 4, 1]][, 1]))), yaxt = "n", ylab = "", xlab = "", xaxs = "i", xlim = c(40, 100))
# # axis(side = 2, at = c(0, 0.5e9, 1e9, 1.5e9, 2e9), labels = c(0, 5000, "10,000", "15,000", "20,000"), las = 1)#c(0, expression(0.5%*%10^9), expression(1%*%10^9), expression(1.5%*%10^9), expression(2%*%10^9)), las = 1)
# for(p in 2:4)	lines(1:y, annualSumm[[i, p, r]][, 1], "o", pch = 19, cex = 0.8, col = cols[p])
# 
# par(new = TRUE)
# plot(1:y, annualSumm[[i, 1, r]][, 2], "n", ylim = c(0, max(c(annualSumm[[2, 1, 1]][, 2], annualSumm[[2, 2, 1]][, 2], annualSumm[[2, 3, 1]][, 2], annualSumm[[2, 4, 1]][, 2]))), yaxt = "n", xaxt = "n", xlab = "", ylab = "", xaxs = "i", xlim = c(40, 100))
# for(p in 1:4)	lines(1:y, annualSumm[[i, p, r]][, 2], "o", pch = 21, bg = "white", cex = 0.8, col = cols[p], lty = 2)
# axis(side = 4, las = 1, at = c(0, 4, 8, 12)*10^8, labels = c(0, 400, 800, 1200))
# 
# 
# #-----------
# # Each arrested scenario in its own plot
# ylimCycles <- array(NA, dim = c(4,2,2))
# for(p in 1:4){
# 	ylimCycles[p, 1, ] <- extendrange(annualSumm[[i, p, r]][60:y, 1], f = 0.15)
# 	ylimCycles[p, 2, ] <- extendrange(annualSumm[[i, p, r]][60:y, 2], f = 0.15)
# }
# # ylimCycles[4,1,] <- c(81.5, 82.5)*10^6 
# # ylimCycles[4,2,] <- c(20.6, 20.9)*10^6
# 
# # quartz(width = 4.5, height = 5, pointsize = 10)
# # pdf(file = "figures/popCycles_GT.pdf", width = 6.8, height = 5, pointsize = 10)
# y0 <- 40
# 
# par(mfrow = c(4,1), mar = c(1,4,1,4), oma = c(3,1,1,1))
# 
# for(p in 1:4){
# 	plot(y0:y, annualSumm[[i, p, r]][y0:y, 1], "o", pch = c(21,22,24,25)[p], col = cols[p], bg = cols[p], cex = 0.8, lwd = 1.2, xaxt = "n", ylab = "", xlab = "", yaxt = "n", ylim = ylimCycles[p,1,])
# 	if(p > 3) axis(side = 1, at = seq(y0, 100, 20), labels = seq(100 + y0, 200, 20)) else axis(side = 1, at = seq(y0, 100, 20), labels = FALSE)
# 	
# 	# axis(side = 2, at = pretty(annualSumm[[i, p, r]][y0:y, 1]), labels = pretty(annualSumm[[i, p, r]][y0:y, 1]*10^-6), las = 1)
# 	axis(side = 2, at = pretty(ylimCycles[p,1,]), labels = pretty(ylimCycles[p,1,]*10^-6), las = 1)
# 	
# 	if(p == 1){
# 		arrows(x0 = 76, x1 = 94, y0 = max(annualSumm[[i, p, r]][y0:y, 1]), y1 = max(annualSumm[[i, p, r]][y0:y, 1]), length = 0.06, cod = 3, lwd = 1.2)
# 		text(85, max(annualSumm[[i, p, r]][y0:y, 1]), pos = 1, "18 years")
# 	
# 		} else if(p == 2){
# 		arrows(x0 = 71, x1 = 81, y0 = max(annualSumm[[i, p, r]][y0:y, 1]), y1 = max(annualSumm[[i, p, r]][y0:y, 1]), length = 0.06, cod = 3, lwd = 1.2)
# 		text(76, max(annualSumm[[i, p, r]][y0:y, 1]), pos = 1, "10 years")
# 	
# 		} else if(p == 3){
# 		arrows(x0 = 73, x1 = 83, y0 = max(annualSumm[[i, p, r]][y0:y, 1]), y1 = max(annualSumm[[i, p, r]][y0:y, 1]), length = 0.06, cod = 3, lwd = 1.2)
# 		text(78, max(annualSumm[[i, p, r]][y0:y, 1]), pos = 1, "10 years")
# 	}
# 	
# 	par(new = TRUE)
# 	plot(y0:y, annualSumm[[i, p, r]][y0:y, 2], "o", pch = c(21,22,24,25)[p], col = cols[p], bg = "white", cex = 0.8, lty = 3, yaxt = "n", xaxt = "n", ylim = ylimCycles[p, 2 ,], ylab = "", xlab = "")
# 	
# 	# axis(side = 4, at = pretty(annualSumm[[i, p, r]][y0:y, 2]), labels = pretty(annualSumm[[i, p, r]][y0:y, 2]*10^-6), las = 1)
# 	axis(side = 4, at = pretty(ylimCycles[p,2,]), labels = pretty(ylimCycles[p,2,]*10^-6), las = 1)
# 	
# 	if(p == 1){
# 		arrows(x0 = 80, x1 = 84, y0 = min(annualSumm[[i, p, r]][y0:y, 2]), y1 = min(annualSumm[[i, p, r]][y0:y, 2]), code = 3, length = 0.06, lwd = 1.2)
# 		text(82, min(annualSumm[[i, p, r]][y0:y, 2]),pos =1, "4 years")
# 	} else if(p == 2){
# 		arrows(x0 = 83, x1 = 86, y0 = min(annualSumm[[i, p, r]][y0:y, 2]), y1 = min(annualSumm[[i, p, r]][y0:y, 2]), code = 3, length = 0.06, lwd = 1.2)
# 		text(84.5, min(annualSumm[[i, p, r]][y0:y, 2]),pos =1, "3 years")
# 	} else if(p == 3){
# 		arrows(x0 = 85, x1 = 87, y0 = min(annualSumm[[i, p, r]][y0:y, 2]), y1 = min(annualSumm[[i, p, r]][y0:y, 2]), code = 3, length = 0.06, lwd = 1.2)
# 		text(86, min(annualSumm[[i, p, r]][y0:y, 2]), pos =1, "2 years")
# 	}
# 	
# if(p == 4) legend("bottomright", pch = c(25, 25), pt.bg = c(cols[4], "white"), lty = c(1, 3), lwd = c(1.2, 1), legend = c("Parasite pressure", "Total hosts"), col = cols[4], bty = "n")
# 
# 	mtext(side =3, adj = 0, line = -1.5, paste(" ", LETTERS[p]))
# 	
# 	
# } # end p
# 
# mtext(side = 2, outer = TRUE, expression(paste("Annual parasite pressure (parasite", {}%*%{}, "days (host)", {}^-1 %*% 10^-6, ")")), line = -1)
# mtext(side =4, outer = TRUE, expression(paste("Total host population size (", {}%*%10^-6, ")")), line = -0.5)
# mtext(side = 1, outer = TRUE, "Year in simulation", line = 2)
# 
# dev.off()


###############################################################################
# For climate scenarios, summarize 
# (1) average parasite burdens in last 80 years
# (2) average host population in the last 80 years.
###############################################################################
boot <- function(x){
	bootX <- apply(matrix(sample(x, length(x) * 1000, replace = TRUE), nrow = 1000), 1, mean)
	return(c(quantile(bootX, 0.025), mean(bootX), quantile(bootX, 0.975)))
}

avg80 <- array(NA, dim = c(3, 4, 3, 2, 3))
for(i in 1:3){
	for(r in 1:3){
		for(p in 1:4){
			for(m in 1:2){
				avg80[i, p, r, m, ] <- boot(annualSumm[[i, p, r]][51:100, m])
			}
		}
	}
}


# Rather than looking at last 80 years, maybe look at the last cycle (18 year, 10 years, 1 year)
cycleYrs <- c(18, 10, 10, 1)

avg80 <- array(NA, dim = c(3, 4, 3, 2, 3))
for(i in 1:3){
	for(r in 1:3){
		for(p in 1:4){
			for(m in 1:2){
				avg80[i, p, r, m, ] <- boot(annualSumm[[i, p, r]][c((100 - cycleYrs[p] + 1):100), m])
			}
		}
	}
}


# Plot as a percentage of the mean in current scenarios
# quartz(width = 4.5, height = 4, pointsize = 10)
# ylims2 <- matrix(c(0.9, 1.15, 0.85, 1.1), 2, 2)
ylims2 <- matrix(c(-0.1, 0.12, -0.15, 0.05)*100, 2, 2)

quartz(width = 4.5, height = 4, pointsize = 10)
for(i in 1:3){
	# if(i == 2) pdf(file = "figures/climateChangeEffects.pdf", width = 3.2, height = 4, pointsize = 10)	
	# if(i == 1) pdf(file = "figures/supplement/climateChangeEffects_lowBeta.pdf", width = 3.2, height = 4, pointsize = 10)	
	# if(i == 3) pdf(file = "figures/supplement/climateChangeEffects_highBeta.pdf", width = 3.2, height = 4, pointsize = 10)	
	
	par(mfrow = c(2,1), mar = c(2,5,2,1), oma = c(2, 0, 0, 6))
	for(m in 1:2){
		plot(1:3, avg80[i, 1, , m, 2], "n", ylim = ylims2[, m], xlim = c(0.5, 3.5), xaxs = "i", xaxt = "n", las = 1, xlab = "", ylab = "", yaxt = "n")
		axis(side = 2, at = seq(-10, 10, 5), labels = paste(seq(-10, 10, 5), "%", sep = ""), las = 1)
		axis(side = 1, at = c(0.5, 1.5, 2.5, 3.5), labels = FALSE)
		axis(side = 1, at = c(1, 2, 3), labels = c("Current", "RCP 2.6", "RCP 8.5"), tck = 0)
		abline(h = 0)
		abline(v = c(1.5, 2.5))
		for(p in 1:4){
			segments(x0 = 1:3 + c(-0.3, -0.1, 0.1, 0.3)[p], x1 = 1:3 + c(-0.3, -0.1, 0.1, 0.3)[p], y0 = (avg80[i, p, , m, 1] - avg80[i, p, 1, m, 2])/avg80[i, p, 1, m, 2]*100, y1 = (avg80[i, p, , m, 3] - avg80[i, p, 1, m, 2])/avg80[i, p, 1, m, 2]*100, col = cols2[p], lwd = 2)
			points(1:3 + c(-0.3, -0.1, 0.1, 0.3)[p], (avg80[i, p, , m, 2] - avg80[i, p, 1, m, 2])/avg80[i, p, 1, m, 2]*100, pch = c(21,22,24,25)[p], col = cols[p], bg = cols2[p])
		}
		mtext(side = 3, adj = 0, line = 0, LETTERS[m])
		mtext(side = 2, c("Change in parasite burdens", "Change in host pop'n")[m], line = 3)
	}
	mtext(side = 1, "Climate scenario", line = 3)
	legend(3.8, 42, pch = c(21,22,24,25), col = cols[1:4], pt.bg = cols2[1:4], title = "Arrested\ndevelopment", legend = c("100%", "50%", "0%", "variable"), bty = "n", cex = 0.8, pt.cex = 1.0, xpd= NA)
	
	dev.off()
}
#------------------------------------------------------------------------------
# Does the timing of infection differ and is that why we see declines in host
# abundance with no significant increase in parasite populations?
#------------------------------------------------------------------------------
# Look at relative summary in last year of simulation

year2plot <- 90 #last 80 years in simulation

annual100 <- array(NA, dim = c(3, 4, 3, 6, 365))
for(i in 1:3){ # Three different transmission parameters
	for(p in 1:4){ # for each level of ppnInhibit (1, 0.5, 0)
		for(r in 1:3){ # for current recp26, rcp 85
			V.ipr <- V[[i,p,r]][, , which(timeDat$year[timeDat$year > 80] == year2plot)]
			for(j in 1:365){
				annual100[i, p, r, 1, j] <- sum(V.ipr[c('P_mov', 'P_stat'), , j]) / sum(V.ipr[c('adult_mov', 'adult_stat'), , j])
				annual100[i, p, r, 2, j] <- sum(V.ipr[c('L4A_mov', 'L4A_stat'), , j]) / sum(V.ipr[c('adult_mov', 'adult_stat'), , j])
				annual100[i, p, r, 3, j] <- sum(V.ipr[c('L4_mov', 'L4_stat'), , j]) / sum(V.ipr[c('adult_mov', 'adult_stat'), , j])
				annual100[i, p, r, 4, j] <- sum(V.ipr['L0', , j])
				annual100[i, p, r, 5, j] <- sum(V.ipr['L3', , j])
				annual100[i, p, r, 6, j] <- sum(V.ipr[c('adult_mov', 'adult_stat', 'yearling_mov', 'yearling_stat', 'calf_mov', "calf_stat"), , j])
			} # end day j
			
			
		}}
} # end m


#-----
# Plot 
#-----
xDate <- as.Date(paste(2100, c(1:365), sep = "-"), format = "%Y-%j") 
stages <- c("P", "L4A", "L4", "L0", "L3", "Hosts")
arrestScenario <- c("100% arrest", "50% arrest", "0% arrest", "variable")
i <- 2
p <- 2
# quartz()
par(mfrow = c(4,1), mar = c(1,5,1,1), oma = c(4,0,2,0))
for(s in c(3,1,5,2)){
	plot(xDate[seq(2, 365, 2)], annual100[i, p, 1, s, seq(2, 365, 2)], "l", col= 1, ylim = extendrange(annual100[i, p, 1:3, s, ], f = 0.05), ylab = stages[s], xlab = "")
	lines(xDate[seq(2, 365, 2)], annual100[i, p, 2, s, seq(2, 365, 2)], col= 4)
	lines(xDate[seq(2, 365, 2)], annual100[i, p, 3, s, seq(2, 365, 2)], col= 2)
	if(s == 3) mtext(side = 3, arrestScenario[p])
	abline(v = as.Date(paste("2100", breedDOY - 240 + 365, sep = "-"), format = "%Y-%j"))
	# points(rep(as.Date(paste("2100", breedDOY - 240 + 365, sep = "-"), format = "%Y-%j"), 2))
}




par(mfrow = c(2,1))
s <- 1
plot(xDate, annual100[i, p, 1, s, ], "l", col= 1, ylim = extendrange(annual100[i, p, 1, s, ], f = 0.05))
par(new = TRUE)
plot(xDate, annual100[i, p, 2, s, ], "l", col= cols[4], ylim = extendrange(annual100[i, p, 2, s, ], f = 0.05))
par(new = TRUE)
plot(xDate, annual100[i, p, 3, s, ], "l", col= cols[1], ylim = extendrange(annual100[i, p, 3, s, ], f = 0.05))
abline(v = as.Date(paste("2100", breedDOY - 240 + 365, sep = "-"), format = "%Y-%j"))


preg <- cbind(seq(0, 60000, 100), 0.8 - 1/(1 + exp(7.025 - 0.000328*seq(0, 60000, 100))))

p <- 2


for(s in c(1,6)){
	plot(xDate, annual100[i, p, 1, s, ], "l", col= 1, ylim = extendrange(annual100[i, p, 1:3, s, ], f = 0.05))
	for(r in 2:3) lines(xDate, annual100[i, p, r, s, ], "l", col= cols[c(4,1)[r-1]])
	abline(v = as.Date(paste("2100", breedDOY - 240 + 365, sep = "-"), format = "%Y-%j"), lty = 3)
	abline(v = as.Date(paste("2100", breedDOY, sep = "-"), format = "%Y-%j"), lty = 3)
	if(s == 1) lines(preg[, 2]/max(preg[,2]) * 180))
}




# Look at annual temp and parameters under those scenarios
temp <- seq(-25, 40, 0.1)
ymax <- c(0.6, 0.03, 0.1)

mu0 <- cbind(predict.mu0(tempDat$current), predict.mu0(tempDat$rcp26), predict.mu0(tempDat$rcp85))
mu3 <- cbind(predict.mu3(tempDat$current), predict.mu3(tempDat$rcp26), predict.mu3(tempDat$rcp85))
rho <- cbind(predict.rho0(tempDat$current), predict.rho0(tempDat$rcp26), predict.rho0(tempDat$rcp85))

par(mfcol =c(3,2), oma = c(3, 4, 1, 0), mar = c(2,2,1,1))
plot(temp, predict.mu0(temp), "l", ylim = c(0, ymax[1]), ylab = "Pre-infective mortality", lwd = 2, las =1, col = grey(0.6))

axis(side = 1, labels = FALSE)
plot(temp, predict.mu3(temp), "l", ylim = c(0, ymax[2]), ylab = "Infective mortality", lwd = 2, las =1, col = grey(0.6))

plot(temp, predict.rho0(temp), "l", ylim = c(0, ymax[3]), ylab = "Development", lwd = 2, las =1, col = grey(0.6))

plot(xDate, mu0[, 1], "l", ylim = c(0, ymax[1]), ylab = "", lwd = 2, las =1)
lines(xDate, mu0[, 2], col = 4, lwd = 2)
lines(xDate, mu0[, 3], col = 2, lwd = 2)

plot(xDate, mu3[, 1], "l", ylim = c(0, ymax[2]), ylab = "", lwd = 2, las =1)
lines(xDate, mu3[, 2], col = 4, lwd = 2)
lines(xDate, mu3[, 3], col = 2, lwd = 2)

plot(xDate, rho[, 1], "l", ylim = c(0, ymax[3]), ylab = "", lwd = 2, las =1)
lines(xDate, rho[, 2], col = 4, lwd = 2)
lines(xDate, rho[, 3], col = 2, lwd = 2)


#################################################################################
# Rather than looking at relative, look at absolute changes
#################################################################################

# Look at change in population from first 20 year to last 20 years for each scenario
i <- 2 # Use base transmission

avgDiff <- array(NA, dim = c(4, 3, 2, 2), dimnames = list(c("100% arrest", "50% arrest", "0% arrest", "variable"), c("current", "rcp26", "rcp85"), c("P", "H"), c("first20", "last20")))
pChange <- array(NA, dim = c(4, 3, 2), dimnames = list(c("100% arrest", "50% arrest", "0% arrest", "variable"), c("current", "rcp26", "rcp85"), c("P", "H")))
for(r in 1:3){
	for(p in 1:4){
		for(m in 1:2){
			avgDiff[p, r, m, 1] <- boot(annualSumm[[i, p, r]][26:45, m])[2]
			avgDiff[p, r, m, 2] <- boot(annualSumm[[i, p, r]][81:100, m])[2]
			pChange[p, r, m] <- (avgDiff[p, r, m, 2] - avgDiff[p, r, m, 1])/avgDiff[p, r, m, 1]
		}
	}
}

# Plot
par(mfrow = c(2,1), mar = c(2,5,2,1), oma = c(2, 0, 0, 6))
for(m in 1:2){
	plot(1:3, pChange[1, , m], ylim = range(pChange[, ,m]), "n", xlim = c(0.5, 3.5), xaxs = "i", xaxt = "n", las = 1, xlab = "", ylab = "")
	axis(side = 1, at = c(0.5, 1.5, 2.5, 3.5), labels = FALSE)
	axis(side = 1, at = c(1, 2, 3), labels = c("Current", "RCP 2.6", "RCP 8.5"), tck = 0)
	abline(h = 0)
	abline(v = c(1.5, 2.5))
	for(p in 1:4){
		points(1:3 + c(-0.3, -0.1, 0.1, 0.3)[p], pChange[p, , m], pch = c(21,22,24,25)[p], col = cols[p], bg = cols2[p])
	}
	mtext(side = 3, adj = 0, line = 0, LETTERS[m])
	mtext(side = 2, c("Change in parasite burdens", "Change in host pop'n")[m], line = 3)
}
mtext(side = 1, "Climate scenario", line = 3)
legend(3.8, 42, pch = c(21,22,24,25), col = cols[1:4], pt.bg = cols2[1:4], title = "Arrested\ndevelopment", legend = c("100%", "50%", "0%", "variable"), bty = "n", cex = 0.8, pt.cex = 1.0, xpd= NA)

#################################################################################
# Uptake rate over space and time in year 100
#################################################################################
x.ind <- seq(300, 600, 1)
t <- seq(120, 300, 1)
t.ind <- which(timeDat$DOY[which(timeDat$year > 80)] %in% t & timeDat$year[which(timeDat$year > 80)] == 100)

library(RColorBrewer)
colUptake <- brewer.ylorbr(100)

V.pr <- array(NA, dim = c(4, 3, length(x.ind), length(t)))
for(p in 1:4){
	for(r in 1:3){
		V.pr[p, r, , ] <- V[[p, r]][15, x.ind, t.ind]
	}
}

tDate <- as.Date(paste("2001", c(6:10), "01", sep = "-"), format = "%Y-%m-%d")
tDOY <- as.numeric(strftime(tDate, format = "%j"))

quartz(width = 6.3, height = 5, pointsize = 10)
par(mfrow = c(4,3), mar = c(1, 0, 0, 0), oma = c(3,5,2,7))

for(p in 1:4){
	uptakeBreaks <- seq(0, max(V.pr[p, , ,])*0.5, length.out = 100)
	# uptakeBreaks <- seq(0, c(100000, 300, 150, 25000)[p], length.out = 100)
	
	for(r in 1:3){
		
		plot(rep(x.ind, length(t)), rep(t, each = length(x.ind)), col = colUptake[findInterval(V.pr[p, r, , ], uptakeBreaks)], pch = 15, xaxs = "i", yaxs = "i", xlab = "", ylab = "", yaxt = "n", xaxt = "n", cex = 0.2)
		# if(r == 1) axis(side = 2, las = 1)
		if(r == 1) axis(side = 2, at = tDOY, labels = month.abb[c(6:10)], las= 1)
		if(p == 4) axis(side = 1, at = seq(350, 550, 100), labels = seq(350, 550, 100)*2)
		
		if(r == 3){
			# polygon(x = c(630, 630, 650, 650), y = c(110, 290, 290, 110), border = NA)
			for(i in seq(1, 100, 2)){
				points(c(631:635), rep(seq(120, 290, length.out = 100)[i], 5), pch = 15, cex = 0.5, col = colUptake[i], xpd = NA)
			}
			for(i in c(25, 50, 75, 100)){
				points(632.5, seq(120, 290, length.out = 100)[i], pch = "-", xpd = NA)
				text(635, seq(120, 290, length.out = 100)[i], pos = 4, round(uptakeBreaks[i]), xpd = NA)
			}
		}
		
		if(r == 1) mtext(side = 3, adj = 0, line = -1.5, paste0("  ", LETTERS[p]))
		
		abline(h = breedDOY, lty = 3)
		
		
	}}

mtext(side = 1, "Distance along migration (km)", outer = TRUE, line = 2)
mtext(side = 3, outer = TRUE, "Parasite uptake rate per host per day", line = 0.5)
