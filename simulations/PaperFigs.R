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
N <- 6
oneYear <- readRDS(paste0("simulations/output/oneYear30_migPath", N, ".rds"))
# dimensions = mig/non-mig (2) * beta (3) * inhibition (4) * variable * day

#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------

# Plot all levels of transmission or just the base level?
allBeta <- TRUE

# Plot resident populations in dashed
plotRes <- TRUE

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


if(allBeta == FALSE){
	par(oma = c(4, 5, 4, 1))
	par(mar = c(0,4,0,1))
} else {
	par(oma = c(4, 9, 4, 1))
	par(mar = c(0,0,0,1))
}

for(i in I){
	for(s in 1:5){ # For each stage
		if(allBeta == FALSE){
			ylim.s <- extendrange(oneYear[1, i, , s, ], f = c(0, yPerc))
		} else {
			ylim.s <- extendrange(oneYear[1, 2, , s, ], f = c(0, yPerc))
		}
		
		if(s == 3|s == 5){
			plot(xDate, oneYear[1, i, 1, 1, ], "n", ylim = ylim.s, yaxt = "n", ylab = "", xaxs = "i")
		}else{
			plot(xDate, oneYear[1, i, 1, 1, ], "n", ylim = ylim.s, yaxt = "n", ylab = "", xaxt = "n", xaxs = "i")
		}
		
		# Plot migration periods as shaded
		u <- par('usr')
		polygon(x = as.Date(paste(2300, c(110, 153, 153, 110), sep = "-"), format = "%Y-%j"), y = u[c(3,3,4,4)], col = grey(0.9), border = NA) # Spring migration
		polygon(x = as.Date(paste(2300, c(169, 180, 180, 169), sep = "-"), format = "%Y-%j"), y = u[c(3,3,4,4)], col = grey(0.9), border = NA) # Post-calving migration
		polygon(x = as.Date(paste(2300, c(251, 355, 355, 251), sep = "-"), format = "%Y-%j"), y = u[c(3,3,4,4)], col = grey(0.9), border = NA) # Fall migration
		
		
		if(s == 1 & i == min(I)) legend("topright", lty = 1:4, legend = paste("Scenario", c(1:4)), bg = "white", xpd = NA, bty = "n")
		
		for(p in 1:4){
			lines(xDate, oneYear[1, i, p, s, ], lty = p)
			if(plotRes == TRUE) lines(xDate, oneYear[2, i, p, s, ], lty = p, col = "#00000050")
		}
		
		if(i == 1 | allBeta == FALSE) mtext(side = 2, line = 1, c("Adult\nparasites", "Arrested\nlarvae", "Developing\nlarvae", "Pre-infective", "Infective")[s])
		if((i == 1 | allBeta == FALSE) & s == 1){
		}
		if(allBeta == FALSE){
			mtext(side = 3, line = -1.5, paste("", LETTERS[s]), adj = 0)
			if(s == 5) legend(21980, 1.6e14, fill = c(1, grey(0.5)), border = NA, xpd = NA, bty = "n", legend = c("Migratory", "Non-migratory"), ncol = 2)
		} else {
			if(i == 1 & s == 5) legend(21980, 1.6e14, fill = c(1, grey(0.5)), border = NA, xpd = NA, bty = "n", legend = c("Migratory", "Non-migratory"), ncol = 2)
		}
		
		
	} #  end s
	
	if(i == min(I)){
		u <- par('usr')
		segments(x0 = u[1] - 0.25*(u[2] - u[1]), x1 = u[1] - 0.25*(u[2] - u[1]), y0 = 5.5*u[4], y1 = 2.3*u[4], xpd = NA)
		text(u[1] - 0.3*(u[2] - u[1]), mean(c(5.5*u[4], 2.3*u[4])), "Within-host parasite stages", xpd = NA, srt = 90, cex = 1.5)
		segments(x0 = u[1] - 0.25*(u[2] - u[1]), x1 = u[1] - 0.25*(u[2] - u[1]), y0 = 2*u[4], y1 = 0, xpd = NA)
		text(u[1] - 0.3*(u[2] - u[1]), u[4], "Free-living parasite stages", xpd = NA, srt = 90, cex = 1.5)
	}
	
} #end i


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



#-----------
# Each arrested scenario in its own plot
period <- array(NA, dim = c(3, 4, 3))

pdf(file = "figures/popCycles_all.pdf", width = 4.5, height = 5, pointsize = 10)
for(i in 1:3){
	for(r in 1:3){
		
		ylimCycles <- array(NA, dim = c(4,2,2))
		
		for(p in 1:4){
			ylimCycles[p, 1, ] <- extendrange(annualSumm[[i, p, r]][40:y, 1], f = 0.15)
			ylimCycles[p, 2, ] <- extendrange(annualSumm[[i, p, r]][40:y, 2], f = 0.15)
		}
		
		y0 <- 40
		
		# quartz(width = 4.5, height = 5, pointsize = 10)
		par(mfrow = c(4,1), mar = c(1,4,1,4), oma = c(3,2,1,2))
		
		for(p in 1:4){
			
			#-------
			# Parasites
			
			y1p <- annualSumm[[i, p, r]][y0:y, 1]
			y2p <- annualSumm[[i, p, r]][y0:y, 2]
			
			plot(y0:y, y1p, "l", lwd = 1.2, xaxt = "n", ylab = "", xlab = "", yaxt = "n", ylim = ylimCycles[p,1,])
			
			if(p > 3) axis(side = 1, at = seq(y0, 100, 20), labels = seq(100 + y0, 200, 20)) else axis(side = 1, at = seq(y0, 100, 20), labels = FALSE)
			
			axis(side = 2, at = pretty(ylimCycles[p,1,]), labels = pretty(ylimCycles[p,1,]*10^-6), las = 1)
		
			
			if(p < 4){
				y1Peaks <- tail(c(2:(length(y1p)-1))[which(y1p[2:(length(y1p)-1)] >  y1p[3:length(y1p)] & y1p[2:(length(y1p)-1)] > y1p[1:(length(y1p)-2)])], 2)
			} else {
				y1Peaks <- head(c(2:(length(y1p)-1))[which(y1p[2:(length(y1p)-1)] >  y1p[3:length(y1p)] & y1p[2:(length(y1p)-1)] > y1p[1:(length(y1p)-2)])], 2)
			}
			segments(x0 = c(y0:y)[y1Peaks], x1 = c(y0:y)[y1Peaks], y0 = y1p[y1Peaks], y1 = max(y1p), lty = 3)
			arrows(x0 = c(y0:y)[y1Peaks[1]], x1 = c(y0:y)[y1Peaks[2]], y0 = max(y1p), y1 = max(y1p), length = 0.06, cod = 3, lwd = 1.2)
			text(mean(c(y0:y)[y1Peaks]), max(annualSumm[[i, p, r]][y0:y, 1]), pos = 1, paste(diff(y1Peaks), "years"))
			
			# Record period for calculating change
			period[i, p, r] <- diff(y1Peaks)
			
			if(p < 4){
				y1Trough <- tail(c(2:(length(y1p)-1))[which(y1p[2:(length(y1p)-1)] <  y1p[3:length(y1p)] & y1p[2:(length(y1p)-1)] < y1p[1:(length(y1p)-2)])], 1)
			} else {
				y1Trough <- c(2:(length(y1p)-1))[which(y1p[2:(length(y1p)-1)] <  y1p[3:length(y1p)] & y1p[2:(length(y1p)-1)] < y1p[1:(length(y1p)-2)])][2]
			}
			segments(x0 = c(y0:y)[y1Trough], x1 = c(y0:y)[y1Trough], y0 = y1p[y1Trough], y1 = min(y1p), lty = 3)
			
			#-------
			# Hosts
			par(new = TRUE)
			
			
			plot(y0:y, y2p, "l", col = "#00000050", yaxt = "n", xaxt = "n", ylim = ylimCycles[p, 2 ,], ylab = "", xlab = "")
			
			axis(side = 4, at = pretty(ylimCycles[p,2,]), labels = pretty(ylimCycles[p,2,]*10^-6), las = 1, col = grey(0.5))
			
			
			y2Trough <- c(2:(length(y2p)-1))[which(y2p[2:(length(y2p)-1)] <  y2p[3:length(y2p)] & y2p[2:(length(y2p)-1)] < y2p[1:(length(y2p)-2)])]
			y2Trough <- y2Trough[which(abs(y2Trough-y1Trough) == min(abs(y2Trough-y1Trough)))]
			
			segments(x0 = c(y0:y)[y2Trough], x1 = c(y0:y)[y2Trough], y0 = y2p[y2Trough], y1 = min(y2p), lty = 3, col = grey(0.5))
			arrows(x0 = c(y0:y)[y1Trough], x1 = c(y0:y)[y2Trough], y0 = min(y2p), y1 = min(y2p), code = 3, length = 0.06, lwd = 1.2)
				
			text(mean(c(y0:y)[c(y1Trough, y2Trough)]), min(annualSumm[[i, p, r]][y0:y, 2]),pos =1, paste(y1Trough[1]- y2Trough[1], " years"))
			
			if(p == 4) legend("bottomright", lwd = c(1.2, 1), legend = c("Parasite pressure", "Total hosts"), col = c(1, grey(0.5)), bty = "n")
			
			mtext(side =3, adj = 0, line = -1.5, paste(" ", LETTERS[p]))
			
			
		} # end p
		
		mtext(side = 2, outer = TRUE, expression(paste("Annual parasite pressure (parasite", {}%*%{}, "days (host)", {}^-1 %*% 10^-6, ")")))
		mtext(side =4, outer = TRUE, expression(paste("Total host population size (", {}%*%10^-6, ")")), col = grey(0.5))
		mtext(side = 1, outer = TRUE, "Year in simulation", line = 2)
		mtext(side = 3, outer = TRUE, paste(c("Low transmission", "Base transmission", "High transmission")[i], c("Current", "RCP 2.6", "RCP 8.5")[r], sep = ": "))

	} # end r
} # end i

dev.off()

###############################################################################
# For climate scenarios, summarize 
# (1) average parasite burdens in last 80 years
# (2) average host population in the last 80 years.
###############################################################################


avg <- array(NA, dim = c(3, 4, 3, 2, 2))
for(i in 1:3){ # For each transmission scenario
	for(r in 1:3){
		for(p in 1:4){
			for(m in 1:2){
				irmp <- annualSumm[[i, p, r]][(100 - period[i, p, r] + 1):100, m]*10^-6
				
				avg[i, p, r, m, 1] <- mean(irmp)
				avg[i, p, r, m, 2] <- sd(irmp)
			}
		}
	}
}

avg[2, ,1, ,]

for(i in 1:3){
	avgTable <- rbind(
		
		scenario1 = c(paste(round(avg[i, 1, c(1:3), 1, 1], 1), " (", round(avg[i, 1, c(1:3), 1, 2], 1), ")", sep = ""), 
									paste(round(avg[i, 1, c(1:3), 2, 1], 1), " (", round(avg[i, 1, c(1:3), 2, 2], 1), ")", sep = "")),
		
		scenario2 = c(paste(round(avg[i, 2, c(1:3), 1, 1], 2), " (", round(avg[i, 2, c(1:3), 1, 2], 2), ")", sep = ""), 
									paste(round(avg[i, 2, c(1:3), 2, 1], 2), " (", round(avg[i, 2, c(1:3), 2, 2], 2), ")", sep = "")),
		
		scenario3 = c(paste(round(avg[i, 3, c(1:3), 1, 1], 2), " (", round(avg[i, 3, c(1:3), 1, 2], 2), ")", sep = ""), 
									paste(round(avg[i, 3, c(1:3), 2, 1], 2), " (", round(avg[i, 3, c(1:3), 2, 2], 2), ")", sep = "")),
		
		scenario4 = c(paste(round(avg[i, 4, c(1:3), 1, 1], 1), " (", round(avg[i, 4, c(1:3), 1, 2], 1), ")", sep = ""), 
									paste(round(avg[i, 4, c(1:3), 2, 1], 2), " (", round(avg[i, 4, c(1:3), 2, 2], 2), ")", sep = ""))
		
	)
	
	write.csv(avgTable, paste0("figures/CCoutput_", c("lowBeta", "baseBeta", "highBeta")[i], ".csv"))
}


# Calculate percent change

percChange <- array(NA, dim = c(3, 4, 2, 2), dimnames = list(c("lowBeta", "baseBeta", "highBeta"), paste0("Scenario", c(1:4)), c("RCP2.6", "RCP8.5"), c("Parasite" ,"Host")))
for(i in 1:3){
	for(p in 1:4){
		for(r in 1:2){
			for(m in 1:2){
				percChange[i, p, r, m] <- (avg[i, p, r + 1, m, 1] - avg[i, p, 1, m, 1])/avg[i, p, 1, m, 1]*100
			}}}}

percChange[2, , , 1]

# Plot as a percentage of the mean in current scenarios
i <- 2
quartz(width = 6.3, height = 2.8, pointsize = 10)
par(mfrow = c(1,2), mar = c(4, 2, 2, 1), oma= c(0, 2, 0, 5))

# Parasite
plot(1:4, percChange[i, , 1, 1], "n", ylim = extendrange(percChange[i, , , 1], f = 0.2), las = 1, xlim = c(0.5, 4.5), ylab = "", xlab = "Scenario")
mtext(side = 2, "Percent change", line = 3)
abline(h = 0)

for(r in 1:2){
	points(1:4, percChange[i, , r, 1], col = c(4,2)[r], pch = 21, bg = c("#2297E620", "#DF536B80")[r])
}
mtext(side = 3, adj = 0, line = 0.5, "A) Annual parasite pressure")

# Host
plot(1:4, percChange[i, , 1, 2], "n", ylim = extendrange(percChange[i, , , 2], f = 0.2), las = 1, xlim = c(0.5, 4.5), ylab = "", xlab = "Scenario")
abline(h = 0)
for(r in 1:2){
	points(1:4, percChange[i, , r, 2], col = c(4,2)[r], pch = 21, bg = c("#2297E620", "#DF536B80")[r])
}
mtext(side = 3, adj = 0, line = 0.5, "B) Total host pop'n")

legend(4.8, 5.5, pch = 21, col = c(4,2), pt.bg = c("#2297E620", "#DF536B80"), c("RCP 2.6", "RCP 8.5"), xpd = NA, bty = "n")
mtext(side = 3, line = 1, outer=TRUE, c("Low transmission", "Base transmission", "High transmission")[i])

#################################################################################
# Uptake rate over space and time in year 100
#################################################################################
V <- readRDS(paste0("simulations/output/V_baseBeta_migPath", N, ".rds"))

x.ind <- seq(300, 600, 1)
t <- seq(120, 300, 1)
t.ind <- which(timeDat$DOY[which(timeDat$year > 80)] %in% t & timeDat$year[which(timeDat$year > 80)] == 100)

colUptake <- colorRampPalette(c("white", 2, 1))(n = 100)

V.pr <- array(NA, dim = c(4, 3, length(x.ind), length(t)))
for(p in 1:4){
	for(r in 1:3){
		V.pr[p, r, , ] <- V[[p, r]][15, x.ind, t.ind]
	}
}

tDate <- as.Date(paste("2001", c(6:10), "01", sep = "-"), format = "%Y-%m-%d")
tDOY <- as.numeric(strftime(tDate, format = "%j"))

yTicks <- list(c(0, 75e4, 150e4, 225e4, 300e4),
								c(0, 75, 150, 300),
								c(0, 35, 70, 105, 140),
								c(0, 7500, 15000, 22500))
	
# quartz(width = 6.3, height = 5.5, pointsize = 10)
png(file = "figures/uptakeRate_baseBeta_migPath6.png", width = 6.3*150, height = 5.5*150, pointsize = 10, res = 150)
par(mfrow = c(4,3), mar = c(1, 0, 0, 0), oma = c(3,5,4,7))

for(p in 1:4){
	uptakeBreaks <- seq(0, max(V.pr[p, , ,])*0.5, length.out = 100)
	# uptakeBreaks <- seq(0, c(100000, 300, 150, 25000)[p], length.out = 100)
	
	for(r in 1:3){
		
		plot(rep(x.ind, length(t)), rep(t, each = length(x.ind)), col = colUptake[findInterval(V.pr[p, r, , ], uptakeBreaks)], pch = 15, xaxs = "i", yaxs = "i", xlab = "", ylab = "", yaxt = "n", xaxt = "n", cex = 0.2)
		
		# abline(h = breedDOY, lty = 3)
		abline(h = tDOY, col = "#00000040", lwd = 0.8)
		abline(v = seq(650, 1200, 50)/2, col = "#00000040", lwd = 0.8)
		
		if(r == 1){
			axis(side = 2, at = tDOY, labels = month.abb[c(6:10)], las= 1)
			mtext(side = 3, adj = 0, line = -1.5, paste0("  ", LETTERS[p]))
		}
		if(p == 4) axis(side = 1, at = seq(350, 550, 100), labels = seq(350, 550, 100)*2)
		
		if(p == 1){
			mtext(side = 3, c("Current", "RCP 2.6", "RCP 8.5")[r], line = 1)
		}
	} # end r
	
	for(i in seq(1, 100, 2)){
		points(c(631:635), rep(seq(120, 290, length.out = 100)[i], 5), pch = 15, cex = 0.5, col = colUptake[i], xpd = NA)
	}
	
	for(j in 1:length(yTicks[[p]])){
		segments(x0 = 631, x1 = 640, y0 = seq(120, 290, length.out = 100)[findInterval(yTicks[[p]][j], uptakeBreaks)] , y1 = seq(120, 290, length.out = 100)[findInterval(yTicks[[p]][j], uptakeBreaks)], xpd = NA)
		text(640, seq(120, 290, length.out = 100)[findInterval(yTicks[[p]][j], uptakeBreaks)], pos = 4, yTicks[[p]][j], xpd = NA)
	}
	 
		
		
	}

mtext(side = 1, "Distance along migration (km)", outer = TRUE, line = 2)
# mtext(side = 3, outer = TRUE, "Parasite uptake rate per host per day", line = 0.5)
text(640, 935, "Parasite\nuptake", xpd = NA)

dev.off()

#------------------------------------------------------------------------------
# Difference in uptake rate?
#------------------------------------------------------------------------------

Vdiff <- array(NA, dim = c(4, 2, length(x.ind), length(t)))
for(p in 1:4){
	for(r in 1:2){
		Vdiff[p, r, , ] <- (V[[p, r + 1]][15, x.ind, t.ind] - V[[p, 1]][15, x.ind, t.ind])
	}
}

yTicks <- list(c(-2e6, -1e6, 0, 1e6, 2e6),
							 c(0, 100, 200, 300),
							 c(-10, 0, 10, 20),
							 c(-4000, 0, 4000, 8000, 12000, 16000))

#quartz(width = 4.5, height = 5.5, pointsize = 10)
png(file = "figures/uptakeRateChange_baseBeta_migPath6.png", width = 4.5*150, height = 5.5*150, pointsize = 10, res = 150)

par(mfrow = c(4,2), mar = c(1, 0, 0, 0), oma = c(3,5,4,7))

for(p in 1:4){
	uptakeBreaks <- seq(min(Vdiff[p, r, ,]), max(Vdiff[p, r, ,]), length.out = 100)
	colUptake <- c(colorRampPalette(c(4,"white"))(n = sum(uptakeBreaks < 0)), colorRampPalette(c("white", 2))(n = sum(uptakeBreaks > 0)))
	for(r in 1:2){
		
		plot(rep(x.ind, length(t)), rep(t, each = length(x.ind)), col = colUptake[findInterval(Vdiff[p, r, , ], uptakeBreaks)], pch = 15, xaxs = "i", yaxs = "i", xlab = "", ylab = "", yaxt = "n", xaxt = "n", cex = 0.3)
		
		# abline(h = breedDOY, lty = 3)
		abline(h = tDOY, col = "#00000040", lwd = 0.8)
		abline(v = seq(650, 1200, 50)/2, col = "#00000040", lwd = 0.8)
		
		if(r == 1){
			axis(side = 2, at = tDOY, labels = month.abb[c(6:10)], las= 1)
			# mtext(side = 3, adj = 0, line = -1.5, paste0("  ", LETTERS[p]))
		}
		if(p == 4) axis(side = 1, at = seq(350, 550, 100), labels = seq(350, 550, 100)*2)
		
		if(p == 1){
			mtext(side = 3, c("RCP 2.6", "RCP 8.5")[r], line = 1)
		}
	} # end r
	
	for(i in seq(1, 100, 2)){
		points(c(631:635), rep(seq(120, 290, length.out = 100)[i], 5), pch = 15, cex = 0.5, col = colUptake[i], xpd = NA)
	}
	
	for(j in 1:length(yTicks[[p]])){
		segments(x0 = 631, x1 = 640, y0 = seq(120, 290, length.out = 100)[findInterval(yTicks[[p]][j], uptakeBreaks)] , y1 = seq(120, 290, length.out = 100)[findInterval(yTicks[[p]][j], uptakeBreaks)], xpd = NA)
		text(640, seq(120, 290, length.out = 100)[findInterval(yTicks[[p]][j], uptakeBreaks)], pos = 4, yTicks[[p]][j], xpd = NA)
	}
	
	
	
}

mtext(side = 1, "Distance along migration (km)", outer = TRUE, line = 2)
# mtext(side = 3, outer = TRUE, "Parasite uptake rate per host per day", line = 0.5)
text(640, 935, "Change in\nuptake", xpd = NA)

dev.off()

#------------------------------------------------------------------------------
# Current uptake only
#------------------------------------------------------------------------------

yTicks <- list(c(0, 75e4, 150e4, 225e4, 300e4),
							 c(0, 75, 150, 300),
							 c(0, 35, 70, 105, 140),
							 c(0, 7500, 15000, 22500))

# quartz(width = 3, height = 5.5, pointsize = 10)
png(file = "figures/uptakeRateCurrent_baseBeta_migPath6.png", width = 3*150, height = 5.5*150, pointsize = 10, res = 150)
par(mfrow = c(4,1), mar = c(1, 0, 0, 0), oma = c(3,5,4,7))
colUptake <- colorRampPalette(c("white", 1))(n = 100)

for(p in 1:4){
	uptakeBreaks <- seq(0, max(V.pr[p, , ,])*0.5, length.out = 100)
	r <- 1
		
		plot(rep(x.ind, length(t)), rep(t, each = length(x.ind)), col = colUptake[findInterval(V.pr[p, r, , ], uptakeBreaks)], pch = 15, xaxs = "i", yaxs = "i", xlab = "", ylab = "", yaxt = "n", xaxt = "n", cex = 0.3)
		
		abline(h = tDOY, col = "#00000040", lwd = 0.8)
		abline(v = seq(650, 1200, 50)/2, col = "#00000040", lwd = 0.8)
		
		axis(side = 2, at = tDOY, labels = month.abb[c(6:10)], las= 1)
		mtext(side = 3, adj = 0, line = -1.5, paste0("  ", LETTERS[p]))
		
		if(p == 4) axis(side = 1, at = seq(350, 550, 100), labels = seq(350, 550, 100)*2)
		
		if(p == 1){
			mtext(side = 3, c("Current", "RCP 2.6", "RCP 8.5")[r], line = 1)
		}
	
	for(i in seq(1, 100, 2)){
		points(c(631:635), rep(seq(120, 290, length.out = 100)[i], 5), pch = 15, cex = 0.5, col = colUptake[i], xpd = NA)
	}
	
	for(j in 1:length(yTicks[[p]])){
		segments(x0 = 631, x1 = 640, y0 = seq(120, 290, length.out = 100)[findInterval(yTicks[[p]][j], uptakeBreaks)] , y1 = seq(120, 290, length.out = 100)[findInterval(yTicks[[p]][j], uptakeBreaks)], xpd = NA)
		text(640, seq(120, 290, length.out = 100)[findInterval(yTicks[[p]][j], uptakeBreaks)], pos = 4, yTicks[[p]][j], xpd = NA)
	}
	
	
	
}

mtext(side = 1, "Distance along migration (km)", outer = TRUE, line = 2)
# mtext(side = 3, outer = TRUE, "Parasite uptake rate per host per day", line = 0.5)
text(640, 935, "Parasite\nuptake", xpd = NA)

dev.off()

#------------------------------------------------------------------------------
# Within-host adult parasites
#------------------------------------------------------------------------------
colUptake <- colorRampPalette(c("white", 1))(n = 100)

V.P <- array(NA, dim = c(4, 3, length(x.ind), length(t)))
for(p in 1:4){
	for(r in 1:3){
		V.P[p, r, , ] <- V[[p, r]]["P_mov", x.ind, t.ind] + V[[p, r]]["P_stat", x.ind, t.ind]
	}
}

yTicks <- list(c(0, 75e4, 150e4, 225e4, 300e4),
							 c(0, 75, 150, 300),
							 c(0, 35, 70, 105, 140),
							 c(0, 7500, 15000, 22500))

# quartz(width = 6.3, height = 5.5, pointsize = 10)
png(file = "figures/uptakeRate_baseBeta_migPath6.png", width = 6.3*150, height = 5.5*150, pointsize = 10, res = 150)
par(mfrow = c(4,3), mar = c(1, 0, 0, 0), oma = c(3,5,4,7))

for(p in 1:4){
	uptakeBreaks <- seq(0, max(V.P[p, , ,])*0.5, length.out = 100)
	# uptakeBreaks <- seq(0, c(100000, 300, 150, 25000)[p], length.out = 100)
	
	for(r in 1:3){
		
		plot(rep(x.ind, length(t)), rep(t, each = length(x.ind)), col = colUptake[findInterval(V.P[p, r, , ], uptakeBreaks)], pch = 15, xaxs = "i", yaxs = "i", xlab = "", ylab = "", yaxt = "n", xaxt = "n", cex = 0.2)
		
		# abline(h = breedDOY, lty = 3)
		abline(h = tDOY, col = "#00000040", lwd = 0.8)
		abline(v = seq(650, 1200, 50)/2, col = "#00000040", lwd = 0.8)
		
		if(r == 1){
			axis(side = 2, at = tDOY, labels = month.abb[c(6:10)], las= 1)
			mtext(side = 3, adj = 0, line = -1.5, paste0("  ", LETTERS[p]))
		}
		if(p == 4) axis(side = 1, at = seq(350, 550, 100), labels = seq(350, 550, 100)*2)
		
		if(p == 1){
			mtext(side = 3, c("Current", "RCP 2.6", "RCP 8.5")[r], line = 1)
		}
		
		abline(h = breedDOY - 240 + 360, col = 2, lty = 3)
	} # end r
	
	for(i in seq(1, 100, 2)){
		points(c(631:635), rep(seq(120, 290, length.out = 100)[i], 5), pch = 15, cex = 0.5, col = colUptake[i], xpd = NA)
	}
	
	# for(j in 1:length(yTicks[[p]])){
	# 	segments(x0 = 631, x1 = 640, y0 = seq(120, 290, length.out = 100)[findInterval(yTicks[[p]][j], uptakeBreaks)] , y1 = seq(120, 290, length.out = 100)[findInterval(yTicks[[p]][j], uptakeBreaks)], xpd = NA)
	# 	text(640, seq(120, 290, length.out = 100)[findInterval(yTicks[[p]][j], uptakeBreaks)], pos = 4, yTicks[[p]][j], xpd = NA)
	# }
	
	
	
}

mtext(side = 1, "Distance along migration (km)", outer = TRUE, line = 2)
# mtext(side = 3, outer = TRUE, "Parasite uptake rate per host per day", line = 0.5)
text(640, 935, "Parasite\nuptake", xpd = NA)

dev.off()

#################################################################################
# Chnage in pre-infective survival
#################################################################################
allTemps <- readRDS("tempData/migRouteOutput/allTemps_migPath6.rds")
t <- 1:365
xSelected <- 1000

params <- array(NA, dim = c(4,3, 365))
for(r in 1:3){
	params[1, r, ] <- predict.mu0(allTemps[[r]][t, xSelected])
	params[2, r, ] <- predict.rho0(allTemps[[r]][t, xSelected])
	params[3, r, ] <- exp(- predict.mu0(allTemps[[r]][t, xSelected])*(1/predict.rho0(allTemps[[r]][t, xSelected])))
	params[4, r, ] <- predict.mu3(allTemps[[r]][t, xSelected])
	
}

quartz(width = 6.3, height = 5, pointsize = 10)
par(mfrow = c(2, 2), mar = c(3,4,2,1), oma = c(0, 2, 0, 0))

for(j in 1:4){
	
	plot(xDate, params[j, 1, ], "n", ylim = range(params[j, , ]), xlab = "", ylab = "", las = 1)
	for(r in 1:3) lines(xDate, params[j, r, ], col = c(1,4,2)[r])

	if(j == 1) legend("topright", lwd = 1, col = c(1,4,2), c("Current", "RCP 2.6", "RCP 8.5"), bty = "n")
	mtext(side = 3, line = 0.5, adj = 0, paste0(LETTERS[j], ") ", c("Pre-infective mortality", "Development", "Proportion surviving to infective", "Infective mortality")[j]))
	}
mtext(side= 2, outer = TRUE, "Parameter value")


# Plot burdens in fall breed
p <- 1
Pfall <- rbind(V[[p,1]]["P_mov", , max(which(timeDat$DOY[timeDat$year > 80] == (breedDOY - 240 + 365)))] + V[[p, 1]]["P_stat", , max(which(timeDat$DOY[timeDat$year > 80] == (breedDOY - 240 + 365)))],
							 V[[p, 2]]["P_mov", , max(which(timeDat$DOY[timeDat$year > 80] == (breedDOY - 240 + 365)))] + V[[p, 2]]["P_stat", , max(which(timeDat$DOY[timeDat$year > 80] == (breedDOY - 240 + 365)))],
							 V[[p, 3]]["P_mov", , max(which(timeDat$DOY[timeDat$year > 80] == (breedDOY - 240 + 365)))] + V[[p, 3]]["P_stat", , max(which(timeDat$DOY[timeDat$year > 80] == (breedDOY - 240 + 365)))])

par(mfrow = c(1,1), mar = c(4,4,2,1), om a= rep(0,4))
plot(x, Pfall[1, ], "l", ylim = range(Pfall))
lines(x, Pfall[2, ],col = 4)
lines(x, Pfall[3, ],col = 2)

# Look during peak

i <- 2
p <- 4
y0 <- 40

affectingCalving <- list(
	Pfall = array(NA, dim = c(3, 1135)),
	Hcalv = array(NA, dim = c(3, 1135)),
	calf= array(NA, dim = c(3, 1135))
)

for(r in 1:3){
	# Peak host population
	y2p <- annualSumm[[i, p, r]][, 2]
	yearPeakH <- tail(c(2:(length(y2p)-1))[which(y2p[2:(length(y2p)-1)] >  y2p[3:length(y2p)] & y2p[2:(length(y2p)-1)] > y2p[1:(length(y2p)-2)])], 1) 
	timeIndex <- which(timeDat$year[timeDat$year > 80] == yearPeakH & timeDat$DOY[timeDat$year > 80] == (breedDOY - 240 + 365))
	
	affectingCalving$Pfall[r, ] <- V[[p,r]]["P_mov", , timeIndex] + V[[p, r]]["P_stat", , timeIndex]
	
	timeIndex <- which(timeDat$year[timeDat$year > 80] == yearPeakH & timeDat$DOY[timeDat$year > 80] == breedDOY)
	affectingCalving$Hcalv[r, ] <- V[[p,r]]["adult_mov", , timeIndex] + V[[p, r]]["adult_stat", , timeIndex]
	affectingCalving$calf[r, ] <- V[[p,r]]["calf_mov", , timeIndex+2] + V[[p, r]]["calf_stat", , timeIndex+2]
}

par(mfrow = c(3,1))
for(m in 1:3){
	plot(x, affectingCalving[[m]][1, ], "n", ylim = range(affectingCalving[[m]]))
	for(r in 1:3){
		lines(x, affectingCalving[[m]][r, ], col = c(1,4,2)[r])
	}
	mtext(side = 3, c("Parasites at breeding", "Hosts at calving", "New calves 2 days after calving")[m])
}

plot(x, affectingCalving[[m]][1, ], "n", ylim = range(affectingCalving[[m]]), xlim = c(400, 800))
for(r in 1:3){
	lines(x, affectingCalving[[3]][r, ], col = c(1,4,2)[r])
}

plot(y0:y, annualSumm[[2, 1, 1]][y0:y, 2], "n")
for(r in 1:3){
	lines(y0:y, annualSumm[[2, 1, r]][y0:y, 2], col = c(1,4,2)[r])
}

# Summarize all space through time
# Peak host population year
Ppop <- array(NA, dim = c(3, 365))
for(r in 1:3){
	# Peak host population
	y2p <- annualSumm[[i, p, r]][, 2]
	yearPeakH <- tail(c(2:(length(y2p)-1))[which(y2p[2:(length(y2p)-1)] >  y2p[3:length(y2p)] & y2p[2:(length(y2p)-1)] > y2p[1:(length(y2p)-2)])], 1) 
	timeIndex <- which(timeDat$year[timeDat$year > 80] == yearPeakH)
	
	Ppop[r, ] <- apply(V[[p,r]]["P_mov", , timeIndex] + V[[p, r]]["P_stat", , timeIndex], 2, sum)
	}

plot(1:365, Ppop[1,], "n", xlim = c(100, 200), ylim = range(Ppop))
for(r in 1:3){
	lines(1:365, Ppop[r, ], col = c(1,4,2)[r])
}
