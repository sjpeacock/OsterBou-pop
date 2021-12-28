# What is the transmission rate


i <- 2 # Transmissionm
p <- 3 # Arrested development scenario
r <- 1 # Current climate

if(i == 2) beta <- 10^-6 else if(i == 3) beta <- 10^-5 else if(i == 1) beta <- 10^-7


# Extract last year of simulation
V.ipr <- V[[i, p, r]][, , which(timeDat$year[timeDat$year > 80] == 100)]

dim(V.ipr)

#------------------------------------------------------------------------------
# Calculate parasite uptake for each DOY
#------------------------------------------------------------------------------
uptakePerHost <- array(NA, dim = c(4, 3, 365), dimnames = list(c("100Arrest", "50Arrest", "0Arrest", "variable"), c("current", "RCP2.6", "RCP8.5"), NULL))

for(p in 1:4){
	for(r in 1:3){
		
		V.ipr <- V[[i, p, r]][, , which(timeDat$year[timeDat$year > 80] == 100)]
	
		for(j in 1:365){
			uptakePerHost[p, r, j] <- sum(beta * V.ipr["L3", ,j] * (V.ipr["adult_mov", ,j] + V.ipr["adult_stat", ,j]))/sum((V.ipr["adult_mov", ,j] + V.ipr["adult_stat", ,j])) * 10^-6
		}
	}}

xDate <- as.Date(paste(2000, c(1:12), 01, sep = "-"))
par(mfrow = c(4,1), mar = rep(0, 4), oma = c(5, 5, 2, 1))
for(p in 1:4){
	plot(DOY, uptakePerHost[p, 1, ], "n", ylim = range(uptakePerHost[p, , seq(1, 365, 2)]), xaxt = "n", xaxs = "i", las = 1, ylab = "")
	if(p == 4) axis(side = 1, at = as.numeric(strftime(xDate, format = "%j")), labels = month.abb)
	for(r in 1:3) lines(DOY[seq(2, 365, 2)], uptakePerHost[p, r, seq(2, 365, 2)],col = c(1, 4, 2)[r])
	mtext(side = 3, adj = 0, line = -2, paste0("   ", letters[p], ") ", c("100% arresting", "50% arresting", "0% arresting", "variable arresting")[p]))
	abline(v = breedDOY+360-240, lty = 3)
}

#------------------------------------------------------------------------------
# Adult parasite burdens per host
#------------------------------------------------------------------------------

pBar <- array(NA, dim = c(4, 3, 365), dimnames = list(c("100Arrest", "50Arrest", "0Arrest", "variable"), c("current", "RCP2.6", "RCP8.5"), NULL))

for(p in 1:4){
	for(r in 1:3){
		
		V.ipr <- V[[i, p, r]][, , which(timeDat$year[timeDat$year > 80] == 100)]
		
		for(j in 1:365){
			mov <- which(V.ipr["adult_mov", ,j] > 10^-10)
			stat <- which(V.ipr["adult_stat", ,j] > 10^-10)
			
			pBar[p, r, j] <- sum(V.ipr["P_mov", mov,j] / (V.ipr["adult_mov", mov,j]) + V.ipr["P_stat", stat,j] / (V.ipr["adult_stat", stat, j])) * 10^-3
		}
	}}

xDate <- as.Date(paste(2000, c(1:12), 01, sep = "-"))
par(mfrow = c(4,1), mar = rep(0, 4), oma = c(5, 5, 2, 1))
for(p in 1:4){
	plot(DOY, pBar[p, 1, ], "n", ylim = range(pBar[p, , ]), xaxt = "n", xaxs = "i", las = 1, ylab = "")
	if(p == 4) axis(side = 1, at = as.numeric(strftime(xDate, format = "%j")), labels = month.abb)
	for(r in 1:3) lines(DOY, pBar[p, r,],col = c(1, 4, 2)[r])
	mtext(side = 3, adj = 0, line = -2, paste0("   ", letters[p], ") ", c("100% arresting", "50% arresting", "0% arresting", "variable arresting")[p]))
	abline(v = breedDOY+360-240, lty = 3)
}

# What the hell is going on?
colBay <- PNWColors::pnw_palette("Bay", n = 100)

# Determine yMax
yMax <- matrix(NA, nrow = 3, ncol = 4)
for(s in 1:3){
	for(p in 1:4){
		for(r in 1:3){
			V.ipr <- V[[i, p, r]][, , which(timeDat$year[timeDat$year > 80] == 100)]
			
			if(s == 1){ #Hosts
				Y <- V.ipr["adult_stat", seq(1, 1135, 10), seq(1, 365, 3)] + V.ipr["adult_mov", seq(1, 1135, 10), seq(1, 365, 3)]
				
			} else if(s == 2){	# All parasites
				Y <- V.ipr["P_stat", seq(1, 1135, 10), seq(1, 365, 3)] + V.ipr["P_mov", seq(1, 1135, 10), seq(1, 365, 3)]
		} else if(s == 3){ # L3 larvae
			Y <- V.ipr["L3", seq(1, 1135, 10), seq(1, 365, 3)] 
		}
		yMax[s, p] <- max(0, max(Y), yMax[s, p], na.rm = TRUE)
		}}}

yMax[2,1] <- 1e11
yMax[2,4] <- 1e11
# quartz(width = 8, height = 8, pointsize = 10)
pdf("figures/gridComparison_migPath6.pdf", width = 8, height = 12, pointsize = 12)
for(s in 1:3){
	par(mfrow = c(4, 3), mar = c(0.2,0.2,0.2,0.2), oma = c(5, 6, 4, 9))
	for(p in 1:4){
		for(r in 1:3){
			V.ipr <- V[[i, p, r]][, , which(timeDat$year[timeDat$year > 80] == 100)]
			
			
			if(s == 1){ #Hosts
				Y <- V.ipr["adult_stat", seq(1, 1135, 10), seq(1, 365, 3)] + V.ipr["adult_mov", seq(1, 1135, 10), seq(1, 365, 3)]
			} else if(s == 2){ # All parasites
				Y <- V.ipr["P_stat", seq(1, 1135, 10), seq(1, 365, 3)] + V.ipr["P_mov", seq(1, 1135, 10), seq(1, 365, 3)]
		} else if(s == 3){
			Y <- V.ipr["L3", seq(1, 1135, 10), seq(1, 365, 3)] 
		}
			yMax.ipr <- max(Y)
			# yMax.ipr <- yMax[s]
			yMax.ipr <- yMax[s, p]
		# par(mar = c(5,4,3,10))
		plot(rep(c(1:1135)[seq(1, 1135, 10)], 122), rep(c(1:365)[seq(1, 365, 3)], each = 114), pch = 15, col = colBay[findInterval(Y, seq(0, yMax.ipr, length.out = 100))], cex = 0.5, xlab = "", xaxt = "n", las = 1, ylab = "", yaxt = "n", xaxs = "i", yaxs = "i")
		if(r == 1) axis(side = 2, at = as.numeric(strftime(xDate, format = "%j")), labels = month.abb, las = 1)
		if(p == 4) axis(side = 1, at = seq(0, 1000, 200), labels = seq(0, 1000, 200)*2)
		abline(v = seq(0, 1135, 100))
		abline(h = c(110, 153, 168, 180, 250, 290, 305, 335))
		if(r == 3) axis(side = 4, at = c(110, 153, 168, 180, 250, 290, 305, 335), label = c("Spring migration", "Calving", "Post-calving migration", "Summer", "Fall migration", "Rut/breeding", "Fall migration", "Winter"), las = 1, cex = 0.8)
		# mtext(side = 3, adj = 0, line = -1.5, paste0(c(" 100% arresting", " 50% arresting", " 0% arresting", " variable arresting")[p], ": ", c("current", "rcp2.6", "rcp8.5")[r]), col = "white", font = 2)
		if(p == 1) mtext(side = 3, line = 2, c("current", " RCP 2.6", "RCP 8.5")[r])
		if(r == 1) mtext(side = 2, line = 3, c("100% arresting", "50% arresting", "0% arresting", "variable arresting")[p])
	}}
mtext(side = 1, outer = TRUE, "Distance along migration (km)", line = 3)
} # end s
dev.off()


#------------
# Look at adult parasite burdens at time of rut
t <- breedDOY + 360 - 240

rut <- array(NA, dim = c(2, 4, 3, 1135))
for(p in 1:4){
	for(r in 1:3){
		V.ipr <- V[[i, p, r]][, , which(timeDat$year[timeDat$year > 80] == 100)]
		rut[1, p, r, ] <- V.ipr["P_stat", , t] + V.ipr["P_mov", s, t]
		rut[2, p, r, ] <- V.ipr["adult_stat", , t] + V.ipr["adult_mov", s, t]
}}

for(p in 1:4){
plot(x, rut[1, p, 1, ], "n", ylim = range( rut[1, p, , ]), xlim = c(1250, 2000), main = c("100% arresting", "50% arresting", "0% arresting", "variable arresting")[p])
for(r in 1:3) lines(x, rut[1, p, r, ], col = c(1,4,2)[r])
}

#------------
# Mean parasite burdens iat time of breeding
#------------
t <- breedDOY + 360 - 240

Pbreed <- array(NA, dim = c(3, 4, 3, 20))
for(i in 1:3){
	for(p in 1:4){
		for(r in 1:3){
			for(y in 1:20){
				V.ipry <- V[[i, p, r]][, , which(timeDat$year[timeDat$year > 80] == 80 + y)]
				
				Pbreed[i, p, r, y] <- sum(c(V.ipry['P_stat', , t], V.ipry['P_mov', , t]))/sum(c(V.ipry['adult_stat', , t],  V.ipry['adult_mov', , t]))
			}
		}
	}
}

par(mfrow = c(2,2), mar = c(3,3,4,1), oma = c(2,2,0,0))
i <- 2
for(p in 1:4){
	plot(c(181:200), Pbreed[i, p, 1, ], "n", ylim = range(Pbreed[i, p, , ]), ylab = "", xlab = "")
	mtext(side = 3, adj = 0, line = 2, c("100% arresting", "50% arresting", "0% arresting", "variable arresting")[p])
	for(r in 1:3){
		lines(c(181:200), Pbreed[i, p, r, ], "o", col = c(1,4,2)[r])
	}
	par(new = TRUE)
	Pp <- seq(min(Pbreed[i, p, , ]), max(Pbreed[i, p, , ]), length.out = 200)
	plot(numCalves(P_mean = Pp, numFemales = 1, pCalf0 = 0.8), Pp, "l", lwd = 3, col = "#00000040", yaxt = "n", xaxt = "n", bty = "n")
	axis(side = 3, col = grey(0.8))
}
mtext(side = 3, line = 0, outer = TRUE, "Mean parasite burdens 8 month prior to calving")
mtext(side = 1, outer= TRUE, "Year")


#------------
# What is the timing of peak parasite burdens

peakDOY <- array(NA, dim = c(3, 4, 3, 20))
for(i in 1:3){
	for(p in 1:4){
		for(r in 1:3){
			for(y in 1:20){
				V.ipry <- V[[i, p, r]][, , which(timeDat$year[timeDat$year > 80] == 80 + y)]
				Pj <- numeric(365)
				for(j in 1:365){
					Pj[j] <- sum(c(V.ipry['P_stat', , j], V.ipry['P_mov', , j]))/sum(c(V.ipry['adult_stat', , j],  V.ipry['adult_mov', , j]))
				}
				peakDOY[i, p, r, y] <- which(Pj == max(Pj))
			}
		}
	}
}

par(mfrow = c(2,2), mar = c(3,3,4,1), oma = c(2,2,0,0))
i <- 2
plot(c(181:200), peakDOY[i, p, 1, ], "n", ylim = c(1, 365), ylab = "", xlab = "")
mtext(side = 3, adj = 0, line = 2, c("100% arresting", "50% arresting", "0% arresting", "variable arresting")[p])
	for(p in 1:4){for(r in 1:3){
		lines(c(181:200), peakDOY[i, p, r, ], "o", col = c(1,4,2)[r])
	}
}
mtext(side = 3, line = 0, outer = TRUE, "Mean parasite burdens 8 month prior to calving")
mtext(side = 1, outer= TRUE, "Year")

#------------
# Mean parasite burdens at calving
#------------
t <- breedDOY
t <- 215 # mid summer

Pcalv <- array(NA, dim = c(3, 4, 3, 20))
for(i in 1:3){
	for(p in 1:4){
		for(r in 1:3){
			for(y in 1:20){
				V.ipry <- V[[i, p, r]][, , which(timeDat$year[timeDat$year > 80] == 80 + y)]
				
				Pcalv[i, p, r, y] <- sum(c(V.ipry['P_stat', , t], V.ipry['P_mov', , t]))/sum(c(V.ipry['adult_stat', , t],  V.ipry['adult_mov', , t]))
			}
		}
	}
}

par(mfrow = c(2,2), mar = c(3,3,4,1), oma = c(2,2,0,0))
i <- 2
for(p in 1:4){
	plot(c(181:200), Pcalv[i, p, 1, ], "n", ylim = range(Pcalv[i, p, , ]), ylab = "", xlab = "")
	mtext(side = 3, adj = 0, line = 2, c("100% arresting", "50% arresting", "0% arresting", "variable arresting")[p])
	for(r in 1:3){
		lines(c(181:200), Pcalv[i, p, r, ], "o", col = c(1,4,2)[r])
	}
}
mtext(side = 3, line = 0, outer = TRUE, "Mean parasite burdens mid-summer")
mtext(side = 1, outer= TRUE, "Year")


#------------------
# Rather than taking mean annual, look at the 95% quanitles of daily parasite burdens

boot <- function(x){
	bootX <- apply(matrix(sample(x, length(x) * 1000, replace = TRUE), nrow = 1000), 1, mean)
	return(c(quantile(bootX, 0.025), mean(bootX), quantile(bootX, 0.975)))
}

avg80 <- array(NA, dim = c(3, 4, 3, 2, 3))
i <- 2
	for(r in 1:3){
		for(p in 1:4){
			V.ipr <- V[[i, p, r]]
			Pj <- numeric(7300)
			Hj <- numeric(7300)
			for(j in 1:7300){
				Pj[j] <- sum(c(V.ipr['P_stat', , j], V.ipr['P_mov', , j]))/sum(c(V.ipr['adult_stat', , j],  V.ipr['adult_mov', , j]))
				Hj[j] <- sum(c(V.ipr['adult_stat', , j],  V.ipr['adult_mov', , j]))
				
			}
			avg80[i, p, r, 1, ] <- boot(Pj)
			
		}
	}

quartz(width = 4.5, height = 4, pointsize = 10)
par(mfrow = c(2,1), mar = c(2,5,2,1), oma = c(2, 0, 0, 6))
for(m in 1:2){
	plot(1:3, avg80[i, 1, , m, 2], "n", ylim = range(avg80[i, 2:4, , m, ]), xlim = c(0.5, 3.5), xaxs = "i", xaxt = "n", las = 1, xlab = "", ylab = "", yaxt = "n")
	axis(side = 1, at = c(0.5, 1.5, 2.5, 3.5), labels = FALSE)
	axis(side = 1, at = c(1, 2, 3), labels = c("Current", "RCP 2.6", "RCP 8.5"), tck = 0)
	abline(v = c(1.5, 2.5))
	for(p in 2:4){
		segments(x0 = 1:3 + c(-0.3, -0.1, 0.1, 0.3)[p], 
						 x1 = 1:3 + c(-0.3, -0.1, 0.1, 0.3)[p], 
						 y0 = avg80[i, p, , m, 1], 
						 y1 = avg80[i, p, , m, 3], 
						 col = cols2[p], lwd = 2)
		points(1:3 + c(-0.3, -0.1, 0.1, 0.3)[p], avg80[i, p, , m, 2], pch = c(21,22,24,25)[p], col = cols[p], bg = cols2[p])
	}
	mtext(side = 3, adj = 0, line = 0, LETTERS[m])
	mtext(side = 2, c("Change in parasite burdens", "Change in host pop'n")[m], line = 3)
}
mtext(side = 1, "Climate scenario", line = 3)
legend(3.8, 42, pch = c(21,22,24,25), col = cols[1:4], pt.bg = cols2[1:4], title = "Arrested\ndevelopment", legend = c("100%", "50%", "0%", "variable"), bty = "n", cex = 0.8, pt.cex = 1.0, xpd= NA)


