###############################################################################
# Host stopping and starting
###############################################################################

source("simulations/bouSetup.R")
source("simulations/popFunctions.R")

y <- 1
dt <- 1
startstop <- "agg"

source("simulations/calcGrid.R")

breaks <- seq(0, 0.99, 0.01)
length(breaks)
startCol <- colorRampPalette(colors = c("white", 3, 1))(n = 100)
stopCol <- colorRampPalette(colors = c("white", 2, 1))(n = 100)

axisDates <- as.numeric(strftime(as.Date(paste(2001, c(1:12), 01, sep = "-"), format = "%Y-%m-%d"), format = "%j"))

# quartz(width = 6.3, height = 3.5, pointsize = 10)
png(filename = "figures/startStopMat.png", width = 975, height =  500, pointsize = 14)
par(mfrow = c(1,2), mar = c(1,1,1,0), oma = c(4,5,2,8))

# Starting rates
plot(rep(seq(1, 1135*2, 2), each = 365), rep(1:365, 1135) , pch = 15, cex = 0.2, col = startCol[findInterval(startMat, breaks)], yaxt = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i")
axis(side = 2, at = axisDates, labels = month.abb, las = 1)
mtext(side = 3, adj = 0, line = 0.5, "A) Starting rate")

# Stopping rates
plot(rep(seq(1, 1135*2, 2), each = 365), rep(1:365, 1135), pch = 15, cex = 0.2, col = stopCol[findInterval(stopMat, breaks)], yaxt = "n", xaxs = "i", yaxs = "i")
axis(side = 2, at = axisDates, labels = FALSE)
abline(v = c(660, 1500, 1650), lty = 2, col = grey(0.8))
mtext(side = 3, adj = 0, line = 0.5, "B) Stopping rate")

# Legend
yBreaks <- seq(1,365,length.out = 100)[round(seq(1, 100, length.out = 11))]
for(i in 1:100){
	polygon(x = c(2400, 2400, 2550, 2550), y = seq(1,365,length.out = 101)[c(i, i+1, i+1, i)], col = startCol[i], border = NA, xpd = NA)
	polygon(x = c(2550, 2550, 2700, 2700), y = seq(1,365,length.out = 101)[c(i, i+1, i+1, i)], col = stopCol[i], border = NA, xpd = NA)
}
for(i in 1:10){
	# polygon(x = c(2400, 2400, 2550, 2550), y = yBreaks[c(i, i+1, i+1, i)], col = startCol[(i-1)*10 + 1], border = NA, xpd = NA)
	# polygon(x = c(2550, 2550, 2700, 2700), y = yBreaks[c(i, i+1, i+1, i)], col = stopCol[(i-1)*10 + 1], border = NA, xpd = NA)
	segments(x0 = 2400, x1 = 2750, y0 = yBreaks[i], y1 = yBreaks[i], xpd = NA)
	text(2750, yBreaks[i], seq(0, 1, 0.1)[i], xpd = NA, pos = 4, cex = 0.8)
}
segments(x0 = 2400, x1 = 2750, y0 = 365, y1 = 365, xpd = NA)
text(2750, 365, "1.0", xpd = NA, pos = 4, cex = 0.8)
polygon(x = c(2400, 2400, 2700, 2700), y = c(1, 365, 365, 1), xpd = NA)

# Axes labels
mtext(side = 1, line = 2, outer = TRUE, "Distance along migration (km)")
mtext(side = 2, line = 3, outer = TRUE, "Time of year")


dev.off()

###############################################################################
# Host density through time
###############################################################################
