# Plot current temps (temp[[1]])
ind.DOY <- seq(1, 365, 7)
n.DOY <- length(ind.DOY)

ind.x <- seq(1, 1135, 10)
n.x <- length(ind.x)

par(mfrow = c(2,2), oma = c(5,5,2,1), mar = rep(0.2, 4))

for(N in c(2, 6, 8, 10)){
temp <- readRDS(paste0("tempData/allTemps_migPath", N, ".rds"))


plot(rep(ind.DOY, n.x), rep(ind.x, each = n.DOY), pch = 15, col = heatCol[findInterval(temp[[1]][ind.DOY, ind.x], seq(-25, 20, length.out = 100))], xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i")
mtext(side = 3, line = -2, font = 2, col = "white", paste0("Migration path #", N))
points(110, 600, pch = 1, col = 1, lwd = 2, cex = 2)
points(153, 600, pch = 1, col = "white", lwd = 2, cex = 2)
}
mtext(side = 2, outer= TRUE, "Spatial grid location (km/2)", line = 2)
mtext(side = 1, outer= TRUE, "DOY", line = 2)
