# Map of the Bathurst caribou migration and ranges

library(PBSmapping)
gshhg <- "~/Google Drive/Mapping/gshhg-bin-2.3.7/"
xlim <- c(-125, -100) + 360
ylim <- c(60, 70)

# Inset
xlim <- c(-150, -80) + 360
ylim <- c(45, 75)

land <- importGSHHS(paste0(gshhg,"gshhs_i.b"), xlim = xlim, ylim = ylim, maxLevel = 2, useWest = TRUE)
rivers <- importGSHHS(paste0(gshhg,"wdb_rivers_i.b"), xlim = xlim, ylim = ylim, useWest = TRUE)
borders <- importGSHHS(paste0(gshhg,"wdb_borders_i.b"), xlim = xlim, ylim = ylim, useWest = TRUE, maxLevel = 1)

#------------------------------------------------------------------------------
# Plot map
quartz(width = 5, height = 6)
plotMap(land, xlim = xlim - 360, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6))

# Inset
plotMap(land, xlim = xlim - 360, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6), xaxt = "n", yaxt = "n")
addLines(borders, col = grey(0.4))
polygon(x = c(-122, -106, -106, -122), y = c(61, 61, 68, 68), border = 1, col = NA, lwd = 2)

#######################################
# Need a migration path that is 2268 km
# too hard to figure it out that way

L <- locateEvents(type = "o", col = 2)
LP <- as.PolySet(data.frame(PID = rep(1, dim(L)[1]), POS = L$EID, X = L$X, Y = L$Y), projection = "LL")

# LP <- LP[which(LP$POS %in% c(54, 55) == FALSE), ]
# LP$POS <- c(1:dim(LP)[1])
calcLength(LP, rollup = 3, close = FALSE)

addPolys(LP, border = 2, lwd = 2)
text(LP$X, LP$Y, LP$POS)
# write.csv(LP, "migPath1.csv")

LP2 <- LP

# Create polygons for different migration seasons
# Look at mean temperature within that range

spring <- locateEvents(type = "o", col = 3, pch = 19, cex = 0.8)
springPoly <- as.PolySet(data.frame(PID = rep(1, dim(spring)[1]), POS = spring$EID, X = spring$X, Y = spring$Y), projection = "LL")
write.csv(springPoly, "springRange.csv")

# Plot lat/lon points

points(lon, lat, cex = 0.3, pch = 19)

grid <- as.EventData(data.frame(
	EID = c(1:(349*277)),
	X = as.numeric(lon),
	Y = as.numeric(lat),
	x = rep(c(1:349), 277),
	y = rep(c(1:277), each = 349)
))

points(lon, lat, cex = 0.3, pch = 19)
addPoints(grid, col = 2, cex = 0.3)

annualRange <- locateEvents(type = "o", col = 2)
annualRangePoly <- as.PolySet(data.frame(PID = rep(1, dim(annualRange)[1]), POS = annualRange$EID, X = annualRange$X, Y = annualRange$Y), projection = "LL")
# write.csv(annualRangePoly, "annualRangePoly.csv")

boop <- findPolys(grid, annualRangePoly)
head(boop)

gridRange <- grid[grid$EID %in% boop$EID, ]
addPoints(gridRange, col = 4, cex = 0.3)

# Winter
# Use average over winter range
winter <- locateEvents(type = "o", col = 3, pch = 19, cex = 0.8)
springPoly <- as.PolySet(data.frame(PID = rep(1, dim(spring)[1]), POS = spring$EID, X = spring$X, Y = spring$Y), projection = "LL")
write.csv(gridRange, "gridRange.csv")

# For the unique (x,y) pairs, how much difference is there in the annual ground temps?

dates <- as.Date(as.Date("2020-01-01"):as.Date("2020-12-31"), origin = "1970-01-01")

# Plot 10 different points
q <- gridRange[sample(1:dim(gridRange)[1], size = 10), c('x', 'y')]
quartz(width = 7, height = 4)
par(mar = c(4, 4, 2, 12))
plot(dates, allTemps[q[1,1], q[1,2], , 3], "l", xlab = y, ylab = "Temperature (*C)", las = 1, ylim = c(-50, 30), col = "#00000050", xaxs = "i")
for(i in 2:10) lines(dates, allTemps[q[i,1], q[i,2], , 3], col = "#00000050")
abline(h = 0)


# Look at annual migration path