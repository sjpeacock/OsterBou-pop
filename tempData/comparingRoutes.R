heatCol <- colorRampPalette(c(4, "gold", 2, 1))(n = 100)

###############################################################################
# Load basemap
###############################################################################

library(PBSmapping)

gshhg <- "~/Google Drive/Mapping/gshhg-bin-2.3.7/"
xlim <- c(-122, -106) + 360
ylim <- c(61, 68)
land <- importGSHHS(paste0(gshhg,"gshhs_i.b"), xlim = xlim, ylim = ylim, maxLevel = 2, useWest = TRUE)
rivers <- importGSHHS(paste0(gshhg,"wdb_rivers_i.b"), xlim = xlim, ylim = ylim, useWest = TRUE)
borders <- importGSHHS(paste0(gshhg,"wdb_borders_i.b"), xlim = xlim, ylim = ylim, useWest = TRUE, maxLevel = 1)

# Color for ranges
colRange <- c(winter = 4, spring = 3, calving = 7, summer = 2, fall = 1)

# Function to plot circle
plotCircle <- function(
	LonDec, #LonDec = longitude in decimal degrees
	LatDec, #LatDec = latitude in decimal degrees of the center of the circle
	Km, #Km = radius of the circle in kilometers
	plot = TRUE # Plot the circle or return coordinates?
) {
	ER <- 6371 #Mean Earth radius in kilometers. C
	AngDeg <- seq(1:360) #angles in degrees 
	Lat1Rad <- LatDec*(pi/180) #Latitude of the center of the circle in radians
	Lon1Rad <- LonDec*(pi/180) #Longitude of the center of the circle in radians
	AngRad <- AngDeg*(pi/180) #angles in radians
	Lat2Rad <-asin(sin(Lat1Rad)*cos(Km/ER)+cos(Lat1Rad)*sin(Km/ER)*cos(AngRad)) #Latitude of each point of the circle rearding to angle in radians
	Lon2Rad <- Lon1Rad+atan2(sin(AngRad)*sin(Km/ER)*cos(Lat1Rad),cos(Km/ER)-sin(Lat1Rad)*sin(Lat2Rad))#Longitude of each point of the circle rearding to angle in radians
	Lat2Deg <- Lat2Rad*(180/pi)#Latitude of each point of the circle rearding to angle in degrees (conversion of radians to degrees deg = rad*(180/pi) )
	Lon2Deg <- Lon2Rad*(180/pi)#Longitude of each point of the circle rearding to angle in degrees (conversion of radians to degrees deg = rad*(180/pi) )
	if(plot == FALSE){
		return(cbind(X = Lon2Deg, Y = Lat2Deg))
	} else {
		polygon(Lon2Deg, Lat2Deg,lty=2)
	}
}
###############################################################################
# CARMA MERRA data on Bathurst ranges
###############################################################################

ranges <- list(
	winter = importShapefile("ranges/BAHwinter_latlong.shp"),
	spring = importShapefile("ranges/BAHspring_latlong.shp"),
	calving = importShapefile("ranges/BAHcalving_latlong.shp"),
	summer = importShapefile("ranges/BAHsummer_latlong.shp"),
	fall = importShapefile("ranges/BAHfall_latlong.shp")
)
annualRange <- importShapefile("ranges/BAHwinter_latlong.shp")

rangeCol <-c()
plotMap(land, xlim = xlim - 360, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6))
for(i in 1:5){
	addPolys(ranges[[i]], border = colRange[i], col = NULL, lwd = 2)
}

rangesCombined <- as.PolySet(data.frame(
	PID = c(rep(1, nrow(ranges$winter)), rep(2, nrow(ranges$spring)), rep(3, nrow(ranges$calving)), rep(4, nrow(ranges$summer)), rep(5, nrow(ranges$fall))),
	SID = c(ranges$winter$SID, ranges$spring$SID, ranges$calving$SID, ranges$summer$SID, ranges$fall$SID),
	POS = c(ranges$winter$POS, ranges$spring$POS, ranges$calving$POS, ranges$summer$POS, ranges$fall$POS),
	X = c(ranges$winter$X, ranges$spring$X, ranges$calving$X, ranges$summer$X, ranges$fall$X),
	Y = c(ranges$winter$Y, ranges$spring$Y, ranges$calving$Y, ranges$summer$Y, ranges$fall$Y)
))

###############################################################################
# Load NARR data
###############################################################################

library(ncdf4) # package for netcdf manipulation
library(raster)

#-----------------------------------------------------------------------------
# Grid
#-----------------------------------------------------------------------------

# Load single month to get grid information

soilTemp <- nc_open("tempData/NARR/tsoil/tsoil.200001.nc")

# Extract variables of interest
lat <- ncvar_get(soilTemp, "lat")
lon <- ncvar_get(soilTemp, "lon")
tsoil <- ncvar_get(soilTemp, "tsoil")
depth <- ncvar_get(soilTemp, "level")

x <- c(1:(dim(tsoil)[1]))
y <- c(1:(dim(tsoil)[2]))

fullGrid <- as.EventData(data.frame(
	EID = c(1:(dim(tsoil)[1]*dim(tsoil)[2])),
	X = c(lon),
	Y = c(lat),
	x = rep(x, length(y)),
	y = rep(y, each = length(x))
))

z <- findPolys(fullGrid, rangesCombined)

clippedGrid <- fullGrid[fullGrid$EID %in% z$EID, ]

plotMap(land, xlim = xlim - 360, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6))
for(i in 1:5){
	addPolys(ranges[[i]], border = colRange[i], col = NULL, lwd = 2)
}
addPoints(clippedGrid, pch = 19, cex = 0.5, col = "#00000040")
text(clippedGrid$X, clippedGrid$Y, c(1:nrow(clippedGrid)), cex = 0.6)
text(clippedGrid$X, clippedGrid$Y, paste(clippedGrid$x, clippedGrid$y, sep = "\n"), cex = 0.4)

clippedGrid$num <- c(1:nrow(clippedGrid))

#-----------------------------------------------------------------------------
# Extract these grid points for all days
#-----------------------------------------------------------------------------

# Create array to store temps for each DOY from 2000-2020 (21 years)
allTemps <- array(NA, dim = c(21, 365, nrow(clippedGrid)))

m <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
startDay <- numeric(12)
startDay[1] <- 1
endDay <- numeric(12)

for(y in 1:21){
	for(i in 1:12){
		zz <- nc_open(paste0("tempData/NARR/tsoil/tsoil.", c(2000:2020)[y], m[i], ".nc"))
		z.tsoil <- ncvar_get(zz, "tsoil")
		
		# # What are the depths?
		# ncvar_get(zz, "level")
		
		# In first year, store number of days; don't need to repeat for all years
		endDay[i] <- min(365, startDay[i] + dim(z.tsoil)[4] - 1) 
		
		DOM <- 1 # Keep track of day of month (DOM)
		for(j in startDay[i]:endDay[i]){ # For each DOY in that month
			for(k in 1:nrow(clippedGrid)){
				allTemps[y, j, k] <- z.tsoil[clippedGrid$x[k], clippedGrid$y[k], 1, DOM] - 273.15
				}
			DOM <- DOM + 1
		}
		if(i < 12) startDay[i+1] <- endDay[i] + 1		
		rm(zz, z.tsoil)
		
	} # end months i
} # end years y

###############################################################################
# Run MCMC for different points
###############################################################################
days_in_each_range <- c(109, 43, 15, 12, 70, 40, 15, 30, 31)
DOY_range <- rep(c(1, 2, 3, 3, 4, 5, 5, 5, 1), times = days_in_each_range)

n.mcmc <- 10
#-----------------------------------------------------------------------------
# For n.mcmc different runs, select a random location within each range polygon
#-----------------------------------------------------------------------------

# # What range is each point in space? When/where are the caribou?
# V <- readRDS("simulations/output/V_host_100.rds")
# xPeak <- numeric(365)
# for(i in 1:365){
# 	xPeak[i] <- which(apply(V[, , i], 2, sum) == max(apply(V[, , i], 2, sum)))[1]
# }
# plot(1:365, xPeak, "l", ylab = "spatial grid space of peak host population", xlab = "Day of year", xaxs = "i", yaxs = "i")
# abline(h = seq(0, 2200, 100), lty = 3)
# abline(v = c(110, 153, 168, 180, 250, 290, 305, 335), col = c(3, 7, 7, 2, 1, 1, 1, 4), lwd = 2)

xRange <- c(
	rep(1, 50), # winter range
	rep(2, 350), # spring migration
	rep(3, 100), #calving and post-calving
	rep(4, 100), # summer
	rep(5, 450), # fall
	rep(1, 85)) # winter

# for(i in 1:5){
# 	segments(
# 	x0 = 400, 
# 	x1 = 400, 
# 	y0 = min(which(xRange == i)), 
# 	y1= max(which(xRange == i)),
# 	col = colRange[i], lwd = 5, xpd = NA)
# }

# Migration route is 2268 km (2 km grid = 1135 grid spaces)
mcmcGridPoints <- array(NA, dim = c(1135, n.mcmc))
mcmcSteps <- list(); length(mcmcSteps) <- n.mcmc

set.seed(398475)
for(n in 1:n.mcmc){
	#---
	# Select a point in each range (*excluding spring range*)
	# (Spring range is pretty broad and made for very circuitous migrations to calving grounds)
	#---
	stopover <- clippedGrid[1:4,]
	for(i in 1:4){
		eid <- sample(findPolys(clippedGrid, ranges[[c(1,3:5)[i]]])$EID, size = 1)
		stopover[i, ] <- clippedGrid[clippedGrid$EID == eid,]
	}
	# # Add winter back to the end of the stopover
	# addPoints(stopover, col = colRange[c(1,3:5)], pch = 8, cex = 1.5, lwd = 2)
	
	#---
	# Create path from those points
	#---
	step <- list(); length(step) <- 4
	# find next grid points in that circle
	for(i in 1:4){ # for each stopover location
		step[[i]] <- stopover[i, ]
		
		# Set next stopover to winter (1) if in fall
		if(i < 4) nextStopover <- stopover[i+1, ] else nextStopover <- stopover[1, ]
		j <- 1
		
		distToNext <- calcGCdist(lon1 = nextStopover$X, lat1 = nextStopover$Y, lon2 = step[[i]]$X[j], step[[i]]$Y[j])$d
		
		while(distToNext > 0){
			j <- j + 1
			# For a current step, create a circle with 50 km radius
			circle <- plotCircle(
				LonDec = tail(step[[i]]$X, 1), 
				LatDec =  tail(step[[i]]$Y, 1),
				Km = 50,
				plot = FALSE)
			circle <- as.PolySet(x = data.frame(PID = rep(1, 360), POS = c(1:360), X = circle[,"X"], Y = circle[, "Y"]), projection = "LL")
			# addPolys(circle, border = colRange[c(2, 3, 4, 5)[i]], col = NULL)
			
			# Identify potential next steps
			zz <- findPolys(clippedGrid, circle)
			
			# Calculate distance between each of those points and the next stopover
			d.zz <- calcGCdist(lon1 = nextStopover$X, lat1 = nextStopover$Y, lon2 = clippedGrid$X[clippedGrid$EID %in% zz$EID], clippedGrid$Y[clippedGrid$EID %in% zz$EID])$d
			
			if(sum(d.zz == 0) == 0){
				# standardize the distances
				n.zz <- (d.zz - mean(d.zz))/sd(d.zz)
				step[[i]][j, ] <- clippedGrid[clippedGrid$EID == zz$EID[sample(1:length(d.zz), prob = 1/(1 + exp(3*n.zz)), size = 1)], ]
			} else {
				step[[i]][j, ] <- clippedGrid[clippedGrid$EID == zz$EID[which(d.zz == 0)], ]
			}
			
			# Calculate distance between that step and the next stopover; if that distance
			# is > 0 then keep going with another step
			distToNext <- calcGCdist(lon1 = nextStopover$X, lat1 = nextStopover$Y, lon2 = step[[i]]$X[j], step[[i]]$Y[j])$d
		} # end steps for that season
	} #end path
	
	mcmcSteps[[n]] <- step
	nPoints <- length(unlist(step))
	# migPathShort <- as.PolySet(data.frame(
	# 	PID = rep(1:4, times = sapply(step, nrow)),
	# 	POS = unlist(sapply(sapply(step, nrow), function(X){c(1:X)})),
	# 	X = c(step[[1]]$X, step[[2]]$X, step[[3]]$X, step[[4]]$X),
	# 	Y = c(step[[1]]$Y, step[[2]]$Y, step[[3]]$Y, step[[4]]$Y)
	# ))
	# 
	# addLines(migPathShort, lwd = 3, col = "#00000080")
	
	nP <- sapply(step, nrow)
	migPathLong <- c(
		rep(step[[1]]$EID, each = ceiling(330/nP[1]))[c(1:330)], # spring migration
		rep(step[[2]]$EID, each = ceiling(170/nP[2]))[c(1:170)],
		step[[3]]$EID,
		rep(step[[4]]$EID, each = ceiling((1135-330-170-nP[3])/nP[4]))[c(1:(1135-330-170-nP[3]))]
	)
	mcmcGridPoints[, n] <- migPathLong
	
}

#-----------------------------------------------------------------------------
# Look at 10 of those migration paths
#-----------------------------------------------------------------------------
par(mfrow = c(1,1), oma = rep(0, 4))
plotMap(land, xlim = xlim - 360, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6))
for(i in 1:5){
 	addPolys(ranges[[i]], border = colRange[i], col = NULL, lwd = 2)
 }
# addPoints(clippedGrid, pch = 19, cex = 0.5, col = "#00000040")
legend("bottomleft", lwd = 2, col = colRange, legend = names(colRange), bg = "white")
# legend("topleft", pch = c(19, NA, 19), pt.cex = c(0.5, NA, 1), col = c("#00000040", "#00000040", 1), lty = c(NA, 2, 1), bg = "white", legend = c("NARR grid", "CIMP5 grid", "Simulated migration path"))

for(i in 1:10){
	dumPoly <- as.PolySet(data.frame(
		PID = rep(1, 1135),
		POS = c(1:1135),
		X = clippedGrid$X[match(mcmcGridPoints[, i], clippedGrid$EID)],
		Y = clippedGrid$Y[match(mcmcGridPoints[, i], clippedGrid$EID)]
	))
	addPolys(dumPoly, lwd = 3, border = "#00000080")
}

i <- 2
for(j in 1:4){
	addPoints(mcmcSteps[[i]][[j]][1, ], col = colRange[c(1,3:5)][j], pch = 8, cex = 1.5, lwd = 2)
	for(k in 2:nrow(mcmcSteps[[i]][[j]])){
		arrows(x0 = mcmcSteps[[i]][[j]][k-1, 'X'], 
					 x1 = mcmcSteps[[i]][[j]][k, 'X'],
					 y0 = mcmcSteps[[i]][[j]][k-1, 'Y'], 
					 y1 = mcmcSteps[[i]][[j]][k, 'Y'],
					 length = 0.06, col = colRange[c(2, 3, 5, 1)][j])
	}
}

#-----------------------------------------------------------------------------
# Extract ground temperature data for 21 years along migration route
#-----------------------------------------------------------------------------

# Period smoother: to make ends match, paste together and then cut out middle
DOY <- c(1:365)
DOY.stitched <- rbind(t(array(rep(DOY, 21), dim = c(365, 21))),t(array(rep(DOY+365, 21), dim = c(365, 21))), t(array(rep(DOY+2*365, 21), dim = c(365, 21)))) 

mcmcTemp <- array(NA, dim = c(n.mcmc, 365, 1135))

for(n in 1:n.mcmc){
	for(i in 1:1135){
		
		# Extract temperatuers for given grid point
		allTemps.ni <- allTemps[, , clippedGrid$num[which(clippedGrid$EID == mcmcGridPoints[i, n])]]
		
		# Apply smoother
		y <- rbind(allTemps.ni, allTemps.ni, allTemps.ni)
		L2 <- lowess(y ~ DOY.stitched, f = 1/12)
		
		# Record smoothed annual curve along migration route
		mcmcTemp[n, , i] <- matrix(L2$y, ncol = 365*3, nrow = 21)[1, 366:(2*365)]
	
	}
}

# write.csv(allTemps[, , clippedGrid$num[which(clippedGrid$EID == mcmcGridPoints[330, n])]], file = "tempData/groundTemps_migPath6_330km.csv", row.names = FALSE)


# Randomly choose 4 of the temp matrices to plot
par(mfrow = c(2,2), mar = c(0,0,0,0), oma=c(5,5,2,2))
for(n in 1:4){
	plot(rep(seq(1, 365, 7), 114), rep(seq(1, 1135, 10), each = 53), pch = 19, cex = 0.8, col = heatCol[findInterval(mcmcTemp[n, seq(1, 365, 7), seq(1, 1135, 10)], seq(-30, 30, length.out = 100))])
}

#---------
# Plot range points and annual temperature curve at those points for four runs
#---------

par(mfrow = c(3,2), oma = rep(0, 4), mar = c(4,4,2,1))
for(n in c(1,3,6,8)){
	plotMap(land, xlim = xlim - 360, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6))
	for(i in 1:5){
		addPolys(ranges[[i]], border = colRange[i], col = NULL, lwd = 0.8)
	}
	
	for(j in 1:4){
		addPoints(mcmcSteps[[n]][[j]][1, ], col = colRange[c(1,3:5)][j], pch = 8, cex = 1.5, lwd = 2)
	}
	
	xPlot <- c(1, 330, 170, 170 + nrow(mcmcSteps[[n]][[3]]))
	plot(DOY, mcmcTemp[n, , xPlot[1]], "n", ylim = c(-22, 20))
	for(j in 1:4){
		lines(DOY, mcmcTemp[n, , xPlot[j]], col = colRange[c(1,3:5)][j])
		for(q in c(1:10)[c(1:10) %in% n == FALSE]){
			xPlot.q <- c(1, 330, 170, 170 + nrow(mcmcSteps[[q]][[3]]))
			lines(DOY, mcmcTemp[q, , xPlot.q[j]], col = paste0(colRange[c(1,3:5)][j], "40"))
		}
	}
}
	


par(mfrow = c(1,1), mar = c(4,4,2,1))
plot(DOY, mcmcTemp[1, , xPlot[i]], "n", ylim = c(-22, 20))
for(n in 1:n.mcmc){
	for(i in 1:4){
		xPlot <- c(1, 330, 170, 170 + nrow(mcmcSteps[[n]][[3]]))
		lines(DOY, mcmcTemp[n, , xPlot[i]], col = paste0(c("#1F98E7", "#F6C70E", "#D05D6D", "#000000")[i], "60"), lwd = 2)
	}}


# saveRDS(mcmcTemp, file = "tempData/NARR/tempGrid_mcmc10.rds")

###############################################################################
# Load Climate Change data
###############################################################################

cc26 <- nc_open("tempData/climateChange/monthly/cmip5_anomaly_tas_monthly_mean_multi-model-ensemble_rcp26_2080-2099.nc")
cc85 <- nc_open("tempData/climateChange/monthly/cmip5_anomaly_tas_monthly_mean_multi-model-ensemble_rcp85_2080-2099.nc")

# Climat change data are a grid rather than points.

lat <- ncvar_get(cc26, "latitude_bounds")
lon <- ncvar_get(cc26, "longitude_bounds")

tas <- list(rcp26 = ncvar_get(cc26, "tas"),
						rcp85 = ncvar_get(cc85, "tas"))

dim(tas[[1]])# longitude, latitude, time (month)

# Find values within the annual range
limsBath <- list(X = range(clippedGrid$X), Y = range(clippedGrid$Y))
lat.in <- unique(which(lat > limsBath$Y[1] & lat < limsBath$Y[2], arr.ind = TRUE)[, 2])
lon.in <- unique(which(lon > limsBath$X[1] & lon < limsBath$X[2], arr.ind = TRUE)[, 2])

#------------------------------------------------------------------------------
# Create polyset out of grid space for climate change projections
#------------------------------------------------------------------------------

nPolys <- length(lon.in)*length(lat.in)

ccGrid <- data.frame(
	PID = rep(c(1:nPolys), each = 4),
	POS = rep(c(1:4), nPolys),
	X = rep(NA, nPolys),
	Y = rep(NA, nPolys),
	lat.in = rep(NA, nPolys),
	lon.in = rep(NA, nPolys)
)

for(i in 1:nPolys){
	ccGrid$X[ccGrid$PID == i] <- lon[c(1, 1, 2, 2), lon.in[rep(1:length(lon.in), each = length(lat.in))[i]]]
	ccGrid$Y[ccGrid$PID == i] <- lat[c(1, 2, 2, 1), lat.in[rep(1:length(lat.in), length(lon.in))[i]]]
	
	ccGrid$lat.in[ccGrid$PID == i] <- lat.in[rep(1:length(lat.in), length(lon.in))[i]]
	ccGrid$lon.in[ccGrid$PID == i] <- lon.in[rep(1:length(lon.in), each = length(lat.in))[i]]
	
}

ccGrid <- as.PolySet(ccGrid, projection = "LL")

###############################################################################
# Extract cc projections along migration path
###############################################################################

#------------------------------------------------------------------------------
# Use migration route n
#------------------------------------------------------------------------------

for(n in 1:n.mcmc){

mcmcPath <- rbind(mcmcSteps[[n]][[1]], mcmcSteps[[n]][[2]], mcmcSteps[[n]][[3]], mcmcSteps[[n]][[4]])

migPoints <- as.EventData(data.frame(
	EID = c(1:nrow(mcmcPath)),
	X = mcmcPath$X,
	Y = mcmcPath$Y
), projection = "LL")

migPathShort <- as.PolySet(data.frame(
	PID = rep(1, nrow(mcmcPath)),
	POS = migPoints$EID,
	X = migPoints$X,
	Y = migPoints$Y
), projection = "LL")

#------------------------------------------------------------------------------
# Plot map
#------------------------------------------------------------------------------
pdf(file = paste0("tempData/map_migPath", n, ".pdf"), width = 5, height = 5, pointsize = 10)
par(mfrow = c(1,1), oma = rep(0, 4))
plotMap(land, xlim = xlim - 360, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6))

addPoints(clippedGrid, pch = 19, cex = 0.5, col = "#00000040")
addPolys(ccGrid, lty = 3, border = "#00000040")


for(i in 1:5){
	addPolys(ranges[[i]], border = colRange[i], col = NULL, lwd = 2)
}
addLines(migPathShort)
addPoints(migPoints, pch = 19, cex = 0.8)

for(j in 1:4){
	addPoints(mcmcSteps[[n]][[j]][1, ], col = colRange[c(1,3:5)][j], pch = 1, cex = 1.5, lwd = 2)
	# for(k in 2:nrow(mcmcSteps[[i]][[j]])){
	# 	arrows(x0 = mcmcSteps[[i]][[j]][k-1, 'X'], 
	# 				 x1 = mcmcSteps[[i]][[j]][k, 'X'],
	# 				 y0 = mcmcSteps[[i]][[j]][k-1, 'Y'], 
	# 				 y1 = mcmcSteps[[i]][[j]][k, 'Y'],
	# 				 length = 0.06, col = colRange[c(2, 3, 5, 1)][j])
	# }
}

legend("bottomleft", lwd = 2, col = colRange, legend = names(colRange), bg = "white")
legend("topleft", pch = c(19, NA, 19), pt.cex = c(0.5, NA, 1), col = c("#00000040", "#00000040", 1), lty = c(NA, 2, 1), bg = "white", legend = c("NARR grid", "CIMP5 grid", "Simulated migration path"))
dev.off()

#------------------------------------------------------------------------------
# Find grid space for climate change projections that corresponds to each DOY
# in space given the migration route
#------------------------------------------------------------------------------

gridMatch <- findPolys(events = migPoints, polys = ccGrid)

PID_step <- gridMatch$PID[order(gridMatch$EID)] #The polys for each "step" along migration

#------------------------------------------------------------------------------
# Extract Climate projections for these grid spaces
#-----------------------------------------------------------------------------

ccAnomalies <- array(data = NA, dim = c(2, length(PID_step), 12))

for(r in 1:2){ # for the low and high emissions scenarios
	for(i in 1:length(PID_step)){ # for each unqiue location along the migration
	ccAnomalies[r, i,] <- tas[[r]][unique(ccGrid$lon.in[ccGrid$PID == PID_step[i]]), unique(ccGrid$lat.in[ccGrid$PID == PID_step[i]]), ]
}}

#------------------------------------------------------------------------------
# Expand to entire spatial grid for each DOY, and smooth
#------------------------------------------------------------------------------
# dim(mcmcTemp)
# [1]   10  365 1135

nP <- sapply(mcmcSteps[[n]], nrow)

month <- as.numeric(strftime(as.Date(paste(2000, c(1:365), sep = "-"), format = "%Y-%j"), format = "%m"))
ccAnomaliesLong <- array(NA, dim = c(2, 365, 1135))
for(r in 1:2){
	for(i in 1:12){ # for each month
		for(j in which(month == i)){
			ccAnomaliesLong[r, j, ] <- c(
			rep(ccAnomalies[r, 1:(nP[1]), i], each = ceiling(330/nP[1]))[c(1:330)], # spring migration
			rep(ccAnomalies[r, (nP[1] + 1):sum(nP[c(1:2)]), i], each = ceiling(170/nP[2]))[c(1:170)],
			ccAnomalies[r, (sum(nP[c(1:2)]) + 1):sum(nP[c(1:3)]), i],
			rep(ccAnomalies[r, (sum(nP[c(1:3)]) + 1):sum(nP[c(1:4)]), i], each = ceiling((1135-330-170-nP[3])/nP[4]))[c(1:(1135-330-170-nP[3]))]
		)
		}
	}}


# Period smoother: to make ends match, paste together and then cut out middle
DOY.stitched.cc <- c(1:(365*3)) 
ccAnomaliesSmooth <- array(NA, dim = c(2, 365, 1135))
for(r in 1:2){
	for(j in 1:1135){
		y1 <- rep(ccAnomaliesLong[r, ,j], 3)
		L1 <- lowess(y1 ~ DOY.stitched.cc, f = 1/12)
		ccAnomaliesSmooth[r, ,j] <- L1$y[366:(2*365)]
	}
}

# write.csv(t(rbind(ccAnomaliesLong[, , 330], ccAnomaliesSmooth[, ,330])), "tempData/ccAnomalies_migPath6_330km.csv", row.names = FALSE)


# Plot check:
xDate <- as.Date(paste("2000", c(1:12), 1, sep = "-"))

par(mfrow = c(1,2), mar = c(5,5,5,1))
for(x in c(1,330)){
plot(DOY, ccAnomaliesSmooth[1, , x], "l", col = 4, ylim = range(ccAnomalies), lwd = 2, ylab = "Projected temperature anomalies", xaxs = "i")
lines(DOY, ccAnomaliesLong[1, , x], col = 4, lty = 3)
lines(DOY, ccAnomaliesSmooth[2, , x], col = 2)
lines(DOY, ccAnomaliesLong[2, , x], col = 2, lty = 3)
axis(side = 3, at = as.numeric(strftime(xDate, format = "%j")), labels = month.abb)
if(x == 330) mtext(side = 3, line = 2.5, "Location: calving grounds (x = 660 km)")
if(x == 1) mtext(side = 3, line = 2.5, "Location: winter range (x = 0 km)")
}

for(t in c(1, 110)){
	plot(c(1:1135), ccAnomaliesSmooth[1, t, ], "l", col = 4, ylim = range(ccAnomalies), lwd = 2, xlab = "Spatial grid location (km/2)", ylab = "Projected temperature anomalies", xaxs = "i")
	lines(c(1:1135), ccAnomaliesLong[1, t, ], col = 4, lty = 3)
	lines(c(1:1135), ccAnomaliesSmooth[2, t, ], col = 2)
	lines(c(1:1135), ccAnomaliesLong[2, t, ], col = 2, lty = 3)
	if(t == 1) mtext(side = 3, line = 2.5, "Time: winter (t = 0 days)")
	if(t == 110) mtext(side = 3, line = 2.5, "Time: calving (t = 110 days)")
}

###############################################################################
# Overall output
###############################################################################

temps <- list(
	current = mcmcTemp[n, , ],
	rcp26 = mcmcTemp[n, , ] + ccAnomaliesSmooth[1, , ],
	rcp85 = mcmcTemp[n, , ] + ccAnomaliesSmooth[2, , ]
)

saveRDS(temps, file = paste0("tempData/allTemps_migPath", n, ".rds"))
}