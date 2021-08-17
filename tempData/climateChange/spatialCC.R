###############################################################################
# Look at output data from climate model
###############################################################################

# Monthly projections from CIMP5 suite of models
# for 2080-2099 under RCP 2.6 and RCP 8.5

# Compare monolevel temperature at surface
# to subsurface "soil temperatures"

library(ncdf4) # package for netcdf manipulation
library(raster)

#-----------------------------------------------------------------------------
# Import RCP 2.6 data
#-----------------------------------------------------------------------------

cc26 <- nc_open("tempData/climateChange/monthly/cmip5_anomaly_tas_monthly_mean_multi-model-ensemble_rcp26_2080-2099.nc")

cc26C_can <- nc_open("tempData/climateChange/monthly/geomet-climate-CMIP5.TT.RCP26.ENS.ABS_PCTL50.nc")
# sink('tempData/climateChange/monthly/geomet-climate-CMIP5.TT.RCP26.ENS.ABS_PCTL50_metadata.txt')
# print(cc26C_can)
# sink()

# # Save the print(nc) dump to a text file
# This tells you what variables are in the netcdf file and their dimensions
# Don't need to run it every time, just when you get the dataset for the first
# time and are wondering what's in there
# sink('ClimateChangeData/cmip5_anomaly_tas_monthly_mean_multi-model-ensemble_rcp26_2080-2099_metadata.txt')
# print(cc26)
# sink()

lat <- ncvar_get(cc26, "latitude_bounds")
lon <- ncvar_get(cc26, "longitude_bounds")

tas <- ncvar_get(cc26, "tas")

dim(tas)# longitude, latitude, time (month)


#-----------------------------------------------------------------------------
# Look at grid
#-----------------------------------------------------------------------------
library(PBSmapping)

gshhg <- "~/Google Drive/Mapping/gshhg-bin-2.3.7/"
xlim <- c(-121.5, -106.5)
ylim <- c(60.5, 68)
land <- importGSHHS(paste0(gshhg,"gshhs_i.b"), xlim = xlim  + 360, ylim = ylim, maxLevel = 2, useWest = TRUE)
rivers <- importGSHHS(paste0(gshhg,"wdb_rivers_i.b"), xlim = xlim  + 360, ylim = ylim, useWest = TRUE)
borders <- importGSHHS(paste0(gshhg,"wdb_borders_i.b"), xlim = xlim  + 360, ylim = ylim, useWest = TRUE, maxLevel = 1)

# Find values within this grid
lat.in <- unique(which(lat > ylim[1] & lat < ylim[2], arr.ind = TRUE)[, 2])
lon.in <- unique(which(lon > xlim[1] & lon < xlim[2], arr.ind = TRUE)[, 2])

# As an example, extract projected change for January for each point
JanAnom <- array(data = NA, dim = c(length(lon.in), length(lat.in)))
gridPoint <- array(data = NA, dim = c(length(lon.in), length(lat.in), 2))
for(i in 1:length(lon.in)){
	for(j in 1:length(lat.in)){
		JanAnom[i, j] <- tas[lon.in[i], lat.in[j], 1]
		gridPoint[i, j, 1] <- mean(lon[, lon.in[i]])
		gridPoint[i, j, 2] <- mean(lat[, lat.in[j]])
	}
}

quartz()
plotMap(land, xlim = xlim, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6))

for(i in lon.in){
	for(j in lat.in){
		polygon(x = lon[c(1, 1, 2, 2), i], y = lat[c(1,2,2,1), j])
	}
}

points(gridPoint[, , 1], gridPoint[, , 2], col = heat.colors(n = 100)[findInterval(JanAnom, seq(2.3, 3.2, length.out = 100))], pch = 19)

###############################################################################
# Match migration to grid
###############################################################################

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

#------------------------------------------------------------------------------
# Create spatial dataset for points for NARR ground temperature data 
#------------------------------------------------------------------------------

gr <- read.csv("gridRange.csv")

# Path over spatial grid for each DOY
# These are the EIDs from gr corresponding to the approx. location of the herd on each DOY
path <- c(rep(78, 110), #winter range
					rep(c(79, 97, 119, 143, 166, 167, 190, 212, 231, 250, 266, 276, 285, 292), each = 3), #spring migration
					rep(299, 14), #calving
					rep(c(298, 297, 290, 283, 275), each = 3), #post-calving migration
					rep(c(265, 249, 230), each = 16), #summer grazing dispersal
					rep(230, 21), #late summer grazing (all stopped)
					rep(c(210, 187, 186, 185, 184, 161, 138, 115, 93, 77), each = 5), # fall migration DOY 251-290
					rep(78, 65)) #winter

migPoints <- as.EventData(data.frame(
	EID = c(1:365),
	X = gr$X[match(path, gr$EID)],
	Y = gr$Y[match(path, gr$EID)]
), projection = "LL")

addPoints(migPoints)

#------------------------------------------------------------------------------
# Find grid space for climate change projections that corresponds to each DOY
# in space given the migration route
#------------------------------------------------------------------------------

gridMatch <- findPolys(events = migPoints, polys = ccGrid)

PID_DOY <- gridMatch$PID[order(gridMatch$EID)] #The polys for each DOY

################################################################################
# Extract Climate projections for these grid spaces
###############################################################################

#-----------------------------------------------------------------------------
# Import data
#-----------------------------------------------------------------------------

cc26 <- nc_open("ClimateChangeData/cmip5_anomaly_tas_monthly_mean_multi-model-ensemble_rcp26_2080-2099.nc")
tas26 <- ncvar_get(cc26, "tas")

cc85 <- nc_open("ClimateChangeData/cmip5_anomaly_tas_monthly_mean_multi-model-ensemble_rcp85_2080-2099.nc")
tas85 <- ncvar_get(cc85, "tas")

#-----------------------------------------------------------------------------
# Create data frame
#-----------------------------------------------------------------------------
month <- as.numeric(strftime(as.Date(paste(2000, c(1:365), sep = "-"), format = "%Y-%j"), format = "%m"))

ccAnomalies <- array(data = NA, dim = c(365, 2))

for(i in 1:365){
	ccAnomalies[i, 1] <- tas26[unique(ccGrid$lon.in[ccGrid$PID == PID_DOY[i]]), unique(ccGrid$lat.in[ccGrid$PID == PID_DOY[i]]), month[i]]
	ccAnomalies[i, 2] <- tas85[unique(ccGrid$lon.in[ccGrid$PID == PID_DOY[i]]), unique(ccGrid$lat.in[ccGrid$PID == PID_DOY[i]]), month[i]]
}

#-----------------------------------------------------------------------------
# Look at climate change anomalies
#-----------------------------------------------------------------------------

plot(DOY, ccAnomalies[, 1], "l", col = 4, ylim = range(ccAnomalies), ylab = "Temperature anomaly (*C)", lwd = 2, las = 1)
lines(DOY, ccAnomalies[, 2], col = 2, lwd = 2)

lines(lowess(ccAnomalies[, 1], f = 1/6), col =4)
lines(lowess(ccAnomalies[, 2], f = 1/6), col =2)

abline(v = c(110, 155), lty = 3)

par(new = TRUE)
plot(DOY, PID_DOY, "l")


################################################################################
# Create master dataframe
###############################################################################

tempOverMigration <- data.frame(
	DOY = c(1:365),
	lon = migPoints$X,
	lat = migPoints$Y,
	temp = tempDOY$temp, 
	rcp26 = ccAnomalies[, 1],
	rcp85 = ccAnomalies[, 2])

write.csv(tempOverMigration, file = "ClimateChangeData/tempOverMigration.csv")

################################################################################
# Create maps illustrating process
###############################################################################

# Map of NARR data and points
# including migration route
# Layer on CC data grid

tempDat <- read.csv("ClimateChangeData/tempOverMigration.csv")
allTemps <- read.csv("LST/NARR/migTemp_allYears.csv")

DOY <- c(1:365)
xDate <- as.Date(paste(2000, DOY, sep = "-"), format = "%Y-%j")

#------------------------------------------------------------------------------
# Smooth climate change predictions as done for annual temps
#------------------------------------------------------------------------------

# Period smoother: to make ends match, paste together and then cut out middle
x <- c(1:(365*3)) 
y1 <- rep(tempDat$rcp26, 3)
L1 <- lowess(y1 ~ x, f = 1/20)
tempDat$rcp26Smoothed <- L1$y[366:(2*365)]

y2 <- rep(tempDat$rcp85, 3)
L2 <- lowess(y2 ~ x, f = 1/20)
tempDat$rcp85Smoothed <- L2$y[366:(2*365)]

tempDat$rcp26added <- tempDat$temp + tempDat$rcp26Smoothed
tempDat$rcp85added <- tempDat$temp + tempDat$rcp85Smoothed

#------------------------------------------------------------------------------
# Plot current, anomaly, and projected
#------------------------------------------------------------------------------
quartz(width = 4, height = 5.5)
par(mfrow = c(3, 1), mar = c(2,6,1,1), oma = c(3,0,0,0))
plot(xDate, tempDat$temp, "n", xlab = "", ylab = "Smoothed historical\nsurface temperature", las = 1, ylim = c(-20, 30))
abline(v = xDate[DOYhighlight], col = 3, lwd = 1.5)
for(i in 1:21) lines(xDate, allTemps[i, ], col = grey(0.8), lwd = 0.8)
lines(xDate, tempDat$temp, lwd = 2)
abline(h = 0, lty = 2)

plot(xDate, tempDat$rcp26, "l", col = 4, ylim = range(tempDat[, c('rcp26', 'rcp85')]), ylab = "Projected\ntemperature anomalies", las = 1)
abline(v = xDate[DOYhighlight], col = 3, lwd = 1.5)
lines(xDate, tempDat$rcp85, col = 2)
lines(xDate, tempDat$rcp26Smoothed, col = 4, lwd = 2)
lines(xDate, tempDat$rcp85Smoothed, col = 2, lwd = 2)

plot(xDate, tempDat$temp, "n", xlab = "", ylab = "Projected\nsurface temperature", las = 1, ylim = c(-20, 30))
abline(v = xDate[DOYhighlight], col = 3, lwd = 1.5)
lines(xDate, tempDat$temp)
lines(xDate, tempDat$rcp26added, col = 4, lwd = 2)
lines(xDate, tempDat$rcp85added, col = 2, lwd = 2)
abline(h = 0, lty = 2)

# tempDat$rcp85added[which(round(tempDat$rcp85added) == 0)]
segments(x0 = xDate[c(147,283)], x1 = xDate[c(147,283)], y0 = 0, y1 = -25, lty = 3)
segments(x0 = xDate[c(139, 292)], x1 = xDate[c(139, 292)], y0 = 0, y1 = -25, col = 4, lty = 3)
segments(x0 = xDate[c(123, 324)], x1 = xDate[c(123, 324)], y0 = 0, y1 = -25, col = 2, lty = 3)


#------------------------------------------------------------------------------
# Plot map of assumed migration route
#------------------------------------------------------------------------------
library(PBSmapping)

gshhg <- "~/Google Drive/Mapping/gshhg-bin-2.3.7/"
xlim <- c(-125, -105) + 360
ylim <- c(61, 69)
land <- importGSHHS(paste0(gshhg,"gshhs_i.b"), xlim = xlim, ylim = ylim, maxLevel = 2, useWest = TRUE)
rivers <- importGSHHS(paste0(gshhg,"wdb_rivers_i.b"), xlim = xlim, ylim = ylim, useWest = TRUE)
borders <- importGSHHS(paste0(gshhg,"wdb_borders_i.b"), xlim = xlim, ylim = ylim, useWest = TRUE, maxLevel = 1)

# Load list of grid points within Bathurst caribou range
gr <- read.csv("gridRange.csv")

# Path over spatial grid for each DOY
# These are the EIDs from gr corresponding to the approx. location of the herd on each DOY
path <- c(rep(78, 110), #winter range
					rep(c(79, 97, 119, 143, 166, 167, 190, 212, 231, 250, 266, 276, 285, 292), each = 3), #spring migration
					rep(299, 14), #calving
					rep(c(298, 297, 290, 283, 275), each = 3), #post-calving migration
					rep(c(265, 249, 230), each = 16), #summer grazing dispersal
					rep(230, 21), #late summer grazing (all stopped)
					rep(c(210, 187, 186, 185, 184, 161, 138, 115, 93, 77), each = 5), # fall migration DOY 251-290
					rep(78, 65)) #winter

DOYhighlight <- as.numeric(strftime(as.Date(paste("2000", c("01-01", "05-15", "06-07", "07-01", "10-01", "12-01"), sep = "-"), format = "%Y-%m-%d"), format = "%j"))

quartz(width = 6.3, height = 5.5)
plotMap(land, xlim = xlim - 360, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6))
addPolys(ccGrid, border = 2, lwd = 0.8)


points(gr$X, gr$Y, pch = 19, cex = 0.3)
lines(gr$X[path], gr$Y[path], lwd = 3, col = "#00000050")
points(gr$X[path[DOYhighlight]], gr$Y[path[DOYhighlight]], col = 3, lwd = 1.5, cex=1.5)

legend("topleft", pch = c(NA, 1, 19, NA), lwd = c(3, NA, NA, 0.8), pt.lwd = c(NA, 1.5, NA, NA), pt.cex = c(NA, 1.5, 0.3, NA), col = c("#00000050", 3, 1, 2), legend = c("Migration route", "Anchor time points", "NARR grid", "Climate projection grid"), bg = "white")
