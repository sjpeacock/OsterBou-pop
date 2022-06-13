###############################################################################
# Look at output data from climate model
###############################################################################

# Climate data from
# https://psl.noaa.gov/data/gridded/data.narr.monolevel.html

# Compare monolevel temperature at surface
# to subsurface "soil temperatures"

library(ncdf4) # package for netcdf manipulation
library(raster)


#-----------------------------------------------------------------------------
# Monolevel: air at surface
#-----------------------------------------------------------------------------

surfTemp <- nc_open("tempData/NARR/otherDat/air.sfc.1979.nc")
			
# # Save the print(nc) dump to a text file
# This tells you what variables are in the netcdf file and their dimensions
# Don't need to run it every time, just when you get the dataset for the first
# time and are wondering what's in there
# sink('LST/air.sfc.1979_metadata.txt')
# print(surfTemp)
# sink()

# Extract variables of interest
lat <- ncvar_get(surfTemp, "lat")
lon <- ncvar_get(surfTemp, "lon")
airTemp <- ncvar_get(surfTemp, "air")
t <- ncvar_get(surfTemp, "time")

# Find which point we want for exploration purposes
# These coordinates are approximately around the Bathurst calving range
ind <- c(which(lon == -110.5, arr.ind = TRUE), findInterval(66.5, lat))
plot(t, airTemp[200,160,] -273.15, "l")
lat[200,160]
lon[200,160]

samplePoints <- cbind(sample(c(1:349), size = 5000, replace = TRUE), sample(c(1:277), size = 5000, replace = TRUE))
plot(c(1,349), c(1, 277), "n")
points(samplePoints[, 1], samplePoints[, 2], col = heat.colors(n = 100)[round(airTemp[cbind(samplePoints, 1)] - 273.15) + 60], pch = 19)
points(170, 200, pch = 8)

lat[170,200]																																				 
lon[170,200]																																				

plot(range(lon), range(lat), "n")
points(lon[samplePoints], lat[samplePoints], col = heat.colors(n = 100)[round(airTemp[cbind(samplePoints, 1)] - 273.15) + 60], pch = 19, cex = 0.5)
points(170, 200, pch = 8)

#-----------------------------------------------------------------------------
# Monolevel: air at 2 m
#-----------------------------------------------------------------------------
twoTemp <- nc_open("tempData/NARR/otherDat/air.2m.1979.nc")

# sink('LST/NARR/air.2m.1979_metadata.txt')
# print(twoTemp)
# sink()

airTemp2m <- ncvar_get(twoTemp, "air")

# Pull out air temps June 1 - 30
ind <- as.numeric(strftime(as.Date("2020-06-01"), format = "%j"))
ind <- c(ind:(ind + 29))
tair2 <- airTemp2m[170, 200, ind]

#-----------------------------------------------------------------------------
# Soil tempertaure 
#-----------------------------------------------------------------------------

soilTemp <- nc_open("tempData/NARR/otherDat/tsoil.197906.nc")

# sink('LST/tsoil.197906_metadata.txt')
# print(soilTemp)
# sink()

lat2 <- ncvar_get(soilTemp, "lat")
lon2 <- ncvar_get(soilTemp, "lon")

lat2[170,200] # Good, they are the same as the other dataset
lon2[170,200]

tsoil <- ncvar_get(soilTemp, "tsoil")

depth <- ncvar_get(soilTemp, "level")

# Look at 0 cm depth
tsoil0 <- tsoil[170, 200, 1, ]

# Pull out air temps June 1 - 30
ind <- as.numeric(strftime(as.Date("2020-06-01"), format = "%j"))
ind <- c(ind:(ind + 29))
tair0 <- airTemp[170, 200, ind]

dim(tsoil0)
dim(tair0)

plot(c(1:30), tair0 - 273.15, "l", xlab = "June 1979", ylab = "Temperature (*C)", las = 1, ylim = c(-20, 25), col = 4)
lines(c(1:30), tair2 - 273.15, col = grey(0.8))
lines(c(1:30), tsoil0 - 273.15,  col = 2)

lines(c(1:30), tsoil[170, 200, 2, ]- 273.15, pch = 18, col = 7)
lines(c(1:30), tsoil[170, 200, 3, ]- 273.15,  pch = 18, col = 7, lty = 2)
lines(c(1:30), tsoil[170, 200, 4, ]- 273.15, pch = 18, col = 7, lty = 4)
lines(c(1:30), tsoil[170, 200, 5, ]- 273.15,  pch = 18, col = 7, lty = 3)

legend("topleft", col = c(grey(0.8), 4, 2, 7, 7, 7, 7),  c("Air (2m)",  "Air (surface)", "Ground (0 cm)", "Ground (10 cm)", "Ground (40 cm)", "Ground (100 cm)", "Ground (800 cm)"), lty = c(rep(1, 4), 2, 4, 3), bty = "n")

###############################################################################
# Try downloading whole year for that single location
###############################################################################

y <- 2020 # Year

allTemps <- array(NA, dim = c(349, 277, 366, 7), dimnames = list(NULL, NULL, NULL, c("air.2m", "air.sfc", "soil0", "soil10", "soil40", "soil100", "soil800")))

# Air 2m
zz <- nc_open(paste0("LtempData/NARR/otherDat/air.2m.", y, ".nc"))
allTemps[, , , 1] <- ncvar_get(zz, "air")
rm(zz)

# Air surface
zz <- nc_open(paste0("tempData/NARR/otherDat/air.sfc.", y, ".nc"))
allTemps[, , , 2] <- ncvar_get(zz, "air")
rm(zz)

# Soil 
m <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
startDay <- numeric(12)
startDay[1] <- 1
endDay <- numeric(12)
for(i in 1:12){
	zz <- nc_open(paste0("tempData/NARR/otherDat/tsoil.", y, m[i], ".nc"))
	z.tsoil <- ncvar_get(zz, "tsoil")
	
	endDay[i] <- startDay[i] + dim(z.tsoil)[4] - 1
	for(j in 1:5){ # for each soil depth
		allTemps[, , startDay[i]:endDay[i], j + 2] <- z.tsoil[, , j, ]
	}
	
	if(i < 12) startDay[i+1] <- endDay[i] + 1
	rm(zz, z.tsoil)
}


allTemps <- allTemps - 273.15

#------------------------------------------------------------------------------
# Plot whole year for 2020
#------------------------------------------------------------------------------
dates <- as.Date(as.Date("2020-01-01"):as.Date("2020-12-31"), origin = "1970-01-01")

quartz(width = 7, height = 4)
par(mar = c(4, 4, 2, 12))
plot(dates, allTemps[170, 200, , 2], "l", xlab = y, ylab = "Temperature (*C)", las = 1, ylim = range(allTemps[170,200,,], na.rm = TRUE), col = 4, xaxs = "i")
abline(h = 0)
lines(dates, allTemps[170, 200,, 1], col = grey(0.8))
lines(dates, allTemps[170, 200,, 3],  col = 2, lwd = 2)

lines(dates, allTemps[170, 200,, 4], pch = 18, col = 7)
lines(dates, allTemps[170, 200,, 5], pch = 18, col = 7, lty = 2)
lines(dates, allTemps[170, 200,, 6], pch = 18, col = 7, lty = 4)
lines(dates, allTemps[170, 200,, 7], pch = 18, col = 7, lty = 3)

u <- par('usr')
legend(u[2] + 0.01 * (u[2] - u[1]), u[4], col = c(grey(0.8), 4, 2, 7, 7, 7, 7),  c("Air (2m)",  "Air (surface)", "Ground (0 cm)", "Ground (10 cm)", "Ground (40 cm)", "Ground (100 cm)", "Ground (800 cm)"), lty = c(rep(1, 4), 2, 4, 3), bty = "n", xpd = NA)

mtext(side = 3, line = 0.5, substitute(paste("Location: Approximate Bathurst calving grounds (", la, degree, "N, ", lo, degree, "E)", sep = ""), list(la = round(lat[170,200],3), lo = round(lon[170,200], 3))))
			
###############################################################################
###############################################################################
# Look at spatial grid
###############################################################################
###############################################################################

# This is where I was trying to figure out what points to use along the "migration
# route" of mh modelled carinbou.

# Might not matter if you're just wanting a single point in space (as above)

# I include it here just FYI, but you'll need additional datasets (shoreline)
# to produce the maps.

library(PBSmapping)

gshhg <- "~/Google Drive/Mapping/gshhg-bin-2.3.7/"
xlim <- c(-121.3, -106.6) + 360
ylim <- c(61, 67.8)
land <- importGSHHS(paste0(gshhg,"gshhs_i.b"), xlim = xlim, ylim = ylim, maxLevel = 2, useWest = TRUE)
rivers <- importGSHHS(paste0(gshhg,"wdb_rivers_i.b"), xlim = xlim, ylim = ylim, useWest = TRUE)
borders <- importGSHHS(paste0(gshhg,"wdb_borders_i.b"), xlim = xlim, ylim = ylim, useWest = TRUE, maxLevel = 1)

# What does the map look like on June 1 (DOY 153)?
temps <- allTemps[cbind(gridRange$x, gridRange$y, 153, 3)]
hc <- colorRampPalette(c("blue", "gold", "orange", "red"))(n = 227)

plotMap(land, xlim = xlim - 360, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6))

points(gridRange$X, gridRange$Y, col = hc[findInterval(temps, seq(-10.5, 12.1, 0.1))], pch = 19)
legend("topright", pch = c(rep(19, 5), NA), col = c(hc[findInterval(c(-10, -5, 0, 5, 10), seq(-10.5, 12.1, 0.1))], 1), legend = c(c(-10, -5, 0, 5, 10), "Range"), lwd = c(rep(NA, 5), 2), bg = "white")
mtext(side = 3, line = 1, "Soil temperature at 0 cm on June 1, 2020")
addPolys(annualRangePoly, border = 1, lwd = 2, col = NULL)


# Add EID to each point
gridRange$EID <- c(1: dim(gridRange)[1])
gridRange <- as.EventData(gridRange, projection = "LL")
plotMap(land, xlim = xlim - 360, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6))
addPoints(gridRange, col = 3)
text(gridRange$X, gridRange$Y, gridRange$EID, cex = 0.3)

############
# What the heck is going on?

gr <- read.csv("gridRange.csv")
# These are all the grid points within the Bathurst annual range


plot(gr$X, gr$Y, pch = 21, bg = grey(seq(0, 0.99, length.out = nrow(gr))))

plotMap(land, xlim = xlim - 360, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6))
text(gr$X, gr$Y, gr$EID, cex = 0.6)

path <- c(rep(78, 110), #winter range
					rep(c(79, 97, 119, 143, 166, 167, 190, 212, 231, 250, 266, 276, 285, 292), each = 3), #spring migration
					rep(299, 14), #calving
					rep(c(298, 297, 290, 283, 275), each = 3), #post-calving migration
					rep(c(265, 249, 230), each = 16), #summer grazing dispersal
					rep(230, 21), #late summer grazing (all stopped)
					rep(c(210, 187, 186, 185, 184, 161, 138, 115, 93, 77), each = 5), # fall migration DOY 251-290
					rep(78, 65)) #winter
					
plotMap(land, xlim = xlim - 360, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6))
points(gr$X, gr$Y, pch = 19, cex = 0.5, col = grey(0.5))
lines(gr$X[path], gr$Y[path])
points(jitter(gr$X[path], factor = 30), jitter(gr$Y[path], factor = 30), col = rainbow(n = 365))

text(gr$X[78], gr$Y[78]-0.1, pos = 1, col = rainbow(n = 365)[1], "winter")
text(gr$X[299]+0.1, gr$Y[299], pos = 4, col = rainbow(n = 365)[which(path == 299)[1]], "calving")
text(gr$X[230]-0.1, gr$Y[230], pos = 2, col = rainbow(n = 365)[tail(which(path == 230), 1)], "summer")

###############################################################################
###############################################################################
# Look at spatial grid
###############################################################################
###############################################################################

# Load list of grid points within Bathurst caribou range
gr <- read.csv("tempData/NARR/gridRange.csv")

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

# Create array to store temps for each DOY from 2000-2020 (21 years)
allTemps <- array(NA, dim = c(21, 365))

#------------------------------------------------------------------------------
# Load NARR sub-surface data
#------------------------------------------------------------------------------

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
			allTemps[y, j] <- z.tsoil[gr$x[path[j]], gr$y[path[j]], 1, DOM] - 273.15
			DOM <- DOM + 1
		}
		if(i < 12) startDay[i+1] <- endDay[i] + 1		
		rm(zz, z.tsoil)

	} # end months i
} # end years y

#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------
colY <- PNWColors::pnw_palette("Bay", n = 21)

plot(DOY, allTemps[1, ], "n", ylim = range(allTemps), ylab = "Temperature (*C)", xlab = "Day of year", las = 1)
abline(h = 0)
for(y in 1:21){
	lines(DOY, allTemps[y, ], col = colY[y])
}

L1 <- lowess(allTemps ~ t(array(rep(DOY, 21), dim = c(365, 21))), f = 1/4)
lines(L1)

# Period smoother: to make ends match, paste together and then cut out middle
x <- rbind(t(array(rep(DOY, 21), dim = c(365, 21))),t(array(rep(DOY+365, 21), dim = c(365, 21))), t(array(rep(DOY+2*365, 21), dim = c(365, 21)))) 
y <- rbind(allTemps, allTemps, allTemps)
L2 <- lowess(y ~ x, f = 1/12)

y.smooth <- matrix(L2$y, ncol = 365*3, nrow = 21)[1, 366:(2*365)]
lines(DOY, y.smooth, lwd = 2) # Matching ends. Good.


# SMoothed to one year, does it look different?
y <- 15
lines(lowess(allTemps[y, ] ~ DOY, f = 1/3), col = colY[y], lwd = 2)
# Sometimes, but on average no.

migTemp <- data.frame(DOY = DOY, temp = y.smooth)
write.csv(migTemp, file = "tempData/NARR/migTemp.csv")
