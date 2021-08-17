DOY <- c(1:365)
# mean temperature
ck <- c(-10.974553,  -8.029248,  -9.808518,  -9.750996, -11.466570)
# half annual temperature range (amplitude)
dk <- c(22.546752,  23.653142,  22.808474,  22.990450, 21.577856)
# DOY corresponding to max temperature
t0 <- c(199.922313, 198.146132, 199.304408, 199.167954, 201.066183) 

j <- 1
tempDOY <- ck[j] + dk[j] * cos((DOY - t0[j])* 2 * pi / 365)


plot(DOY, tempDOY, "l", xaxs = "i")
abline(h = ck[1])
mean(tempDOY)

# Climate scenarios low = RCP2.6 and high = RCP8.5
# Projections are for 2081 - 2100
# Temp increases for Dec (DOY 336) - Feb (DOY 59): + 3*C (low) and + 12*C (high)
# Temp increases for June (153) - Aug (216): + 1.5*C (low) and + 6*C (high)

# Where are these points
strftime(as.Date("2020-08-031"), format = "%j")
u <- par('usr')
polygon(x = c(153, 216, 216, 153), y = u[c(3, 3, 4, 4)], col = "#FF000030", border = NA)
polygon(x = c(1, 59, 59, 1), y = u[c(3, 3, 4, 4)], col = "#0000FF30", border = NA)
polygon(x = c(336, 365, 365, 336), y = u[c(3, 3, 4, 4)], col = "#0000FF30", border = NA)

# smooth increase across the year, and then add?

###############################################################################
# Look at output data from climate model
###############################################################################

# Cliamte data from 
# https://climate-scenarios.canada.ca/?lang=en&page=cccr-data

library(ncdf4) # package for netcdf manipulation

tChange <- data.frame(
	model = rep(c('rcp26', 'rcp85'), each = 4),
	season = rep(c("DJF", "MAM", "JJA", "SON"), 2),
	tas25 = rep(NA, 8),
	tas50 = rep(NA, 8),
	tas75 = rep(NA, 8)
)

for(p in 1:3){ # for 3 percentiles (25, 50, 75)
	for(i in 1:2){ # for two climate model scenarios
	for(j in 1:4){  #for each of four seasons + annual
		nc_data <- nc_open(paste("ClimateChangeData/tas_Amon_", c("rcp26", "rcp85")[i], "_canada_", c("DJF", "MAM", "JJA", "SON", "ANN")[j], "_2080_enspctl", c(25, 50, 75)[p], ".nc", sep = ""))
		
		# # Save the print(nc) dump to a text file
		# {
		# 	sink('tas_Amon_rcp26_canada_ANN_2080_enspctl50_metadata.txt')
		# 	print(nc_data)
		# 	sink()
		# }

		lat <- ncvar_get(nc_data, "lat")
		lon <- ncvar_get(nc_data, "lon")
		tasArray <- ncvar_get(nc_data, "tas")

		# r <- raster(t(tasArray), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
		# r <- flip(r, direction='y')
		# 
		# plot(r, xlab = "Lon", ylab = "Lat")
		# points(360-110.5, 66.5)

		# Centroid of calving grounds
		# 66.5 lat, -110.2 lon
		tChange[tChange$model == c("rcp26", "rcp85")[i] & tChange$season == c("DJF", "MAM", "JJA", "SON")[j], p + 2] <- tasArray[which(lon == 360-110.5), which(lat == 66.5)]
		
		rm(tasArray, nc_data, lat, lon)
	}
	}
}

write.csv(tChange, file = "ClimateChangeData/tasSummaryCalving.csv")

xDate <- c(as.Date("2100-01-15"), as.Date("2100-04-15"), as.Date("2100-07-15"), as.Date("2100-10-15"))
x <- as.numeric(strftime(xDate, format = "%j"))
y_rcp26 <- tChange[1:4, 'tas50']
y_rcp85 <- tChange[5:8, 'tas50']

fit26 <- nls(
	y ~ ck + dk * cos((x - t0)* 2 * pi / 365), 
	data = data.frame(x = x, y = y_rcp26), 
	start = list(ck = 2, dk = 2, t0 = 200))
fit85 <- nls(
	y ~ ck + dk * cos((x - t0)* 2 * pi / 365), 
	data = data.frame(x = x, y = y_rcp85), 
	start = list(ck = 2, dk = 2, t0 = 200))


parChange <- cbind(rcp26 = summary(fit26)$coefficients[, 1], rcp85 = summary(fit85)$coefficients[, 1])
write.csv(parChange, file = "ClimateChangeData/parChange.csv")


xDate <- as.Date(paste("2100", DOY, sep = "-"), format = "%Y-%j")
plot(xDate, rep(1, 365), "n", col = 4, xlab = "", ylab = expression(paste("Temperature change (", degree, "C)", sep = "")), las = 1, ylim = c(-2, 12), xaxs = "i")
u <- par('usr')
polygon(x = c(as.Date("2100-03-01"), as.Date("2100-05-31"), as.Date("2100-05-31"), as.Date("2100-03-01")), y = u[c(3,3,4,4)], col= grey(0.8), border = NA)

polygon(x = c(as.Date("2100-09-01"), as.Date("2100-11-30"), as.Date("2100-11-30"), as.Date("2100-09-01")), y = u[c(3,3,4,4)], col= grey(0.8), border = NA)

abline(h = 0)

# Winter
segments(x0 = as.Date("2100-01-01"), x1 = as.Date("2100-02-28"), y0 = y_rcp26[1], y1 = y_rcp26[1], col = 4, lwd =2)
segments(x0 = as.Date("2100-01-01"), x1 = as.Date("2100-02-28"), y0 = y_rcp85[1], y1 = y_rcp85[1], col = 2, lwd =2)
segments(x0 = as.Date("2100-12-01"), x1 = as.Date("2100-12-31"), y0 = y_rcp26[1], y1 = y_rcp26[1], col = 4, lwd =2)
segments(x0 = as.Date("2100-12-01"), x1 = as.Date("2100-12-31"), y0 = y_rcp85[1], y1 = y_rcp85[1], col = 2, lwd =2)

# Spring
segments(x0 = as.Date("2100-03-01"), x1 = as.Date("2100-05-31"), y0 = y_rcp26[2], y1 = y_rcp26[2], col = 4, lwd =2)
segments(x0 = as.Date("2100-03-01"), x1 = as.Date("2100-05-31"), y0 = y_rcp85[2], y1 = y_rcp85[2], col = 2, lwd =2)

# Summer
segments(x0 = as.Date("2100-06-01"), x1 = as.Date("2100-08-31"), y0 = y_rcp26[3], y1 = y_rcp26[3], col = 4, lwd =2)
segments(x0 = as.Date("2100-06-01"), x1 = as.Date("2100-08-31"), y0 = y_rcp85[3], y1 = y_rcp85[3], col = 2, lwd =2)

# Autumn
segments(x0 = as.Date("2100-09-01"), x1 = as.Date("2100-11-30"), y0 = y_rcp26[2], y1 = y_rcp26[2], col = 4, lwd =2)
segments(x0 = as.Date("2100-09-01"), x1 = as.Date("2100-11-30"), y0 = y_rcp85[2], y1 = y_rcp85[2], col = 2, lwd =2)



lines(as.Date(paste("2100", DOY, sep = "-"), format = "%Y-%j"), summary(fit26)$coefficients['ck', 1] + summary(fit26)$coefficients['dk', 1] * cos((DOY - summary(fit26)$coefficients['t0', 1]) * 2 * pi / 365), col = 4, lty = 2)
lines(as.Date(paste("2100", DOY, sep = "-"), format = "%Y-%j"), summary(fit85)$coefficients['ck', 1] + summary(fit85)$coefficients['dk', 1] * cos((DOY - summary(fit85)$coefficients['t0', 1]) * 2 * pi / 365), col = 2, lty = 2)

# fit85 <- nls(
# 	y ~ ck + dk * cos((x - t0)* 2 * pi / 365), 
# 	data = data.frame(
# 		x = DOY, 
# 		y = c(
# 			rep(y_rcp85[1], 59), 
# 			rep(y_rcp85[2], 92), 
# 			rep(y_rcp85[3], 92), 
# 			rep(y_rcp85[4], 91), 
# 			rep(y_rcp85[1], 31))), 
# 	start = list(ck = -10, dk = 2, t0 = 200))
# 
# 
