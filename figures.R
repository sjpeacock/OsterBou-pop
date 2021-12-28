library(gplots)
library(ggsci)


library(pals)
colPal <- pal_locuszoom()#scico(n = 7, palette = "berlin")
cols <- colPal(7)
colPal2 <- pal_locuszoom(alpha = 0.5)#scico(n = 7, palette = "berlin")
cols2 <- colPal2(7)

names(cols) <- c("Calving", "Summer", "Spring", "Fall", "Winter", "purple", "grey")
names(cols2) <- c("Calving", "Summer", "Spring", "Fall", "Winter", "purple", "grey")

###############################################################################
# Annual cycle
###############################################################################
library(REdaS)

# Color for ranges
colRange <- c(winter = 4, spring = 3, calving = 7, summer = 2, fall = 1)

# There are in reality 150 animals in our database. 
# I added two for Summer 2015 just so that the season appears in the boxplot

nv <- 365
angle.inc <- 2 * pi/nv
angles <- rev(seq(0, 2 * pi - angle.inc, by = angle.inc)+pi/2)
xv <- cos(angles)
yv <- sin(angles)
m <- as.numeric(strftime(as.Date(paste("2014", c(1:12), "01", sep="-"), format="%Y-%m-%d"), format="%j"))


winter <- c(1:109)
spring <- c(110:152)
calving <- c(153:167)
postCalving <- c(168:179)
summer <- c(180:249)
fall <- c(250:334)
rut <- c(290:304)
earlyWinter <- c(335:365)

# Breeding date, when animals move up a class, is June 7 (DOY = 158)
breedDOY <- as.numeric(strftime(as.Date("1985-06-07"), "%j"))

# L4 resume development at the start of spring migration
L4startDOY <- as.numeric(strftime(as.Date("1985-04-20"), "%j"))

# jpeg(filename="Annual_cycle.jpg", width=480, height=480, pointsize=12, quality=600)
quartz(width = 3, height = 3.5, pointsize = 8)
par(mar=c(2,1,4,1))
plot(c(-1.5,1.5), c(-1.5, 1.5), "n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
polygon(xv, yv, border = 1)
# points(c(0, -1, 0, 1), c(1, 0, -1, 0), pch=19)
text((cos(angles)*1.2)[m], (sin(angles)*1.2)[m], c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

lines(xv[earlyWinter], yv[earlyWinter], col = colRange[winter], lwd = 7)
lines(xv[winter], yv[winter], col = colRange["winter"], lwd = 7)

lines(xv[spring], yv[spring], col = colRange["spring"], lwd = 3)

lines(xv[calving], yv[calving], col = colRange["calving"], lwd = 7)
lines(xv[postCalving], yv[postCalving], col = colRange["calving"], lwd = 3)

lines(xv[summer], yv[summer], col = colRange["summer"], lwd = 7)

lines(xv[fall], yv[fall], col = colRange["fall"], lwd = 3)


legend(-1.7, 2.2, lwd = c(7, 3, 7, 3, 7, 3), col = colRange[c("winter", "spring", "calving", "calving", "summer", "fall")], legend = c("winter", "spring", "calving", "post-calving", "summer", "fall"), bty="n", title = "Season", xpd = NA, ncol = 3)

points(xv[spring[1]], yv[spring[1]])
text(xv[spring[1]], yv[spring[1]], paste0(strftime(as.Date(paste("2001", spring[1], sep = "-"), format = "%Y-%j"), format = "%b %d"), "   "), srt = rad2deg(angles[spring[1]]), adj = 1)

points(xv[breedDOY], yv[breedDOY])
text(xv[breedDOY], yv[breedDOY], "Calving - Jun 7   ", srt = rad2deg(angles[breedDOY]), adj = 1)

points(xv[summer[1]], yv[summer[1]])
text(xv[summer[1]], yv[summer[1]], paste0(strftime(as.Date(paste("2001", summer[1], sep = "-"), format = "%Y-%j"), format = "%b %d"), "   "), srt = rad2deg(angles[summer[1]]), adj = 1)

points(xv[fall[1]], yv[fall[1]])
text(xv[fall[1]], yv[fall[1]], "   Sept 7", srt = 180+rad2deg(angles[fall[1]]), adj = 0)

points(xv[earlyWinter[1]], yv[earlyWinter[1]])
text(xv[earlyWinter[1]], yv[earlyWinter[1]], "   Dec 1", srt = 180+rad2deg(angles[earlyWinter[1]]), adj = 0)

lines((cos(angles)*0.9)[rut], (sin(angles)*0.9)[rut], lwd = 1)
text(xv[mean(rut)], yv[mean(rut)], paste0("    ", "Rut"), srt = 180+rad2deg(angles[mean(rut)]), adj = 0)

# points((cos(angles)*0.9)[breedDOY - 240 + 365], (sin(angles)*0.9)[breedDOY - 240 + 365], pch = 8)
points(xv[breedDOY - 240 + 365], yv[breedDOY - 240 + 365], pch = 8, cex = 1.5)
text(xv[breedDOY - 240 + 365], yv[breedDOY - 240 + 365], paste0("   ", strftime(as.Date(paste("2001", breedDOY - 240 + 365, sep = "-"), format = "%Y-%j"), format = "%b %d")), srt = 180+rad2deg(angles[breedDOY - 240 + 365]), adj = 0)


legend("bottomleft", pch = 8, pt.cex = 1.5, bty = "n", "Date when parasite burdens influence pregnancy")
###############################################################################
# Conceptual illustration of MTE relationships
###############################################################################

Kelvin <- function(temp){ return(temp + 273.15)}

E <- 1
Eh <- 3.25
El <- 3.25
Tl <- Kelvin(0)
Th <- Kelvin(30)
T0 <- Kelvin(15)
a <- 0.1
k <- 8.62 * 10^-5

Temp <- seq(-10, 40, 0.1)
nT <- length(Temp)

param <- cbind(
	
	constant = rep(a, nT),
	
	A =        a * exp(-E/k * (1/Kelvin(Temp) - 1/T0)),
	
	SSUL =     a * exp(-E/k * (1/Kelvin(Temp) - 1/T0)) * (1 + exp(El/k * (1/Kelvin(Temp) - 1/Tl)) + exp(Eh/k*(-1/Kelvin(Temp) + 1 / Th))),
	
	SSUL.rho = a * exp(-E/k * (1/Kelvin(Temp) - 1/T0)) * (1 + exp(El/k * (1/Kelvin(Temp) - 1/Tl)) + exp(Eh/k*(-1/Kelvin(Temp) + 1 / Th)))^(-1)

	)

quartz(width = 6.3, height = 3, pointsize = 10)
layout(matrix(c(1,1,2), nrow = 1))
 par(mar = c(4,3,1,1))
plot(Temp, param[, 1], "l", ylim = c(0, 1), col= cols[1], las = 1, bty="l", yaxt="n", xaxt="n", xlab = "", ylab = "", lwd = 1.5)
u <- par('usr')
abline(v = c(Tl, Th)-273.15, col = grey(0.8), lwd = 2)

arrows(x0 = u[1], x1 = u[2], y0 = u[3], y1 = u[3], length = 0.1, xpd = NA)
arrows(x0 = u[1], x1 = u[1], y0 = u[3], y1 = u[4], length = 0.1, xpd = NA)
lines(Temp, param[, 2], col = cols[2], lwd = 1.5)
lines(Temp, param[, 3], col = cols[4], lty = 3, lwd = 1.5)
lines(Temp, param[, 4], col = cols[3], lty = 2, lwd = 1.5)
axis(side = 1, at = c(Tl, Th)-273.15, labels = c(expression(paste(italic(T[L]))), expression(paste(italic(T[H])))))
mtext(side = 1, line = 3, "Temperature")
mtext(side = 2, line = 1, "Rate (y)")

segments(x0=45, x1=50, y0=u[3] + 0.8*(u[4]-u[3]), y1=u[3] + 0.8*(u[4]-u[3]), xpd=NA, col = cols[1], lwd = 1.5)
text(50, u[3] + 0.8*(u[4]-u[3]), pos = 4,"Constant rate", xpd = NA, cex = 1.2, col = cols[1])

segments(x0=45, x1=50, y0=u[3] + 0.7*(u[4]-u[3]), y1=u[3] + 0.7*(u[4]-u[3]), col = cols[2], xpd=NA, lwd = 1.5)
text(50, u[3] + 0.7*(u[4]-u[3]), pos = 4,"Arrhenius", xpd = NA, cex = 1.2, col = cols[2])
# text(45, u[3] + 0.6*(u[4]-u[3]), expression(y = y[0]*exp(-E/k*(1/T - 1/T[0]))), xpd = NA)

segments(x0=45, x1=50, y0=u[3] + 0.45*(u[4]-u[3]), y1=u[3] + 0.45*(u[4]-u[3]), col = cols[3], lty = 2, xpd=NA, lwd = 1.5)
text(50, u[3] + 0.45*(u[4]-u[3]), pos = 4,"Sharpe-Schoolfield\n(development)", xpd = NA, cex = 1.2, col = cols[3])

segments(x0=45, x1=50, y0=u[3] + 0.2*(u[4]-u[3]), y1=u[3] + 0.2*(u[4]-u[3]), col = cols[4], lty = 3, xpd=NA, lwd = 1.5)
text(50, u[3] + 0.2*(u[4]-u[3]), pos = 4,"Inverse Sharpe-Schoolfield\n(mortality)", xpd = NA, cex = 1.2, col = cols[4])

# # CHECK EQUATIONS
# T <- seq(-10, 50, 0.1)
# rates <- cbind(
# 	Molar2013 = a * exp(-E/k * (1/Kelvin(T) - 1/T0)) * (1 + exp(El/k * (1/Kelvin(T) - 1/Tl)) + exp(Eh/k*(-1/Kelvin(T) + 1 / Th)))^(-1),
# 	Molnar2017 = a * exp(-E/k * (1/Kelvin(T) - 1/T0)) * (1 + exp(El/k * (1/Tl - 1/Kelvin(T))) + exp(Eh/k*(-1/Kelvin(T) + 1 / Th)) ^ (-1))
# )
# 
# plot(T, rates[, 2], "l", col = 2, las = 1, ylab = "r(T)", ylim = c(0, 200))
# lines(T, rates[, 1], col = 4)
# abline(v = c(Tl, Th)-273.15, col = grey(0.8), lwd = 2)
# legend("bottomright", lwd = 1, col = c(2,4), c("Molnar2013", "Schoolfield1981/Molnar2017"))
# 
###############################################################################
# Inhibition scenarios
###############################################################################
z <- as.Date(paste("1985", c(1:365), sep = "-"), format = "%Y-%j")
DOY <- c(1:365)

pI4 <- c(rep(1, 109), 1/(1 + exp(-0.08*(c(110:365) - 172))))

quartz(width = 3.2, height = 2, pointsize = 10)
par(mfrow = c(1,1), mar = c(3,4,1,1))

plot(z, DOY, "n", ylim = c(0,1), ylab  = "Proportion arresting", las = 1, xlab = "", xaxs = "i", bty = "l")
lines(z, pI4, lty = 4, lwd = 1.5)
abline(h = 1, lty = 1)
abline(h = 0.5, lty = 2)
abline(h = 0, lty = 3)
legend(as.Date("1985-08-15"), 0.55, lwd = c(1,1,1,1.5), lty = c(1:4), c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4"), xpd = NA, bg = "white", cex = 0.8)

###############################################################################
# Temperature data vs model
###############################################################################

allTemps <- readRDS("tempData/migRouteOutput/allTemps_migPath6.rds")
narr21 <- read.csv("tempData/groundTemps_migPath6_330km.csv")
ccAnomalies <- read.csv("tempData/ccAnomalies_migPath6_330km.csv")

x.ind <- 330
xDate <- as.Date(paste(2001, c(1:365), sep = "-"), format = "%Y-%j")
quartz(width = 3.2, height = 6, pointsize = 10)
par(mfrow = c(3,1), mar = c(2,5,1,1), oma = c(2,0,0,0))

# Ground temps current
plot(xDate, narr21[1,], "n", ylim = c(-30, 30), xaxs = "i", ylab = "", xlab = "", las = 1)
for(i in 1:21) lines(xDate, narr21[i,], col = grey(0.8), lwd = 0.8)
abline(h = 0,lty = 2)
lines(xDate, allTemps[[1]][, x.ind], lwd = 1.5)
mtext(side = 3, line = -1.5, "  A", adj = 0)

# Anomalies
plot(xDate, ccAnomalies[, 1], "n", ylim = c(0,15), xaxs = "i", ylab = "", xlab = "", las = 1)
for(r in 1:2){
	lines(xDate, ccAnomalies[, r], col = c(4,2)[r], lwd = 0.8)
	lines(xDate, ccAnomalies[, r+2], col = c(4,2)[r], lwd = 1.5)
}
mtext(side = 3, line = -1.5, "  B", adj = 0)

# Together
plot(xDate, allTemps[[1]][, x.ind], "n", ylim = c(-30, 30), xaxs = "i", ylab = "", xlab = "", las = 1)
abline(h = 0,lty = 2)

for(r in 1:3){
	lines(xDate, allTemps[[r]][, x.ind], col = c(1,4,2)[r], lwd = 1.5)
}
mtext(side = 3, line = -1.5, "  C", adj = 0)

mtext(side = 2, outer = TRUE, expression(paste("Temperature (", degree, "C)")), line = -2)
###############################################################################
# Larvae parameters and MTE relationships
###############################################################################

par1 <- read.csv("MTE_dataAnalysis/output/bestModelParameters.csv")

# Function to convert character numbers from .csv above to numeric for plotting
convNum <- function(x, type = "m"){
	x.out <- length(x)
	for(i in 1:length(x)){
		dum <- strsplit(x[i], split = "()")[[1]]
		if(type=="m"){
			x.out[i] <- as.numeric(paste(as.numeric(dum[1]), ".", as.numeric(dum[3]), as.numeric(dum[4]), as.numeric(dum[5]), sep = ""))
		} else if(type == "l"){
			x.out[i] <- as.numeric(paste(as.numeric(dum[8]), ".", as.numeric(dum[10]), as.numeric(dum[11]), as.numeric(dum[12]), sep = ""))
		} else if(type == "u"){
			x.out[i] <- as.numeric(paste(as.numeric(dum[15]), ".", as.numeric(dum[17]), as.numeric(dum[18]), as.numeric(dum[19]), sep = ""))
		}
	} # end i
	
	return(x.out)
}

# Function to predict param over temp from MTE hyperparameters
predict.MTE <- function(
	params = c(a = 0.068, E = 0.884, Eh = NA, Th = NA, El = NA, Tl = NA, z = 1),
	temp = 15){
	# Arrhenius	
	if(is.na(params['Eh']) == TRUE & is.na(params['El']) == TRUE){
		return(params['a'] * exp(-params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))))
		
	}else{
		# SS upper bound
		if(is.na(params['El']) == TRUE){
			return(params['a'] * exp(-params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['Eh'] / (8.62*10^-5)*(-1/(temp+273.15)+1/(params['Th']+273.15))))^params['z'])
		}
		# SS lower bound
		if(is.na(params['Eh']) == TRUE){
			return(params['a'] * exp(-params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['El'] / (8.62*10^-5)*(1/(temp+273.15)-1/(params['Tl']+273.15))))^params['z'])
		}
		# SS upper AND lower bound
		if(sum(is.na(params)) == 0){
			return(params['a'] * exp(-params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['El'] / (8.62*10^-5)*(1/(temp+273.15)-1/(params['Tl']+273.15))) + exp(params['Eh'] / (8.62*10^-5)*(-1/(temp+273.15)+1/(params['Th']+273.15))))^params['z'])
		}
	}
} # end function

# Read in temperature data including CC scenarios
# Can choose different temperature profiles depending on migration route
tempDat <- readRDS(paste0("tempData/migRouteOutput/allTemps_migPath", 6, ".rds"))
tempRanges <- rbind(range(tempDat[[1]][120:300, 330]), range(tempDat[[2]][120:300, 330]), range(tempDat[[3]][120:300, 330]))

source("simulations/popFunctions.R")
temp <- seq(5, 35, 5)
temp.all <- seq(-30, 40, 0.2)
quartz(width = 8.3, height = 2.8, pointsize = 10)
par(mfrow = c(1,3), mar = c(3,3,2,1), oma = c(1, 2, 0, 8))

#------------------------------------------------------------------------------
#mu0
#------------------------------------------------------------------------------
plot(temp.all, rep(0, length(temp.all)), "n", las = 1, ylab = "", xlab = "", ylim = c(0, 1), xaxs = "i")

for(r in 1:3){
	# temptemp <- seq(tempRanges[r, 1], tempRanges[r,2], length.out = 50)
	# lines(temptemp,  predict.mu0(temptemp), col = c("#00000060", "#4995E090", "#CF5D6D90")[r], lwd = c(8,6,4)[r])
	temptemp <- seq(tempRanges[r, 1], tempRanges[r,2], 0.2)
	if(r == 1){
		polygon(x = c(rep(tempRanges[r,1], 2), rep(-50, 2), temptemp, tempRanges[r,2]), 
					y = c(-1, rep(min(predict.mu0(temptemp)), 2), max(predict.mu0(temptemp)), predict.mu0(temptemp), -1), col = c("#00000020", "#4995E040", "#CF5D6D60")[r], border = NA)
	} else {
		polygon(x = rep(c(tempRanges[r,1], -50, tempRanges[r,2]), each = 2), y = c(-1, rep(predict.mu0(tempRanges[r,1]), 2), rep(predict.mu0(tempRanges[r,2]), 2), -1), col = c("#00000020", "#4995E040", "#CF5D6D60")[r], border = NA)
	}
}

lines(temp.all, predict.MTE(params = c(a = convNum(par1[9,3]), E = convNum(par1[10,3]), Eh = NA, Th = NA, z = 1), temp = temp.all), lty = 3)
lines(temp.all, predict.mu0(temp.all))

plotCI(temp, convNum(par1[1:7, 'mu0'], type = "m"), li = convNum(par1[1:7, 'mu0'], type = "l"), ui = convNum(par1[1:7, 'mu0'], type = "u"), gap = 0,  sfrac= 0.008, add = TRUE, pch = 21, pt.bg = "white", cex = 1.5)

# Points from Aleuy et al. (2020) Marshallagia study
# freezing survival of L1
points(-9, 0.0855, pch = 8, xpd = NA)
points(-20, 0.9857, pch = 8, xpd = NA)
mtext(side = 3, line = -1.5, paste(" ", LETTERS[1]), adj = 0)

#------------------------------------------------------------------------------
# rho
#------------------------------------------------------------------------------

plot(temp.all, rep(0, length(temp.all)), "n", las = 1, ylab = "", xlab = "", ylim = c(0, 0.1), xaxs = "i")

# Polygon
for(r in 1:3){
	# temptemp <- seq(tempRanges[r, 1], tempRanges[r,2], length.out = 50)
	# lines(temptemp,  predict.rho0(temptemp), col = c("#00000060", "#4995E090", "#CF5D6D90")[r], lwd = c(8,6,4)[r])
	polygon(x = rep(c(tempRanges[r,1], -50, tempRanges[r,2]), each = 2), y = c(-1, rep(predict.rho0(tempRanges[r,1]), 2), rep(predict.rho0(tempRanges[r,2]), 2), -1), col = c("#00000020", "#4995E040", "#CF5D6D60")[r], border = NA)
}

lines(temp.all, predict.MTE(params = c(a = convNum(par1[9, 'rho']), E = convNum(par1[10, 'rho']), Eh = convNum(par1[11, 'rho']), Th = 30.568, z = -1), temp = temp.all), lty = 3, lwd = 1.5)
lines(temp.all, predict.rho0(temp.all))

plotCI(temp[1:7], convNum(par1[1:7, 'rho'], type = "m"), li = convNum(par1[1:7, 'rho'], type = "l"), ui = convNum(par1[1:7, 'rho'], type = "u"), gap = 0,  sfrac= 0.008, add = TRUE, pch = 21, pt.bg = "white", cex = 1.5)

mtext(side = 3, line = -1.5, paste(" ", LETTERS[2]), adj = 0)
#------------------------------------------------------------------------------
# mu3
#------------------------------------------------------------------------------

plot(temp.all, rep(0, length(temp.all)), "n", las = 1, ylab = "", xlab = "", ylim = c(0, 0.1), xaxs = "i")
for(r in 1:3){
	# temptemp <- seq(tempRanges[r, 1], tempRanges[r,2], length.out = 50)
	# lines(temptemp,  predict.mu3(temptemp), col = c("#00000060", "#4995E090", "#CF5D6D90")[r], lwd = c(8, 6, 4)[r])
	polygon(x = rep(c(tempRanges[r,1], -50, tempRanges[r,2]), each = 2), y = c(-1, rep(predict.mu3(tempRanges[r,1]), 2), rep(predict.mu3(tempRanges[r,2]), 2), -1), col = c("#00000020", "#4995E040", "#CF5D6D60")[r], border = NA)
}
lines(temp.all, rep(convNum(par1[9, 'mu1']), length(temp.all)), lty = 3)
lines(temp.all, predict.mu3(temp.all))

plotCI(temp[1:6], convNum(par1[1:6, 'mu1'], type = "m"), li = convNum(par1[1:6, 'mu1'], type = "l"), ui = convNum(par1[1:6, 'mu1'], type = "u"), gap = 0,  sfrac= 0.008, add = TRUE, pch = 21, pt.bg = "white", cex = 1.5)

# Points from Aleuy et al. (2020) Marshallagia study
# freezing survival of L3
points(-9, 0.002344, pch = 8, xpd = NA)
points(-20, 0.01762, pch = 8, xpd = NA)
# points(-35, 3, pch = 8, xpd = NA)
mtext(side = 3, line = -1.5, paste(" ", LETTERS[3]), adj = 0)


#------------------------------------------------------------------------------
legend(42, 0.1, pch = c(21, 8, NA, NA), lty = c(NA, NA, 3, 1), lwd = c(NA, NA, 1.5, 1), legend = c("Exp't estimates", "Marshallagia", "Fitted MTE", "Assumed MTE"), xpd = NA, bty = "n")

# legend(42, 0.06, lty = 1, lwd = c(8,6,4), col = c("#00000060", "#4995E090", "#CF5D6D90"), legend = c("current", "RCP 2.6", "RCP 8.5"), xpd = NA, bty = "n", title = "Summer\ntemperatures")

legend(46, 0.06, fill = c("#00000020", "#4995E040", "#CF5D6D60"), border = NA, legend = c("current", "RCP 2.6", "RCP 8.5"), xpd = NA, bty = "n", title = "Summer\ntemperatures")


mtext(side = 1, outer = TRUE, expression(paste("Temperature (", degree, "C)")))
mtext(side = 2, outer = TRUE, "Parameter value")

