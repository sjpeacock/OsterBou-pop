library(gplots)
library(ggsci)
library(pals)
colPal <- pal_locuszoom()#scico(n = 7, palette = "berlin")
cols <- colPal(7)
colPal2 <- pal_locuszoom(alpha = 0.5)#scico(n = 7, palette = "berlin")
cols2 <- colPal2(7)


# From Stien et al. 2002, relationship between mean Parasite burdens in the previous October and calf production:

# pCalf = pCalf0 * (1 - 1/(1 + exp(7.025 - 0.000328 * parasites))))

# In Stien et al. (2002), they report the "difference" which is the proportional difference. I.e., 
# (pCalf - pCalf0)/pCalf0 = 1/(1 + exp(7.025 - 0.000328 * parasites))))
AlbonData <- read.csv("parameterization/Albon2002_1a.csv",header = TRUE)
OG <- seq(0, 16000, 50)

# But if this is on an individual level, how can we adjust for the distribution of parasites among
# hosts, which is assumed to be binomial?

# Plot to compare stochastic and deterministic approaches
P_mean.all <- seq(0, 16000, 100)
calves.all <- matrix(NA, ncol = 3, nrow = length(P_mean.all))
calves.all[, 1] <- 0.8 * (1  - 1/(1 + exp(7.025 - 0.000328 * P_mean.all)))
for(i in 1:length(P_mean.all)){
	numParasitesPerFemale <- rnbinom(n = 10^6, size = 0.9940684, mu = P_mean.all[i])
	calf <- rbinom(
		n = 10^6,
		size = 1,
		prob = 0.8 * (1 - 1/(1 + exp(7.025 - 0.000328 * numParasitesPerFemale))))
	
	calves.all[i, 2] <- sum(calf) * 10^-6
}

for(i in 1:length(P_mean.all)){
	numParasitesPerFemale <- rnbinom(n = 10^3, size = 0.9940684, mu = P_mean.all[i])
	calf <- rbinom(
		n = 10^3,
		size = 1,
		prob = 0.8 * (1 - 1/(1 + exp(7.025 - 0.000328 * numParasitesPerFemale))))
	
	calves.all[i, 3] <- sum(calf) * 10^-3
}


#------------------------------------------------------------------------------
# Plot showing Albon data and neg binomial difference
#------------------------------------------------------------------------------
quartz(width = 4, height = 4.5, pointsize = 10)
par(mfrow = c(2,1), mar = c(1,6,2,2), oma = c(4, 0, 0, 0))

# Albon data
plot(AlbonData$intensity, AlbonData$calfDiff, xlim=c(0, 16000), ylim = c(0, 0.2), las = 1, xaxt = "n", ylab = expression((p[0] - p(x))/p[0]), col = cols[1])
axis(side = 1, labels = FALSE)
mtext(side = 3, line = 0.5, "A", adj = 0)
lines(OG, 1/(1 + exp(7.025 - 0.000328 * OG)), col = cols[1], lwd = 1.5)
mtext(side = 2, line = 4.5, "Difference in calving rates")

# This study
plot(P_mean.all, calves.all[, 1], "l", ylim = range(calves.all, na.rm = TRUE), xlab = "", ylab = "", las = 1, xlim=c(0, 16000), lwd = 1.5, col = cols[1])
mtext(side = 2, line = 3.5, "Calving rate [p(x)]")
mtext(side = 1, line = 3, "Mean parasite burden [x]")
abline(h = 0.8, lty = 3)
lines(P_mean.all, calves.all[,2], lwd = 1.5)
lines(P_mean.all, calves.all[, 3], lwd = 0.8)
axis(side = 4, at = 0.8, labels = expression(p[0]), las = 1)
lines(x, 0.8*(1 - paramsT3['a']*x^paramsT3['m']/(1 + paramsT3['b']*x^paramsT3['m'])), col = cols[4], lwd = 1.5)
mtext(side = 3, line = 0.5, "B", adj = 0)

legend("bottomleft", lwd = c(1.5, 1.5, 0.8, 1.5), col = c(cols[1], 1, 1, cols[4]), c("Albon et al. (2002)", expression(paste("NB sim (n = ", 10^6, ")")), expression(paste("NB sim (n = ", 10^3, ")")), "Fit used"), cex = 0.8, bty = "n")
#------------------------------------------------------------------------------
# Fit model to estimate new parameters for calving rate of 10^6 females
#------------------------------------------------------------------------------
x <- seq(0, 100000, 500)
y <- numeric(length(x))
for(i in 1:length(x)){
	numParasitesPerFemale <- rnbinom(n = 10^6, size = 0.9940684, mu = x[i])
	calf <- rbinom(
		n = 10^6,
		size = 1,
		prob = 0.8 * (1 - 1/(1 + exp(7.025 - 0.000328 * numParasitesPerFemale))))
	
	y[i] <- sum(calf) * 10^-6
}

fit <- nls(
	y ~ 0.8 * (1 - 1/(1 + exp(a - b * x))), 
	data = data.frame(y, x),
	start = c(a = 5, b = 0.0003))
params <- coefficients(fit)

par(mfrow = c(2,1), mar = c(1,6,2,2), oma = c(4, 0, 0, 0))

# Difference in calving rates over many many parasites
plot(x,(0.8 - y)/0.8, "l", ylim = c(0, 1), las = 1, ylab = expression((p[0] - p(x))/p[0]))
# abline(h = 0.8, lty = 2)
# lines(x, 0.8 * (1 - 1/(1 + exp( params['a'] -  params['b'] * x))), col = cols[4])
lines(x, 1/(1 + exp(7.025 - 0.000328 * x)), col = cols[1], lwd = 1.5)
lines(x, paramsT3['a']*x^paramsT3['m']/(1 + paramsT3['b']*x^paramsT3['m']), col = cols[4], lwd = 1.5)
mtext(side = 2, line = 4.5, "Difference in calving rates")
mtext(side = 1, line = 3, "Mean parasite burden [x]")

# Different form: Type 3 functional response??
plot(x, (0.8 - y)/0.8, "l")

diff <- (0.8 - y)/0.8
fitT3 <- nls(
	y ~ a*x^m/(1 + b*x^m), 
	data = data.frame(y = diff, x),
	start = c(a = 3e-9, b = 4e-9, m = 2))

paramsT3 <- coefficients(fitT3)

lines(x, paramsT3['a']*x^paramsT3['m']/(1 + paramsT3['b']*x^paramsT3['m']), col = cols[4])
