# Compare monolevel temperature at surface
# to subsurface "soil temperatures"

library(ncdf4) # package for netcdf manipulation
library(raster)

#-----------------------------------------------------------------------------
# Import RCP 2.6 data
#-----------------------------------------------------------------------------

cc26 <- nc_open("tempData/climateChange/monthly/geomet-climate-CMIP5.TT.RCP26.ENS.ABS_PCTL50.nc")
# sink('tempData/climateChange/monthly/geomet-climate-CMIP5.TT.RCP26.ENS.ABS_PCTL50_metadata.txt')
# print(cc26)
# sink()

lat <- ncvar_get(cc26, "lat")
lon <- ncvar_get(cc26, "lon")

tas <- ncvar_get(cc26, "Band1")

dim(tas)# longitude, latitude, time (month)

points(rep(lon, length(lat)), rep(lat, each = length(lon)), pch = 19, cex = 0.3)
abline(h = lat)
abline(v = lon)

points(rep(lon, each = length(lat)), rep(lat, length(lon)), pch = 19, cex = 1.5, col = heat.colors(n = 100)[findInterval(tas, seq(-20, 15, length.out = 100))])

for(i in which(lon > c(-120) & lon < c(-105))){
	for(j in which(lat>60 & lat < 68)){
		# points(lon[i], lat[j])
		text(lon[i], lat[j], round(tas[i,j], 1), cex = 0.8)
}}