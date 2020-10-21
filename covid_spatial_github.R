#################
# COVID-19 spatial maps
# Citation: Goldstein ND, Wheeler DC, Gustafson P, Burstyn I. A Bayesian Approach to Improving Spatial Estimates of Prevalence of COVID-19 After Accounting for Misclassification Bias in Surveillance Data in Philadelphia, PA. Manuscript in preparation.
# 4/10/20 -- Neal Goldstein
#################


### FUNCTIONS ###

library("tidycensus") #retrieve ACS data, note if error installing on MacOS see: https://github.com/r-quantities/units/issues/1
library("rgdal") #read shapefile
library("psych") #PCA


### READ DATA ###

#NYC data retreived Apr 18 2020 from: https://github.com/nychealth/coronavirus-data
nyc_zcta = read.csv("NYC/tests-by-zcta.csv", as.is=T, stringsAsFactors=F)

#Philly data retreived Jun 10 2020 from: https://www.opendataphilly.org/dataset/covid-cases
philly_zip = read.csv("Philly/covid_cases_by_zip.csv", as.is=T, stringsAsFactors=F)

#zctas for crosswalk: from UDS mapper https://www.udsmapper.org/zcta-crosswalk.cfm
zcta = read.csv("zip_to_zcta_2019.csv", as.is=T, stringsAsFactors=F)

#retrieve census variables of interest: using tidycensus but could manually obtain from FactFinder
#2018 is most recent year data are available per package
census_api_key("paste api key here")
population = get_acs(geography="zcta", table="B01003", year=2018, output="wide")

#census-based deprivation index: details provided here https://towardsdatascience.com/a-census-based-deprivation-index-using-r-7aa738da697c, https://www.ncbi.nlm.nih.gov/pubmed/17031568
deprivation = get_acs(geography="zcta", variables=c("B17001_002", "B17001_001", "B06009_002" , "B06009_001","B09008_011","B09008_001","B08124_002", "B08124_001", "B25014_005","B25014_006",  "B25014_007","B25014_011", "B25014_012", "B25014_013","B25014_001", "B19058_002", "B19058_001","C23002C_021", "C23002D_008","C23002C_017", "C23002D_003","B19001_002", "B19001_003", "B19001_004","B19001_005", "B19001_006", "B19001_001"), output="wide", year=2018)

#ZCTA shapefile national data from Census: https://www.census.gov/cgi-bin/geo/shapefiles/index.php?year=2019&layergroup=ZIP+Code+Tabulation+Areas
us_zcta = readOGR("tl_2019_us_zcta510/", "tl_2019_us_zcta510")


### CLEAN and JOIN DATA ###

#NEW YORK

#join population to NYC zcta
nyc_data = merge(nyc_zcta,population,by.x="MODZCTA",by.y="GEOID")

#recode and clean
nyc_data$Negative = nyc_data$Total - nyc_data$Positive
nyc_data$Population = nyc_data$B01003_001E
nyc_data$ZCTA = nyc_data$MODZCTA
nyc_data$zcta_cum.perc_pos = NULL
nyc_data$NAME = NULL
nyc_data$B01003_001E = NULL
nyc_data$B01003_001M = NULL
nyc_data$MODZCTA = NULL

#PHILLY

#crosswalk Philly zip to zctas
philly_zip$zcta = NA
for (i in 1:nrow(philly_zip)) {
  z = zcta$ZCTA[which(zcta$ZIP_CODE==as.numeric(philly_zip$zip_code[i]))]
  philly_zip$zcta[i] = ifelse(length(z)>0, z, NA)
}
rm(i,z,zcta)

philly_data = data.frame("Positive"=NA, "Total"=NA, "Negative"=NA, "Population"=NA, "ZCTA"=unique(philly_zip$zcta), stringsAsFactors=F)

for (i in 1:nrow(philly_data)) {
  philly_data$Positive[i] = sum(philly_zip$count[philly_zip$zcta==philly_data$ZCTA[i] & philly_zip$covid_status=="POS"], na.rm=T)
  philly_data$Negative[i] = sum(philly_zip$count[philly_zip$zcta==philly_data$ZCTA[i] & philly_zip$covid_status=="NEG"], na.rm=T)
}
rm(i)

#join population to Philly zcta
philly_data = merge(philly_data,population,by.x="ZCTA",by.y="GEOID")

#recode and clean
philly_data$Total = philly_data$Positive + philly_data$Negative
philly_data$Population = philly_data$B01003_001E
philly_data$NAME = NULL
philly_data$B01003_001E = NULL
philly_data$B01003_001M = NULL

rm(nyc_zcta, philly_zip, population)


### PREDICTORS of SURVEILLANCE ###

#create deprivation index: https://towardsdatascience.com/a-census-based-deprivation-index-using-r-7aa738da697c, https://www.ncbi.nlm.nih.gov/pubmed/17031568
deprivation$pct_poverty = deprivation$B17001_002E / deprivation$B17001_001E
deprivation$pct_noHS = deprivation$B06009_002E / deprivation$B06009_001E
deprivation$pct_FHH = deprivation$B09008_011E / deprivation$B09008_001E
deprivation$pct_mgmt = deprivation$B08124_002E / deprivation$B08124_001E 
deprivation$pct_crowd = (deprivation$B25014_005E + deprivation$B25014_006E + deprivation$B25014_007E + deprivation$B25014_011E + deprivation$B25014_012E + deprivation$B25014_013E) / deprivation$B25014_001E
deprivation$pct_pubassist = deprivation$B19058_002E / deprivation$B19058_001E
deprivation$pct_unempl = (deprivation$C23002C_021E + deprivation$C23002D_008E) / (deprivation$C23002C_017E + deprivation$C23002D_003E)
deprivation$pct_under30K = ((deprivation$B19001_002E + deprivation$B19001_003E + deprivation$B19001_004E + deprivation$B19001_005E + deprivation$B19001_006E) / deprivation$B19001_001E)
deprivation_matrix = as.matrix(deprivation[, c("pct_poverty","pct_noHS","pct_FHH","pct_mgmt","pct_crowd","pct_pubassist", "pct_unempl","pct_under30K")])
deprivation_matrix[is.nan(deprivation_matrix)] = 0
deprivation$census_NDI = principal(deprivation_matrix,nfactors = 1)$scores 

philly_data = merge(philly_data, deprivation[,c("GEOID","census_NDI")], by.x="ZCTA", by.y="GEOID", all.x=T, all.y=F)
nyc_data = merge(nyc_data, deprivation[,c("GEOID","census_NDI")], by.x="ZCTA", by.y="GEOID", all.x=T, all.y=F)

rm(deprivation, deprivation_matrix)


### CREATE ZCTA SHAPEFILE DATA for LOCAL JURISDICTIONS ###

#subset shapefile based on covid data
nyc_sf = us_zcta[us_zcta$ZCTA5CE10 %in% unique(nyc_data$ZCTA), ]
philly_sf = us_zcta[us_zcta$ZCTA5CE10 %in% unique(philly_data$ZCTA), ]
rm(us_zcta)


### SAVE DATA ###

save.image("covid_spatial.RData")


### FUNCTIONS ###

library("MASS") #negative binomial
library("sp") #shapefile
library("plotrix") #plotting functions
#library("ggmap") #plotting functions
library("cartogram") #mapping cartograms 
library("sf") #spatial data
library("spdep") #spatial modeling
library("RColorBrewer") #color palette
library("maptools") #shapefile manipulations
library("raster") #area function


### READ DATA ###

load("covid_spatial.RData")


### SET SPATIAL PARAMETERS ###

#set CRS; https://epsg.io/3652
philly_sf_proj = st_transform(st_as_sf(philly_sf), crs=3652)

#merge case data with shapefile
philly_sf_joined = merge(x=philly_sf_proj, y=philly_data, by.x="ZCTA5CE10", by.y="ZCTA", all.x=T, duplicateGeoms=T)

#create a population density variable
philly_sf_joined$Density = philly_sf_joined$Population / area(philly_sf) 


# ### EXPLORATORY PREDICTIONS ###
# 
# #Total = number tested
# #Total_notest = number not tested
# #Positive_true = true number infected among those test
# #Positive = observed number infected among those tested
# #Positive_missed = measurement error due to imperfect diagnostic test
# #Negative = observed number not infected among those tested
# #Infect_notest = measurement error due to lack of testing
# #Infect_missed = total measurement error as function of imperfect test and lack of testing
# #Prev_zcta = True community prevalance
# 
# #estimates of prev(disease|asymptomatic) as of Apr 22 2020
# #COVID-19 Antibody Seroprevalence in Santa Clara County, California: https://doi.org/10.1101/2020.04.14.20062463
# #Estimating the number of SARS-CoV-2 infections in the United States: https://doi.org/10.1101/2020.04.13.20064519
# prev_disase_nottest = 0.025
# 
# #estimate of PCR accuracy
# SN = 0.75
# SP = 0.98
# 
# #hypothetical
# philly_data$Total_notest = philly_data$Population - philly_data$Total
# philly_data$Positive_true = round((((philly_data$Positive/philly_data$Total) - 1 + SP) / (SN + SP - 1)) * philly_data$Total)
# philly_data$Positive_missed = philly_data$Positive_true - philly_data$Positive
# philly_data$Infect_notest = round(prev_disase_nottest * philly_data$Total_notest)
# philly_data$Infect_missed = philly_data$Positive_missed + philly_data$Infect_notest
# philly_data$Prev_zcta = (philly_data$Infect_missed + philly_data$Positive_true)/philly_data$Population
# 
# #per 10k capita adjustment
# philly_data$Infect_missed_percap = round(philly_data$Infect_missed / philly_data$Population * 10000)
# 
# #misclassificatiom burden: aiming for 50 per 100k for 2 weeks from https://www.dailyitem.com/news/wolf-state-will-start-easing-restrictions-in-northcentral-northwestern-pa/article_044673c9-81ea-5ae6-8960-9634d97b687b.amp.html
# #philly_data$Observed = round((philly_data$Positive)/philly_data$Population * 100000) / 6 * 2
# #philly_data$True = round((philly_data$Infect_missed + philly_data$Positive_true)/philly_data$Population * 100000) / 6 * 2


### DEPRIVATION and TESTING ###

par(mfrow=c(3,1))

plot(philly_data$census_NDI, (philly_data$Positive/philly_data$Population*10000), xlab="ZIP deprivation index", ylab="Positive tests per 10,000 population")
abline(lm((philly_data$Positive/philly_data$Population*10000) ~ philly_data$census_NDI), col="red")

plot(philly_data$census_NDI, (philly_data$Positive/philly_data$Total*10000), xlab="ZIP deprivation index", ylab="Positive tests per 10,000 tests")
abline(lm((philly_data$Positive/philly_data$Total*10000) ~ philly_data$census_NDI), col="red")

plot(philly_data$census_NDI, (philly_data$Total/philly_data$Population*10000), xlab="ZIP deprivation index", ylab="Total tests per 10,000 population")
abline(lm((philly_data$Total/philly_data$Population*10000) ~ philly_data$census_NDI), col="red")

#model of observed prevalence by ADI
summary(glm.nb((Positive/Population) ~ census_NDI, data=philly_data))
summary(glm(Positive ~ census_NDI, offset=log(Population), family=quasipoisson(), data=philly_data))

#expected change in prevalance by unit change in ADI
x=exp(rnorm(100000,-3.5,0.5)) - (exp(rnorm(100000,-3.5,0.5) + rnorm(100000,0.04,0.04)))
mean(x*100)
quantile(x*100,probs=c(0.025,0.975))
rm(x)


### EXPLORATORY MAP ###

par(mar=rep(0.1,4))
plot(philly_sf_joined$geometry)

#create a scalar between 0.01 and 1.01 based on missed cases per capita
case_scale = (philly_data$Infect_missed_percap - min(philly_data$Infect_missed_percap)) / (max(philly_data$Infect_missed_percap) - min(philly_data$Infect_missed_percap)) + 0.01

#place a circle sized relative to the number of missed infections in that zip code
centroids = st_centroid(philly_sf_joined$geometry)

for (i in 1:length(centroids)) {
  floating.pie(xpos=centroids[[i]][1], ypos=centroids[[i]][2], x=c((philly_data$Positive_missed[i] / philly_data$Infect_missed[i]), (1-(philly_data$Positive_missed[i] / philly_data$Infect_missed[i]))), col=c("red","blue"), radius=2500*case_scale[i])
  #draw.circle(x=centroids[i,1], y=centroids[i,2], radius=0.01*case_scale[i])
}
rm(i)


### EXPLORATORY CARTOGRAM ###

#cartogram distorted by population
carto = cartogram_cont(philly_sf_joined, "Population", itermax=4)

par(mar=rep(0.1,4))
plot(carto$geometry)

#create a scalar between 0.01 and 1.01 based on missed cases
case_scale = (philly_data$Infect_missed - min(philly_data$Infect_missed)) / (max(philly_data$Infect_missed) - min(philly_data$Infect_missed)) + 0.01

#place a circle sized relative to the number of missed infections in that zip code
centroids = st_centroid(carto$geometry)

for (i in 1:length(centroids)) {
  floating.pie(xpos=centroids[[i]][1], ypos=centroids[[i]][2], x=c((philly_data$Positive_missed[i] / philly_data$Infect_missed[i]), (1-(philly_data$Positive_missed[i] / philly_data$Infect_missed[i]))), col=c("red","blue"), radius=2500*case_scale[i])
}
rm(i)


### EXPLORATORY SPATIAL ANALYSIS ###

#define neighbors based on shared boundary
sa.nb = poly2nb(philly_sf_joined, queen=T)

#create a representation of the binary weights
sa.wt = nb2listw(neighbours=sa.nb, style="B")

#moran's spatial autocorrelation on missed infections per capita
moran.mc(philly_sf_joined$Infect_missed_percap, listw=sa.wt, nsim=1000)

#spatial correlogram, shows autocorrelation by order of neighbor (1=neighbor, 2=neighbor's neighbor, etc)
plot(sp.correlogram(neighbours=sa.nb,var=philly_sf_joined$Infect_missed_percap,order=4,method="I",style="B",zero.policy=T), main="")

#select eigenvectors, see: https://www.ncbi.nlm.nih.gov/pubmed/24571639 
spatial_eigen = SpatialFiltering(Infect_missed_percap ~ 1, data=philly_sf_joined, nb=sa.nb, style="W", ExactEV=T)
ncol(fitted(spatial_eigen)) #4 eigenvectors selected

#merge eigengvector data to mapping data
philly_sf_joined = cbind(philly_sf_joined,fitted(spatial_eigen))

#choropleth plots by 4-levels of shading
spplot(as_Spatial(philly_sf_joined), "vec1", cuts=3, col.regions=brewer.pal(4, "Spectral"))
spplot(as_Spatial(philly_sf_joined), "vec3", cuts=3, col.regions=brewer.pal(4, "Spectral"))
spplot(as_Spatial(philly_sf_joined), "vec5", cuts=3, col.regions=brewer.pal(4, "Spectral"))
spplot(as_Spatial(philly_sf_joined), "vec6", cuts=3, col.regions=brewer.pal(4, "Spectral"))

#create subsets of spatial data by eigenvectors (using eigenvector #1 with manual adjustments to ensure contiguous boundaries)
vec1_zcta1 = c(19116,19154,19115,19114,19111,19152,19136,19149,19135,19120,19124)
vec1_zcta2 = c(19118,19128,19127,19150,19119,19138,19126,19141,19129,19140,19137,19134,19144,19133)
vec1_zcta3 = c(19151,19131,19132,19121,19122,19123,19125,19106)
vec1_zcta4 = c(19139,19143,19142,19153,19145,19148,19146,19147,19104,19130,19103,19102,19107)

#reference map
par(mar=rep(0.1,4))
plot(philly_sf_joined$geometry, border="gray")
plot(unionSpatialPolygons(as_Spatial(philly_sf_joined)[philly_sf_joined$ZCTA5CE10 %in% vec1_zcta1, ], rep(T,length(vec1_zcta1))), add=T, border=brewer.pal(4, "Spectral")[4], lwd=3, col=adjustcolor(brewer.pal(4, "Spectral")[4],0.2))
plot(unionSpatialPolygons(as_Spatial(philly_sf_joined)[philly_sf_joined$ZCTA5CE10 %in% vec1_zcta2, ], rep(T,length(vec1_zcta2))), add=T, border=brewer.pal(4, "Spectral")[3], lwd=3, col=adjustcolor(brewer.pal(4, "Spectral")[3],0.2))
plot(unionSpatialPolygons(as_Spatial(philly_sf_joined)[philly_sf_joined$ZCTA5CE10 %in% vec1_zcta3, ], rep(T,length(vec1_zcta3))), add=T, border=brewer.pal(4, "Spectral")[2], lwd=3, col=adjustcolor(brewer.pal(4, "Spectral")[2],0.2))
plot(unionSpatialPolygons(as_Spatial(philly_sf_joined)[philly_sf_joined$ZCTA5CE10 %in% vec1_zcta4, ], rep(T,length(vec1_zcta4))), add=T, border=brewer.pal(4, "Spectral")[1], lwd=3, col=adjustcolor(brewer.pal(4, "Spectral")[1],0.2))
text(t(sapply(slot(as_Spatial(philly_sf_joined), "polygons"), function(i) slot(i, "labpt"))), cex=0.6, labels=philly_sf_joined$ZCTA5CE10)

#cartogram distorted by population density
carto = cartogram_cont(philly_sf_joined, "Density", itermax=3)

par(mar=rep(0.1,4))
plot(carto$geometry, border="gray")

#draw eigenvector boundaries
plot(unionSpatialPolygons(as_Spatial(carto)[carto$ZCTA5CE10 %in% vec1_zcta1, ], rep(T,length(vec1_zcta1))), add=T, border=brewer.pal(4, "Spectral")[4], lwd=3, col=adjustcolor(brewer.pal(4, "Spectral")[4],0.2))
plot(unionSpatialPolygons(as_Spatial(carto)[carto$ZCTA5CE10 %in% vec1_zcta2, ], rep(T,length(vec1_zcta2))), add=T, border=brewer.pal(4, "Spectral")[3], lwd=3, col=adjustcolor(brewer.pal(4, "Spectral")[3],0.2))
plot(unionSpatialPolygons(as_Spatial(carto)[carto$ZCTA5CE10 %in% vec1_zcta3, ], rep(T,length(vec1_zcta3))), add=T, border=brewer.pal(4, "Spectral")[2], lwd=3, col=adjustcolor(brewer.pal(4, "Spectral")[2],0.2))
plot(unionSpatialPolygons(as_Spatial(carto)[carto$ZCTA5CE10 %in% vec1_zcta4, ], rep(T,length(vec1_zcta4))), add=T, border=brewer.pal(4, "Spectral")[1], lwd=3, col=adjustcolor(brewer.pal(4, "Spectral")[1],0.2))

#place a pie chart sized relative to the number of missed infections in that boundary
#the scalar here is just proportion of the largest #
floating.pie(xpos=getSpPPolygonsLabptSlots(unionSpatialPolygons(as_Spatial(carto)[carto$ZCTA5CE10 %in% vec1_zcta1, ], rep(T,length(vec1_zcta1))))[1], 
             ypos=getSpPPolygonsLabptSlots(unionSpatialPolygons(as_Spatial(carto)[carto$ZCTA5CE10 %in% vec1_zcta1, ], rep(T,length(vec1_zcta1))))[2], 
             x=c(mean(carto$Positive_missed[carto$ZCTA5CE10 %in% vec1_zcta1] / carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta1]), mean(1-(carto$Positive_missed[carto$ZCTA5CE10 %in% vec1_zcta1] / carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta1]))), 
             col=c("yellow","black"), radius=5000*(sum(carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta1])/sum(carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta1])))


floating.pie(xpos=getSpPPolygonsLabptSlots(unionSpatialPolygons(as_Spatial(carto)[carto$ZCTA5CE10 %in% vec1_zcta2, ], rep(T,length(vec1_zcta2))))[1], 
             ypos=getSpPPolygonsLabptSlots(unionSpatialPolygons(as_Spatial(carto)[carto$ZCTA5CE10 %in% vec1_zcta2, ], rep(T,length(vec1_zcta2))))[2], 
             x=c(mean(carto$Positive_missed[carto$ZCTA5CE10 %in% vec1_zcta2] / carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta2]), mean(1-(carto$Positive_missed[carto$ZCTA5CE10 %in% vec1_zcta2] / carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta2]))), 
             col=c("yellow","black"), radius=5000*(sum(carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta2])/sum(carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta1])))

floating.pie(xpos=10000+getSpPPolygonsLabptSlots(unionSpatialPolygons(as_Spatial(carto)[carto$ZCTA5CE10 %in% vec1_zcta3, ], rep(T,length(vec1_zcta3))))[1], 
             ypos=getSpPPolygonsLabptSlots(unionSpatialPolygons(as_Spatial(carto)[carto$ZCTA5CE10 %in% vec1_zcta3, ], rep(T,length(vec1_zcta3))))[2], 
             x=c(mean(carto$Positive_missed[carto$ZCTA5CE10 %in% vec1_zcta3] / carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta3]), mean(1-(carto$Positive_missed[carto$ZCTA5CE10 %in% vec1_zcta3] / carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta3]))), 
             col=c("yellow","black"), radius=5000*(sum(carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta3])/sum(carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta1])))

floating.pie(xpos=getSpPPolygonsLabptSlots(unionSpatialPolygons(as_Spatial(carto)[carto$ZCTA5CE10 %in% vec1_zcta4, ], rep(T,length(vec1_zcta4))))[1], 
             ypos=getSpPPolygonsLabptSlots(unionSpatialPolygons(as_Spatial(carto)[carto$ZCTA5CE10 %in% vec1_zcta4, ], rep(T,length(vec1_zcta4))))[2], 
             x=c(mean(carto$Positive_missed[carto$ZCTA5CE10 %in% vec1_zcta4] / carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta4]), mean(1-(carto$Positive_missed[carto$ZCTA5CE10 %in% vec1_zcta4] / carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta4]))), 
             col=c("yellow","black"), radius=5000*(sum(carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta4])/sum(carto$Infect_missed[carto$ZCTA5CE10 %in% vec1_zcta1])))


### NAIVE ANALYSIS ###

sum(philly_data$Total)
sum(philly_data$Positive)
sum(philly_data$Positive)/sum(philly_data$Total)*100
sum(philly_data$Population)
sum(philly_data$Total)/sum(philly_data$Population)*1000
sum(philly_data$Positive)/sum(philly_data$Population)*100

zcta = unique(philly_data$ZCTA)
zcta_test = NA
zcta_prev = NA
for (i in 1:length(zcta)) {
  zcta_test[i] = sum(philly_data$Total[philly_data$ZCTA==zcta[i]])/sum(philly_data$Population[philly_data$ZCTA==zcta[i]])*1000
  zcta_prev[i] = sum(philly_data$Positive[philly_data$ZCTA==zcta[i]])/sum(philly_data$Population[philly_data$ZCTA==zcta[i]])*100
}
rm(i)

summary(zcta_test)
summary(zcta_prev)

zcta_prev[which(zcta==19120)]

par(mar=c(4, 4, 4, 5) + 0.1)
plot(x=1:47, y=zcta_prev, pch=1, xaxt="n", yaxt="n", ylab="Prevalence (%)", xlab="", ylim=c(0,20))
points(x=1:47, y=zcta_test/10, pch=8)
axis(side=1, at=1:47, labels=zcta, las=2, cex=0.5)
axis(side=2, at=seq(0, 20, by=1))
axis(side=4, at=seq(0, 20, by=1), labels=seq(0, 200, by=10))
mtext("Tests per 1,000", side=4, line=3)
legend("topright",legend=c("Tests","Prevalence"), pch=c(8,1))

#define neighbors based on shared boundary
sa.nb = poly2nb(philly_sf_joined, queen=T)

#create a representation of the binary weights
sa.wt = nb2listw(neighbours=sa.nb, style="B")

#moran's spatial autocorrelation on tests per capita using rate test
EBImoran.mc(n=philly_sf_joined$Total, x=philly_sf_joined$Population, listw=sa.wt, nsim=10000)
EBImoran.mc(n=philly_sf_joined$Positive, x=philly_sf_joined$Population, listw=sa.wt, nsim=10000)

#spatial correlogram, shows autocorrelation by order of neighbor (1=neighbor, 2=neighbor's neighbor, etc)
par(mfrow=c(1,2))
plot(sp.correlogram(neighbours=sa.nb,var=(philly_sf_joined$Total/philly_sf_joined$Population),order=4,method="I",style="B",zero.policy=T), main="")
plot(sp.correlogram(neighbours=sa.nb,var=(philly_sf_joined$Positive/philly_sf_joined$Population),order=4,method="I",style="B",zero.policy=T), main="")

#choropleth map: see https://edzer.github.io/sp/ for helpful tips
philly_sf_joined$Naive_prev = philly_sf_joined$Positive/philly_sf_joined$Population
sp_arrow = list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(2730000,220000), scale = 10000)
sp_scale = list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(2725000,215000), scale = 21120, fill=c("transparent","black"))
sp_scale_text1 = list("sp.text", c(2726000,211000), "0")
sp_scale_text2 = list("sp.text", c(2742000,211000), "4 mi")
spplot(as_Spatial(philly_sf_joined), "Naive_prev", cuts=7, col.regions=brewer.pal(8, "Reds"), sp.layout=list(sp_arrow,sp_scale,sp_scale_text1,sp_scale_text2))
rm(sp_arrow,sp_scale,sp_scale_text1,sp_scale_text2)
#text(t(sapply(slot(as_Spatial(philly_sf_joined), "polygons"), function(i) slot(i, "labpt"))), cex=0.6, labels=philly_sf_joined$ZCTA5CE10)

head(philly_sf_joined[order(philly_sf_joined$Naive_prev), c("ZCTA5CE10","Naive_prev")])
tail(philly_sf_joined[order(philly_sf_joined$Naive_prev), c("ZCTA5CE10","Naive_prev")])
