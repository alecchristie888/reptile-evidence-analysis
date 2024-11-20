library(data.table)
library(tidyverse)
library(DHARMa)
library(glmmTMB)
library(car)
library(emmeans)
library(ggmap)
library(mapdata)
library(viridis)
library(ggplot2)
library(grid)
library(sf)
library(ggalt)
library(ggthemes)
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(devtools)
library(data.table)
library(RColorBrewer)
library(MuMIn)

setwd("~/reptile-evidence-analysis")
rep.dat <- fread("Reptile-geog-data-revised.csv")

#########Geography GLM
rep.dat$Continent <- factor(rep.dat$continent, levels=c("Asia", "North America", "Oceania", "Europe", "Africa", "South America"))

### geographic model if poisson - results show geog model cannot use poisson distribution
Geogmodel1 <- glm(num.studies ~ gdp.per.cap + Number.of.Species*prop.thr + Continent, family="poisson", data=rep.dat)
testDispersion(Geogmodel1)
simulationOutput <- simulateResiduals(fittedModel = Geogmodel1, plot = F)
residuals(simulationOutput)
plot(simulationOutput)
testZeroInflation(simulationOutput)

### ZERO-INFLATED POISSON REGRESSION suggests that the excess zeros are generated from a separate process from the count values - not sure if this is best option for this data
### Tests for the quasipoisson distribution for geog model - not symmetrical
Geogmodelqp <- glm(num.studies ~ gdp.per.cap + Number.of.Species*prop.thr  + Continent, family="quasipoisson", data=rep.dat)
dev_residuals <- residuals(Geogmodelqp, type = "deviance")
# Plot deviance residuals against predicted values
plot(fitted(Geogmodelqp), dev_residuals, main = "Deviance Residuals vs. Fitted Values")
#patterns in residuals with some large outliers.

### poisson - results show tax model cannot use poisson distribution
scale_0.5SD <- function(x){
  (x - mean(x)) / (2*sd(x))
}

rep.dat$gdp.per.cap_scaled <- scale_0.5SD(rep.dat$gdp.per.cap)
rep.dat$Number.of.Species_scaled <- scale_0.5SD(rep.dat$Number.of.Species)
rep.dat$prop.thr_scaled <- scale_0.5SD(rep.dat$prop.thr)

cor.test(rep.dat$gdp.per.cap_scaled, rep.dat$prop.thr_scaled)
cor.test(rep.dat$Number.of.Species_scaled, rep.dat$prop.thr_scaled)
cor.test(rep.dat$Number.of.Species_scaled, rep.dat$gdp.per.cap_scaled)

Geogmodelnb <- glmmTMB(num.studies ~ gdp.per.cap_scaled + Number.of.Species_scaled*prop.thr_scaled + Continent, family=nbinom2, ziformula = ~0, data=rep.dat)

options(na.action = "na.fail") #Must run this code once to use dredge
all_models_geog <- dredge(Geogmodelnb)
m.top.models.2aic <- get.models(all_models_geog, subset = delta <2)
length(m.top.models.2aic)

#Geogmodelnb2 without continent has lower AIC by >2 and not significant difference with log likelihood test
Geogmodelnb2 <- glmmTMB(num.studies ~ gdp.per.cap_scaled + Number.of.Species_scaled+prop.thr_scaled, family=nbinom2, ziformula = ~0, data=rep.dat)

summary(get.models(all_models_geog, subset = delta <2)[[1]])[6][[1]][1]
confint(get.models(all_models_geog, subset = delta <2)[[1]],full=TRUE)

#write.csv(summary(get.models(all_models_geog, subset = delta <2)[[1]])[6][[1]][1],"geog-modelsummary-update.csv",row.names=FALSE)
#write.csv(confint(get.models(all_models_geog, subset = delta <2)[[1]],full=TRUE),"geog-modelsummary-updateconfint.csv",row.names=FALSE)


################################################## map plotting

#devtools::install_github("eliocamp/ggalt@new-coord-proj",force=TRUE)

### This first section of code downloads a standard world map and makes it look professional

URL <- "http://naciscdn.org/naturalearth/10m/physical/ne_10m_land.zip"
fil <- basename(URL)
if (!file.exists(fil)) download.file(URL, fil)

shp <- grep("shp$", unzip(fil), value=TRUE)

world <- st_read(shp, stringsAsFactors=FALSE)


WorldData <- data.table(map_data('world'))
WorldData$country <- WorldData$region
setdiff(rep.dat$country,WorldData$country)

sort(unique(WorldData$country))

WorldData[country=="Antigua",country:="Antigua and Barbuda"]
WorldData[country=="Bahamas",country:="Bahamas, The"]
WorldData[country=="Turkey",country:="Turkiye"]
WorldData[country=="UK",country:="United Kingdom"]
WorldData[country=="USA",country:="United States"]
WorldData[country=="Brunei",country:="Brunei Darussalam"]
WorldData[country=="Democratic Republic of the Congo",country:="Congo, Dem. Rep."]
WorldData[country=="Republic of Congo",country:="Congo, Rep."]
WorldData[country=="Ivory Coast",country:="Cote d'Ivoire"]
WorldData[country=="Egypt",country:="Egypt, Arab Rep."]
WorldData[country=="Czech Republic",country:="Czechia"]
WorldData[country=="Swaziland",country:="Eswatini"]
WorldData[country=="Gambia",country:="Gambia, The"]
WorldData[country=="Kyrgyzstan",country:="Kyrgyz Republic"]
WorldData[country=="Laos",country:="Lao PDR"]
WorldData[country=="Micronesia",country:="Micronesia, Fed. Sts."]
WorldData[country=="Slovakia",country:="Slovak Republic"]
WorldData[country=="Saint Kitts",country:="St. Kitts and Nevis"]
WorldData[country=="Saint Lucia",country:="St. Lucia"]
WorldData[country=="Saint Martin",country:="St. Martin (French part)"]
WorldData[country=="Saint Vincent",country:="St. Vincent and the Grenadines"]
WorldData[country=="Syria",country:="Syrian Arab Republic"]
WorldData[country=="Trinidad",country:="Trinidad and Tobago"]
WorldData[country=="Virgin Islands",country:="Virgin Islands (U.S.)"]
WorldData[country=="Yemen",country:="Yemen, Rep."]
WorldData[country=="Vietnam",country:="Viet Nam"]
WorldData[country=="Cape Verde",country:="Cabo Verde"]
WorldData[country=="Iran",country:="Iran, Islamic Rep."]
WorldData[country=="Russia",country:="Russian Federation"]
WorldData[country=="Sint Maarten",country:="Sint Maarten (Dutch part)"]
WorldData[country=="Venezuela",country:="Venezuela, RB"]
WorldData[country=="South Korea",country:="Korea, Rep."]

rep.dat1 <- rep.dat

#combine these for the purposes of plotting
rep.dat1[country=="Hong Kong SAR, China",country:="China"]
rep.dat1[country=="Macao SAR, China",country:="China"]

rep.dat2 <- data.table(rep.dat1 %>% 
  group_by(country) %>%
  summarise(num.studies=sum(num.studies)))

#setdiff(rep.dat2$country,WorldData$country)
#setdiff(WorldData$country,rep.dat2$country)

map.df <- merge(rep.dat2, WorldData, by="country")
str(map.df)

map.df2 <- unique(map.df[,list(country,long,lat,group,order,region,num.studies)])

unique(map.df[order(num.studies),list(country,num.studies)])
#setdiff(map.df$country,rep.dat2$country)

cleanup_light <- 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'grey90', colour = 'white'), 
        legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size=25),
        legend.title = element_text(size=30),
        axis.line = element_line(colour = "white"), 
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())

ggplot() +
  geom_map(
    data = map.df2, map = map.df2,
    colour="grey80",
    aes(long, lat, map_id = region,
        fill=log10(num.studies))
  ) +
  scale_fill_viridis_c(name="Number\nof studies",breaks=log10(c(1,3,10,25,100,250)),labels=c(1,3,10,25,100,250), limits=c(0,log10(300)),option="plasma")+
  cleanup_light +
  xlab("") + ylab("")

#ggsave("Reptilemap.png",units="cm",dpi=600,width=50,height=27.5,device="png")

