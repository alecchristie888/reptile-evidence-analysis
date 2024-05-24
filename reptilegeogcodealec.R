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

setwd("~/reptile-evidence-analysis")
rep.dat <- fread("~/Reptile-geog-data.csv")

#########Geography GLM

### geographic model Poisson - results show geog model cannot use poisson distribution
Geogmodel1 <- glm(num.studies ~ gdp.per.cap + Number.of.Reptile.Species + Continent, family="poisson", data=rep.dat)
testDispersion(Geogmodel1)
simulationOutput <- simulateResiduals(fittedModel = Geogmodel1, plot = F)
residuals(simulationOutput)
plot(simulationOutput)
testZeroInflation(simulationOutput)

### ZERO-INFLATED POISSON REGRESSION suggests that the excess zeros are generated from a separate process from the count values
### Tests for the quasipoisson distribution for geog model - not symmetrical
Geogmodelqp <- glm(num.studies ~ gdp.per.cap + Number.of.Reptile.Species + Continent, family="quasipoisson", data=rep.dat)
dev_residuals <- residuals(Geogmodelqp, type = "deviance")
# Plot deviance residuals against predicted values
plot(fitted(Geogmodelqp), dev_residuals, main = "Deviance Residuals vs. Fitted Values")
#patterns in residuals with some large outliers.

### Using glmmTMB negative binomial for geog - passes Dharma tests
Geogmodelnb <- glmmTMB(num.studies ~ gdp.per.cap + Number.of.Reptile.Species + Continent, family=nbinom2, ziformula = ~0 , data=rep.dat)

Geogmodelnb2 <- glmmTMB(num.studies ~ gdp.per.cap + Number.of.Reptile.Species, family=nbinom2, ziformula = ~0 , data=rep.dat)
Geogmodelnb3 <- glmmTMB(num.studies ~ gdp.per.cap + Continent, family=nbinom2, ziformula = ~0, data=rep.dat)
Geogmodelnb4 <- glmmTMB(num.studies ~ Number.of.Reptile.Species + Continent, family=nbinom2, ziformula = ~0, data=rep.dat)
Geogmodelnb5 <- glmmTMB(num.studies ~ gdp.per.cap, family=nbinom2, ziformula = ~0, data=rep.dat)
Geogmodelnb6 <- glmmTMB(num.studies ~ Number.of.Reptile.Species, family=nbinom2, ziformula = ~0, data=rep.dat)
Geogmodelnb7 <- glmmTMB(num.studies ~ Continent, family=nbinom2, ziformula = ~0, data=rep.dat)

anova(Geogmodelnb, Geogmodelnb2) #Geogmodelnb2 without continent has lower AIC by >2 and not significant difference with log likelihood test
anova(Geogmodelnb, Geogmodelnb3) 
anova(Geogmodelnb, Geogmodelnb4) 
anova(Geogmodelnb, Geogmodelnb5)
anova(Geogmodelnb, Geogmodelnb6) 
anova(Geogmodelnb, Geogmodelnb7) 
## Best model is one that does not include continent - which is Geogmodelnb2

#loglikelihood ratio test to check assumptions of neg.binomial - passes - sig. better than poisson
m3 <- glm(num.studies ~ gdp.per.cap + Continent, family = "poisson", data = rep.dat)
pchisq(2 * (logLik(Geogmodelnb2) - logLik(m3)), df = 1, lower.tail = FALSE)
2 * (logLik(Geogmodelnb2) - logLik(m3))

#Dharma plots look fine
simulationOutput <- simulateResiduals(fittedModel = Geogmodelnb2)
plot(simulationOutput)
testZeroInflation(simulationOutput)

#ANOVA analysis of deviance tests and model summary results
Anova(Geogmodelnb2, type="II")
summary(Geogmodelnb2)

#number of studies by country
rep.dat[rev(order(num.studies)),list(country,num.studies)]

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
WorldData[country=="China",country:="China, P.R."]
WorldData[country=="Netherlands",country:="The Netherlands"]
WorldData[country=="Turkey",country:="Turkiye"]
WorldData[country=="UK",country:="United Kingdom"]
WorldData[country=="USA",country:="United States"]
WorldData[country=="Brunei",country:="Brunei Darussalam"]
WorldData[country=="Democratic Republic of the Congo",country:="Congo, Dem. Rep."]
WorldData[country=="Republic of Congo",country:="Congo, Rep."]
WorldData[country=="Ivory Coast",country:="Cote d'Ivoire"]
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

rep.dat2 <- rep.dat
rep.dat2[country=="Macao SAR, China",country:="China, P.R."]
rep.dat2[country=="Hong Kong SAR, China",country:="China, P.R."]

rep.dat2 <- data.table(rep.dat2 %>% 
  group_by(country) %>%
  summarise(num.studies=sum(num.studies)))

setdiff(rep.dat2$country,WorldData$country)
setdiff(WorldData$country,rep.dat2$country)

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

