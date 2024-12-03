# if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
# BiocManager::install("ggtree")
library(ggtree)
library(dplyr)
library(tidytree)
library(treeio)
library(dplyr)
library(data.table)
library(dplyr)
library(tidyr)
library(glm2)
library(emmeans)
library(car)
library(DHARMa)
library(glmmTMB)
library(MuMIn)
library(ape)
library(cluster)
library(RColorBrewer)
library(ggplot2)


setwd("~/reptile-evidence-analysis")
MasterTax <- fread("Reptile-taxa-data-revised.csv")

## Taxonomy GLM
str(MasterTax)
unique(MasterTax$redlistCategory)
MasterTax$redlistCategory <- factor(MasterTax$redlistCategory, levels=c( "Least Concern","Data Deficient", "Not Evaluated", "Near Threatened","Vulnerable", "Endangered", "Critically Endangered"))
unique(MasterTax$redlistCategory)
unique(MasterTax$orderName)
MasterTax$orderName <- factor(MasterTax$orderName, levels=c("Squamata","Testudines","Crocodilia"))
MasterTax$page_views.x <- as.integer(MasterTax$page_views.x)
MasterTax$endemic <- factor(MasterTax$`Insular/endemic (yes or no)`,levels=c("No","Yes"))
MasterTax$venom  <- factor(MasterTax$`Venomous (yes or no)`,levels=c("No","Yes","Unknown"))
MasterTax$bodymass <- MasterTax$`Maximum body mass (g)`

MasterTax <- na.omit(MasterTax)

cor.test(MasterTax$ED,MasterTax$page_views.x)
cor.test(MasterTax$ED,MasterTax$`Maximum body mass (g)`)
cor.test(MasterTax$page_views.x,MasterTax$`Maximum body mass (g)`)

### taxonomic model if poisson - results show tax model cannot use poisson distribution
Taxmodel1 <- glm(num.studies ~ redlistCategory + ED + page_views.x +orderName + endemic + venom + bodymass, data=MasterTax, family=poisson)
testDispersion(Taxmodel1)
simulationOutput <- simulateResiduals(fittedModel = Taxmodel1, plot = F)
residuals(simulationOutput)
plot(simulationOutput)

### ZERO-INFLATED POISSON REGRESSION suggests that the excess zeros are generated from a separate process from the count values - not sure if this is best option for this data

### Tests for the quasipoisson distribution for Tax model - not symmetrical
Taxmodelqp <- glm(num.studies ~ redlistCategory + ED + page_views.x +orderName + endemic + venom + bodymass, data=MasterTax, family=quasipoisson)
dev_residuals <- residuals(Taxmodelqp, type = "deviance")
# Plot deviance residuals against predicted values
plot(fitted(Taxmodelqp), dev_residuals, main = "Deviance Residuals vs. Fitted Values")

### poisson - results show tax model cannot use poisson distribution
scale_0.5SD <- function(x){
  (x - mean(x)) / (2*sd(x))
}

MasterTax$page_views.x_scaled <- scale_0.5SD(MasterTax$page_views.x)
MasterTax$ED_scaled <- scale_0.5SD(MasterTax$ED)
MasterTax$bodymass_scaled <- scale_0.5SD(MasterTax$bodymass)

Taxmodelnb <- glmmTMB(num.studies ~ redlistCategory + ED_scaled + page_views.x_scaled +orderName  + endemic + venom + bodymass_scaled, family=nbinom2, ziformula = ~0, data=MasterTax)

options(na.action = "na.fail") #Must run this code once to use dredge
all_models_tax <- dredge(Taxmodelnb)
#write.csv(data.frame(all_models_tax),"all_models_tax_new.csv")

all_models_tax
nrow(all_models_tax)

#several models with deltaAICc <2
m.top.models.2aic <- get.models(all_models_tax, subset = delta <2)
length(m.top.models.2aic)

all_models_tax_avg <- summary(model.avg(m.top.models.2aic))

all_models_tax_avg
sw(all_models_tax_avg)
confint(all_models_tax_avg,full=TRUE)

#averaged model summary
#write.csv(summary(all_models_tax_avg)[9],"glm_tax_summaryupdate_new.csv")
#write.csv(sw(all_models_tax_avg),"glm_tax_relimp_new.csv")
#write.csv(confint(all_models_tax_avg,full=TRUE),"glm_tax_confint_new.csv")

#top model summaries
#write.csv(summary(get.models(all_models_tax, subset = delta <2)[[1]])[6][[1]][1],"glm_tax_summaryupdatetop3_1_new.csv")
#write.csv(summary(get.models(all_models_tax, subset = delta <2)[[2]])[6][[1]][1],"glm_tax_summaryupdatetop3_2_new.csv")
#write.csv(summary(get.models(all_models_tax, subset = delta <2)[[3]])[6][[1]][1],"glm_tax_summaryupdatetop3_3_new.csv")

#write.csv(confint(get.models(all_models_tax, subset = delta <2)[[1]],full=TRUE),"glm_tax_confint_top3_1new.csv")
#write.csv(confint(get.models(all_models_tax, subset = delta <2)[[2]],full=TRUE),"glm_tax_confint_top3_2new.csv")
#write.csv(confint(get.models(all_models_tax, subset = delta <2)[[3]],full=TRUE),"glm_tax_confint_top3_3new.csv")

#top 100 studied species
#write.csv(MasterTax[rev(order(page_views.x)),][1:100],"top100species_new.csv")

#all species with at least one study
#write.csv(MasterTax[num.studies>0,][rev(order(num.studies))],"allstudiedspecies_new.csv")

#test assumption of neg.binom.
# m3 <- glm(num.studies ~ redlistCategory + ED_scaled + page_views.x_scaled +orderName  + endemic + venom + bodymass_scaled, family = "poisson", data = MasterTax)
# pchisq(2 * (logLik(Taxmodelnb) - logLik(m3)), df = 1, lower.tail = FALSE)
# 2 * (logLik(Taxmodelnb) - logLik(m3))

####################################################
## Tree of families

head(MasterTax)

rt <- "(((Chelidae, (Pelomedusidae, Podocnemididae)), ((Carettochelyidae, Trionychidae), 
((((Dermochelyidae, Cheloniidae), (Chelydridae, (Dermatemydidae, (Kinosternidae)))), ((Platysternidae, Emydidae), (Testudinidae, Geoemydidae)))))), 
(((Sphenodontidae, (Dibamidae, ((((Carphodactylidae, Pygopodidae), Diplodactylidae), ((Eublepharidae, ( Sphaerodactylidae, (Phyllodactylidae, Gekkonidae)))))
                       ,(((((Xantusiidae, (Gerrhosauridae, Cordylidae)), Scincidae), 
((((Teiidae, Gymnophthalmidae), Alopoglossidae, (Lacertidae, (Rhineuridae, (Bipedidae, (Blanidae, (Cadeidae, (Trogonophidae, Amphisbaenidae))))))))
, (((((Helodermatidae,(Xenosauridae,(Diploglossidae, (Anguidae)))), 
(Shinisauridae, (Lanthanotidae, Varanidae)))), 
  ((Chamaeleonidae, Agamidae), 
  (Leiocephalidae,(Iguanidae, (Tropiduridae, (((Hoplocercidae,(Crotaphytidae, Corytophanidae), 
  (Phrynosomatidae,(Polychrotidae, Dactyloidae), ((Liolaemidae, (Opluridae, Leiosauridae))))))))))))
, (Leptotyphlopidae, (Gerrhopilidae, (Xenotyphlopidae, Typhlopidae))), 
(Anomalepididae, ((Aniliidae, Tropidophiidae), (Boidae,(Xenophidiidae, Bolyeriidae)),
(Xenopeltidae, (Loxocemidae, Pythonidae))),  ((Uropeltidae, (Cylindrophiidae, Anomochilidae)))))
, ((Acrochordidae, (Xenodermidae, (Pareidae, (Viperidae, ((Homalopsidae, (Colubridae, 
((Prosymnidae, Psammophiidae), ((Atractaspididae), 
(Pseudaspididae, (Pseudoxyrhophiidae, Elapidae, Lamprophiidae, Cyclocoridae))))))))))))))))))))), 
((Alligatoridae, (Crocodylidae, Gavialidae))));"


#place holder structure of the tree

rttree <- read.tree(text = rt)

# Check the structure of the phylogenetic tree
plot(rttree)

suborder_labs <- rttree$tip.label

x = as_tibble(rttree)

# get species data
treeddat <- MasterTax[,list(familyName,num.studies,binom)]
treeddat[familyName=="Elapoidea",familyName:="Lamprophiidae"] #move Buhoma genus here, despite its uncertainty

d <- treeddat %>%
  group_by(familyName) %>%
  summarise(avg.num.studies=mean(num.studies))

d$label <- d$familyName
d <- rbind(d,tibble(familyName="Sphenodontidae",avg.num.studies=6,label="Sphenodontidae"))

setdiff(d$label,x$label)
setdiff(x$label,d$label)


# Join tree data with species data
treecon <- full_join(x, d, by = "label")
treecon$sqrt.avg.num.studies <- sqrt(treecon$avg.num.studies) 
#make zeros NAs so we can colour them easily with grey
treecon$sqrt.avg.num.studies[which(treecon$sqrt.avg.num.studies==0)] <- NA

treecon1 <-  as.treedata(treecon) 

str(treecon1)


breaks <- c(0.01,0.25,1,3,5,10,20,30)
data.limits <- c(0,31)

# max(treecon$avg.num.studies,na.rm=TRUE)
# min(treecon$avg.num.studies,na.rm=TRUE)
# max(treecon$sqrt.avg.num.studies,na.rm=TRUE)
# min(treecon$sqrt.avg.num.studies,na.rm=TRUE)**2

ggtree(treecon1, layout="fan", open.angle=45) + 
  geom_tiplab2(size=6.75, offset=0.75) + 
  geom_tippoint(aes(fill=sqrt.avg.num.studies),shape=21,size=10)+
  scale_fill_gradientn(name="Studies per species",
                       breaks=sqrt(breaks),
                       limits=sqrt(data.limits),
                       labels=breaks,
                       colours=c(brewer.pal(n=9, name="OrRd"))[c(1,2,5,7,9)],
                       na.value="grey50")+
  theme(legend.position=c(0.90,0.35),
        legend.key.size = unit(2, 'cm'),
        axis.text = element_text(size=15),
        legend.text = element_text(size=20),
        legend.title = element_text(size=22))+xlim(0,28)

#ggsave("reptiletree.png",height=40,width=40,units="cm",dpi=600)
#ggsave("reptiletree.svg",height=40,width=40,units="cm",dpi=600,device="svg")



