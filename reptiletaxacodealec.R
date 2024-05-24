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

library(ape)
library(cluster)
library(RColorBrewer)
library(ggplot2)


setwd("~/reptile-evidence-analysis")
MasterTax <- fread("~/Reptile-taxa-data.csv")

## Taxonomy GLM
unique(MasterTax$redlistCategory)
MasterTax$redlistCategory <- factor(MasterTax$redlistCategory, levels=c("Data Deficient", "Not Assessed", "Least Concern", "Near Threatened", "Vulnerable", "Endangered", "Critically Endangered"))
unique(MasterTax$redlistCategory)
unique(MasterTax$orderName)
MasterTax$orderName <- factor(MasterTax$orderName, levels=c("Squamata","Testudines","Crocodilia"))
MasterTax$page_views.x <- as.integer(MasterTax$page_views.x)

### taxonomic model if poisson - results show tax model cannot use poisson distribution
Taxmodel1 <- glm(num.studies ~ redlistCategory + ED + page_views.x +orderName, data=MasterTax, family=poisson)
testDispersion(Taxmodel1)
simulationOutput <- simulateResiduals(fittedModel = Taxmodel1, plot = F)
residuals(simulationOutput)
plot(simulationOutput)

### ZERO-INFLATED POISSON REGRESSION suggests that the excess zeros are generated from a separate process from the count values
### Tests for the quasipoisson distribution for Tax model - not symmetrical
Taxmodelqp <- glm(num.studies ~ redlistCategory + ED + page_views.x +orderName, data=MasterTax, family=quasipoisson)
dev_residuals <- residuals(Taxmodelqp, type = "deviance")
# Plot deviance residuals against predicted values
plot(fitted(Taxmodelqp), dev_residuals, main = "Deviance Residuals vs. Fitted Values")

MasterTax$page_views.x_scaled <- scale(MasterTax$page_views.x)
MasterTax$ED_scaled <- scale(MasterTax$ED)

#model selection
Taxmodelnb <- glmmTMB(num.studies ~ redlistCategory + ED_scaled + page_views.x_scaled +orderName, family=nbinom2, ziformula = ~0, data=MasterTax)

Taxmodelnb3 <- glmmTMB(num.studies ~ redlistCategory + ED_scaled + page_views.x_scaled, family=nbinom2, ziformula = ~0 ,data=MasterTax)
Taxmodelnb4 <- glmmTMB(num.studies ~ redlistCategory + ED_scaled + orderName, family=nbinom2, ziformula = ~0 ,data=MasterTax)
Taxmodelnb5 <- glmmTMB(num.studies ~ redlistCategory + page_views.x_scaled + orderName, family=nbinom2, ziformula = ~0 ,data=MasterTax)
Taxmodelnb6 <- glmmTMB(num.studies ~ ED_scaled + page_views.x_scaled +orderName, family=nbinom2, ziformula = ~0 ,data=MasterTax)

Taxmodelnb7 <- glmmTMB(num.studies ~ ED_scaled + page_views.x_scaled, family=nbinom2, ziformula = ~0 ,data=MasterTax)
Taxmodelnb8 <- glmmTMB(num.studies ~ ED_scaled + orderName, family=nbinom2, ziformula = ~0 ,data=MasterTax)
Taxmodelnb9 <- glmmTMB(num.studies ~ ED_scaled + redlistCategory, family=nbinom2, ziformula = ~0 ,data=MasterTax)
Taxmodelnb10 <- glmmTMB(num.studies ~ page_views.x_scaled + orderName, family=nbinom2, ziformula = ~0 ,data=MasterTax)
Taxmodelnb11 <- glmmTMB(num.studies ~ page_views.x_scaled + redlistCategory, family=nbinom2, ziformula = ~0 ,data=MasterTax)
Taxmodelnb12 <- glmmTMB(num.studies ~ redlistCategory + orderName, family=nbinom2, ziformula = ~0 ,data=MasterTax)

Taxmodelnb13 <- glmmTMB(num.studies ~ page_views.x_scaled, family=nbinom2, ziformula = ~0 ,data=MasterTax)
Taxmodelnb14 <- glmmTMB(num.studies ~ ED_scaled , family=nbinom2, ziformula = ~0 ,data=MasterTax)
Taxmodelnb15 <- glmmTMB(num.studies ~ orderName, family=nbinom2, ziformula = ~0 ,data=MasterTax)
Taxmodelnb16 <- glmmTMB(num.studies ~ redlistCategory, family=nbinom2, ziformula = ~0 ,data=MasterTax)

anova(Taxmodelnb, Taxmodelnb3) 
anova(Taxmodelnb, Taxmodelnb4) 
anova(Taxmodelnb, Taxmodelnb5)##nb5 and original model are close in AIC but within 2, so keep original model
anova(Taxmodelnb, Taxmodelnb6) 
anova(Taxmodelnb, Taxmodelnb7) 
anova(Taxmodelnb, Taxmodelnb8) 
anova(Taxmodelnb, Taxmodelnb9) 
anova(Taxmodelnb, Taxmodelnb10) 
anova(Taxmodelnb, Taxmodelnb11) 
anova(Taxmodelnb, Taxmodelnb12) 
anova(Taxmodelnb, Taxmodelnb13) 
anova(Taxmodelnb, Taxmodelnb14) 
anova(Taxmodelnb, Taxmodelnb15) 
anova(Taxmodelnb, Taxmodelnb16) 


simulationOutput <- simulateResiduals(fittedModel = Taxmodelnb)
plot(simulationOutput) ## deviation not significant
testZeroInflation(simulationOutput)

## Taxmodelnb is the best model that passes DHARMa
#look at model summary and analysis of deviance test
summary(Taxmodelnb)
Anova(Taxmodelnb, type="II")

emmeans(Taxmodelnb, pairwise ~ redlistCategory, type="response")
write.csv(emmeans(Taxmodelnb, pairwise ~ redlistCategory, type="response")$contrasts,"emmeanstaxredlist.csv")
emmeans(Taxmodelnb, pairwise ~ orderName, type="response")
write.csv(emmeans(Taxmodelnb, pairwise ~ orderName, type="response")$contrasts,"emmeanstaxorder.csv")
mean(MasterTax[orderName=="Squamata",num.studies])
mean(MasterTax[orderName=="Testudines",num.studies])

#test assumption of neg. binom model
m3 <- glm(num.studies ~ redlistCategory + ED_scaled + page_views.x_scaled +orderName, family = "poisson", data = MasterTax)
pchisq(2 * (logLik(Taxmodelnb) - logLik(m3)), df = 1, lower.tail = FALSE)
2 * (logLik(Taxmodelnb) - logLik(m3))


#most viewed species
MasterTax[rev(order(page_views.x))][1:100]
#most studied species
MasterTax[num.studies>0][rev(order(num.studies))]

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



