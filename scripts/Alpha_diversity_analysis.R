#This script considers the alpha diversity index (rarefaction curve, richness, Shannon diversity) on the filtered OTUs abundance.
#The dataset can be obtained on request writing to giorgio.giraffa@crea.gov.it

library(vegan)
library(agricola)

#rarefaction curve
aprile18_filtrato_finale_definitivo_senza.c1660 <- read.csv("data/aprile18_filtrato_finale_definitivo_senza c1660.csv", sep=";", row.names = 1)

raremax <- min(rowSums(aprile18_filtrato_finale_definitivo_senza.c1660))

data <- rarecurve(aprile18_filtrato_finale_definitivo_senza.c1660, step = 1, col = "blue", label = TRUE, cex = 0.6)

abline(v = raremax)

#Richness and Shannon diversity
data <- read.csv("/aprile18_filtrato_finale_nopiemonteec1660.csv", sep=";")

data$Regione <- data[order(data[,1]),]

abundance.matrix <- data[,2:254]

indices <- data[,c("Regione")]

indices$Richness <- rowSums(abundance.matrix>0)

indices$Shannon <- diversity(abundance.matrix) 

par(family = "serif")

colors = rainbow(11)[10:1]

boxplot(Richness~Region, data=indices, boxwex= 0.90 , col=colors, 
        cex.axis= 0.90 , ylab="Richness")

boxplot(Shannon~Region, data=indices, boxwex=0.90, col=colors, 
        cex.axis=0.90, ylab="Shannon diversity")

mod.richness <- lm(Richness~Region, data=indices)
mod.shannon <- lm(Shannon~Region, data=indices)

anova(mod.richness)

summary(mod.richness)

plant.av <- aov(mod.richness)

summary(plant.av)

tukey.test2 <- HSD.test(plant.av, trt = 'Region')
tukey.test2

anova(mod.shannon)

summary(mod.shannon)

prova.av <- aov(mod.shannon)

summary(prova.av)

tukey.test2 <- HSD.test(prova.av, trt = 'Region')
tukey.test2