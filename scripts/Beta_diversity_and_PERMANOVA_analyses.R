#This script considers the beta diversity and PERMANOVA analyses on the filtered OTUs abundance.
#The dataset can be obtained on request writing to giorgio.giraffa@crea.gov.it
#Beta diversity analysis#

library(vegan)
library(ggplot2)

aprile18_filtrato_finale <- read.csv2("data/aprile18_filtrato_finale.csv", sep=";", header = TRUE)

bio <- ((aprile18_filtrato_finale)[-1])

prova <- metaMDS(bio)

NMDS_bio <- data.frame(NMDS1=prova$points[,1], NMDS2=prova$points[,2])

complete <- as.data.frame(cbind(aprile18_filtrato_finale[c(1)],NMDS_bio))

NMDS <- ggplot(complete, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=Region), size=5 ) +
  theme_bw()
NMDS 

#PERMANOVA analysis

BC.dist=vegdist(bio, distance="bray")

adonis(BC.dist ~ Region, data = aprile18_filtrato_finale, permutations = 1000)