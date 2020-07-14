#This script considers how we performed the taxonomic analysis on the relative abundance of the dominant species (> 1% total reads).
#The dataset can be obtained on request writing to giorgio.giraffa@crea.gov.it

library(vegan)
library(reshape2)
library(ggplot2)

dominant.species <- read.csv2("data/Apr18_specie_magg_unopercento_per manuscript.csv", header = TRUE, stringsAsFactors = FALSE) 

specie.data <- dominant.species[order(dominant.species[,1]),] 
specie.data

dataspecie.dominanti <- ((specie.data)[,-2])
dataspecie.dominanti

data_long.dominanti <- melt(dataspecie.dominanti, id.vars = "Campioni", variable.name = "Species")

species.dominanti <- ggplot(data_long.dominanti, aes(x = Campioni, y = value, fill = Species)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#FF0000", "#0000FF", "#006600", "#666666", "#6600FF", "#FF6633", "#663300", "#FFCC00", "#FF6699", "#000000", "#330066", "#CC99FF", "#FF99FF", "#CC99CC", "#FF9999", "#FFFF99", "#CCCC99", "#FFFFCC", "#99FF99", "#99FFCC", "#66CCCC", "#66FFFF", "#66CCFF", "#00FF99", "#33FF99", "#00FFFF", "#66CCFF", "#6666FF", "#00CC99", "#00FFCC", "#00FF99", "#00FF66")) +
  xlab ("Samples") + ylab ("Abundance") +
  theme_classic()+ theme(legend.text = element_text(face = "italic"), axis.text.x=element_text(hjust=1, vjust=0.8, angle = 90)) + 
  coord_cartesian(ylim=c(0.0,1.0))
species.dominanti