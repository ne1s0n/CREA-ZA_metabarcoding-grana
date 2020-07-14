#This script does on the heatmap of relative abundance of the dominant (> 1% total reads) and subdominant (0.1-1% total reads) species.
#The dataset can be obtained on request writing to giorgio.giraffa@crea.gov.it

library(reshape2)
library(plyr)
library(scales)
library(ggplot2)

#Dominant species (> 1% total reads)
dominanti_grana <- read.csv2("data/Apr18_specie_magg_unopercento_prime15 specie_regione_per manuscript.csv", header = TRUE, stringsAsFactors = FALSE)

meta_grana <- ((dominanti_grana)[-1])

meta_grana.ordered <- meta_grana [order(meta_grana[,1]),] 
meta_grana.ordered

df_grana<-melt(meta_grana.ordered)

colnames(df_grana)<-c("Regions","Species","Value")
colnames(df_grana)

df_grana<-ddply(df_grana,.(Regions),transform,rescale=sqrt(Value))

grana <- ggplot(df_grana, aes(Species, Regions)) + 
  geom_tile(aes(fill = Value),colour = "gray") + 
  scale_fill_gradient(low = "white",high = "blue") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + theme(legend.spacing.y = unit(1,"lines"),axis.ticks = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1,size=11, face = "italic"),axis.text.y = element_text(hjust = 0.5,size=11,vjust = 0.8))
grana

#subdominant species (0.1-1% total reads)

secondarie_grana <- read.csv2("data/Specie-secondarie_prime15 specie_regione", header = TRUE, stringsAsFactors = FALSE)

meta_grana <- ((secondarie_grana)[-1])

meta_grana.ordered <- meta_grana [order(meta_grana[,1]),] 

df_grana<-melt(meta_grana.ordered)

colnames(df_grana)<-c("Regions","Species","Value")

df_grana<-ddply(df_grana,.(Regions),transform,rescale=sqrt(Value))

grana <- ggplot(df_grana, aes(Species, Regions)) + 
  geom_tile(aes(fill = Value),colour = "gray") + 
  scale_fill_gradient(low = "white",high = "red")+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + theme(legend.spacing.y = unit(1,"lines"),axis.ticks = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1,size=11, face = "italic"),axis.text.y = element_text(hjust = 1,size=11,vjust = 0.8))
grana

pdf("Heatmap.pdf")
print(p)
dev.off()