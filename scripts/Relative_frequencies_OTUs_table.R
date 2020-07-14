#This script considers how we determined the relative frequencies on the filtered OTUs abundance.
#The dataset can be obtained on request writing to giorgio.giraffa@crea.gov.it

library(vegan)

Aprile18 <-read.csv2("data/aprile18_filtrato_finale.csv", header=TRUE)

df_Aprile18 <- ((Aprile18)[-1:-2])

sample.tot = rowSums(df_Aprile18)
datarow = df_Aprile18
for(i in 1:nrow(df_Aprile18)){
  datarow[i,] = df_Aprile18[i,]/sample.tot[i]
}

datarow = cbind(Aprile18[,1:2], datarow)

rowSums(datarow[,-1:-2])

write.table(x = datarow, file = "data/Apr18_freq relative.csv", sep=";", dec=",", quote = TRUE)