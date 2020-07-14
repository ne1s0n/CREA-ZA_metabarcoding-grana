#This script does OTUs filtering.
#It starts from data/OTUs_Aprile18.csv containing the original OTUs table from the Base-Clear company.
#The dataset can be obtained on request writing to giorgio.giraffa@crea.gov.it

# OTUs less to 5 reads were removed

OTUs_Apr18_Filt1 <- replace (OTUs_Aprile18, OTUs_Aprile18 == 1, 0)

OTUs_Apr18_Filt2<- replace (OTUs_Apr18_Filt1, OTUs_Apr18_Filt1 == 2, 0)

OTUs_Apr18_Filt3 <- replace (OTUs_Apr18_Filt2, OTUs_Apr18_Filt2 == 3, 0)

OTUs_Apr18_Filt4 <- replace (OTUs_Apr18_Filt3, OTUs_Apr18_Filt3 == 4, 0)

otu2018 <- ((OTUs_Apr18_Filt4)[-1:-2])

# OTUs more or equal to 5 reads were considered in our study

otuApr18clean <- as.data.frame(otu2018[,colSums(otu2018) >= 5])

APR18finale <- as.data.frame(cbind(infoOTUapr18, otuApr18clean))
 
write.table(x = APR18finale, file = "/data/aprile18_filtrato_finale.csv", sep = ";", dec = ",", quote = TRUE)