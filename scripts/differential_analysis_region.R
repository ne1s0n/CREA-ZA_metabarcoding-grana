#This script does differential analysis on all possible pair of regions.
#It starts from a data/aprile18_filtrato_finale.csv containing the filtered
#OTUs and that can be obtained on request writing to giorgio.giraffa@crea.gov.it

library(plyr)
library(phyloseq)
library(DESeq2)
library(ggplot2)

#this should be changed accordingly
setwd('~/research/crea-za_newtech_microbioma/')

#chosen significance level
alpha = 0.01

#outfile name root
outfile_root = paste(sep='', 'results/differential_analysis_region_alpha', alpha)

# OTU and sample data -----------------------------------------------------
OTU = read.table(
  'data/aprile18_filtrato_finale.csv', 
  stringsAsFactors = FALSE, sep=';', header = TRUE, row.names = 1)

#removing bad sample
OTU = OTU[rownames(OTU) != 'C1660_sp_BS601', ]

#break into OTU and meta data, describing region and province
samples = OTU[,1:2]
OTU = OTU[,-1:-2]

#removing data on unclassified samples
OTU$Unclassified = NULL

#moving to phyloseq format
OTU.phyloseq = otu_table(OTU, taxa_are_rows = FALSE)

#converting sample data, with some attention to string/factors
samples$Provincia = factor(samples$Provincia)
samples$Region = gsub(samples$Region, pattern = '-', replacement = '_')
samples$Region = gsub(samples$Region, pattern = ' ', replacement = '_')
samples$Region = factor(samples$Region)
samples.phyloseq = sample_data(samples)

# PHYLOGENETIC INFO -------------------------------------------------------
#extracting taxonomy info from the colnames in the OTU table
tax = data.frame(
  blob = colnames(OTU)
)
tax = ddply(tax, .(blob), function(x){
  #genus.species need to be split
  pieces = strsplit(x$blob, split = '.', fixed = TRUE)[[1]]
  
  #when the genus is "unclassified" we have something like "unclassified.Lactobacillaceae"
  #and we swap the terms (e.g. genus=Lactobacillaceae species=unclassified)
  if (pieces[1] == 'unclassified'){
    tmp = pieces[1]
    pieces[1] = pieces[2]
    pieces[2] = tmp
  }
  
  return(data.frame(
    Genus = pieces[1],
    Species = pieces[2]
  ))
})

#blob need to be used as rownames at this point
rownames(tax) = tax$blob
tax$blob = NULL

#moving to phyloseq format
tax.phyloseq = tax_table(as.matrix(tax))

# DIFFERENTIAL ANALYSIS ---------------------------------------------------
#creating a single phyloseq object
grana.phyloseq = phyloseq(OTU.phyloseq, samples.phyloseq, tax.phyloseq)

#DeSeq2 format conversion
grana.deseq2 = phyloseq_to_deseq2(grana.phyloseq, ~ Region)

#actual differential analysis
grana.deseq2 = DESeq(grana.deseq2, test="Wald", fitType="parametric")

#for each pair or region, we extract results
res = NULL
pairs = t(combn(unique(samples$Region), m=2))

for (i in 1:nrow(pairs)){
  region1 = as.character(pairs[i, 1])
  region2 = as.character(pairs[i, 2])
  
  #extracting the results, with shrinkage
  res.now = lfcShrink(grana.deseq2, contrast = c('Region', region1, region2), type = 'ashr')

  #subsetting on alpha level (posthoc)
  res.now = res.now[which(res.now$padj < alpha), ]
  
  if (nrow(res.now) == 0){
    writeLines(paste('Nothing significant for', region1, '-', region2))
    #no significantly discriminant species found
    next
  }
  
  #if we get here, something is significant
  res.now = cbind(as(res.now, "data.frame"), as(tax_table(grana.phyloseq)[rownames(res.now), ], "matrix"))
  res.now$region1 = region1
  res.now$region2 = region2

  #putting all together
  res = rbind(res, res.now)
}

#screen names
res$region1 = gsub(res$region1, pattern = '_', replacement = '-')
res$region2 = gsub(res$region2, pattern = '_', replacement = '-')

#saving results
write.csv(res, paste(sep='', outfile_root, '.csv'))

# MAKING A PLOT -----------------------------------------------------------
#sorting by region pairs, then by effect
res2 = res[order(res$region1, res$region2, res$log2FoldChange, decreasing = TRUE),]

#showable labels
res2$regions = paste(res2$region1, 'vs.', res2$region2)
res2$organisms = paste(res2$Genus, res2$Species)

#a plot, by organism
p = ggplot(res2, aes(y=log2FoldChange, color=regions, x=organisms, group=organisms)) + 
  geom_point(size = 4, position = position_dodge2(width=0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, color = 'black')) +
  theme(axis.text.y = element_text(color = 'black')) +
  scale_x_discrete('') +
  scale_y_continuous('log2 Fold Change') +
  scale_color_brewer(palette="Dark2", 'Regions') +
  guides(colour = guide_legend(title.position = "top", ncol = 2)) + 
  theme(legend.position = 'bottom')

p

ggsave(paste(sep='', outfile_root, '.png'), p, height = 5)

#a second plot, by region pair
res2$regions = paste(sep='', res2$region1, ' vs.\n', res2$region2)
p = ggplot(res2, aes(y=log2FoldChange, color=organisms, x=regions, group=regions)) + 
  geom_point(size = 4, position = position_dodge2(width=0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, color = 'black')) +
  theme(axis.text.y = element_text(color = 'black')) +
  scale_x_discrete('') +
  scale_y_continuous('log2 Fold Change') +
  scale_color_brewer(palette="Dark2", 'Organism') +
  guides(colour = guide_legend(title.position = "top", ncol = 2)) + 
  theme(legend.position = 'bottom')

p

ggsave(paste(sep='', outfile_root, '.by_region.png'), p, height = 5)
