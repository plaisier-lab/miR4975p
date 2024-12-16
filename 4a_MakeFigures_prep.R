###########
## Genes ##
###########
## Run in R ##
#setwd('D:/Dropbox (ASU)/Bi_hsa_miR_497/')
#setwd('C:/Users/wilfe/Dropbox (ASU)/biMIR497_shared/Data_and_Code/0_Code')
setwd('D:/Dropbox (ASU)/biMIR497_shared/Data_and_Code/0_Code')


# Apoptosis: GO:0042981
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id'), filters = 'go', values = 'GO:0042981', mart = ensembl)
write.csv(gene.data,paste('../5_Figures/NeededInfo/gene.data.csv'),sep='')

# Cell cycle: GO:0007049
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id'), filters = 'go', values = 'GO:0007049', mart = ensembl)
write.csv(gene.data,paste('../5_Figures/NeededInfo/cellCycle_GO0007049.csv'),sep='')

# apoptosis symbol converter - did not convert all using translator below
library('SeqGSEA')
listR = c('AIFM1', 'AKT1', 'AKT2', 'AKT3', 'APAF1', 'ATM', 'BAD', 'BAX', 'BCL2', 'BCL2L1', 'BID', 'BIRC2', 'BIRC3', 'CAPN1', 'CAPN2', 'CASP10', 'CASP3', 'CASP6', 'CASP7', 'CASP8', 'CASP9', 'CFLAR', 'CHP1', 'CHP2', 'CHUK', 'CSF2RB', 'CYCS', 'DFFA', 'DFFB', 'ENDOD1', 'ENDOG', 'EXOG', 'FADD', 'FAS', 'FASLG', 'IKBKB', 'IKBKG', 'IL1A', 'IL1B', 'IL1R1', 'IL1RAP', 'IL3', 'IL3RA', 'IRAK1', 'IRAK2', 'IRAK3', 'IRAK4', 'MAP3K14', 'MYD88', 'NFKB1', 'NFKBIA', 'NGF', 'NTRK1', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIK3R5', 'PPP3CA', 'PPP3CB', 'PPP3CC', 'PPP3R1', 'PPP3R2', 'PRKACA', 'PRKACB', 'PRKACG', 'PRKAR1A', 'PRKAR1B', 'PRKAR2A', 'PRKAR2B', 'PRKX', 'RELA', 'RIPK1', 'TNF', 'TNFRSF10A', 'TNFRSF10B', 'TNFRSF10C', 'TNFRSF10D', 'TNFRSF1A', 'TNFSF10', 'TP53', 'TRADD', 'TRAF2', 'XIAP')

yesman = convertSymbol2Ensembl(listR)
write.csv(yesman,paste('../5_Figures/NeededInfo/apoptosisPathwayGenes.csv'),sep='')

# cell cycle pathway genes
library(KEGGREST)
library(org.Hs.eg.db)
library(tidyverse)

hsa_path_eg  <- keggLink("pathway", "hsa") %>% tibble(pathway = ., eg = sub("hsa:", "", names(.)))
hsa_kegg_anno <- hsa_path_eg %>% mutate(symbol = mapIds(org.Hs.eg.db, eg, "SYMBOL", "ENTREZID"),ensembl = mapIds(org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID"))

write.csv(hsa_kegg_anno,paste('../5_Figures/NeededInfo/hsa_kegg_anno.csv'),sep='')

# fin
