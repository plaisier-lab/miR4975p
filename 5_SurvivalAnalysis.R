###########
## Bueno ##
###########

#setwd('D:/Dropbox (ASU)/bi_hsa_miR497_PlaisierAnalysis/Manuscript/Survival_SW/')
setwd('D:/Dropbox (ASU)/biMIR497_shared/Data_and_Code/0_Code')

# Gene counts and RPKMs
load('../4_SurvivalAnalyses/mesothelioma_genentech/meso_genecounts_18Nov2016.rda')

# Normalize
library(DESeq2)
p1 = read.csv('../4_SurvivalAnalyses/mesothelioma_genentech/phenotypes_meso.csv',header=T,row.names=1)
conds = data.frame(patients=p1[colnames(gene),'consensus_cluster'])
rownames(conds) = colnames(gene)
dds = DESeqDataSetFromMatrix(countData = gene, colData = conds, design= ~ 1)
vst1 = vst(dds)
o1 = assay(vst1)


# Genes to map
#genes1 = intersect(read.csv('330genesInOverlap.csv', header=T, row.names=1)[,1],rownames(o1))
genes1 = "CCND2"

library(survival)

pdf('../4_SurvivalAnalyses/results/genesSurvial_Bueno_CCDN2.pdf')
cph1 = matrix(ncol=2, nrow=length(genes1))
rownames(cph1) = genes1
colnames(cph1) = c('Bueno.coef','Bueno.p')
for(gene1 in genes1) {
    cox.obj <- coxph(Surv(p1[colnames(o1),'survival'],p1[colnames(o1),'status']==1)~o1[gene1,]+p1[colnames(o1),'age_at_surgery'])
    cph1[gene1,] = c(summary(cox.obj)$coefficients[1,'coef'], summary(cox.obj)$coefficients[1,'Pr(>|z|)'])
    sfit <- survfit(Surv(p1[colnames(o1),'survival'],p1[colnames(o1),'status']==1)~o1[gene1,]>median(as.numeric(o1[gene1,])))
    sd1 <- survdiff(Surv(p1[colnames(o1),'survival'],p1[colnames(o1),'status']==1)~o1[gene1,]>median(as.numeric(o1[gene1,])))
    sd1
    plot(sfit,col=c('dodgerblue2','orchid2'),lty=c(1,2),ylab='Survival Probability (%)',xlab='Time',main=paste(gene1, ': Cox PH p-value',signif(summary(cox.obj)$coefficients[1,'Pr(>|z|)'],2),'; Diff p-value = ',signif(1-pchisq(sd1$chisq,1),2)))
    legend("topright", legend=c('Less than median','Greater than median'), col=c('dodgerblue2','orchid2'), lty=c(1,2))
}
dev.off()

eg_genes2 = read.csv('../5_Figures/results/330genesInOverlap.csv', header=T)[rownames(cph1),]
cph1 = cbind(cph1, Bueno.p.adjust=p.adjust(cph1[,'Bueno.p'], method='BH'), eg_genes2[,1])
write.csv(cph1, '../4_SurvivalAnalyses/results/survival_Bueno_330genesOverlap.csv')

eg_genes2 = "CCND2"
cph1 = cbind(cph1, Bueno.p.adjust=p.adjust(cph1[,'Bueno.p'], method='BH'), eg_genes2)
write.csv(cph1, '../4_SurvivalAnalyses/results/survival_Bueno_CCND2.csv')


##########
## TCGA ##
##########
#setwd('D:/Dropbox (ASU)/bi_hsa_miR497_PlaisierAnalysis/Manuscript/Survival_SW/')

# Load RNA-seq
d1 = read.csv(gzfile('../4_SurvivalAnalyses/mesothelioma_TCGA/mesoTCGA_ratios.csv.gz'),header=T,row.names=1)
rownames(d1) = paste('egID',rownames(d1),sep='_')

# Load miRNA-seq
miRseq = read.csv('../4_SurvivalAnalyses/mesothelioma_TCGA/MESO_miRNAseq.csv',header=T,row.names=1)

# Load miRNA-seq
p1 = read.csv('../4_SurvivalAnalyses/mesothelioma_TCGA/phenotypes_mesoTCGA.csv',header=T,row.names=1)
rownames(p1) = gsub('-','.',rownames(p1))
p1[,'OS.time'] = (p1[,'OS.time']/30.42)/12

# Need to translate symbol to entrez and combine with egID
library('org.Hs.eg.db')
symbols = read.csv('../5_Figures/results/330genesInOverlap.csv', header=T, row.names=1)[,1]
mappedAnswers = as.data.frame(as.matrix(mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')))

eg_genes1 = intersect(paste0('egID_',mappedAnswers[,1]),rownames(d1))

pdf('../4_SurvivalAnalyses/results/genesSurvial_TCGA_330genesInOverlap.pdf')
cph1 = matrix(ncol=2, nrow=length(eg_genes1))
rownames(cph1) = eg_genes1
colnames(cph1) = c('TCGA.coef','TCGA.p')
for(gene1 in eg_genes1) {
    cox.obj <- coxph(Surv(p1[colnames(d1),'OS.time'],p1[colnames(d1),'OS']==1)~as.numeric(d1[gene1,])+p1[colnames(d1),'age_at_initial_pathologic_diagnosis'])
    cph1[gene1,] = c(summary(cox.obj)$coefficients[1,'coef'], summary(cox.obj)$coefficients[1,'Pr(>|z|)'])
    sfit <- survfit(Surv(p1[colnames(d1),'OS.time'],p1[colnames(d1),'OS']==1)~as.numeric(d1[gene1,])>median(as.numeric(d1[gene1,])))
    sd1 <- survdiff(Surv(p1[colnames(d1),'OS.time'],p1[colnames(d1),'OS']==1)~as.numeric(d1[gene1,])>median(as.numeric(d1[gene1,])))
    sd1
    plot(sfit,col=c('dodgerblue2','orchid2'),lty=c(1,2),ylab='Survival Probability (%)',xlab='Time',main=paste(gene1, ': Cox PH p-value',signif(summary(cox.obj)$logtest['pvalue'],2),'; Diff p-value = ',signif(1-pchisq(sd1$chisq,1),2)),xlim=c(0,11.5))
    legend("topright", legend=c('Less than median','Greater than median'), col=c('dodgerblue2','orchid2'), lty=c(1,2))
}
dev.off()

eg_genes2 = read.csv('../5_Figures/results/330genesInOverlap.csv')[unlist(sapply(rownames(cph1), function(x) { strsplit(x,'_')[[1]][2] } )),]

#mappedAnswers1 = as.data.frame(as.matrix(mapIds(org.Hs.eg.db, gsub('egID_', '', rownames(cph1)), 'SYMBOL', 'ENTREZID')))
cph1 = cbind(cph1, TCGA.p.adjust=p.adjust(cph1[,'TCGA.p'], method='BH'))

write.csv(cph1, '../4_SurvivalAnalyses/results/survival_TCGA_330genesInOverlap.csv')
