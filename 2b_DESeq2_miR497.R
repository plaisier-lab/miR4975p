#setwd('D:/Dropbox (ASU)/Bi_hsa_miR_497/')
#setwd('/Users/cplaisie/Dropbox (ASU)/Bi_hsa_miR_497/')
#setwd('C:/Users/wilfe/Dropbox (ASU)/biMIR497_shared/Data_and_Code/0_Code')
setwd('D:/Dropbox (ASU)/biMIR497_shared/Data_and_Code/0_Code')

# Load counts
d0 = read.csv('../2_DESeq2/gexp_counts.csv',header=T,row.names=1)
library(DESeq2)
conds = data.frame(
    id_time=sapply(colnames(d0), function(x) { paste0(substr(strsplit(x,'_')[[1]][1],1,2),'_',strsplit(x,'_')[[1]][2]) }))
    #time=sapply(colnames(d0), function(x) { as.numeric(strsplit(x,'_')[[1]][2]) }))
    #,
    #subtype=factor(c(rep('Epi',6), rep('Sarc',6),rep('Epi',6))))

d0 = d0[which(apply(d0,1,sum)!=0),]
dds = DESeqDataSetFromMatrix(countData = d0, colData = conds, design= ~ id_time) # subtype
#rld = rlog(dds,blind=T)
rld = vst(dds)
write.table(assay(rld),'../2_DESeq2/results/gexp_norm_all.tsv',sep='\t')
write.table(assay(rld),'../2_DESeq2/results/gexp_norm_all.csv',sep=',')

# Cluster samples
library(pheatmap)
pdf('../2_DESeq2/results/correlationReplicates.pdf')
#c1 = cor(assay(rld)[names(sort(apply(assay(rld), 1, mad),decreasing=T))[1:10000],],method='spearman')
c1 = cor(assay(rld),method='pearson')
diag(c1) = NA
pheatmap(c1,color=colorRampPalette(c('blue4','white','gold'))(32),clustering_method='single') # Looks much bett
dev.off()

plotPCA_v2 <- function(x, intgroup = "condition", ntop = 500, col)
{
    rv = rowVars(assay(x))
    select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(x)[select, ]))
    fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop = FALSE]),
        1, paste, collapse = " : "))
    if (missing(col)) {
        col = if (nlevels(fac) >= 3)
            brewer.pal(nlevels(fac), "Paired")
        else c("lightgreen", "dodgerblue")
    }
    imp = summary(pca)$importance[2,]
    par(mar=c(1,1,1,1))
    xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x),
        pch = 16, cex = 2, aspect = "iso", col = col, main = draw.key(key = list(rect = list(col = col),
            text = list(levels(fac)), rep = FALSE)), xlab = paste('PC1 (',signif(imp[[1]],2)*100,'%)',sep=''), ylab = paste('PC2 (',signif(imp[[2]],2)*100,'%)',sep=''))
}

library(matrixStats)
library(RColorBrewer)
library(lattice)

pdf('../2_DESeq2/results/pca_rlog_counts.pdf')
plotPCA_v2(rld, intgroup=c('id_time'),ntop=2000)
#plotPCA_v2(rld, intgroup=c('time'),ntop=2000)
dev.off()

# Do differential expression analyses
dseq1 = DESeq(dds, betaPrior=TRUE)

# Cnvert to entrz
# Biomart gene id conversioin
library(AnnotationDbi)
library(org.Hs.eg.db)
symb = select(org.Hs.eg.db, keys = rownames(dseq1), keytype = 'ENSEMBL', columns = 'SYMBOL')

# ei vs. ci 12H
res.ei_vs_ci_12H = results(dseq1, contrast=list('id_timeei_12', c('id_timeci_12')), listValues=c(1, -1/5))
symbs <- symb$SYMBOL[match(rownames(res.ei_vs_ci_12H), symb$ENSEMBL, nomatch = NA)]
res.ei_vs_ci_12H$symbols = symbs

# ep vs. cp 12H
res.ep_vs_cp_12H = results(dseq1, contrast=list('id_timeep_12', c('id_timecp_12')), listValues=c(1, -1/5))
symbs <- symb$SYMBOL[match(rownames(res.ep_vs_cp_12H), symb$ENSEMBL, nomatch = NA)]
res.ep_vs_cp_12H$symbols = symbs

# ei vs. ci 72H
res.ei_vs_ci_72H = results(dseq1, contrast=list('id_timeei_72', c('id_timeci_72')), listValues=c(1, -1/5))
symbs <- symb$SYMBOL[match(rownames(res.ei_vs_ci_72H), symb$ENSEMBL, nomatch = NA)]
res.ei_vs_ci_72H$symbols = symbs

# ep vs. cp 72H
res.ep_vs_cp_72H = results(dseq1, contrast=list('id_timeep_72', c('id_timecp_72')), listValues=c(1, -1/5))
symbs <- symb$SYMBOL[match(rownames(res.ep_vs_cp_72H), symb$ENSEMBL, nomatch = NA)]
res.ep_vs_cp_72H$symbols = symbs


# A list of all DE genes
allDEUp = c() #
allDEDown = c() #
allDEUpList = list()
allDEDownList = list()

deTmp = list(res.ei_vs_ci_12H, res.ep_vs_cp_12H, res.ei_vs_ci_72H, res.ep_vs_cp_72H)
tmp2 = c('ei_vs_ci_12H', 'ep_vs_cp_12H','ei_vs_ci_72H', 'ep_vs_cp_72H')
tmp3 = c('ei_vs_ci_12H', 'ep_vs_cp_12H','ei_vs_ci_72H', 'ep_vs_cp_72H')

# Plot volcano plots
pdf('../2_DESeq2/results/volcano_plots_DESeq2.pdf')
#par(mfrow=c(2,2))
for(i in 1:length(deTmp)) {
    tmp = deTmp[[i]]
    plot(tmp$log2FoldChange, -log10(tmp$padj),pch=20,col=rgb(0,0,1,0.25),ylab='-log10(Adjusted P-Value)',xlab='log2(Fold-Change)')
    # Add markers
    # Add lines
    abline(h=-log10(0.05),col=rgb(1,0,0,0.5),lwd=2,lty=2)
    abline(v=c(-1,1),col=rgb(1,0,0,0.5),lwd=2,lty=2)
    up = intersect(rownames(tmp)[which(tmp$log2FoldChange>=1)],rownames(tmp)[which(tmp$padj<=0.05)])
    down = intersect(rownames(tmp)[which(tmp$log2FoldChange<=-1)],rownames(tmp)[which(tmp$padj<=0.05)])
    # Add title
    title(paste(tmp3[[i]],'\nDown: ',length(down),', Up: ',length(up)))
    allDEUp = c(allDEUp, up)
    allDEDown = c(allDEDown, down)
    allDEUpList[[tmp2[i]]] = up
    allDEDownList[[tmp2[i]]] = down
    # add significance
    for(geneName in 1:length(rownames(deTmp[[i]]))) {
      if(rownames(deTmp[[i]])[geneName] %in% up) {
        deTmp[[i]][geneName,]$DifferentialExpression = 'UP'
      }
      else if(rownames(deTmp[[i]])[geneName] %in% down) {
        deTmp[[i]][geneName,]$DifferentialExpression = 'DOWN'
      }
      else {
        deTmp[[i]][geneName,]$DifferentialExpression = 'NA'
      }
    }

    # Write out differentially expressed transcripts
    write.csv(cbind(data.frame(tmp)[which(rownames(tmp) %in% up),]),paste('../2_DESeq2/results/res.',tmp2[[i]],'_UP.csv',sep=''))
    write.csv(cbind(data.frame(tmp)[which(rownames(tmp) %in% down),]),paste('../2_DESeq2/results/res.',tmp2[[i]],'_DOWN.csv',sep=''))
    # save out all data
    col_order <- c("symbols", "DifferentialExpression", "baseMean",	"log2FoldChange",	"lfcSE",	"stat",	"pvalue",	"padj")
    deTmp2 <- deTmp[[i]][, col_order]
    write.csv(deTmp2, paste0('../2_DESeq2/results/results_',tmp2[[i]],'.csv'))
}
dev.off()



'''Basically you need to test ei vs. ci and take the down-regulated genes which are the miR-497 targets (and the seconday/tertiary/etc.)
Then test ep vs. cp and take the up-regualted genes which are the mRNAs pulled down by bi-miR-497
Then intersect ei_vs_ci_down and ep_vs_cp_up
'''
