import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import hypergeom
from matplotlib_venn import venn2, venn3
import seaborn as sns
import mygene as mg
mj = mg.MyGeneInfo()
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42

######################
## Build translator ##
######################
# Use one results with all genes - common across each results file
data = pd.read_csv('../2_DESeq2/results/results_ei_vs_ci_12H.csv', sep = ',', index_col = 0)
genelist = [i.split('.')[0] for i in list(data.index)]

# Translate
geneSymbols = mj.querymany(genelist, scopes='ensemblgene', fields='symbol', species='human')

# Build translator for easy use
translator = pd.DataFrame(index = genelist, columns = ['symbols'])
for i in geneSymbols:
    for j in i.keys():
        if 'symbol' in i.keys():
            translator.loc[i['query'], 'symbols'] = i['symbol']
        #print(i[j])

translator = translator.dropna()

#######################
## Genes for heatmap ##
#######################
# everything is in ensembl - translate to symbols before you put into figure
hsa_KEGG_genes = pd.read_csv('../5_Figures/NeededInfo/hsa_kegg_anno.csv', sep = ',', index_col = 'ensembl')
#apoptosisGenes_GOBP = list(set(list(pd.read_csv('gene.data.csv', sep = ',', index_col = 'ensembl_gene_id').index)))
apoptosisGenes_GOBP = list(set(list(hsa_KEGG_genes.loc[hsa_KEGG_genes['pathway'] == 'path:hsa04210'].index)))
apoptosisGenes_pathway = list(set(list(pd.read_csv('../5_Figures/NeededInfo/apoptosisPathwayGenes.csv', sep = ',', index_col = 'ensembl_gene_id').index)))

cellCycleGenes_GOBP = list(set(list(pd.read_csv('../5_Figures/NeededInfo/cellCycle_GO0007049.csv', sep = ',', index_col = 'ensembl_gene_id').index)))
cellCycleGenes_pathway = list(set(list(hsa_KEGG_genes.loc[hsa_KEGG_genes['pathway'] == 'path:hsa04110'].index)))

## data
pulldown_12h_DExp_UP = pd.read_csv('../2_DESeq2/results/res.ep_vs_cp_12H_UP.csv', sep = ',', index_col = 0)
pulldown_72h_DExp_UP = pd.read_csv('../2_DESeq2/results/res.ep_vs_cp_72H_UP.csv', sep = ',', index_col = 0)
pulldown_12h_DExp_DOWN = pd.read_csv('../2_DESeq2/results/res.ep_vs_cp_12H_DOWN.csv', sep = ',', index_col = 0)
pulldown_72h_DExp_DOWN = pd.read_csv('../2_DESeq2/results/res.ep_vs_cp_72H_DOWN.csv', sep = ',', index_col = 0)
input_12h_DExp_UP = pd.read_csv('../2_DESeq2/results/res.ei_vs_ci_12H_UP.csv', sep = ',', index_col = 0)
input_72h_DExp_UP = pd.read_csv('../2_DESeq2/results/res.ei_vs_ci_72H_UP.csv', sep = ',', index_col = 0)
input_12h_DExp_DOWN = pd.read_csv('../2_DESeq2/results/res.ei_vs_ci_12H_DOWN.csv', sep = ',', index_col = 0)
input_72h_DExp_DOWN = pd.read_csv('../2_DESeq2/results/res.ei_vs_ci_72H_DOWN.csv', sep = ',', index_col = 0)

enrichmentDictionary = {'pulldown_12h_DExp_UP': pulldown_12h_DExp_UP,'pulldown_12h_DExp_DOWN':pulldown_12h_DExp_DOWN,'pulldown_72h_DExp_UP':pulldown_72h_DExp_UP,'pulldown_72h_DExp_DOWN':pulldown_72h_DExp_DOWN,'input_12h_DExp_UP':input_12h_DExp_UP,'input_12h_DExp_DOWN':input_12h_DExp_DOWN,'input_72h_DExp_UP':input_72h_DExp_UP,'input_72h_DExp_DOWN':input_72h_DExp_DOWN}

#########################################
## Gene expression from prerank inputs ##
#########################################
input_12h = pd.read_csv('../2_DESeq2/results/results_ei_vs_ci_12H.csv', sep = ',', index_col = 0)
input_12h_padj = input_12h.loc[input_12h['padj']<=0.05]
pulldown_12h = pd.read_csv('../2_DESeq2/results/results_ep_vs_cp_12H.csv', sep = ',', index_col = 0)
input_72h = pd.read_csv('../2_DESeq2/results/results_ei_vs_ci_72H.csv', sep = ',', index_col = 0)
input_72h_padj = input_72h.loc[input_72h['padj']<=0.05]
pulldown_72h = pd.read_csv('../2_DESeq2/results/results_ep_vs_cp_72H.csv', sep = ',', index_col = 0)

##########################################
## Build dataframe to turn into heatmap ##
##########################################
## APOPTOSIS
apoptosisHeatmap_GOBP = pd.DataFrame(index = apoptosisGenes_GOBP, columns = ['input_12h', 'input_72h', 'pulldown_12h', 'pulldown_72h'])

for gene in apoptosisHeatmap_GOBP.index:
    if gene in set(list(input_12h_padj.index)+list(input_72h_padj.index)):
        apoptosisHeatmap_GOBP.loc[gene, 'input_12h'] = input_12h.loc[gene, 'log2FoldChange']
        apoptosisHeatmap_GOBP.loc[gene, 'input_72h'] = input_72h.loc[gene, 'log2FoldChange']
        apoptosisHeatmap_GOBP.loc[gene, 'pulldown_12h'] = pulldown_12h.loc[gene, 'log2FoldChange']
        apoptosisHeatmap_GOBP.loc[gene, 'pulldown_72h'] = pulldown_72h.loc[gene, 'log2FoldChange']

apoptosisHeatmap_GOBP.dropna()


# Pathway - reasonable
#apoptosisHeatmap_pathway = pd.DataFrame(index = apoptosisGenes_pathway, columns = ['input_12h', 'pulldown_12h', 'input_72h', 'pulldown_72h'])
#apoptosisHeatmap_pathway = pd.DataFrame(index = apoptosisGenes_pathway, columns = ['input_12h', 'pulldown_12h', 'input_72h', 'pulldown_72h'])
#apoptosisHeatmap_pathway = pd.DataFrame(index = apoptosisGenes_pathway, columns = ['pulldown_12h','pulldown_72h'])
apoptosisHeatmap_pathway = pd.DataFrame(index = apoptosisGenes_pathway, columns = ['input_12h','input_72h'])

for gene in apoptosisHeatmap_pathway.index:
    if gene in input_12h.index:
        apoptosisHeatmap_pathway.loc[gene, 'input_12h'] = input_12h.loc[gene, 'log2FoldChange']
        apoptosisHeatmap_pathway.loc[gene, 'input_72h'] = input_72h.loc[gene, 'log2FoldChange']
        #apoptosisHeatmap_pathway.loc[gene, 'pulldown_12h'] = pulldown_12h.loc[gene, 'log2FoldChange']
        #apoptosisHeatmap_pathway.loc[gene, 'pulldown_72h'] = pulldown_72h.loc[gene, 'log2FoldChange']

apoptosisHeatmap_pathway_trans = apoptosisHeatmap_pathway.dropna()
back2symbols = list(translator.loc[[i for i in apoptosisHeatmap_pathway_trans.index if i in translator.index], 'symbols'])
apoptosisHeatmap_pathway_trans.index = back2symbols


## Heatmap
pp3 = PdfPages('../5_Figures/results/Apoptosis_KEGGPathway.pdf')
plt.figure(figsize=(10, 10))
plt.title('KEGG apoptosis pathway')
sns.clustermap(apoptosisHeatmap_pathway.sort_index().dropna().astype(float), square = True, annot = False, cmap = 'vlag', xticklabels = True, yticklabels = True, vmin = -3, vmax = 3,col_cluster=False)
#sns.heatmap(apoptosisHeatmap_pathway_trans.sort_index().astype(float), square = True, annot = False, cmap = 'vlag', xticklabels = True, yticklabels = True, vmin = -3, vmax = 3)
#plt.tight_layout()
pp3.savefig()
pp3.close()
plt.clf()
plt.close()

# Using rld counts
counts_data = pd.read_csv('../2_DESeq2/results/gexp_norm_all.csv', sep = ',')

# plot apoptosis genes from normalized data
apoptosis_normalizedData = counts_data.loc[list(set(apoptosisGenes_pathway).intersection(counts_data.index))]
back2symbols = list(translator.loc[[i for i in apoptosis_normalizedData.index if i in translator.index], 'symbols'])
apoptosis_normalizedData.index = back2symbols

hours12 = ['ci1_12_hours', 'ci2_12_hours', 'ci3_12_hours', 'cp1_12_hours', 'cp2_12_hours', 'cp3_12_hours', 'ei1_12_hours', 'ei2_12_hours', 'ei3_12_hours', 'ep1_12_hours', 'ep2_12_hours', 'ep3_12_hours']
#hours_input = ['ci1_12_hours', 'ci2_12_hours', 'ci3_12_hours', 'ei1_12_hours', 'ei2_12_hours', 'ei3_12_hours','ci1_72_hours', 'ci2_72_hours', 'ci3_72_hours', 'ei1_72_hours', 'ei2_72_hours', 'ei3_72_hours']
hours_input = ['ei1_12_hours', 'ei2_12_hours', 'ei3_12_hours','ei1_72_hours', 'ei2_72_hours', 'ei3_72_hours']
hours_pulldown = ['cp1_12_hours', 'cp2_12_hours', 'cp3_12_hours', 'ep1_12_hours', 'ep2_12_hours', 'ep3_12_hours','cp1_72_hours', 'cp2_72_hours', 'cp3_72_hours', 'ep1_72_hours', 'ep2_72_hours', 'ep3_72_hours']

indexLabels = apoptosisHeatmap_pathway_trans.loc[apoptosisHeatmap_pathway_trans.abs().max(axis = 1) > 1].sort_index().index

# significance
ugh12 = input_12h.loc[[translator.loc[translator['symbols'] == i].index[0] for i in indexLabels], 'padj'] <= 0.05
translator.loc[ugh12[ugh12 == False].index]

ugh72 = input_72h.loc[[translator.loc[translator['symbols'] == i].index[0] for i in indexLabels], 'padj'] <= 0.05
translator.loc[ugh72[ugh72 == False].index]

## Selected apoptosis genes
pp3 = PdfPages('../5_Figures/results/Apoptosis_KEGGPathway_rld.pdf')
plt.figure(figsize=(10, 40))
plt.title('KEGG apoptosis pathway')
s1 = sns.clustermap(apoptosis_normalizedData.loc[indexLabels, hours_input].dropna().astype(float), annot = False, cmap = 'Greens', xticklabels = True, yticklabels = True, col_cluster=False, row_cluster=True) #,standard_scale = 1)
#sns.clustermap(apoptosis_normalizedData.loc[indexLabels, hours_pulldown].sort_index().dropna().astype(float), annot = False, cmap = 'vlag', xticklabels = True, yticklabels = True, col_cluster=False, row_cluster=False) #,standard_scale = 1)
#sns.heatmap(apoptosisHeatmap_pathway_trans.sort_index().astype(float), square = True, annot = False, cmap = 'vlag', xticklabels = True, yticklabels = True, vmin = -3, vmax = 3)
#plt.tight_layout()
pp3.savefig()
pp3.close()
plt.clf()
plt.close()

# L2FC figure
pp3 = PdfPages('../5_Figures/results/Apoptosis_heatmap.pdf')
#plt.figure(figsize=(10, 1000))
#plt.title('KEGG Cell Cycle pathway')
sns.clustermap(apoptosisHeatmap_pathway_trans.loc[s1.data2d.index].dropna().astype(float), square = True, annot = False, cmap = 'vlag', xticklabels = True, yticklabels = True, figsize = (10,40), col_cluster=False, row_cluster=False, vmin=-2.5, vmax=2.5)# , standard_scale = 1)
#sns.heatmap(apoptosisHeatmap_pathway_trans.sort_index().astype(float), square = True, annot = False, cmap = 'vlag', xticklabels = True, yticklabels = True, vmin = -3, vmax = 3)
#plt.tight_layout()
pp3.savefig()
pp3.close()
plt.clf()
plt.close()



## CELL CYCLE
cellCycleHeatmap_GOBP = pd.DataFrame(index = cellCycleGenes_GOBP, columns = ['input_12h', 'input_72h', 'pulldown_12h', 'pulldown_72h'])

for gene in cellCycleHeatmap_GOBP.index:
    if gene in input_12h.index:
        cellCycleHeatmap_GOBP.loc[gene, 'input_12h'] = input_12h.loc[gene, 'log2FoldChange']
        cellCycleHeatmap_GOBP.loc[gene, 'input_72h'] = input_72h.loc[gene, 'log2FoldChange']
        cellCycleHeatmap_GOBP.loc[gene, 'pulldown_12h'] = pulldown_12h.loc[gene, 'log2FoldChange']
        cellCycleHeatmap_GOBP.loc[gene, 'pulldown_72h'] = pulldown_72h.loc[gene, 'log2FoldChange']

cellCycleHeatmap_GOBP.dropna()

# Pathway - reasonable
cellCycleHeatmap_pathway = pd.DataFrame(index = cellCycleGenes_pathway, columns = ['input_12h','input_72h'])

for gene in cellCycleHeatmap_pathway.index:
    if gene in input_12h.index:
        cellCycleHeatmap_pathway.loc[gene, 'input_12h'] = input_12h.loc[gene, 'log2FoldChange']
        cellCycleHeatmap_pathway.loc[gene, 'input_72h'] = input_72h.loc[gene, 'log2FoldChange']

cellCycleHeatmap_pathway_trans = cellCycleHeatmap_pathway.dropna()
back2symbols = list(translator.loc[[i for i in cellCycleHeatmap_pathway_trans.index if i in translator.index], 'symbols'])
cellCycleHeatmap_pathway_trans.index = back2symbols


## Heatmap
pp3 = PdfPages('../5_Figures/results/cellCycle_heatmap.pdf')
#plt.figure(figsize=(10, 1000))
#plt.title('KEGG Cell Cycle pathway')
sns.clustermap(cellCycleHeatmap_pathway_trans.sort_index().dropna().astype(float), square = True, annot = False, cmap = 'vlag', xticklabels = True, yticklabels = True, vmin = -3, vmax = 3,figsize = (10,40), col_cluster=False)
#sns.heatmap(apoptosisHeatmap_pathway_trans.sort_index().astype(float), square = True, annot = False, cmap = 'vlag', xticklabels = True, yticklabels = True, vmin = -3, vmax = 3)
#plt.tight_layout()
pp3.savefig()
pp3.close()
plt.clf()
plt.close()

pp3 = PdfPages('../5_Figures/results/cellCycle_Pathway.pdf')
#plt.figure(figsize=(10, 1000))
#plt.title('KEGG Cell Cycle pathway')
sns.clustermap(cellCycleHeatmap_pathway_trans.loc[cellCycleHeatmap_pathway_trans.abs().max(axis = 1) > 1].sort_index().dropna().astype(float), square = True, annot = False, cmap = 'vlag', xticklabels = True, yticklabels = True, figsize = (10,40), col_cluster=False) #, standard_scale = 1)
#sns.heatmap(apoptosisHeatmap_pathway_trans.sort_index().astype(float), square = True, annot = False, cmap = 'vlag', xticklabels = True, yticklabels = True, vmin = -3, vmax = 3)
#plt.tight_layout()
pp3.savefig()
pp3.close()
plt.clf()
plt.close()

# significance
interestingGenes_apoptosis = counts_data.loc[list(set(apoptosisGenes_pathway).intersection(counts_data.index))].index
answers_12h_apoptosis = pulldown_12h_DExp_UP.loc[[i for i in interestingGenes_apoptosis if i in pulldown_12h_DExp_UP.index]]
symbol_12h_apoptosis = list(translator.loc[[i for i in answers_12h_apoptosis.index if i in translator.index], 'symbols'])

pulldown_72h_DExp_UP.loc[[i for i in interestingGenes_apoptosis if i in pulldown_72h_DExp_UP.index]]
answers_72h_apoptosis = pulldown_72h_DExp_UP.loc[[i for i in interestingGenes_apoptosis if i in pulldown_72h_DExp_UP.index]]
symbol_72h_apoptosis = list(translator.loc[[i for i in answers_72h_apoptosis.index if i in translator.index], 'symbols'])

interestingGenes_CC = counts_data.loc[list(set(cellCycleGenes_GOBP).intersection(counts_data.index))].index
answers_12h_CC = pulldown_12h_DExp_UP.loc[[i for i in interestingGenes_CC if i in pulldown_12h_DExp_UP.index]]
symbol_12h_CC = list(translator.loc[[i for i in answers_12h_CC.index if i in translator.index], 'symbols'])

pulldown_72h_DExp_UP.loc[[i for i in interestingGenes_apoptosis if i in pulldown_72h_DExp_UP.index]]
answers_72h_CC = pulldown_72h_DExp_UP.loc[[i for i in interestingGenes_CC if i in pulldown_72h_DExp_UP.index]]
symbol_72h_CC = list(translator.loc[[i for i in answers_72h_CC.index if i in translator.index], 'symbols'])



#apoptosisGenes_GOBP
#cellCycleGenes_GOBP

# Using rld counts
# plot apoptosis genes from normalized data
cellCycle_normalizedData = counts_data.loc[list(set(cellCycleGenes_pathway).intersection(counts_data.index))]
back2symbols = list(translator.loc[[i for i in cellCycle_normalizedData.index if i in translator.index], 'symbols'])
cellCycle_normalizedData.index = back2symbols

hours12 = ['ci1_12_hours', 'ci2_12_hours', 'ci3_12_hours', 'cp1_12_hours', 'cp2_12_hours', 'cp3_12_hours', 'ei1_12_hours', 'ei2_12_hours', 'ei3_12_hours', 'ep1_12_hours', 'ep2_12_hours', 'ep3_12_hours']
hours_pulldown = ['cp1_12_hours', 'cp2_12_hours', 'cp3_12_hours', 'ep1_12_hours', 'ep2_12_hours', 'ep3_12_hours','cp1_72_hours', 'cp2_72_hours', 'cp3_72_hours', 'ep1_72_hours', 'ep2_72_hours', 'ep3_72_hours']
hours_input = ['ci1_12_hours', 'ci2_12_hours', 'ci3_12_hours', 'ei1_12_hours', 'ei2_12_hours', 'ei3_12_hours','ci1_72_hours', 'ci2_72_hours', 'ci3_72_hours', 'ei1_72_hours', 'ei2_72_hours', 'ei3_72_hours']

indexLabels1 = cellCycleHeatmap_pathway_trans.loc[cellCycleHeatmap_pathway_trans.abs().max(axis = 1) > 1].sort_index().index

ugh12 = input_12h.loc[[translator.loc[translator['symbols'] == i].index[0] for i in indexLabels1], 'padj'] <= 0.05
translator.loc[ugh12[ugh12 == False].index]

ugh72 = input_72h.loc[[translator.loc[translator['symbols'] == i].index[0] for i in indexLabels1], 'padj'] <= 0.05
translator.loc[ugh72[ugh72 == False].index]


pp3 = PdfPages('../5_Figures/results/cellCycle_KEGGPathway_rld.pdf')
plt.figure(figsize=(10, 40))
plt.title('KEGG Cell Cycle pathway')
sns.clustermap(cellCycle_normalizedData.loc[indexLabels1, hours_input].sort_index().dropna().astype(float), annot = False, cmap = 'vlag', xticklabels = True, yticklabels = True, col_cluster=False, row_cluster=False, figsize=(10,40)) #, standard_scale = 1)
#sns.heatmap(apoptosisHeatmap_pathway_trans.sort_index().astype(float), square = True, annot = False, cmap = 'vlag', xticklabels = True, yticklabels = True, vmin = -3, vmax = 3)
#plt.tight_layout()
pp3.savefig()
pp3.close()
plt.clf()
plt.close()


# apoptosis
pdf = PdfPages('../5_Figures/results/vennDiagram_apoptosisGO0042981.pdf')
venn3([set(pulldown_12h_DExp_UP.index), set(pulldown_72h_DExp_UP.index), set(apoptosisGenes_GOBP)], ('12 hour pulldown up regulated', '72 hour pulldown up regulated', 'GO:0042981'))
pdf.savefig()
plt.clf()

venn3([set(pulldown_12h_DExp_UP.index), set(input_72h_DExp_DOWN.index), set(apoptosisGenes_GOBP)], ('12 hour pulldown up regulated', '72 hour input down regulated', 'GO:0042981'))
pdf.savefig()
plt.clf()
pdf.close()

# cell cycle
pdf = PdfPages('../5_Figures/results/vennDiagram_cellcycleGO0007049.pdf')
venn3([set(pulldown_72h_DExp_UP.index), set(input_72h_DExp_DOWN.index), set(cellCycleGenes_GOBP)], ('72 hour pulldown up regulated', '72 hour input down regulated', 'GO:0007049'))
pdf.savefig()
plt.clf()

venn3([set(pulldown_12h_DExp_UP.index), set(pulldown_72h_DExp_UP.index), set(cellCycleGenes_GOBP)], ('12 hour pulldown up regulated', '72 hour pulldown up regulated', 'GO:0007049'))
pdf.savefig()
plt.clf()
pdf.close()

## TargetScan venn diagrams
miR_497_targetscan = list(pd.read_csv('../1_CollectTargetScan8Info/miR_497_5p_targetscan_V8_0.csv', index_col =0)['Gene_ID_noVer'])
enrichmentScores = {}

pdf = PdfPages('../5_Figures/results/vennDiagram_TargetScan.pdf')
venn3([set(pulldown_12h_DExp_UP.index), set(input_72h_DExp_DOWN.index), set(miR_497_targetscan)], ('12 hour pulldown up regulated', '72 hour input down regulated', 'TargetScan'))
pdf.savefig()
plt.clf()

venn3([set(pulldown_12h_DExp_UP.index), set(pulldown_72h_DExp_UP.index), set(miR_497_targetscan)], ('12 hour pulldown', '72 hour pulldown', 'targetscan_targets'))
pdf.savefig()
plt.clf()

# 12 hour pulldown
venn2([set(pulldown_12h_DExp_UP.index), set(miR_497_targetscan)], ('12 hour pulldown UP', 'TargetScan'))
pdf.savefig()
plt.clf()

# 72 hour pulldown
venn2([set(pulldown_72h_DExp_UP.index), set(miR_497_targetscan)], ('72 hour pulldown UP', 'TargetScan'))
pdf.savefig()
plt.clf()

pdf.close()

for dataSet in enrichmentDictionary:
    M = len(input_12h.index)
    N = len(set(input_12h.index).intersection(enrichmentDictionary[dataSet].index))
    n = len(set(input_12h.index).intersection(miR_497_targetscan))
    x = len(set(input_12h.index).intersection(set(enrichmentDictionary[dataSet].index).intersection(miR_497_targetscan)))
    freq = float(x)/float(N)
    pv = 1-hypergeom.cdf(x,M,n,N)
    enrichmentScores[dataSet + ' X ' + 'TargetScan'] = {'Overlap (x)':x, 'All genes (M)':M, 'set 1 (n)':n, 'set 2 (N)':N, 'freq (x/N)':freq, 'p-value':pv}

pd.DataFrame(enrichmentScores).T.to_csv('../5_Figures/results/enrichmentScores_TargetScan.csv')


# Comparisons with hypergeometrics
enrichmentScores = {}

pdf = PdfPages('../5_Figures/results/vennDiagrams.pdf')
plt.figure(figsize=(10, 10))
# 12 hours
venn2([set(pulldown_12h_DExp_UP.index), set(input_12h_DExp_DOWN.index)], ('12 hour pulldown UP', '12 hour input DOWN'))
pdf.savefig()
plt.clf()

M = len(input_12h.index)
N = len(set(input_12h.index).intersection(pulldown_12h_DExp_UP.index))
n = len(set(input_12h.index).intersection(input_12h_DExp_DOWN.index))
x = len(set(input_12h.index).intersection(set(pulldown_12h_DExp_UP.index).intersection(input_12h_DExp_DOWN.index)))
freq = float(x)/float(N)
pv = 1-hypergeom.cdf(x,M,n,N)
enrichmentScores['12 hour pulldown UP' + ' X ' + '12 hour input DOWN'] = {'Overlap (x)':x, 'All genes (M)':M, 'set 1 (n)':n, 'set 2 (N)':N, 'freq (x/N)':freq, 'p-value':pv}


# 72 hours
venn2([set(pulldown_72h_DExp_UP.index), set(input_72h_DExp_DOWN.index)], ('72 hour pulldown UP', '72 hour input DOWN'))
pdf.savefig()
plt.clf()

M = len(input_12h.index)
N = len(set(input_12h.index).intersection(pulldown_72h_DExp_UP.index))
n = len(set(input_12h.index).intersection(input_72h_DExp_DOWN.index))
x = len(set(input_12h.index).intersection(set(pulldown_72h_DExp_UP.index).intersection(input_72h_DExp_DOWN.index)))
freq = float(x)/float(N)
pv = 1-hypergeom.cdf(x,M,n,N)
enrichmentScores['72 hour pulldown UP' + ' X ' + '72 hour input DOWN'] = {'Overlap (x)':x, 'All genes (M)':M, 'set 1 (n)':n, 'set 2 (N)':N, 'freq (x/N)':freq, 'p-value':pv}


# 12 x 72 hours
venn2([set(pulldown_12h_DExp_UP.index), set(input_72h_DExp_DOWN.index)], ('12 hour pulldown UP', '72 hour input DOWN'))
pdf.savefig()
plt.clf()

M = len(input_12h.index)
N = len(set(input_12h.index).intersection(pulldown_12h_DExp_UP.index))
n = len(set(input_12h.index).intersection(input_72h_DExp_DOWN.index))
x = len(set(input_12h.index).intersection(set(pulldown_12h_DExp_UP.index).intersection(input_72h_DExp_DOWN.index)))
freq = float(x)/float(N)
pv = 1-hypergeom.cdf(x,M,n,N)
enrichmentScores['12 hour pulldown UP' + ' X ' + '72 hour input DOWN'] = {'Overlap (x)':x, 'All genes (M)':M, 'set 1 (n)':n, 'set 2 (N)':N, 'freq (x/N)':freq, 'p-value':pv}

# inputs
venn2([set(input_12h_DExp_DOWN.index), set(input_72h_DExp_DOWN.index)], ('12 hour input DOWN', '72 hour input DOWN'))
pdf.savefig()
plt.clf()

M = len(input_12h.index)
N = len(set(input_12h.index).intersection(input_12h_DExp_DOWN.index))
n = len(set(input_12h.index).intersection(input_72h_DExp_DOWN.index))
x = len(set(input_12h.index).intersection(set(input_12h_DExp_DOWN.index).intersection(input_72h_DExp_DOWN.index)))
freq = float(x)/float(N)
pv = 1-hypergeom.cdf(x,M,n,N)
enrichmentScores['12 hour input DOWN' + ' X ' + '72 hour input DOWN'] = {'Overlap (x)':x, 'All genes (M)':M, 'set 1 (n)':n, 'set 2 (N)':N, 'freq (x/N)':freq, 'p-value':pv}

# pulldowns
venn2([set(pulldown_12h_DExp_UP.index), set(pulldown_72h_DExp_UP.index)], ('12 hour pulldown UP', '72 hour pulldown UP'))
pdf.savefig()
plt.clf()

M = len(input_12h.index)
N = len(set(input_12h.index).intersection(pulldown_12h_DExp_UP.index))
n = len(set(input_12h.index).intersection(pulldown_72h_DExp_UP.index))
x = len(set(input_12h.index).intersection(set(pulldown_12h_DExp_UP.index).intersection(pulldown_72h_DExp_UP.index)))
freq = float(x)/float(N)
pv = 1-hypergeom.cdf(x,M,n,N)
enrichmentScores['12 hour pulldown UP' + ' X ' + '72 hour pulldown UP'] = {'Overlap (x)':x, 'All genes (M)':M, 'set 1 (n)':n, 'set 2 (N)':N, 'freq (x/N)':freq, 'p-value':pv}

pdf.close()
pd.DataFrame(enrichmentScores).T.to_csv('../5_Figures/results/enrichmentScores.csv')

## overlaps with pathways

geneOverlap = {}
geneOverlap['12H_pulldown X cell cycle'] = list(translator.loc[list(set(pulldown_12h_DExp_UP.index).intersection(set(cellCycleGenes_pathway))), 'symbols'])
geneOverlap['72H_pulldown X cell cycle'] = list(translator.loc[list(set(pulldown_72h_DExp_UP.index).intersection(set(cellCycleGenes_pathway))), 'symbols'])
geneOverlap['12H_pulldown X apoptosis'] = list(translator.loc[list(set(pulldown_12h_DExp_UP.index).intersection(set(apoptosisGenes_pathway))), 'symbols'])
geneOverlap['72H_pulldown X apoptosis'] = list(translator.loc[list(set(pulldown_72h_DExp_UP.index).intersection(set(apoptosisGenes_pathway))), 'symbols'])

pd.DataFrame({'top':{i:' '.join(geneOverlap[i]) for i in geneOverlap}}).to_csv('../5_Figures/results/geneOverlap_pathwayTerms.csv')

# Venn Diagrams
enrichmentScores = {}
pdf = PdfPages('../5_Figures/results/vennDiagrams_PathwayOverlap.pdf')
plt.figure(figsize=(5, 5))

venn2([set(pulldown_12h_DExp_UP.index), set(cellCycleGenes_pathway)], ('12 hour pulldown UP', 'Cell cycle pathway'))
pdf.savefig()
plt.clf()
M = len(input_12h.index)
N = len(set(input_12h.index).intersection(pulldown_12h_DExp_UP.index))
n = len(set(input_12h.index).intersection(cellCycleGenes_pathway))
x = len(set(input_12h.index).intersection(set(pulldown_12h_DExp_UP.index).intersection(cellCycleGenes_pathway)))
freq = float(x)/float(N)
pv = 1-hypergeom.cdf(x,M,n,N)
enrichmentScores['12 hour pulldown UP' + ' X ' + 'CellCycle pathway'] = {'Overlap (x)':x, 'All genes (M)':M, 'Pathway (n)':n, 'set 2 (N)':N, 'freq (x/N)':freq, 'p-value':pv}


venn2([set(pulldown_72h_DExp_UP.index), set(cellCycleGenes_pathway)], ('72 hour pulldown UP', 'Cell cycle pathway'))
pdf.savefig()
plt.clf()
M = len(input_12h.index)
N = len(set(input_12h.index).intersection(pulldown_72h_DExp_UP.index))
n = len(set(input_12h.index).intersection(cellCycleGenes_pathway))
x = len(set(input_12h.index).intersection(set(pulldown_72h_DExp_UP.index).intersection(cellCycleGenes_pathway)))
freq = float(x)/float(N)
pv = 1-hypergeom.cdf(x,M,n,N)
enrichmentScores['72 hour pulldown UP' + ' X ' + 'CellCycle pathway'] = {'Overlap (x)':x, 'All genes (M)':M, 'Pathway (n)':n, 'set 2 (N)':N, 'freq (x/N)':freq, 'p-value':pv}


venn2([set(pulldown_12h_DExp_UP.index), set(apoptosisGenes_pathway)], ('12 hour pulldown UP', 'Apoptosis pathway'))
pdf.savefig()
plt.clf()
M = len(input_12h.index)
N = len(set(input_12h.index).intersection(pulldown_12h_DExp_UP.index))
n = len(set(input_12h.index).intersection(apoptosisGenes_pathway))
x = len(set(input_12h.index).intersection(set(pulldown_12h_DExp_UP.index).intersection(apoptosisGenes_pathway)))
freq = float(x)/float(N)
pv = 1-hypergeom.cdf(x,M,n,N)
enrichmentScores['12 hour pulldown UP' + ' X ' + 'Apoptosis pathway'] = {'Overlap (x)':x, 'All genes (M)':M, 'Pathway (n)':n, 'set 2 (N)':N, 'freq (x/N)':freq, 'p-value':pv}


venn2([set(pulldown_72h_DExp_UP.index), set(apoptosisGenes_pathway)], ('72 hour pulldown UP', 'Apoptosis pathway'))
pdf.savefig()
plt.clf()
M = len(input_12h.index)
N = len(set(input_12h.index).intersection(pulldown_72h_DExp_UP.index))
n = len(set(input_12h.index).intersection(apoptosisGenes_pathway))
x = len(set(input_12h.index).intersection(set(pulldown_72h_DExp_UP.index).intersection(apoptosisGenes_pathway)))
freq = float(x)/float(N)
pv = 1-hypergeom.cdf(x,M,n,N)
enrichmentScores['72 hour pulldown UP' + ' X ' + 'Apoptosis pathway'] = {'Overlap (x)':x, 'All genes (M)':M, 'Pathway (n)':n, 'set 2 (N)':N, 'freq (x/N)':freq, 'p-value':pv}

pdf.close()
pd.DataFrame(enrichmentScores).T.to_csv('../5_Figures/results/enrichmentScores_pathways.csv')


## save out overlaps
# DisGenet Genes
disgenet = list(pd.read_csv('../2_DESeq2/DisGeNET_databases/DisGeNET_genes.csv', index_col = 0, header = 0)['Genes'])

# 12 hour bound targets x 72 hour down regulated genes
save_12hBound_x_72hRegulated = pd.DataFrame(translator.loc[[i for i in set(pulldown_12h_DExp_UP.index).intersection(input_72h_DExp_DOWN.index) if i in translator.index]])
for i in save_12hBound_x_72hRegulated.index:
    symbol = save_12hBound_x_72hRegulated.loc[i, 'symbols']
    if symbol in disgenet:
        save_12hBound_x_72hRegulated.loc[i, 'DisGeNET'] = 'yes'
    else:
        save_12hBound_x_72hRegulated.loc[i, 'DisGeNET'] = 'no'
    if i in miR_497_targetscan:
        save_12hBound_x_72hRegulated.loc[i, 'TargetScan'] = 'yes'
    else:
        save_12hBound_x_72hRegulated.loc[i, 'TargetScan'] = 'no'

save_12hBound_x_72hRegulated.to_csv('../5_Figures/results/12hBound_x_72hRegulated.csv')

# 12 hour bound targets x targetscan
save_12hBound_x_targetscan = pd.DataFrame(translator.loc[[i for i in set(pulldown_12h_DExp_UP.index).intersection(miR_497_targetscan) if i in translator.index]])
for i in save_12hBound_x_targetscan.index:
    symbol = save_12hBound_x_targetscan.loc[i, 'symbols']
    if symbol in disgenet:
        save_12hBound_x_targetscan.loc[i, 'DisGeNET'] = 'yes'
    else:
        save_12hBound_x_targetscan.loc[i, 'DisGeNET'] = 'no'
    if i in miR_497_targetscan:
        save_12hBound_x_targetscan.loc[i, 'TargetScan'] = 'yes'
    else:
        save_12hBound_x_targetscan.loc[i, 'TargetScan'] = 'no'

save_12hBound_x_targetscan.to_csv('../5_Figures/results/12hBound_x_targetscan.csv')

# 72 hour down regulated genes x targetscan
save_72hRegulated_x_targetscan = pd.DataFrame(translator.loc[[i for i in set(input_72h_DExp_DOWN.index).intersection(miR_497_targetscan) if i in translator.index]])
for i in save_72hRegulated_x_targetscan.index:
    symbol = save_72hRegulated_x_targetscan.loc[i, 'symbols']
    if symbol in disgenet:
        save_72hRegulated_x_targetscan.loc[i, 'DisGeNET'] = 'yes'
    else:
        save_72hRegulated_x_targetscan.loc[i, 'DisGeNET'] = 'no'
    if i in miR_497_targetscan:
        save_72hRegulated_x_targetscan.loc[i, 'TargetScan'] = 'yes'
    else:
        save_72hRegulated_x_targetscan.loc[i, 'TargetScan'] = 'no'

save_72hRegulated_x_targetscan.to_csv('../5_Figures/results/72hRegulated_x_targetscan.csv')

# 12 hour bound targets x 72 hour down regulated genes x targetscan
save_12hBound_x_72hRegulated_x_targetscan = pd.DataFrame(translator.loc[[i for i in set(set(pulldown_12h_DExp_UP.index).intersection(input_72h_DExp_DOWN.index)).intersection(miR_497_targetscan) if i in translator.index]])
for i in save_12hBound_x_72hRegulated_x_targetscan.index:
    symbol = save_12hBound_x_72hRegulated_x_targetscan.loc[i, 'symbols']
    if symbol in disgenet:
        save_12hBound_x_72hRegulated_x_targetscan.loc[i, 'DisGeNET'] = 'yes'
    else:
        save_12hBound_x_72hRegulated_x_targetscan.loc[i, 'DisGeNET'] = 'no'
    if i in miR_497_targetscan:
        save_12hBound_x_72hRegulated_x_targetscan.loc[i, 'TargetScan'] = 'yes'
    else:
        save_12hBound_x_72hRegulated_x_targetscan.loc[i, 'TargetScan'] = 'no'

save_12hBound_x_72hRegulated_x_targetscan.to_csv('../5_Figures/results/12hBound_x_72hRegulated_x_targetscan.csv')


## 330 overlap
overlap330_33 = list(set(pulldown_12h_DExp_UP.index).intersection(set(input_72h_DExp_DOWN.index).intersection(set(miR_497_targetscan))))
overlap330_330 = list(set(pulldown_12h_DExp_UP.index).intersection(set(input_72h_DExp_DOWN.index)))

# add some genes
translator.loc['ENSG00000281383','symbols'] = 'JD028642'
translator.loc['ENSG00000280351','symbols'] = 'CTD-2561B21.8'
translator.loc['ENSG00000260121','symbols'] = 'RP5-1142A6.9'
translator.loc['ENSG00000196656','symbols'] = 'RPS26L1'

save_overlap330 = translator.loc[list(set(overlap330_330).intersection(translator.index))]
for i in save_overlap330.index:
    symbol = save_overlap330.loc[i, 'symbols']
    if symbol in disgenet:
        save_overlap330.loc[i, 'DisGeNET'] = 'yes'
    else:
        save_overlap330.loc[i, 'DisGeNET'] = 'no'
    if i in miR_497_targetscan:
        save_overlap330.loc[i, 'TargetScan'] = 'yes'
    else:
        save_overlap330.loc[i, 'TargetScan'] = 'no'

save_overlap330.to_csv('../5_Figures/results/330genesInOverlap.csv')

#finalResults330 = overlap330.copy()

pulldown_12h.loc[overlap330_330.index].to_csv('../5_Figures/results/12hPulldown_330genesToSend.csv')
input_72h.loc[overlap330_330.index].to_csv('../5_Figures/results/72hinput_330genesToSend.csv')

'''
enrichr330 = gp.enrichr(gene_list=list(overlap330['SYMBOLS']), gene_sets=['KEGG_2016'], organism='Human', outdir='test/enrichr_kegg', cutoff=0.05, no_plot=False)
enr_results = enrichr330.results[enrichr330.results['Adjusted P-value']<=0.05]
enr_results_nonsig = enrichr330.results
enr_results.to_csv('FEA_330overlap.csv')
enr_results_nonsig.to_csv('FEA_330overlap_nonsig.csv')
'''


## make table of survival from paper data
# index are SYMBOLS
survivalData = pd.read_csv('../5_Figures/NeededInfo/mRNA_survival_PMID35354049.csv', index_col = 0)

# Make table
describe330genes = pd.DataFrame(index = list(input_72h.index), columns = ['SYMBOL ID', '12h pulldown L2FC', '12h pulldown adj pval', '72h input L2FC', '72h input adj pval', 'Z.uni', 'p.uni', 'Z.multi', 'p.multi', 'TargetScan miRNA?'])

# sub dataframes you want
symbols = translator.loc[overlap330_330, 'symbols']

L2FC_12hPulldown = pulldown_12h.loc[overlap330_330, 'log2FoldChange']
adjpval_12hPulldown = pulldown_12h.loc[overlap330_330, 'padj']

L2FC_72hInput = input_72h.loc[overlap330_330, 'log2FoldChange']
adjpval_72hInput = input_72h.loc[overlap330_330, 'padj']

table_iter1 = pd.concat([symbols, L2FC_12hPulldown, adjpval_12hPulldown, L2FC_72hInput, adjpval_72hInput], axis = 1)
table_iter2 = table_iter1.reset_index()
table_iter2.columns = ['EMSEMBL ID', 'SYMBOL ID', '12h pulldown up L2FC', '12h pulldown up adj pval', '72h input down L2FC', '72h input down adj pval']
table_iter3 = table_iter2.set_index('SYMBOL ID')

survivalOfThe330 = survivalData.loc[[i for i in list(symbols) if i in survivalData.index]]

table_iter4 = pd.concat([survivalOfThe330, table_iter3], axis = 1)

table_iter4['TargetScan miRNA?'] = np.nan
for i in table_iter4['EMSEMBL ID']:
    if i in list(miR_497_targetscan):
        table_iter4.loc[translator.loc[i, 'symbols'], 'TargetScan miRNA?'] = 'yes'
    else:
        table_iter4.loc[translator.loc[i, 'symbols'], 'TargetScan miRNA?'] = 'no'

table_iter5 = table_iter4.loc[:,['EMSEMBL ID','12h pulldown up L2FC', '12h pulldown up adj pval', '72h input down L2FC', '72h input down adj pval','UV COXPH', 'UV COXPH PV', 'MV COXPH', 'MV XOXPH PV', 'TargetScan miRNA?']]
table_iter5.to_csv('../5_Figures/results/TableOf330Genes.csv')

# fin
