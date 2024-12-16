#############
## Imports ##
#############

import pandas as pd
import json
import numpy as np
import mygene as mg
mj = mg.MyGeneInfo()
from scipy import stats
from scipy.stats import hypergeom
import statsmodels.formula.api as sm
import statsmodels
import math
import re
import pickle
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
from matplotlib_venn import venn2, venn3
import seaborn as sns
import gseapy as gp
import gzip
import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.4.0" # change as needed
#os.environ["R_HOME"] = r"C:\Users\swilferd.ASURITE\AppData\Local\Programs\R\R-4.2.1" # for D drive
import rpy2
import rpy2.robjects as robj
from rpy2.robjects import FloatVector, IntVector, StrVector
from rpy2 import rinterface
from matplotlib_venn import venn2, venn3

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

#####################
## import GO terms ##
#####################
gene2GOfile = '../3_FunctionalEnrichmentAnalyses/gene2go.gz'
gene2go = {}
with gzip.open(gene2GOfile, 'rt') as inFile:
    #gene2go = inFile.read()
    inFile.readline()
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split('\t')
        if splitUp[0] == '9606':
            if not splitUp[2] in gene2go:
                gene2go[splitUp[2]] = []
            gene2go[splitUp[2]].append(splitUp[1])

##################
## pull in data ##
##################
pulldown_12hour_UP = pd.read_csv('../2_DESeq2/results/res.ep_vs_cp_12H_UP.csv', index_col = 0) # 12 hour pulldown is up
pulldown_12hour_DOWN = pd.read_csv('../2_DESeq2/results/res.ep_vs_cp_12H_DOWN.csv', index_col = 0)

pulldown_72hour_UP = pd.read_csv('../2_DESeq2/results/res.ep_vs_cp_72H_UP.csv', index_col = 0) # 72 hour pulldown is up
pulldown_72hour_DOWN = pd.read_csv('../2_DESeq2/results/res.ep_vs_cp_72H_DOWN.csv', index_col = 0)

input_12hour_UP = pd.read_csv('../2_DESeq2/results/res.ei_vs_ci_12H_UP.csv', index_col = 0)
input_12hour_DOWN = pd.read_csv('../2_DESeq2/results/res.ei_vs_ci_12H_DOWN.csv', index_col = 0) # 12 hour input is down

input_72hour_UP = pd.read_csv('../2_DESeq2/results/res.ei_vs_ci_72H_UP.csv', index_col = 0)
input_72hour_DOWN = pd.read_csv('../2_DESeq2/results/res.ei_vs_ci_72H_DOWN.csv', index_col = 0) # 72 hour input is down

#############################
## Format data for enrichr ##
#############################
## fill out directed DESeq2 results
directed_miR497data = {'ei_vs_ci_12H':{'down':[], 'up':[]}, 'ep_vs_cp_12H':{'down':[], 'up':[]}, 'ei_vs_ci_72H':{'down':[], 'up':[]}, 'ep_vs_cp_72H':{'down':[], 'up':[]}}

for sampleID in directed_miR497data.keys():
    directed_miR497data[sampleID]['up'] = list(translator.loc[[i for i in pd.read_csv('../2_DESeq2/results/res.' + sampleID + '_UP.csv', sep = ',', index_col = 0).index if i in translator.index], 'symbols'])
    directed_miR497data[sampleID]['down'] = list(translator.loc[[i for i in pd.read_csv('../2_DESeq2/results/res.' + sampleID + '_DOWN.csv', sep = ',', index_col = 0).index if i in translator.index], 'symbols'])

## fill out for the union and intersections
additional_directed_miR497data_unions = {'12h_inputXpulldown':{'down':[], 'up':[]}, '72h_inputXpulldown':{'down':[], 'up':[]}, 'input_12hX72h':{'down':[], 'up':[]}, 'pulldown_12hX72h':{'down':[], 'up':[]}}
additional_directed_miR497data_intersections = {'12h_inputXpulldown':{'down':[], 'up':[]}, '72h_inputXpulldown':{'down':[], 'up':[]}, 'input_12hX72h':{'down':[], 'up':[]}, 'pulldown_12hX72h':{'down':[], 'up':[]}}

for direction in ['up','down']:
    # unions
    additional_directed_miR497data_unions['12h_inputXpulldown'][direction] = list(set(directed_miR497data['ei_vs_ci_12H'][direction]).union(directed_miR497data['ep_vs_cp_12H'][direction]))
    additional_directed_miR497data_unions['72h_inputXpulldown'][direction] = list(set(directed_miR497data['ei_vs_ci_72H'][direction]).union(directed_miR497data['ep_vs_cp_72H'][direction]))
    additional_directed_miR497data_unions['input_12hX72h'][direction] = list(set(directed_miR497data['ei_vs_ci_12H'][direction]).union(directed_miR497data['ei_vs_ci_72H'][direction]))
    additional_directed_miR497data_unions['pulldown_12hX72h'][direction] = list(set(directed_miR497data['ep_vs_cp_12H'][direction]).union(directed_miR497data['ep_vs_cp_72H'][direction]))
    # intersections
    additional_directed_miR497data_intersections['12h_inputXpulldown'][direction] = list(set(directed_miR497data['ei_vs_ci_12H'][direction]).intersection(directed_miR497data['ep_vs_cp_12H'][direction]))
    additional_directed_miR497data_intersections['72h_inputXpulldown'][direction] = list(set(directed_miR497data['ei_vs_ci_72H'][direction]).intersection(directed_miR497data['ep_vs_cp_72H'][direction]))
    additional_directed_miR497data_intersections['input_12hX72h'][direction] = list(set(directed_miR497data['ei_vs_ci_12H'][direction]).intersection(directed_miR497data['ei_vs_ci_72H'][direction]))
    additional_directed_miR497data_intersections['pulldown_12hX72h'][direction] = list(set(directed_miR497data['ep_vs_cp_12H'][direction]).intersection(directed_miR497data['ep_vs_cp_72H'][direction]))

## Dictionary of dictionaries of all the information
allInformation = {'individual' : directed_miR497data, 'unions' : additional_directed_miR497data_unions, 'intersections' : additional_directed_miR497data_intersections}

#############################
## Format data for enrichr ##
#############################
## fill out directed DESeq2 results
directed_miR497data = {'ei_vs_ci_12H':{'down':[], 'up':[]}, 'ep_vs_cp_12H':{'down':[], 'up':[]}, 'ei_vs_ci_72H':{'down':[], 'up':[]}, 'ep_vs_cp_72H':{'down':[], 'up':[]}}

for sampleID in directed_miR497data.keys():
    directed_miR497data[sampleID]['up'] = list(translator.loc[[i for i in pd.read_csv('../2_DESeq2/results/res.' + sampleID + '_UP.csv', sep = ',', index_col = 0).index if i in translator.index], 'symbols'])
    directed_miR497data[sampleID]['down'] = list(translator.loc[[i for i in pd.read_csv('../2_DESeq2/results/res.' + sampleID + '_DOWN.csv', sep = ',', index_col = 0).index if i in translator.index], 'symbols'])

## fill out for the union and intersections
additional_directed_miR497data_unions = {'12h_inputXpulldown':{'down':[], 'up':[]}, '72h_inputXpulldown':{'down':[], 'up':[]}, 'input_12hX72h':{'down':[], 'up':[]}, 'pulldown_12hX72h':{'down':[], 'up':[]}}
additional_directed_miR497data_intersections = {'12h_inputXpulldown':{'down':[], 'up':[]}, '72h_inputXpulldown':{'down':[], 'up':[]}, 'input_12hX72h':{'down':[], 'up':[]}, 'pulldown_12hX72h':{'down':[], 'up':[]}}

for direction in ['up','down']:
    # unions
    additional_directed_miR497data_unions['12h_inputXpulldown'][direction] = list(set(directed_miR497data['ei_vs_ci_12H'][direction]).union(directed_miR497data['ep_vs_cp_12H'][direction]))
    additional_directed_miR497data_unions['72h_inputXpulldown'][direction] = list(set(directed_miR497data['ei_vs_ci_72H'][direction]).union(directed_miR497data['ep_vs_cp_72H'][direction]))
    additional_directed_miR497data_unions['input_12hX72h'][direction] = list(set(directed_miR497data['ei_vs_ci_12H'][direction]).union(directed_miR497data['ei_vs_ci_72H'][direction]))
    additional_directed_miR497data_unions['pulldown_12hX72h'][direction] = list(set(directed_miR497data['ep_vs_cp_12H'][direction]).union(directed_miR497data['ep_vs_cp_72H'][direction]))
    # intersections
    additional_directed_miR497data_intersections['12h_inputXpulldown'][direction] = list(set(directed_miR497data['ei_vs_ci_12H'][direction]).intersection(directed_miR497data['ep_vs_cp_12H'][direction]))
    additional_directed_miR497data_intersections['72h_inputXpulldown'][direction] = list(set(directed_miR497data['ei_vs_ci_72H'][direction]).intersection(directed_miR497data['ep_vs_cp_72H'][direction]))
    additional_directed_miR497data_intersections['input_12hX72h'][direction] = list(set(directed_miR497data['ei_vs_ci_12H'][direction]).intersection(directed_miR497data['ei_vs_ci_72H'][direction]))
    additional_directed_miR497data_intersections['pulldown_12hX72h'][direction] = list(set(directed_miR497data['ep_vs_cp_12H'][direction]).intersection(directed_miR497data['ep_vs_cp_72H'][direction]))

## Dictionary of dictionaries of all the information
allInformation = {'individual' : directed_miR497data, 'unions' : additional_directed_miR497data_unions, 'intersections' : additional_directed_miR497data_intersections}

#################
## Run enrichr ##
#################
for informationType in allInformation.keys():
    print(informationType)
    for sampleID in allInformation[informationType]:
        print(sampleID)
        for direction in ['up', 'down']:
            print(direction)
            #for geneSet in ['GO_Biological_Process_2021', 'KEGG_2021_Human', 'Elsevier_Pathway_Collection', 'WikiPathway_2021_Human', 'MSigDB_Hallmark_2020', 'MSigDB_Oncogenic_Signatures']:
            for geneSet in ['GO_Biological_Process_2021', 'KEGG_2021_Human', 'WikiPathway_2021_Human']:
                print(geneSet)
                GO_results = gp.enrichr(gene_list = allInformation[informationType][sampleID][direction], gene_sets=geneSet, organism='Human', outdir='../3_FunctionalEnrichmentAnalyses/results/enrichr/test/enrichr_' + geneSet + '_' + sampleID + '_' + informationType, cutoff=0.05, no_plot=False)
                GO_results.results.to_csv('../3_FunctionalEnrichmentAnalyses/results/enrichr/'+ sampleID + '_' + direction + '_' + geneSet + '_' + informationType + '.csv')

## pull out significant pathways
stuffWithPathways = {}
stuffWithPathwaysKeys = []
for informationType in allInformation.keys():
    print(informationType)
    for sampleID in allInformation[informationType]:
        print(sampleID)
        for direction in ['up', 'down']:
            print(direction)
            #for geneSet in ['GO_Biological_Process_2021', 'KEGG_2021_Human', 'Elsevier_Pathway_Collection', 'WikiPathway_2021_Human', 'MSigDB_Hallmark_2020', 'MSigDB_Oncogenic_Signatures']:
            for geneSet in ['GO_Biological_Process_2021', 'KEGG_2021_Human', 'WikiPathway_2021_Human']:
                print(geneSet)
                dataSet = pd.read_csv('../3_FunctionalEnrichmentAnalyses/results/enrichr/'+ sampleID + '_' + direction + '_' + geneSet + '_' + informationType + '.csv', index_col ='Term')
                #stuffWithPathways[sampleID + '_' + direction + '_' + geneSet + '_' + informationType] = []
                for term in dataSet.index:
                    if dataSet.loc[term, 'Adjusted P-value'] <= 0.05:
                        if sampleID + '_' + direction + '_' + geneSet + '_' + informationType not in stuffWithPathwaysKeys:
                            stuffWithPathwaysKeys.append(sampleID + '_' + direction + '_' + geneSet + '_' + informationType)
                        #stuffWithPathways[sampleID + '_' + direction + '_' + geneSet + '_' + informationType] += (term)

for filename in stuffWithPathwaysKeys:
    dataSet = pd.read_csv('../3_FunctionalEnrichmentAnalyses/results/enrichr/'+ filename + '.csv', index_col ='Term')
    dataSet.to_csv('../3_FunctionalEnrichmentAnalyses/results/enrichr/significant/'+ filename + '.csv')

#############################
## Format data for prerank ##
#############################
# declare samples
miR497data = ['results_ei_vs_ci_12H', 'results_ep_vs_cp_12H', 'results_ei_vs_ci_72H', 'results_ep_vs_cp_72H']

for sampleID in miR497data:
    prerank_format = pd.read_csv('../2_DESeq2/results/' + sampleID + '.csv', sep = ',', index_col = 0)
    prerank_saveout = {}
    for ensegene in prerank_format.index:
        if ensegene in translator.index:
            prerank_saveout[translator.loc[ensegene.split('.')[0],'symbols']] = float(prerank_format.loc[ensegene, 'log2FoldChange'])
    pd.DataFrame.from_dict(prerank_saveout, orient = 'index').sort_values(by = 0, ascending = False).dropna().to_csv('../3_FunctionalEnrichmentAnalyses/prerank_inputs/' + sampleID + '_PRERANK.csv', sep = ',')

#################
## Run prerank ##
#################
for sampleID in miR497data:
    testRun = pd.read_csv('../3_FunctionalEnrichmentAnalyses/prerank_inputs/' + sampleID + '_PRERANK.csv', sep = ',')
    preranked_results = gp.prerank(rnk=testRun, gene_sets='GO_Biological_Process_2021', min_size=3, max_size=1000, permutation_num=1000, threads = 20, outdir= '../3_FunctionalEnrichmentAnalyses/results/prerank/' + sampleID, seed=6, verbose=True)

#######################################
## Import and format prerank results ##
#######################################
miR497data = ['results_ei_vs_ci_12H', 'results_ep_vs_cp_12H', 'results_ei_vs_ci_72H', 'results_ep_vs_cp_72H']

## Format FDR correction
## Format prerank results - make your own FDR
for sampleID in miR497data:
    if os.path.exists('../3_FunctionalEnrichmentAnalyses/results/prerank/' + sampleID):
        resultsReport = pd.read_csv('../3_FunctionalEnrichmentAnalyses/results/prerank/' + sampleID + '/gseapy.gene_set.prerank.report.csv', index_col = 'Term')
        adjustedPvals = statsmodels.stats.multitest.multipletests(resultsReport['NOM p-val'], alpha=0.05, method='fdr_bh')
        resultsReport['FDR:BH Adjusted p-values'] = adjustedPvals[1]
        resultsReport.to_csv('../3_FunctionalEnrichmentAnalyses/results/prerank/' + sampleID + '/' + sampleID + '_prerankOutput_FDRcorrected.csv', sep = ',')

# pull out prerank results with FDR <= 0.05
for sampleID in miR497data:
    if os.path.exists('../3_FunctionalEnrichmentAnalyses/results/prerank/' + sampleID):
        resultsReport = pd.read_csv('../3_FunctionalEnrichmentAnalyses/results/prerank/' + sampleID + '/' + sampleID + '_prerankOutput_FDRcorrected.csv', index_col = 'Term')
        passData = resultsReport.loc[resultsReport['FDR:BH Adjusted p-values']<=0.05]
        passData.to_csv('../3_FunctionalEnrichmentAnalyses/results/prerank/' + sampleID + '/' + sampleID + '_prerankOutput_FDRcorrected_passed005.csv', sep = ',')
        print(len(passData.index))

#######################################
## Run prerank cell death cell cycle ##
#######################################
## comment out if runing whole script
robj.r("""
    library('GOSemSim')
    hsGO = godata('org.Hs.eg.db',ont='BP')
    """)

## Define the semantic similarity function
def semSim_Mechanisms(term1):
    semSimR = robj.r("""
    semSimR <- function(term1) {
        pyroptosis = mgoSim(c(term1), c('GO:0070269'), semData=hsGO, measure='Jiang', combine='max')
        mitoptosis = mgoSim(c(term1), c('GO:0000423'), semData=hsGO, measure='Jiang', combine='max')
        necrosis = mgoSim(c(term1), c('GO:0060544'), semData=hsGO, measure='Jiang', combine='max')
        autophagy = mgoSim(c(term1), c('GO:1903146', 'GO:0044804', 'GO:0016239', 'GO:0030242', 'GO:1904714', 'GO:0016242', 'GO:1904925', 'GO:0061912'), semData=hsGO, measure='Jiang', combine='max')
        cellDeath = mgoSim(c(term1), c('GO:0010941'), semData=hsGO, measure='Jiang', combine='max')
        oxidativeStress = mgoSim(c(term1), c('GO:1903206', 'GO:1902175'), semData=hsGO, measure='Jiang', combine='max')
        cellCycle = mgoSim(c(term1), c('GO:0045930', 'GO:0045787', 'GO:0033262', 'GO:1902749', 'GO:1902806', 'GO:0090266', 'GO:1901990', 'GO:1901992', 'GO:1902751', 'GO:0051726', 'GO:0000083'), semData=hsGO, measure='Jiang', combine='max')
        apoptosis = mgoSim(c(term1), c('GO:0008637', 'GO:1902165', 'GO:0043280', 'GO:1904035', 'GO:2001236', 'GO:1900119', 'GO:2001268', 'GO:0001844', 'GO:1902230', 'GO:2001239', 'GO:0043066', 'GO:1902175', 'GO:0006309', 'GO:1902172', 'GO:2001243', 'GO:1902235', 'GO:1902237', 'GO:0097194', 'GO:0043653', 'GO:1902043', 'GO:1901029', 'GO:0036462', 'GO:1902041', 'GO:0008635', 'GO:1902110', 'GO:2001235', 'GO:0043154', 'GO:2000353', 'GO:0043065', 'GO:1902255', 'GO:2001238', 'GO:0030263', 'GO:0006915', 'GO:1900118', 'GO:2001269', 'GO:1901028'), semData=hsGO, measure='Jiang', combine='max')
        return(c(pyroptosis, mitoptosis, necrosis, autophagy, cellDeath, oxidativeStress, cellCycle, apoptosis))
    }
    """)
    try:
        ss1 = semSimR(term1)
        return dict(zip(['pyroptosis', 'mitoptosis', 'necrosis', 'autophagy', 'cellDeath', 'oxidativeStress', 'cellCycle', 'apoptosis'], ss1))
    except:
        return 'NA'

# pull out GO terms
prerank_GOterms = []
for sampleID in miR497data:
    if os.path.exists('../3_FunctionalEnrichmentAnalyses/results/prerank/' + sampleID):
        resultsReport = pd.read_csv('../3_FunctionalEnrichmentAnalyses/results/prerank/' + sampleID + '/' + sampleID + '_prerankOutput_FDRcorrected_passed005.csv', index_col = 'Term')
        prerank_GOterms += [i.split('(')[1].rstrip(')') for i in resultsReport.index]

prerank_GOterms1 = list(set(prerank_GOterms))

#####################
## import GO terms ##
#####################
gene2GOfile = '../3_FunctionalEnrichmentAnalyses/gene2go.gz'
gene2go = {}
with gzip.open(gene2GOfile, 'rt') as inFile:
    #gene2go = inFile.read()
    inFile.readline()
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split('\t')
        if splitUp[0] == '9606':
            if not splitUp[2] in gene2go:
                gene2go[splitUp[2]] = []
            gene2go[splitUp[2]].append(splitUp[1])

# Run semantic similarity
if os.path.isfile('../3_FunctionalEnrichmentAnalyses/results/cellCycleCellDeath/semanticSimilarity_data_prerankFinal_systematicallyChosen.pkl'):
    simSem_Final = pd.read_pickle('../3_FunctionalEnrichmentAnalyses/results/cellCycleCellDeath/semanticSimilarity_data_prerankFinal_systematicallyChosen.pkl')
    goData_Final = pd.DataFrame.to_dict(simSem_Final, 'index')
    print('goTerms semantic similarity already calculated! :) Importing pickle...')
else:
    print('Finding GO term semantic similarity in network...')
    allGO = []
    for go in prerank_GOterms1:
        if not go in allGO:
            allGO.append(go)
    goSemSim_Final = {}
    for go1 in allGO:
        # Check to see if in gene2go
        if go1 in gene2go:
            print('Working on term ' + go1)
            goSemSim_Final[go1] = semSim_Mechanisms(go1)
    simSem_Final = pd.DataFrame.from_dict(goSemSim_Final, orient = 'index')
    simSem_Final.to_pickle('../3_FunctionalEnrichmentAnalyses/results/cellCycleCellDeath/semanticSimilarity_data_prerankFinal_systematicallyChosen.pkl')
    goData_Final = copy.copy(goSemSim_Final) # would not deepcopy

# upload prerank data
prerank_GOterms = {}
for sampleID in miR497data:
    if os.path.exists('../3_FunctionalEnrichmentAnalyses/results/prerank/' + sampleID):
        allData = pd.read_csv('../3_FunctionalEnrichmentAnalyses/results/prerank/' + sampleID + '/' + sampleID + '_prerankOutput_FDRcorrected_passed005.csv', index_col = 'Term').dropna()
        prerank_GOterms[sampleID] = list(set(['GO:' + i.split(':')[1].strip(')') for i in list(allData.index)]))

mechanisms = ['apoptosis', 'autophagy', 'mitoptosis', 'necrosis', 'oxidativeStress', 'pyroptosis', 'cellDeath', 'cellCycle']
## Build dataframe structure for results
prerank_index = []
for drugTerm in prerank_GOterms.keys():
    for direction in ['up', 'down']:
        prerank_index.append(drugTerm + ' ' + direction)

prerank_finalResults = pd.DataFrame(index = prerank_index, columns = mechanisms)

# Add a direction
directed_prerank_GOterms = {}
for sampleID in miR497data:
    if os.path.exists('../3_FunctionalEnrichmentAnalyses/results/prerank/' + sampleID):
        allData = pd.read_csv('../3_FunctionalEnrichmentAnalyses/results/prerank/' + sampleID + '/' + sampleID + '_prerankOutput_FDRcorrected_passed005.csv', index_col = 'Term').dropna()
        directed_prerank_GOterms[sampleID + ' up'] = ['GO:' + i.split(':')[1].strip(')') for i in list(allData.index) if allData.loc[i,'NES'] >= 0]
        directed_prerank_GOterms[sampleID + ' down'] = ['GO:' + i.split(':')[1].strip(')') for i in list(allData.index) if allData.loc[i,'NES'] < 0]

# Fill out results
for details in directed_prerank_GOterms:
    print(details)
    for GOterm_in_drug in directed_prerank_GOterms[details]:
        if GOterm_in_drug in goData_Final.keys():
            for semsimResult in goData_Final[GOterm_in_drug]:
                if str(prerank_finalResults.loc[details, semsimResult]) == str(float('nan')):
                    prerank_finalResults.loc[details, semsimResult] = goData_Final[GOterm_in_drug][semsimResult]
                elif prerank_finalResults.loc[details, semsimResult] <= goData_Final[GOterm_in_drug][semsimResult]:
                    prerank_finalResults.loc[details, semsimResult] = goData_Final[GOterm_in_drug][semsimResult]

prerank_finalResults.to_csv('../3_FunctionalEnrichmentAnalyses/results/prerank/prerank_finalResults_systematicallyChosen.csv', sep = ',')

###############################################
## Generate heatmaps from sematic similarity ##
###############################################
# structure the data
heatmapStructure_up = prerank_finalResults.loc[[i for i in prerank_finalResults.index if 'up' in i], [i for i in prerank_finalResults.columns if i != 'cellCycle']]
heatmapStructure_down = prerank_finalResults.loc[[i for i in prerank_finalResults.index if 'down' in i], 'cellCycle']


# plot the data - M12T
pp2 = PdfPages('../3_FunctionalEnrichmentAnalyses/results/prerank/cellDeath.pdf')
plt.title('Methods of cellular death: up regulated processes')
sns.heatmap(heatmapStructure_up.astype(float), square = True, cmap = 'gray_r', xticklabels = True, yticklabels = True, vmin = 0, vmax = 1)
plt.tight_layout()
pp2.savefig()
pp2.close()
plt.clf()
plt.close()

pp3 = PdfPages('../3_FunctionalEnrichmentAnalyses/results/prerank/cellCycle.pdf')
plt.title('Cell cycle: down regulated processes')
sns.heatmap(heatmapStructure_down.to_frame().astype(float), square = True, cmap = 'gray_r', xticklabels = True, yticklabels = True, vmin = 0, vmax = 1)
plt.tight_layout()
pp3.savefig()
pp3.close()
plt.clf()
plt.close()









# fin
