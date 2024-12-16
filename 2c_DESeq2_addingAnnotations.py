import pandas as pd
import scipy.stats

# DEG files
deg_files = ['results_ei_vs_ci_12H.csv', 'results_ei_vs_ci_72H.csv', 'results_ep_vs_cp_12H.csv', 'results_ep_vs_cp_72H.csv']

# Read in DEG files
deg_pds = {}
for deg_file in deg_files:
    name = deg_file.strip('results_').strip('.csv')
    deg_pds[name] = pd.read_csv('../2_DESeq2/results/'+deg_file, index_col = 0, header = 0)

# DisGenet Genes
disgenet = list(pd.read_csv('../2_DESeq2/DisGeNET_databases/DisGeNET_genes.csv', index_col = 0, header = 0)['Genes'])

# Read in DEG files
surv_pds = {}
for surv_file in ['univariate','multivariate']:
    surv_pds[surv_file] = pd.read_excel('../4_SurvivalAnalyses/PMID_35354049/mRNA/mRNA_'+surv_file+'_CPH_Z_PMID3534049.xlsx', index_col = 0, header = 0)

# Add columns
for name1 in deg_pds:
    # Add multivariate
    deg_pds[name1]['CPH uni. Z'] = 'NA'
    deg_pds[name1]['CPH uni. p-value'] = 'NA'
    deg_pds[name1]['CPH multi. Z'] = 'NA'
    deg_pds[name1]['CPH multi. p-value'] = 'NA'
    deg_pds[name1]['DisGeNET'] = 'NA'

    for ens1 in deg_pds[name1].index:
        symbol1 = deg_pds[name1].loc[ens1,'symbols']
        # Univariate
        if symbol1 in surv_pds['univariate'].index:
            deg_pds[name1].loc[ens1,'CPH uni. Z'] = surv_pds['univariate'].loc[symbol1, 'MESO']
            if surv_pds['univariate'].loc[symbol1, 'MESO'] > 0:
                deg_pds[name1].loc[ens1,'CPH uni. p-value'] = scipy.stats.norm.sf(surv_pds['univariate'].loc[symbol1, 'MESO'])
            else:
                deg_pds[name1].loc[ens1,'CPH uni. p-value'] = scipy.stats.norm.cdf(surv_pds['univariate'].loc[symbol1, 'MESO'])

        # Multivariate
        if symbol1 in surv_pds['multivariate'].index:
            deg_pds[name1].loc[ens1,'CPH multi. Z'] = surv_pds['multivariate'].loc[symbol1, 'MESO']
            if surv_pds['multivariate'].loc[symbol1, 'MESO'] > 0:
                deg_pds[name1].loc[ens1,'CPH multi. p-value'] = scipy.stats.norm.sf(surv_pds['multivariate'].loc[symbol1, 'MESO'])
            else:
                deg_pds[name1].loc[ens1,'CPH multi. p-value'] = scipy.stats.norm.cdf(surv_pds['multivariate'].loc[symbol1, 'MESO'])

        # DisGeNET
        if symbol1 in disgenet:
            deg_pds[name1].loc[ens1,'DisGeNET'] = 'Yes'

# Write out as a multitabbed excel file
with pd.ExcelWriter('../6_Supplementary/results/differentially_expressed_genes_plus_6_5_2024.xlsx', engine='openpyxl') as writer:
    for name1 in deg_pds:
        deg_pds[name1].to_excel(writer, sheet_name=name1)

