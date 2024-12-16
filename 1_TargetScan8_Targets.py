import pandas as pd

# Load data
df1 = pd.read_csv('../1_CollectTargetScan8Info/05222024/Predicted_Targets_Info.default_predictions.txt/Predicted_Targets_Info.default_predictions.txt', sep='\t')

# Filter for miR-497-5p
miR_sites = df1.loc[df1['miR Family'].str.contains('497-5p')]

# Dump out file
mir_sites_v2 = miR_sites[['Gene ID','Gene Symbol']].drop_duplicates()
mir_sites_v2['Gene_ID_noVer'] = [i.split('.')[0] for i in mir_sites_v2['Gene ID']]
mir_sites_v2.to_csv('../1_CollectTargetScan8Info/miR_497_5p_targetscan_V8_0.csv')
