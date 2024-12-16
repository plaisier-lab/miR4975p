import pandas as pd
import numpy as np

file_names = ['ci1', 'ci2', 'ci3', 'cp1', 'cp2', 'cp3', 'ei1', 'ei2', 'ei3', 'ep1', 'ep2', 'ep3']

## Verify that the indices are the same
# All lists are true, not the best way of doing this but kind of effective 12/14/2021

data = []
for file in file_names:
    for time in ['12_hours','72_hours']:
        column_name = file
        if time=='12_hours':
            column_name = column_name.capitalize()
        tmp = pd.read_csv('../2_DESeq2/hit-counts/'+time+'/'+file+'.counts.txt', sep = '\t', index_col = 'Geneid', usecols = ['Geneid', column_name])
        tmp.columns = [file+'_'+time]
        data.append(tmp)


## Combine the counts (12 files, tab separated)
# Reference file: Dropbox (ASU)/Mesothelioma/mesoCellLines/mRNA/gexp_counts

gexp_counts = pd.concat(data, axis = 1)


gexp_counts.to_csv('../2_DESeq2/gexp_counts.csv', sep = ',')
