#read in data file
import pandas as pd
import sys

df = pd.read_csv(sys.argv[1])

#select columns by name
columns_to_keep = ['Chromosome', 'Feature',	'Start',	'End',	'Start_UTRref',	'End_UTRref',	'reference_id',	'ref_gene_id', 'Strand']
filtered_df = df[columns_to_keep]

# save to a new file
filtered_df.to_csv('SRR11463564_official.csv', index=False)