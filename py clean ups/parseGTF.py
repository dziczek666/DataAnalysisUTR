import pyranges as pr
import sys 

gtf = pr.read_gtf(sys.argv[1], as_df=True)

print(gtf)

gtf.to_csv(sys.argv[2], index=False)

#read in data file
import pandas as pd
df = pd.read_csv('filtered.csv')

#select the desired columns by name
columns_to_keep = ['Start', 'End', 'transcript_id']
filtered_df = df[columns_to_keep]

# save the selected columns to a new file
filtered_df.to_csv('filteredUTR_SRR11463564.csv', index=False)
