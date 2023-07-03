import os
import pandas as pd

# specify the dir
directory = './'

newColNames = ['Chromosome', 'ExonStart', 'UTREnd', 'Strand', 'Covrage', 'RefID', 'RefGeneID', 'RefUTRStart', 'RefUTREnd', 'RefUTRLength', 'UTRLength', 'UTRDiff', 'UTRType']

for file in os.listdir(directory):
    if file.endswith('.csv'):
        # read the CSV file into a DataFrame
        df = pd.read_csv(os.path.join(directory, file))
        # rename the columns
        df.columns = newColNames
        # change UTRDiff value to absolute numbers
        df['UTRDiff'] = df['UTRDiff'].abs()
        # sort from the bigger UTRDiff
        df = df.sort_values('UTRDiff', ascending=False)
        # drop rows if coverage < 1
        df = df[df['Covrage'] > 1]
        # Save first 1000 rows
        df_head_1000 = df.head(1000)
        newNames = f'{file}_{1000}.csv'
        # save to file    
        df.to_csv(file, index=False)
        df_head_1000.to_csv(newNames, index=False)

