import os
import pandas as pd

# Define the path
directory = './'

# Create an empty list to hold the df
dfs = []

# Loop through the files
for filename in os.listdir(directory):
    if filename.endswith('.csv'):
        df = pd.read_csv(os.path.join(directory, filename))
        # Add a column
        df['filename'] = filename
        # Append the data frame to the list
        dfs.append(df)

# Merge the data frames into one
merged_df = pd.concat(dfs, ignore_index=True)

# Write to a new file
merged_df.to_csv("merged.csv", index=False)