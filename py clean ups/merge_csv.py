import pandas as pd
import sys

# Load the first CSV file into a DataFrame
df1 = pd.read_csv(sys.argv[1])

# Load the second CSV file into a DataFrame
df2 = pd.read_csv(sys.argv[2])

# Merge the two DataFrames on the 'transcript_id' column
merged_df = pd.merge(df1, df2, on='transcript_id2')

# Save the merged DataFrame to a new CSV file
merged_df.to_csv('merged_file.csv', index=False)