import os
import re

# specify the directory containing the CSV files
directory = './'

# specify the new name for the files
new_name = 'filtered'

# loop through all files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.csv'):
        # extract the ID from the filename using a regular expression
        match = re.search(r'SRR\d+', filename)
        if match:
            id = match.group()
            
            # construct the new file name with the ID and the new name
            new_filename = f'{id}_{new_name}.csv'
            
            # rename the file
            os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))