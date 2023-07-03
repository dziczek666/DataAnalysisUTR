import csv

# Input and output file paths
csv_file = 'SRR2148236_filtered.csv_1000.csv'
gtf_file = 'SRR2148236.gtf'

# Open CSV file for reading and GTF file for writing
with open(csv_file, 'r') as csv_in, open(gtf_file, 'w') as gtf_out:
    # Create a CSV reader object
    reader = csv.reader(csv_in)
    
    # Loop over rows in CSV file
    for row in reader:
        # Extract required information from CSV row
        chrom = row[0]
        source = 'MyPipeline'
        feature_type = 'exon'
        start = row[1]
        end = row[2]
        score = '.'
        strand = row[3]
        frame = '.'
        gene_id = row[6]
        transcript_id = row[5]
        
        # Write GTF row to output file
        gtf_out.write(f'{chrom}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tgene_id "{gene_id}"; transcript_id "{transcript_id}";\n')

