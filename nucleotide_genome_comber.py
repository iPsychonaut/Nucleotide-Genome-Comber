# -*- coding: utf-8 -*-
"""
Created on Tue May 16 22:48:58 2023

@author: ian.michael.bollinger@gmail.com
"""
# Import required libraries
import random
import pandas as pd
import time
import os
import sys
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
import glob
import gzip
import subprocess
import tempfile
from tqdm import tqdm

def create_BLAST_db(database_dir):
    # Find all .gbff.gz and .cds files in the directory
    gbff_files = glob.glob(database_dir + "/*.gbff.gz")
    cds_files = glob.glob(database_dir + "/*.cds")

    if not gbff_files and not cds_files:
        print("\nNo .gbff.gz or .cds files found in the directory")

        # Check if there are already BLAST databases in the directory
        files_in_directory = os.listdir()
        filtered_files = [file for file in files_in_directory if "blast_db" in file]
        if filtered_files:
            print("\nBLAST databases already exist in the directory:")
            for file in filtered_files:
                print(file)
            return

    fasta_files = []
    # Convert GBFF to FASTA
    for input_filename in tqdm(gbff_files, desc="Converting GBFF files to FASTA"):
        output_filename = input_filename.replace('.gbff.gz','.fasta')
        with gzip.open(input_filename, "rt") as handle:
            count = SeqIO.write(SeqIO.parse(handle, "genbank"), output_filename, "fasta")
        print(f"\nConverted {count} records from {input_filename}")
        fasta_files.append(output_filename)

    # Add .cds files to the list of FASTA files
    fasta_files.extend(cds_files)

    # Create a BLAST database from each FASTA file
    for fasta_file in tqdm(fasta_files, desc="Creating BLAST databases"):
        db_name = 'blast_db'
        cmd = f"makeblastdb -in {fasta_file} -dbtype nucl -out {db_name}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        print(result.stdout)
        print(result.stderr)


def find_fasta_file(directory):
    """Function to find the .fasta file in the given directory."""
    # Initialize an empty list to hold the names of .fasta files
    fasta_files = []

    # Iterate over all files in the directory
    for file in os.listdir(directory):
        # If the file is a .fasta file, add its name to the list
        if file.endswith(".fasta"):
            fasta_files.append(file)
    
    # If there is exactly one .fasta file, return its name
    if len(fasta_files) == 1:
        return fasta_files[0]
    # If there are no .fasta files, return None
    elif len(fasta_files) == 0:
        return None
    else:
        # If there is more than one .fasta file, raise an error
        raise ValueError("Multiple .fasta files found in the directory.")


def chunk_seq(record, name):
    """Function to divide a sequence into chunks."""
    # Calculate the length of the sequence and the size of each chunk
    length = len(record)

    # If the sequence is shorter than 1000 bases, use the whole sequence for each chunk
    if length < 1000:
        return record.seq, record.seq, record.seq

    # Divide the sequence into three parts
    one_end = int(length/3)
    start = record.seq[0:one_end]
    middle = record.seq[one_end:2*one_end]
    end = record.seq[2*one_end:]

    def handle_chunk(chunk):
        """Function to handle a chunk of a sequence."""
        if len(chunk) > 2500:
            # If the chunk is longer than 2500 bases, randomly select a 2500 base segment
            start = random.randint(0, len(chunk) - 2500)  # ensure there's at least 2500 bases left after start
            return chunk[start:start+2500]
        else:
            # If the chunk is shorter than 1000 bases, return the whole chunk
            return chunk

    # Handle each chunk
    start = handle_chunk(start)
    middle = handle_chunk(middle)
    end = handle_chunk(end)
    
    return start, middle, end


def run_blastn(sequence, db_name):
    # Create a temporary file and write the sequence to it
    with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_file:
        temp_file.write(f">query\n{sequence}")
        temp_file_name = temp_file.name

    # Define the blastn command
    blastn_cline = NcbiblastnCommandline(query=temp_file_name, db=db_name, evalue=0.001, outfmt="6 qseqid sseqid pident qlen qstart qend length gaps", out="blastn_out.txt")
    
    # Run the command
    stdout, stderr = blastn_cline()
    
    # Read the output file into a pandas DataFrame
    blast_results = pd.read_csv("blastn_out.txt", sep="\t", names=["qseqid", "sseqid", "pident", "qlen", "qstart", "qend", "length", "gaps"])
    
    # Clean up the temporary file
    os.unlink(temp_file_name)
    
    # Check if the DataFrame is empty
    if blast_results.empty:
        return None, None, None
        
    # Pull out data based on column names
    alignment_length = blast_results['length'].iloc[0]
    percent_id = blast_results['pident'].iloc[0]
    query_len = blast_results['qlen'].iloc[0]
    qstart = blast_results['qstart'].iloc[0]
    qend = blast_results['qend'].iloc[0]
    gaps = blast_results['gaps'].iloc[0]
    
    # Calculate query coverage
    percent_qmatch = (alignment_length / query_len) * 100

    return alignment_length, percent_id, percent_qmatch
    

if __name__ == '__main__':
    # Set the default values
    DATABASE_DIRECTORY = 'C:/Users/theda/OneDrive/Documents/Python/Entheome/Nucleotide-Genome-Comber/Resources/Oryza_sativa/Nipponbare' # Replace with a folder containing a cds or gbff.gz file to use as search filter
    MIXED_GENOME_DIRECTORY = 'C:/Users/theda/OneDrive/Documents/Python/Entheome/Nucleotide-Genome-Comber/Resources/E1' # Replace with a folder containing only the fasta file for contigs that need filtering

    # Check command line arguments
    if len(sys.argv) == 3:
        # Get the command line arguments if provided
        DATABASE_DIRECTORY = sys.argv[1]
        MIXED_GENOME_DIRECTORY = sys.argv[2]
    else:
        print("No command line arguments provided, using default directories.")
        print(f"DATABASE_DIRECTORY: {DATABASE_DIRECTORY}")
        print(f"MIXED_GENOME_DIRECTORY: {MIXED_GENOME_DIRECTORY}")

    #
    # Step 1 Prepare Filter Database
    #    
    # Load Database(s) to use as filters
    os.chdir(DATABASE_DIRECTORY)
    create_BLAST_db(DATABASE_DIRECTORY)

    # NCBI_database_dir = 'C:/Users/theda/OneDrive/Documents/Python/Entheome/Resources/Oryza_sativa/NCBI'
    # os.chdir(NCBI_database_dir)
    # create_BLAST_db(NCBI_database_dir)


    # 
    # Step 2 Load Mixed Genome Fasta
    #   
    # Define directory where the mixed genome is located
      # Change the current working directory to the mixed genome directory
    os.chdir(MIXED_GENOME_DIRECTORY)

    # Define the columns of the DataFrame
    columns = ['contig_name', 'start_chunk', 'mid_chunk', 'end_chunk', 'start_percent_id', 'mid_percent_id', 'end_percent_id', 'start_percent_qmatch', 'mid_percent_qmatch', 'end_percent_qmatch', 'start_alignment_length', 'mid_alignment_length', 'end_alignment_length']

    # Initialize an empty DataFrame with the specified columns
    df = pd.DataFrame(columns=columns)

    # Use the function to find the .fasta file in the mixed genome directory
    INPUT_FASTA = find_fasta_file(MIXED_GENOME_DIRECTORY)


    # 
    # Step 3 Prep the Fasta file into 1000bp contigs 
    # 
    # Check if a .fasta file was found
    if INPUT_FASTA is not None:
        print("\nFound .fasta file:", INPUT_FASTA)
    else:
        print("\nNo .fasta file found in the directory.")

    # Parse the .fasta file into sequences
    fasta_sequences = SeqIO.parse(open(INPUT_FASTA),'fasta')

    # Initialize a counter for the contig number
    contig_number = 0

    # For each sequence in the .fasta file
    for fasta in tqdm(fasta_sequences, desc="Preparing contigs"):
        # Get the name and sequence
        name, sequence = fasta.id, str(fasta.seq)
       
        # Try to divide the sequence into chunks
        try:
            start_chunk, mid_chunk, end_chunk = chunk_seq(fasta, name)
            # If successful, add the chunks to the DataFrame
            df.loc[contig_number] = [name, start_chunk, mid_chunk, end_chunk, None, None, None, None, None, None, None, None, None]
            # Increment the contig number
            contig_number += 1

        # If an error occurs, print a message and skip to the next sequence
        except ValueError as e:
            print("\nSkipping contig {} due to error: {}".format(name, e))

    # Save the DataFrame as a CSV file
    df.to_csv(INPUT_FASTA.replace('.fasta','.csv'), index=False)


    # Step 4 Run Local BLAST for each chunk
    # 
    # Start the timer
    start_time = time.time()

    # Iterate over each row in the DataFrame and each database
    for database_dir, csv_suffix, database_cmd in [(DATABASE_DIRECTORY, '_blast_REMOVED.csv', 'blast_db')]:
        os.chdir(database_dir)
        for index, row in tqdm(df.iterrows(), total=df.shape[0], desc=f"Running local BLAST on {database_dir.split('/')[-1]}"):
            # For each chunk, run local BLAST and update the percent_id in the DataFrame
            for chunk in ['start_chunk', 'mid_chunk', 'end_chunk']:
                chunk_sequence = row[chunk]
                        
                alignment_length, percent_id, percent_qmatch = run_blastn(chunk_sequence, database_cmd)
                df.at[index, chunk.replace('chunk', 'percent_id')] = percent_id
                df.at[index, chunk.replace('chunk', 'percent_qmatch')] = percent_qmatch
                df.at[index, chunk.replace('chunk', 'alignment_length')] = alignment_length

        # Convert percent_id columns to numeric
        for col in ['start_percent_id', 'mid_percent_id', 'end_percent_id', 'start_percent_qmatch', 'mid_percent_qmatch', 'end_percent_qmatch', 'start_alignment_length', 'mid_alignment_length', 'end_alignment_length']:
            df[col] = pd.to_numeric(df[col], errors='coerce')

        # Filter the DataFrame to remove rows with a percent_id greater than 80 and percent_qmatch greater than 50
        filtered_df = df[~((df['start_percent_id'] > 80.000) & (df['start_percent_qmatch'] > 50.000) |
                            (df['mid_percent_id'] > 80.000) & (df['mid_percent_qmatch'] > 50.000) |
                            (df['end_percent_id'] > 80.000) & (df['end_percent_qmatch'] > 50.000) |
                            ((df['start_percent_id'] == 100.000) & (df['start_alignment_length'] > 100)) |
                            ((df['mid_percent_id'] == 100.000) & (df['mid_alignment_length'] > 100)) |
                            ((df['end_percent_id'] == 100.000) & (df['end_alignment_length'] > 100)))]
    
        removed_df = df[~df.index.isin(filtered_df.index)]
    
        # Save the removed rows as a separate CSV file
        os.chdir(MIXED_GENOME_DIRECTORY)
        removed_df.to_csv(INPUT_FASTA.replace('.fasta', csv_suffix), index=False)
    
        # Update the DataFrame to only include the filtered rows
        df = filtered_df

    # Calculate the percentage of contigs kept
    contig_list = df['contig_name'].tolist()
    print(f'\nKeeping {round((len(contig_list)/contig_number)*100,0)}% of Contigs')

    # Stop the timer and print the duration
    end_time = time.time()
    duration = end_time - start_time
    duration_minutes = round(duration/60, 0)
    print(f"\nStep 4 took {duration_minutes} minutes to complete")

    # Save the remaining contigs as a new .fasta file
    OUTPUT_FASTA = INPUT_FASTA.replace('.fasta', '_filtered.fasta')

    with open(INPUT_FASTA, 'r') as original, open(OUTPUT_FASTA, 'w') as filtered:
        fasta_sequences = SeqIO.parse(original, 'fasta')
        for fasta in fasta_sequences:
            if fasta.id in contig_list:
                SeqIO.write(fasta, filtered, 'fasta')

    print(f"\nFiltered fasta file saved as {OUTPUT_FASTA}")
