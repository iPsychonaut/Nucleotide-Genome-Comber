# Nucleotide Genome Comber

This repository contains a Python script `nucleotide_genome_comber.py` designed for processing genomes. The pipeline involves various steps including preparing a database for Basic Local Alignment Search Tool (BLAST), loading a mixed genome fasta, preparing the fasta file into chunks, and finally running a local BLAST for each chunk, trimming away any matches.

## Requirements

The script is written in Python and requires the following Python libraries:

- BioPython
- glob
- gzip
- os
- sys
- subprocess
- tqdm
- random
- pandas
- time

In addition, it requires a local installation of BLAST+. You can download and install BLAST+ from the [NCBI website](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

The example data uses the [Nipponbare Genome for Oryza sativa](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4958085/) as its filter database.

## Steps

Here's a brief overview of the main steps performed by the script:

1. **Prepare MM Database**: This step involves creating a BLAST database from the GenBank (.gbff.gz) and/or coding sequence (.cds) file found in a specified directory.

2. **Load Mixed Genome Fasta**: The script then loads a mixed genome fasta from the specified directory.

3. **Prepare the Fasta file into chunks**: The script divides the genome into 1000 base pair contigs. Each contig is further divided into three equal chunks - start, middle, and end. If a contig is smaller than 1000 base pairs, it is used as is.

4. **Run Local BLAST for each chunk and Trim Matches**: In this step, the script performs a local BLAST for each chunk against a specified database. It then records the percentage of identity for each chunk. Any contigs that have a match in the database are trimmed away.

Finally, the script filters the contigs based on the percent identity. Contigs with more than 80% identity are removed, and the remaining contigs are saved into a new fasta file.

## Usage
NOTE: The database can be built from a '.gbff.gz' file or a '.cds' file; both the database & mixed genome folders should only have one (1) file in it each

Run the script from the command line and provide the two required directories as arguments, like this:

python nucleotide_genome_comber.py C:/path/to/database/directory C:/path/to/mixed/genome/directory

OR

The main script `nucleotide_genome_comber.py` needs to have the following lines altered to suit your needs:

- Line 142: DATABASE_DIRECTORY = 'C:/Users/theda/OneDrive/Documents/Python/Entheome/Resources/Oryza_sativa/Nipponbare' # Replace with a folder containing a cds or gbff.gz file to use as search filter
- Line 143: MIXED_GENOME_DIRECTORY = 'C:/Users/theda/OneDrive/Documents/Python/Entheome/Resources/E1' # Replace with a folder containing only the fasta file for contigs that need filtering

## Contact

Please feel free to contact me at ian.michael.bollinger@gmail.com if you have any questions or issues.