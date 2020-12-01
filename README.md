# ViPhy
Virus phylogeny system.


## Preprocessing process
Translates a nucleotide sequences(DNA) into a sequence of proteins.

The steps in the script are the following:

1. Opens and reads ``settings.json``

2. Downloads the genbank genomes especified in the Json file from NCBI 

3. If the analysis_type indicated in the Json file is 'nucleotide':

3.1. Navigate into input files directory until there are no more files to read.

3.1.1. When it finds a fasta file, it checks if the file contains a nucleotide sequence. If not, then returns to the step 3.1.. 

3.1.2. When it finds a genbank file, it extracts the nucleotide sequence from the file. 

3.1.3. Translates the nucleotide sequence into a protein sequence.

3.1.4. Stores the combined protein sequence in a new fasta file


4. If the analysis_type indicated in the Json file is 'protein':

4.1. Navigates into the input files directory until there are no more files to read.

4.1.1. When it finds a fasta file, it checks if the file contains a nucleotide sequence. If so, stores the sequence into a new fasta file. 

4.1.2. When it finds a genbank file, it extracts the protein sequence from the file. 

4.1.2.1. Stores the combined protein sequence in a new fasta file


## Blast
Creates a Blast database with the fasta files stored in ``WorkingFolder``. Then compares each protein fasta file to each created database.


## Viphy_env.yml

Conda environment created using Phyton38 that contains a specific collection of conda packages installed to ensure the smooth running of the program.

To activate the `viphy_env.yml` environment:

	$ conda activate virtual-env

To deactivate the conda environment:

	$ source deactivate


### Dependencies

- Python 3.8.3
- Biopython 1.78


## Settings.json

This file contains important information for the correct operation of the code. Please, don't change this file path!

Structure of the JSON file:


{

	"user_email": "A.N.Other@example.com",

	"input_folder" : "genome_data",

	"genome_accessions": [["GQ919031.1", "JX182370.1"], ["NC_015464"], ["NC_042011.1"]],

	"output_folder" : "results",

	"working_folder" : "tmp",

	"analysis_type": "nucleotide",

	"e_value": 0.0001

}



### Details of the JSON file: 

- "user_email": Email that will be used to identify the user in NCBI
- "analysis_type": Can be "nucleotide" or "protein"
- "input_folder": Folder where you can find the fasta or genbank files that the program will read
- "working_folder": Path for the directory in which all the other files are located.
- "output_folder": Folder in which the outgoing files will be located.
- "genome_accessions": List of lists that contains the identifier of the files you want to download from Genbank database. 
- "e_value": Number of expected hits of similar score that could be found just by chance 



## Other files

- Fasta files: Inside `Inputs` folder. Documents that contain a nucleotide or amino acid sequence that the program will read. If the content is a nucleotide sequence, it will be translates into proteins.

- Genbank files: Inside `Inputs` folder. Documents where you can obtained a nucleotide or amino acid sequence, together with other relevant information. 



## Folders

- ``Inputs``: Folder to stores the set of files that the program will read 

- ``WorkingFolder``: Folder to stores the resulting sequences after the translation ends. 

- ``Outputs``: Folder that store the final results


