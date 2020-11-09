# ViPhy
Virus phylogeny system.


## The script
The script ``main.py`` transcribe and translate a nucleotide sequence(DNA) into a sequence of proteins.

The steps in the script are the following:

1. Opens and reads ``settings.json``

2. If the analysis_type indicated in the Json file is 'nucleotide':

2.1. Traverses directories and subdirectories until there are no more files to read.

2.1.1. When it finds a fasta file, it checks if the file contains a nucleotide sequence. If not, then returns to the step 2.1.. 

2.1.2. When it finds a genbank file, it extracts the nucleotide sequence from the file. 

2.1.3. Transcribes and translates the nucleotide sequence into a protein sequence.

2.1.4. Stores the combined protein sequence in a new fasta file

3. If the analysis_type indicated in the Json file is 'protein':

3.1. Traverses directories and subdirectories until there are no more files to read.

3.1.1. When it finds a fasta file, it checks if the file contains a nucleotide sequence. If so, stores the sequence into a new fasta file. 

3.1.2. When it finds a genbank file, it extracts the protein sequence from the file. 

3.1.2.1. Stores the combined protein sequence in a new fasta file


## Viphy-env.yml

Conda environment created using Phyton38 that contains a specific collection of conda packages installed to ensure the smooth running of the program.

To activate the `viphy-env.yml` environment:

	$ conda activate viphy-env

To deactivate the conda environment:

	$ source deactivate


### Dependencies

- Python 3.8.3
- Biopython 1.78


## Settings.json

This file contains important information for the correct operation of the code. Please, don't change this file path!

Structure of the JSON file:


{

	"input_folder" : "genome_data",

	"genome_accessions": [["GQ919031.1", "JX182370.1"], ["NC_015464"], ["NC_042011.1"]],

	"output_folder" : "results",

	"working_folder" : "tmp",

	"analysis_type": "nucleotide",
}



### Details of the JSON file: 

- "input_folder": Folder where you can find the fasta or genbank files that the program will read
- "genome_accessions": List of lists that contains the identifier of the files you want to download from Genbank database. 
- "working_folder": Path for the directory in which all the other files are located.
- "output_folder": Folder in which the outgoing files will be located.
- "analysis_type": Can be "nucleotide" or "protein"


## Other files

- Fasta files: Inside `Inputs` folder. Documents that contain a nucleotide or amino acid sequence that the program will read. If the content is a nucleotide sequence, it will be translates into proteins.

- Genbank files: Inside `Inputs` folder. Documents where you can obtained a nucleotide or amino acid sequence, together with other relevant information. 

- ``expectedResult.fasta``: Contains a expected protein sequence. It can be used for testing the translation process. 

- ``proteinSeq.fasta``: Inside `Results` folder. Contains the sequence of amino acids obtained after running this program.


## Folders

- ``Inputs``: Folder to stores the set of files that the program will read 

- ``Results``: Folder to stores the resulting sequences after the translation ends. 


