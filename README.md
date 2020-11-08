# ViPhy
Virus phylogeny system.


## The script
The script ``main.py`` transcribe and translate a nucleotide sequence(DNA) into a sequence of proteins.

The steps in the script are the following:

1. Opens and reads ``settings.txt``
2. Reads a sequence(DNA) from a fasta or a genbank file
3. Tranforms the sequences into a string
4.1 If the input file contained a nucleotide sequence, then we start a translation process.(steps 5 to 9)
4.2 If the input file contained a protein sequence, then we save it into a new fasta file
5. Generates the three possible frames for each sequence (+1, +2, +3)
6. Reverse the nucleotide sequence (string) and generates the last three possible frames for each sequence (-1, -2, -3)
7. Swaps the DNA sequences for protein sequences
8. Join the six frames in the same sequence
9. Stores the combined protein sequence in a new fasta file


## Viphy-env.yml

Conda environment created using Phyton38 that contains a specific collection of conda packages installed to ensure the smooth running of the program.

To activate the `viphy-env.yml` environment:

	``$ conda activate viphy-env``

To deactivate the conda environment:

	``$ source deactivate``


## Settings.json

This file contains important information for the correct operation of the code. Please, don't change this file path!

Structure of the JSON file:

``
{
"input_folder" : "genome_data",
"genome_accessions": [["GQ919031.1", "JX182370.1"], ["NC_015464"], ["NC_042011.1"]],
"output_folder" : "results",
"working_folder" : "tmp",
"analysis_type": "nucleotide",
}
``


### Features: 

- "input_folder": Folder where you can find the fasta or genbank files that the program will read
- "genome_accessions": List of lists that contains the identifier of the files you want to download from Genbank database. 
- "working_folder": Path for the directory in which all the other files are located.
- "output_folder": Folder in which the outgoing files will be located.
- "analysis_type": Can be "nucleotide" or "protein"


## Other files

- ``sequence.fasta``: Fasta file containing a nucleotide or amino acid sequence, obtained from Genbank database. 

- ``sequence.gb``: Genbank file where you can obtained a nucleotido or amino acid sequence. 

- ``expectedResult.fasta``: Contains a expected protein sequence. It can be used for testing the translation process. 

- ``proteinSeq.fasta``: Inside `Results` folder. Contains the sequence of amino acids obtained after running this program.


## Folders

- ``Results``: Folder to stores the resulting sequences after the translation ends. 


## Dependencies

- Python 3.8.3
- Biopython 1.78
