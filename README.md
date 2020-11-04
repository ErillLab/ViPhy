# ViPhy
Virus phylogeny system.


## The script
The script ``main.py`` translate a nucleotide sequence(DNA) into a sequence of amino acids.

The steps in the script are the following:

1. Opens and reads ``settings.txt``
2. Reads a sequence(DNA) from a fasta or a genbank file
3. Tranforms the sequences into a string
4.1 If the input file contained a nucleotide sequence, then we start a translation process.(steps 5 to 9)
4.2 If the input file contained a protein sequence, then we save it into a new fasta file
5. Generates the three possible frames for each sequence (+1, +2, +3)
<<<<<<< HEAD
6. Reverse the nucleotide sequence (string) and generates the last three possible frames for each sequence (-1, -2, -3)
7. Swaps the DNA sequences for protein sequences
8. Join the six frames in the same sequence
=======
6. Reverses the nucleotide sequence (string) and generates the last three possible frames for each sequence (-1, -2, -3)
7. Swaps the DNA sequences for protein sequences
8. Joins the six frames in the same sequence
>>>>>>> origin/main
9. Stores the combined protein sequence in a new fasta file


## Settings.txt
<<<<<<< HEAD

This file contains important information for the correct operation of the code. Please, don't change this file path!

=======
>>>>>>> origin/main
Structure of the JSON file:

``
{
"input_file": "fasta",
"working_folder": "",
"output_folder": "Results/",
"analysis_type": "nucleotide"
}
``
### Features: 

- "input_file": Can be "fasta" (fasta file) or "gb" (genbank file).
- "working_folder":  Path for the directory in which the inbound files are located.
- "output_folder": Path for the directory in which the outgoing files are located.
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
<<<<<<< HEAD
=======


>>>>>>> origin/main
