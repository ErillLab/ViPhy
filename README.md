# ViPhy
Virus phylogeny system.

## The script
The script ``main.py`` translate a nucleotide sequence(DNA) into a sequence of amino acids.

The steps in the script are the following:

1. Reads the nucleotide sequence(DNA) from the fasta file
2. Tranforms the sequences into a string
3. Generates the three possible frames for each sequence (+1, +2, +3)
4. Reverse the nucleotide sequence (string) and generates the last three possible frames for each sequence (-1, -2, -3)
5. Swaps the DNA sequences for protein sequences
6. Join the six frames in the same sequence
7. Stores the combined protein sequence in a new fasta file

## Other files

``Sequence.fasta``: Contains a DNA sequence example obtained from Genbank database. The sequences should be stored in this file format.
``proteinSeq.fasta``: Contains the sequence of amino acids obtained after running this program.

## Dependencies

- Python 3.8.3
- Biopython 1.78

