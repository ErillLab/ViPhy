# This is a sample Python script.

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
import os.path as path


def nucleotideToAminoacid(DNAseq):
    """
    Translate a nucleotide sequence into six amino acid sequences.
    DNAseq: nucleotide sequence, should be string
    Returns a single amino acid sequence that come from combining the six amino acid sequences, should be string
    """

    protein_seq1 = ""
    protein_seq2 = ""
    protein_seq3 = ""

    # First direction
    for i in range(0, len(DNAseq), 3):
        codon1 = DNAseq[i:i + 3]
        codon2 = DNAseq[i + 1:i + 4]
        codon3 = DNAseq[i + 2:i + 5]
        if (len(codon1)==3):
            codon1 = Seq(codon1)
            protein_seq1 += codon1.translate('1')
        if (len(codon2)==3):
            codon2 = Seq(codon2)
            protein_seq2 += codon2.translate('1')
        if (len(codon3)==3):
            codon3 = Seq(codon3)
            protein_seq3 += codon3.translate('1')

    protein_seq = protein_seq1 + protein_seq2 + protein_seq3

    # Second direction
    DNAseqR = reverse_complement(DNAseq) # Return the reverse complement sequence
    DNAseqR = str(DNAseqR)

    protein_seq1 = ""
    protein_seq2 = ""
    protein_seq3 = ""

    for i in range(0, len(DNAseqR), 3):
        codon1 = DNAseqR[i:i + 3]
        codon2 = DNAseqR[i + 1:i + 4]
        codon3 = DNAseqR[i + 2:i + 5]
        if (len(codon1)==3):
            codon1 = Seq(codon1)
            protein_seq1 += codon1.translate('1')
        if (len(codon2)==3):
            codon2 = Seq(codon2)
            protein_seq2 += codon2.translate('1')
        if (len(codon3)==3):
            codon3 = Seq(codon3)
            protein_seq3 += codon3.translate('1')

    protein_seq += protein_seq1 + protein_seq2 + protein_seq3
    protein_seq = str(protein_seq)

    return protein_seq


def fastaToString(DNAfile):
    """
    Takes a FASTA file with a nucleotide sequence and transforms it to a string.
    DNAfile: file containing nucleotide sequence, should be FASTA format
    Returns the nucleotide sequence obtained from the fasta file and the first line of the fasta file.
    """
    file = open(DNAfile)
    fastaContent = file.read()
    fastaList = fastaContent.split("\n")
    fastaName = fastaList[0]
    fastaList = fastaList[1:]

    DNAseq = ""
    for line in fastaList:
        DNAseq += line
    file.close()
    return DNAseq, fastaName

def exportFasta(proteinSeq, fastaName):
    """
    Saves the protein sequence into a fasta file
    proteinSeq: combination of the six possible protein sequences. Obtained from "nucleotideToAminoacid" function
    fastaName: first line from the original fasta file
    """
    fileContent = ""
    file = open("proteinSeq.fasta", "w")
    for i in proteinSeq:
        fileContent += i
        #if (len(fileContent) % 60) == 0:
           # fileContent += '\n'
    fileContent = fastaName + '\n' + fileContent
    file.write(fileContent)
    file.close()

def compareFiles(inputFile, outputFile):
    """
        Compare the fasta file obtained after executing the code against the fasta file containing the expected results
        inputFile: resulting fasta file
        fastaName: fasta file containing the results that you hope to achieve
        """
    file1 = open(inputFile).readlines()
    file1_line = ""

    for lines in file1:
        lines = lines.rstrip('\n')
        file1_line += lines

    file2 = open(outputFile).readlines()
    file2_line = ""

    for lines in file2:
        lines = lines.rstrip('\n')
        file2_line += lines

    if file1_line == file2_line:
        print("Match. The result is correct")
    else:
        if len(file1_line) != len(file2_line):
            print("No Match. The length of both files is not the same")
        else:
            count = 0
            for a, b in zip(file1_line, file2_line):
                if a != b:
                    count += 1
            print("Not Match. Number of differences:", count)


# Main function
if __name__ == '__main__':

    fasta_file = 'sequence.fasta'
    testFile = 'expectedResult.fasta'

    DNAseq, fastaName = fastaToString(fasta_file)
    proteinSeq = nucleotideToAminoacid(DNAseq)

    exportFasta(proteinSeq, fastaName) # Fasta file containing the resulting protein sequence

    if path.exists(testFile):
        compareFiles('proteinSeq.fasta', testFile)






