# 2020

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
import os.path as path
import json


def nucleotideToAminoacid(DNAseq):
    '''
    Translates a nucleotide sequence into six amino acid sequences.
    :param DNAseq: nucleotide sequence, should be string
    :return: Returns a single amino acid sequence that come from combining the six possible amino acid sequences, should be string
    '''
    directions = 2
    final_protein_seq = ""

    # First direction: j = 0 / Second direction: j = 1
    for j in range(0, directions):
        protein_seq1 = ""
        protein_seq2 = ""
        protein_seq3 = ""
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
        final_protein_seq += protein_seq1 + protein_seq2 + protein_seq3

        if j == 0:
            DNAseq = reverse_complement(DNAseq) # Return the reverse complement sequence
            DNAseq = str(DNAseq)

    final_protein_seq = str(final_protein_seq)

    return final_protein_seq

def isNucleotide(nucleotide_line):
    '''
    Indicates if the input value is nucleotide sequence by returning a boolean
    :param nucleotide_line: line from a fasta file that contains a nucleotide or protein sequence
    :return: Returns a boolean, whose value is True when the input value is a nucleotide sequence
    '''
    is_nucleotide = True

    for nucleobase in nucleotide_line:
        if nucleobase != 'A' and nucleobase != 'C' and nucleobase != 'G' and nucleobase != 'T':
            if nucleobase == 'N':
                nucleobase.replace('N', '')
            else:
                is_nucleotide = False

    return is_nucleotide


def fastaToString(DNAfile):
    '''
    Takes a FASTA file with a nucleotide or protein sequence and transforms it to a string.
    :param DNAfile: file that contains a nucleotide or protein sequence, should be FASTA format
    :return: Returns the sequence obtained from the fasta file as a string, the first line of the fasta file and a boolean value
    that indicates if "DNAfile" contains a nucleotide sequence.
    '''
    DNAseq = ""
    is_nucleotide = True

    file = open(DNAfile)
    fasta_content = file.read()
    fasta_list = fasta_content.split("\n")
    fasta_name_line = fasta_list[0]
    fasta_list = fasta_list[1:]

    for line in fasta_list:
        DNAseq += line
        is_nucleotide = isNucleotide(line)

    file.close()
    return DNAseq, fasta_name_line, is_nucleotide


def gbToString(gb_file):
    '''
    Extracts a nucleotide sequence from a Genbank file and transforms it to a string.
    :param gb_file: file that can store several sequences and extra information.
    :return: Returns the sequence obtained from the genbank file as a string and the nucleotide identifier
    '''
    gb_name = ""
    DNAseq = ""

    try:
        for seq_record in SeqIO.parse(gb_file, "genbank"):
            gb_name = str(seq_record.id)
            DNAseq = str(seq_record.seq)
        return DNAseq, gb_name
    except:
        print("Genbank file is incomplete")


def exportFasta(protein_seq, fasta_name, output_folder):
    '''
    Saves the protein sequence into a fasta file
    :param protein_seq: combination of the six possible protein sequences. Obtained from "nucleotideToAminoacid" function
    :param fasta_name: first line from the original fasta file
    :param output_folder: directory where fasta file will be saved
    :return: Returns the directory where the fasta file is stored
    '''
    file_content = ""
    output_file = output_folder + 'proteinSeq.fasta'
    file = open(output_file, "w")
    for i in protein_seq:
        file_content += i
        #if (len(fileContent) % 60) == 0:
           # fileContent += '\n'
    file_content = fasta_name + '\n' + file_content
    file.write(file_content)
    file.close()
    return output_file


def compareFiles(input_file, output_file):
    '''
    Compares the resulting fasta file against the fasta file containing the expected results
    :param input_file: resulting fasta file
    :param output_file: fasta file containing the results that you hope to achieve
    '''
    # Input file
    file1 = open(input_file).readlines()
    file1_line = ""
    for lines in file1:
        lines = lines.rstrip('\n')
        file1_line += lines

    # Output file
    file2 = open(output_file).readlines()
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

    try:
        setting_file = open('settings.json', )
        json_file = json.load(setting_file)

        fasta_file = json_file["working_folder"] + 'sequence.fasta'
        gb_file = json_file["working_folder"] + 'sequence.gb'
        test_file = json_file["working_folder"] + 'expectedResult.fasta'

        if json_file["analysis_type"] == 'nucleotide':
            if json_file["input_file"] == 'fasta':
                if path.exists(fasta_file):
                    DNAseq, identifier_seq, nucleotide_type = fastaToString(fasta_file)
                    if nucleotide_type:
                        # translate nucleotide to protein
                        protein_seq = nucleotideToAminoacid(DNAseq)

                        # Create the fasta file that will contain resulting protein sequence
                        output_file = exportFasta(protein_seq, identifier_seq, json_file["output_folder"])

                        if path.exists(test_file):
                            compareFiles(output_file, test_file)
                        else:
                            print('Could not find expectedResult.fasta')
                    else:
                        print('Error: fasta file does not contain a nucleotide sequence')
                else:
                    print('Could not find fasta file')
            else:
                if json_file["input_file"] == 'gb':
                    # Check that genbank file do contain sequence
                    if path.exists(gb_file):
                        DNAseq, identifier_seq = gbToString(gb_file)
                        protein_seq = nucleotideToAminoacid(DNAseq)
                        output_file = exportFasta(protein_seq, identifier_seq, json_file["output_folder"])
                    else:
                        print('Could not find the genbank file')
                else:
                    print('Incorrect input_file value in json file')

        else:
            if json_file["analysis_type"] == 'protein':
                protein_list = []
                if path.exists(gb_file):
                    try:
                        for seq_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
                            for seq_feature in seq_record.features:
                                if seq_feature.type == "CDS":
                                    protein_list.append(seq_feature.qualifiers['translation'][0])  # Saves protein sequences
                    except:
                        print("Genbank file is incomplete")
                else:
                    print('Could not find the genbank file')
            else:
                print('Incorrect analysis_type value in json file')
        setting_file.close()

    except:
        print('Could not open settings.txt')





