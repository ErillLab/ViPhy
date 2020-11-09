# 2020

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
import os.path as path
import json
import os



def nucleotideToAminoacid(RNAseq):
    '''
    Translates a nucleotide sequence into six amino acid sequences.
    :param RNAseq: nucleotide sequence, should be string
    :return: Returns a single amino acid sequence that come from combining the six possible amino acid sequences, should be string
    '''
    directions = 2
    final_protein_seq = ""

    # First direction: j = 0 / Second direction: j = 1
    for j in range(0, directions):
        protein_seq1 = ""
        protein_seq2 = ""
        protein_seq3 = ""
        for i in range(0, len(RNAseq), 3):
            codon1 = RNAseq[i:i + 3]
            codon2 = RNAseq[i + 1:i + 4]
            codon3 = RNAseq[i + 2:i + 5]
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
            RNAseq = reverse_complement(RNAseq) # Return the reverse complement sequence
            RNAseq = str(RNAseq)

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
    Takes a FASTA file with a nucleotide or protein sequence, transcribe it into RNA and transforms it to a string.
    :param DNAfile: file that contains a nucleotide or protein sequence, should be FASTA format
    :return: Returns the transcribed RNA sequence obtained from the fasta file as a string, the identifier of the fasta
    file and a boolean value that indicates if "DNAfile" contains a nucleotide sequence.
    '''
    DNAseq = ""
    RNAseq = ""
    fasta_seq_name = ""

    fasta_sequences = SeqIO.parse(open(DNAfile), 'fasta')

    for fasta in fasta_sequences:
        fasta_seq_name, DNAseq = fasta.id, fasta.seq
    is_nucleotide = isNucleotide(str(DNAseq))
    if is_nucleotide:
        RNAseq = str(DNAseq.transcribe())
    else:
        is_nucleotide = False

    return RNAseq, fasta_seq_name, is_nucleotide


def gbToString(gb_file):
    '''
    Extracts a nucleotide sequence from a Genbank file, transcribe it into RNA and transforms it to a string.
    :param gb_file: file that can store several sequences and extra information.
    :return: Returns the transcribed sequence obtained from the genbank file as a string and its identifier
    '''
    gb_name = ""
    RNAseq = ""

    try:
        for seq_record in SeqIO.parse(gb_file, "genbank"):
            gb_name = str(seq_record.id)
            RNAseq = seq_record.seq.transcribe()
            RNAseq = str(RNAseq)

        return RNAseq, gb_name
    except:
        print("Genbank file is incomplete")


def exportFasta(protein_seq, fasta_name, output_folder, file_name):
    '''
    Saves the protein sequence into a fasta file
    :param protein_seq: combination of the six possible protein sequences. Obtained from "nucleotideToAminoacid" function
    :param fasta_name: first line from the original fasta file
    :param output_folder: directory where fasta file will be saved
    :param file_name: name of the input file
    '''
    file_content = ""
    output_file = output_folder + file_name

    file = open(output_file, "w")
    for i in protein_seq:
        file_content += i

    file_content = fasta_name + '\n' + file_content
    file.write(file_content)
    file.close()



# Main function
if __name__ == '__main__':

    # configuration file
    setting_file = ""
    try:
        setting_file = open('settings.json', )
    except:
        print('Could not open settings.txt')

    json_file = json.load(setting_file)  # Read json file

    test_file = 'expectedResult.fasta'

    # Read directories and subdirectories
    input_path = "."

    for dir_path, subdir_list, file_list in os.walk(input_path):
        for fname in file_list:
            full_path = os.path.join(dir_path, fname)

            with open(full_path, 'r') as reader:
                if dir_path == input_path + json_file["working_folder"] + "\\" + json_file["input_folder"]:
                    final_path = dir_path + '\\' + fname
                    final_output_path = input_path + json_file["working_folder"] + "\\" + json_file["output_folder"] + '\\'
                    extension = path.splitext(fname)[1]

                    if json_file["analysis_type"] == 'nucleotide':
                        if extension == '.fasta':
                            if path.exists(final_path):
                                RNAseq, identifier_seq, nucleotide_type = fastaToString(final_path)
                            else:
                                print('Could not find', fname)

                        if extension == '.gb' or extension == '.gbk':
                            if path.exists(final_path):
                                RNAseq, identifier_seq = gbToString(final_path)
                                nucleotide_type = True
                            else:
                                print('Could not find '+ fname)

                        if nucleotide_type:
                            # translate nucleotide to protein
                            protein_seq = nucleotideToAminoacid(RNAseq)

                            # Create the fasta file that will contain the resulting protein sequence
                            exportFasta(protein_seq, identifier_seq, final_output_path, fname)
                        else:
                            print(fname + ' does not contain a nucleotide sequence')

                    # Analysis_type = 'protein'
                    else:
                        if json_file["analysis_type"] == 'protein':
                            protein_seq = ""
                            fasta_seq_name = ""
                            if extension == '.fasta':
                                if path.exists(full_path):
                                    fasta_sequences = SeqIO.parse(open(full_path), 'fasta')
                                    for fasta in fasta_sequences:
                                        fasta_seq_name += fasta.id
                                        protein_seq += fasta.seq # concatenate all protein sequences

                                    is_nucleotide = isNucleotide(protein_seq)
                                else:
                                    print('Could not find', fname)
                            if extension == '.gb' or extension == '.gbk':
                                protein_list = []
                                is_nucleotide == False
                                if path.exists(full_path):
                                    try:
                                        for seq_record in SeqIO.parse(open(full_path, "r"), "genbank"):
                                            for seq_feature in seq_record.features:
                                                if seq_feature.type == "CDS":
                                                    protein_list.append(seq_feature.qualifiers['translation'][0])  # Saves protein sequences
                                        for i in protein_list:
                                                protein_seq += i

                                    except:
                                        print("Genbank file is incomplete")
                                else:
                                    print('Could not find the genbank file')

                            if is_nucleotide == False:
                                exportFasta(protein_seq, fasta_seq_name, final_output_path, fname)
                        else:
                            print('Incorrect analysis_type value in json file')
                    setting_file.close()











