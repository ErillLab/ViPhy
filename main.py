# 2020

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement, UnknownSeq
from Bio import Entrez
import json
import os
import os.path as path
import re
import sys


def accessing_ncbi(accessing_list, user_email, input_folder):
    '''
    Goes through the accessing list to find the files that the user whats to download from NCBI and call
    another function to do so.
    :param accessing_list: list of lists that contains the elements we want to download
    :param user_email: email that will be used to identify the user in NCBI
    :param input_folder: folder where files will be stored
    '''
    for list in accessing_list:
        list_length = len(list)
        if list_length == 1:
            file_name = list[0] + ".gb"
            file_path = input_folder + '/' + file_name

            if not os.path.exists(file_path):  # Checks if the file was downloaded before
                print('Downloading ' + file_name)
                download_file(user_email, file_path, list)

        if list_length > 1:
            folder_name = list[0].split('.')  # Deletes possible dots in the folder name
            folder_path = input_folder + '/' + folder_name[0]

            if not os.path.exists(folder_path):  # Checks if there is a folder with the same name
                os.mkdir(folder_path)  # Creates a folder
            for position in range(list_length):
                file_name = list[position] + ".gb"
                file_path = folder_path + '/' + file_name
                if not os.path.exists(file_path):  # Checks if the file was downloaded before
                    print('Downloading ' + file_name)
                    download_file(user_email, file_path, list)


def download_file(user_email, file_path, list):
    '''
    Downloads files form NCBI in GenBank plain text format and saves it into the input folder
    :param user_email: email that will be used to identify the user in NCBI
    :param file_path: input folder path or subdirectory where the files will be saved
    :param list: list that contains the files to download
    '''
    match = re.search(r'[\w.-]+@[\w.-]+.\w+', user_email)
    if match:
        Entrez.email = user_email  # Always tell NCBI who you are
        try:
            handle = Entrez.efetch(db="nucleotide", id=list[0], rettype="gbwithparts", retmode="text")
            data = handle.read()
            f = open(file_path, "w")
            f.write(data)
            f.close()
        except Exception:
            print("Error! Cannot fetch " + list[0])
    else:
        print("Please, introduce a valid email address format")


def read_fasta_file(input_folder, fasta_file):
    '''
    Reads a fasta file to obtain a nucleotide or protein sequence.
    :param input_folder: folder from where we obtain the file to read
    :param fasta_file: file that contains a nucleotide or protein sequence, should be FASTA format
    :return: Returns a list with the identifier and another list with the nucleotide sequence
    '''
    record_seq = []
    record_id = []
    print('Reading ' + fasta_file)

    # This loop allows us to concatenate multiple sequences
    for record in SeqIO.parse(input_folder + '/' + fasta_file, "fasta"):
        record_id.append(record.id)
        record_seq.append(record.seq)

    return record_id, record_seq


def read_gb_file(input_folder, gb_file):
    '''
    Reads and extracts a nucleotide sequence from a Genbank file
    :param input_folder: folder from where we obtain the genbank file to read
    :param gb_file: file that can store several sequences and extra information
    :return: Returns a list with the identifier and another list with the nucleotide sequence
    '''
    record_id = []
    record_seq = []
    print('Reading ' + gb_file)
    for records in SeqIO.parse(input_folder + '/' + gb_file, "genbank"):
        sequence = records.seq
        if isinstance(sequence, UnknownSeq):
            sys.exit("There seems to be no sequence in your GenBank file!")
        else:
            record_id.append(records.id)
            record_seq.append(records.seq)

    return record_id, record_seq


def get_protein(input_folder, gb_file):
    '''
    Reads and extracts protein sequences from a Genbank file
    :param input_folder: folder from where we obtain the genbank file to read
    :param gb_file: file that can store several sequences and extra information, including protein sequences
    :return: Returns a list of identifiers and another of protein sequences obtained from the genbank file
    '''
    id_list = []
    protein_list = []
    for seq_record in SeqIO.parse(input_folder + '/' + gb_file, "genbank"):
        for seq_feature in seq_record.features:
            if seq_feature.type == "CDS":
                id_list.append(seq_feature.qualifiers['protein_id'][0])
                protein_list.append(
                    seq_feature.qualifiers['translation'][0])  # Saves protein sequences

    return id_list, protein_list


def is_nucleotide_true(fasta_content):
    '''
    Indicates if the input value is nucleotide sequence by returning a boolean
    :param fasta_content: content from a fasta file that contains a nucleotide or protein sequence
    :return: Returns a boolean, whose value is True when the input value is a nucleotide sequence
    '''
    is_nucleotide = True

    for sequence in fasta_content:
        for element in sequence:
            if element != 'A' and element != 'C' and element != 'G' and element != 'T':
                is_nucleotide = False

    return is_nucleotide


def nucleotide_to_protein(file_content):
    '''
    Translates a nucleotide sequence into six amino acid sequences and then concatenate them into one longer sequence
    :param file_content: list of nucleotide sequences
    :return: Returns a single amino acid sequence that come from combining the six possible amino acid sequences
    '''
    final_protein_seq = []

    for sequence in file_content:
        directions = 2
        protein_seq = ""
        nucleotide_seq = str(sequence)

        # First direction: j = 0 / Second direction: j = 1
        for j in range(0, directions):
            protein_seq1 = ""
            protein_seq2 = ""
            protein_seq3 = ""
            for i in range(0, len(nucleotide_seq), 3):
                codon1 = nucleotide_seq[i:i + 3]
                codon2 = nucleotide_seq[i + 1:i + 4]
                codon3 = nucleotide_seq[i + 2:i + 5]
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

            if j == 0:
                nucleotide_seq = reverse_complement(nucleotide_seq)  # Return the reverse complement sequence
                nucleotide_seq = str(nucleotide_seq)

        final_protein_seq.append(protein_seq)
    return final_protein_seq


def export_fasta(id_lists, protein_seq_lists, working_folder):
    '''
    Saves the protein sequence into a fasta file
    :param id_lists: list of lists. Contains the sequence identifiers
    :param protein_seq_lists: list of list. Contains the combination of the six possible protein sequences.
    Obtained using "nucleotide_to_protein" function
    :param working_folder: directory where fasta file will be saved
    '''
    id_lists = correct_format(id_lists)
    working_file = working_folder + '/' + id_lists[0] + '.fasta'

    f = open(working_file, "w")
    position_id_list = 0
    for sequence_list in protein_seq_lists:
        file_content = ""

        for sequence in sequence_list:
            for element in sequence:
                file_content += str(element)

        f.write(">" + str(id_lists[position_id_list]) + '\n' + str(file_content) + "\n")
        position_id_list += 1

    f.close()


def correct_format(name_list):
    '''
    Transforms identifier names to ensure that follows the standard format. Should be a list.
    :param name_list: file names list to check and correct if necessary
    :return: correct list with file names
    '''
    correct_name = []
    for name in name_list:
        new_name= str(name).replace("['", "")
        new_name = new_name.replace("']", "")
        correct_name.append(new_name)
    return correct_name


if __name__ == '__main__':

    # Configuration file
    setting_file = ""
    try:
        setting_file = open('settings.json')
    except IOError:
        print('Could not open settings.json')
    else:
        # Reads json file
        json_file = json.load(setting_file)
        user_email = json_file["user_email"]
        input_folder = json_file["input_folder"]
        if not os.path.exists(input_folder):
            os.mkdir(input_folder)
        working_folder = json_file["working_folder"]
        if not os.path.exists(working_folder):
            os.mkdir(working_folder)
        accessing_list = json_file["genome_accessions"]

        # Get NCBI files
        accessing_ncbi(accessing_list, user_email, input_folder)

        # Deletes all content from working folder
        for file in os.listdir(working_folder):
            os.remove(working_folder + '/' + file)

        # Navigates into the input_folder
        for file in os.listdir(input_folder):
            file_extension = path.splitext(file)[1]

            # Analysis_type = 'nucleotide'
            if json_file["analysis_type"] == 'nucleotide':
                if file_extension == '.fasta':
                    sequence_id, sequence = read_fasta_file(input_folder, file)
                    is_nucleotide = is_nucleotide_true(sequence)  # Checks if the fasta file contains a nucleotide sequence
                    if is_nucleotide:
                        protein_seq = nucleotide_to_protein(sequence)  # Translate to proteins
                        export_fasta(sequence_id, protein_seq,
                                     working_folder)  # Creates a fasta file with the translated sequence
                    else:
                        print('--> ' + file + ' does not contain a nucleotide sequence')
                else:
                    if file_extension == '.gb' or file_extension == '.gbk':
                        sequence_id, sequence = read_gb_file(input_folder, file)
                        protein_seq = nucleotide_to_protein(sequence)
                        export_fasta(sequence_id, protein_seq, working_folder)
                    else:
                        if file_extension == '':  # If it finds a folder rather than a file
                            subdirectory = input_folder + '/' + file
                            sequence_id = []
                            protein_seq = []
                            sequence_list = []
                            for subdirectory_file in os.listdir(subdirectory):
                                id, seq = read_gb_file(subdirectory, subdirectory_file)
                                sequence_id.append(id)
                                sequence_list.append(seq)  # Concatenates sequences from all the files in the subdirectory
                            for list in sequence_list:
                                protein_seq.append(nucleotide_to_protein(list))

                            export_fasta(sequence_id, protein_seq,
                                         working_folder)  # Creates a fasta file with the final translated sequence
                        else:
                            print(file + " is not a genbank or a fasta file")

            # Analysis_type = 'protein'
            else:
                if json_file["analysis_type"] == 'protein':
                    if file_extension == '.fasta':
                        sequence_id, sequence = read_fasta_file(input_folder, file)
                        is_nucleotide = is_nucleotide_true(sequence)
                        if not is_nucleotide:
                            export_fasta(sequence_id, sequence,
                                         working_folder)  # Creates a fasta file with the original protein sequence
                        else:
                            print('--> ' + file + ' does not contain a protein sequence')
                    else:
                        if file_extension == '.gb' or file_extension == '.gbk':
                            print('Reading ' + file)
                            sequence_id, sequence = get_protein(input_folder, file)
                            export_fasta(sequence_id, sequence,
                                         working_folder)
                        else:
                            if file_extension == '':  # If it finds a folder rather than a file
                                subdirectory = input_folder + '/' + file
                                sequence_id = []
                                protein_seq = []
                                sequence_list = []
                                for subdirectory_file in os.listdir(subdirectory):
                                    id, seq = read_gb_file(subdirectory, subdirectory_file)
                                    sequence_id.append(id)
                                    sequence_list.append(seq)

                                for list in sequence_list:
                                    protein_seq.append(nucleotide_to_protein(list))

                                export_fasta(sequence_id, protein_seq,working_folder)
                            else:
                                print(file + " is not a genbank or a fasta file")
                else:
                    print("Please, indicate a correct analysis type. It should be `nucleotide´ or `protein´")

        setting_file.close()









