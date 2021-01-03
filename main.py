# 2020

from Bio import SeqIO, Entrez, Phylo
from Bio.Seq import Seq, reverse_complement, UnknownSeq
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastpCommandline
from Bio.Blast import NCBIXML
from ete3 import Tree
from io import StringIO

import json
import os
import re
import sys
import dendropy
import csv
import six

import warnings

warnings.simplefilter("ignore")


def access_ncbi(accessing_list, user_email, input_folder):
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
            folder_path = input_folder + '/' + list[0]
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
            print("Download unsuccessful! Cannot fetch " + list[0])
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


def read_gb_as_nucleotide(input_folder, gb_file, error_list):
    '''
    Reads and extracts a nucleotide sequence from a Genbank file
    :param input_folder: folder from where we obtain the genbank file to read
    :param gb_file: file that can store several sequences and extra information
    :param error_list: list of errors that have arisen during the sequence preprocessing
    :return: Returns a list with the identifier, another list with the nucleotide sequence and the list of errors
    '''
    record_id = []
    record_seq = []
    print('Reading ' + gb_file)

    for records in SeqIO.parse(input_folder + '/' + gb_file, "genbank"):
        sequence = records.seq
        if isinstance(sequence, UnknownSeq):
            error_list.append([gb_file, "There seems to be no sequence in this GenBank file"])
        else:
            record_id.append(records.id)
            record_seq.append(records.seq)

    return record_id, record_seq, error_list


def read_gb_file_as_protein(input_folder, gb_file):
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
                return False

    return is_nucleotide


def translate_to_protein(file_content):
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
                if (len(codon1) == 3):
                    codon1 = Seq(codon1)
                    protein_seq1 += codon1.translate('1')
                if (len(codon2) == 3):
                    codon2 = Seq(codon2)
                    protein_seq2 += codon2.translate('1')
                if (len(codon3) == 3):
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
    if len(id_lists) > 0:
        str_list = str(id_lists[0])
        if str_list.count("[") > 0:
            id_name = id_lists[0][0]
            working_file = working_folder + '/' + id_name + '.fasta'
        else:
            id_name = id_lists[0]
            working_file = working_folder + '/' + id_name + '.fasta'

        file_content = ""
        for sequence_list in protein_seq_lists:
            for sequence in sequence_list:
                for element in sequence:
                    file_content += str(element)

        f = open(working_file, "w")
        # Writes identifiers and sequences into the fasta file
        f.write(">" + str(id_name) + '\n' + str(file_content) + "\n")
        f.close()


def export_fasta_as_protein_folder(id_lists, protein_seq_lists, working_folder):
    '''
    Saves protein sequences from different files into a single fasta file. This file shall contain each proteins
    separated by its own identifier
    :param id_lists: list of lists. Contains the sequence identifiers
    :param protein_seq_lists: list of list of lists. Contains a group of proteins.
    :param working_folder: directory where fasta file will be saved
    '''
    element_count = 0

    working_file = working_folder + '/' + id_lists[0][0] + '.fasta'
    f = open(working_file, "w")
    for protein_double_list in protein_seq_lists:
        list_count = len(protein_double_list)
        for protein_list in protein_double_list:
            element_count = len(protein_list)
        for i in range(0, list_count):
            for j in range(0, element_count):
                # Writes identifiers and sequences into the fasta file
                f.write(">" + str(id_lists[i][j]) + '\n' + str(protein_double_list[i][j]) + "\n")
    f.close()


def preprocessing_phase(file_name, error_list):
    '''
    Divides the preprocessing phase depending on the analysis type.
    :param file_name: file or folder found inside input folder. It should be a fasta or a genbank file
    :param error_list: list of errors that have arisen during the sequence preprocessing
    :return: Returns the error_list with the detected errors
    '''
    file_extension = os.path.splitext(file_name)[1]  # Gets the file extension

    # Checks if the file is empty and it is not a folder
    if os.stat(input_folder + '/' + file_name).st_size == 0 and not os.path.isdir(input_folder + '/' + file_name):
        error_list.append([file, "Its an empty file"])
    else:
        # Checks if the analysis type indicated is "nucleotide"
        if analysis_type == 'nucleotide':
            error_list = preprocessing_as_nucleotide(file_name, file_extension, error_list)
        # Checks if the analysis type indicated is "protein"
        elif analysis_type == 'protein':
            error_list = preprocessing_as_protein(file_name, file_extension, error_list)
        else:
            print("Please, indicate a correct analysis type. It should be 'nucleotide' or 'protein'")
    return error_list


def preprocessing_as_nucleotide(file, file_extension, error_list):
    '''
    Gets nucleotide sequences from fasta and genbank files, transforms them into protein sequences and saves the final
    obtained sequence into a fasta file
    :param file: file or folder found inside input folder. It should be a fasta or a genbank file
    :param file_extension: file extension. If it is a folder it should be fasta or genbank but it can also be a folder
    :param error_list: list of errors that have arisen during the sequence preprocessing
    :return: Returns the error_list with the detected errors
    '''
    if file_extension in ['.fasta', '.fas', '.fa']:
        sequence_id, sequence = read_fasta_file(input_folder, file)
        if is_nucleotide_true(sequence):  # Checks if the fasta file contains a nucleotide sequence
            protein_seq = translate_to_protein(sequence)  # Translate to proteins
            export_fasta(sequence_id, protein_seq,
                         working_folder)  # Creates a fasta file with the translated sequence
        else:
            error_list.append([file, "Does not contain a nucleotide sequence."])
    else:
        if file_extension in ['.gb', '.gbk', '.genbank']:
            sequence_id, sequence, error_list = read_gb_as_nucleotide(input_folder, file, error_list)
            protein_seq = translate_to_protein(sequence)
            export_fasta(sequence_id, protein_seq, working_folder)
        elif os.path.isdir(input_folder + '/' + file):  # If it finds a folder rather than a file
            subdirectory = input_folder + '/' + file
            sequence_id = []
            protein_seq = []
            sequence_list = []
            for subdirectory_file in os.listdir(subdirectory):
                extension = os.path.splitext(subdirectory_file)[1]
                if os.path.isdir(subdirectory + '/' + subdirectory_file):  # If there is a folder
                    error_list.append([subdirectory_file, "There should not be a folder inside a subdirectory"])
                # Checks file extension:
                elif extension not in ['.fasta', '.fas', '.fa', '.gb', '.gbk', '.genbank']:
                    error_list.append([subdirectory_file, "It is not a genbank or fasta file"])
                else:
                    if os.stat(subdirectory + '/' + subdirectory_file).st_size == 0:  # If there is a empty file
                        error_list.append([subdirectory_file, "It is a empty file"])
                    else:
                        id, seq, error_list = read_gb_as_nucleotide(subdirectory, subdirectory_file, error_list)
                        if is_nucleotide_true(seq):
                            sequence_id.append(id)
                            sequence_list.append(seq)  # Concatenates all sequences from the subdirectory files
            for list in sequence_list:
                protein_seq.append(translate_to_protein(list))
            protein_seq_lists = [protein_seq]
            export_fasta_as_protein_folder(sequence_id, protein_seq_lists,
                                           working_folder)  # Creates a fasta file with the final translated sequence
        else:
            error_list.append([file, "It is not a genbank or a fasta file."])
    return error_list


def preprocessing_as_protein(file, file_extension, error_list):
    '''
    Gets protein sequences from fasta and genbank files and saves them into a fasta file
    :param file_name: file or folder found inside input folder. It should be a fasta or a genbank file
    :param file_extension: file extension. If it is a folder it should be fasta or genbank but it can also be a folder
    :param error_list: list of errors that have arisen during the sequence preprocessing
    :return: Returns the error_list with the detected errors
    '''
    if file_extension in ['.fasta', '.fas', '.fa']:
        sequence_id, sequence = read_fasta_file(input_folder, file)
        if not is_nucleotide_true(sequence):
            export_fasta(sequence_id, sequence,
                         working_folder)  # Creates a fasta file with the original protein sequence
        else:
            error_list.append([file, "Does not contain a protein sequence. "])
    else:
        if file_extension in ['.gb', '.gbk', '.genbank']:
            print('Reading ' + file)
            sequence_id, sequence = read_gb_file_as_protein(input_folder, file)
            export_fasta(sequence_id, sequence,
                         working_folder)
        else:
            if os.path.isdir(input_folder + '/' + file):  # If it finds a folder rather than a file
                subdirectory = input_folder + '/' + file
                sequence_id = []
                protein_seq = []
                sequence_list = []
                for subdirectory_file in os.listdir(subdirectory):
                    extension = os.path.splitext(subdirectory_file)[1]
                    if os.path.isdir(subdirectory + '/' + subdirectory_file):  # If there is a folder
                        error_list.append([subdirectory_file, "There should not be a folder inside a subdirectory"])
                    # Checks file extension:
                    elif extension not in ['.fasta', '.fas', '.fa', '.gb', '.gbk', '.genbank']:
                        error_list.append([subdirectory_file, "It is not a genbank or fasta file"])
                    else:
                        if os.stat(subdirectory + '/' + subdirectory_file).st_size == 0:  # If there is a empty file
                            error_list.append([subdirectory_file, "It is a empty file"])
                        else:
                            id, seq = read_gb_file_as_protein(subdirectory, subdirectory_file)
                            if not is_nucleotide_true(seq):
                                sequence_id.append(id)
                                sequence_list.append(seq)
                protein_seq.append(sequence_list)
                export_fasta_as_protein_folder(sequence_id, protein_seq, working_folder)
            else:
                error_list.append([file, "It is not a genbank or a fasta file."])
    return error_list


def display_error_messages(error_list):
    '''
    Shows the error messages collected during the preprocessing phase
    :param error_list: list of errors that have arisen during the sequence preprocessing
    '''
    if len(error_list) > 0:
        print('\n', "The following files has been ignored:")
        for f in error_list:
            print("-", f[0], ": ", f[1])  # File name + error
    else:
        print("No errors were found")


def make_blast_database(working_folder):
    '''
    Concatenats all the files extracted from the working folder and creates a BLAST database
    :param working_folder: directory from where we get the fasta files with protein sequences
    :return: Returns two created dictionaries
    '''
    db_folder = "dbFolder"
    db_file = "DataBase.fasta"
    id_list = []
    length_list = []

    if os.path.exists(db_folder):
        # Deletes old databases before creating a new one
        for file in os.listdir(db_folder):
            os.remove(db_folder + '/' + file)

    # Concatenates all the working_files
    id_list, length_list = concatenate_files(working_folder, db_file, id_list, length_list)

    # Create a dictionary to store the HSPs and another dictionary that will store the distance between sequences
    dict, distance_dictionary = dictionary_creation(id_list, length_list)

    db_name = db_file.split(".")  # Separates the file name and the file extension
    # Creates a blast database
    cmd = NcbimakeblastdbCommandline(input_file=db_file, out=db_folder + '/' + db_name[0], dbtype="prot")
    os.system(str(cmd))

    return dict, distance_dictionary


def concatenate_files(directory, db_file, id_list, length_list):
    '''
    Concatenates files extracted from a specific folder and saves important information that will be used to create a
    dictionary.
    :param directory: folder from where we get the files we want to concatenate
    :param db_file: name of the file that will contain the information to create the database
    :param id_list: list where the protein names will be stored
    :param length_list: list where the length of the protein sequence will be stored
    :return: Returns two lists. One with the protein names and another with the sequence length
    '''
    db_content = ""
    for working_file in os.listdir(directory):
        with open(directory + '/' + working_file) as f:
            db_content += f.read()
            db_content += '\n'
        for record in SeqIO.parse(directory + '/' + working_file, "fasta"):
            # Saves information that will be used to create a dictionary
            id_list.append(record.id)
            length_list.append(len(record.seq))

    # Writes the content in a separated file that will be used to create the database
    with open(db_file, 'w') as fW:
        fW.write(db_content)

    return id_list, length_list


def dictionary_creation(id_list, length_list):
    '''
    Creates a dictionary. It writes a pair key:value with the information extracted from a several fasta files. It
    will have a composite key and a vector of zeros as a internal value
    :param id_list: list where the protein names are stored
    :param length_list: list where the length of the protein sequence are stored
    :return: Returns the current dictionary
    '''
    dictionary = {}
    distance_dictionary = {}
    count = 0
    for i in id_list:
        l = [0] * length_list[count]  # List of zeros
        for j in id_list:
            key = i + '-' + j  # Dictionary key creation
            dictionary[key] = l
            # The initial distance value will be 1 to indicates that two sequences are completely different
            distance_dictionary[key] = 1.0
        count += 1
    return dictionary, distance_dictionary


def blastp(working_folder, e_value, dict):
    '''
    Compares a protein sequence from a fasta file to a protein sequence from the database
    :param working_folder: directory from where we get the fasta files
    :param e_value: number of expected hits of similar quality (score) that could be found just by chance
    :param dict: dictionary where we will save a vector with the hits position and its scores
    :return: Returns a dictionary that contains the coverage vectors
    '''
    for file in os.listdir(working_folder):
        seq_file_path = working_folder + '/' + file
        seq_file2 = seq_file_path.split(".")  # Removes the file extension from the file name
        seq_file_name = seq_file2[0].split("/")  # Removes the directory from the file name
        alignment_file = 'dbFolder/' + seq_file_name[1] + "_to_DB.xml"  # Output file name

        # Blatp of the file against the database
        cmd = NcbiblastpCommandline(query=seq_file_path, db="dbFolder/DataBase", evalue=e_value, outfmt=5,
                                    out=alignment_file)
        os.system(str(cmd))

        # Gets important information from the database to create a list
        dict.update(parse_xml_file(alignment_file, dict))

    return dict
    # print("Dictionary: ", dict)


def parse_xml_file(alignment_file, dictionary):
    '''
    Parses a blast xml file to get important information we will require to make a list with the hits after blast
    :param alignment_file: file obtained after Blast
    :param dictionary: current dictionary with the hits positions
    :return: Returns the current dictionary with the hits positions
    '''
    result_handle = open(alignment_file, "r")
    blast_records = NCBIXML.parse(result_handle)  # Parsing blastp object

    for rec in blast_records:
        query = rec.query
        query_length = rec.query_length
        for alignment in rec.alignments:
            hit = alignment.hit_def
            list = [0] * query_length  # List of zeros that will be used to replace an specific value in the dictionary

            for hsp in alignment.hsps:
                hit = hit.replace(' ', '')  # Deletes extra spaces in the sequence identifiers
                identity = hsp.identities  # Identical base pairs
                query_start = hsp.query_start
                query_end = hsp.query_end

                # Coverage process
                dictionary = coverage(list, query, hit, query_start, query_end, identity, dictionary)

    result_handle.close()

    return dictionary


def coverage(list, query, hit, query_start, query_end, identity, dict):
    '''
    Writes in a dictionary where the hits are indicated. If two hit are mapped to the same position, we will only write
    the biggest score
    :param list: list of zeros that will be modified with score values
    :param query: sequence that is compared with the database
    :param hit: sequence from the database
    :param query_start: position where the hit start
    :param query_end: position where the hit ends
    :param identity: value of the comparison of both sequences
    :param dict: current dictionary
    :return: Returns the updated dictionary and a boolean that indicates whether or not there is an overlap
    '''
    overlap = False
    # Fills the list in the positions where a hits happens with the score value
    for i in range(query_start, query_end + 1):
        if dict[query + '-' + hit][i - 1] == 0:
            list[i - 1] = identity
        elif dict[query + '-' + hit][i - 1] < identity:
            list[i - 1] = identity
        else:
            list[i - 1] = dict[query + '-' + hit][i - 1]
            overlap = True

    dict[query + '-' + hit] = list
    return dict


def get_distance(dictionary, distance_dictionary, distance_function):
    '''
    Calculates the distance between two sequences with the formula indicated by the user
    :param dictionary: dictionary that includes the coverage vectors
    :param distance_dictionary: dictionary where the distances will be saved
    :param distance_function: indicates the formula that will be used to get the distance
    '''
    if distance_function == 'd0':
        distance_dictionary = d0(dictionary, distance_dictionary)
    else:
        if distance_function == 'd4':
            distance_dictionary = d4(dictionary, distance_dictionary)
        elif distance_function == 'd6':
            distance_dictionary = d0(dictionary, distance_dictionary)
            print("D0: ", distance_dictionary)
            distance_dictionary = d4(dictionary, distance_dictionary)
            print("D4: ", distance_dictionary)
            distance_dictionary = d6(dictionary, distance_dictionary)
            print("D6: ", distance_dictionary)
        else:
            print("Please, introduce a valid distance function. It can be 'd0', 'd4' or 'd6'.")

    d_matrix, key_list = distance_matrix(distance_dictionary)
    # print("Matrix: ", d_matrix)
    return d_matrix, distance_dictionary, key_list


def d0(dict, distance_dictionary):
    '''
    Calculates the distance between two sequences using the HSPs lenght between them and their total length
    :param dict: Dictionary that contains the coverage vectors
    :param distance_dictionary: Dictionary that will contain the calculated distances
    :return: Returns the current distance dictionary
    '''

    for key in dict.keys():
        total_length = 0
        hit_length = 0
        # Sequences
        coverage_vector = dict[key]

        # Gets inverted key
        inverted_key = get_inverted_key(key)
        inverted_coverage_vector = dict[inverted_key]

        # Hits length
        hit_length += vector_no_zeros(coverage_vector)
        hit_length += vector_no_zeros(inverted_coverage_vector)

        # Total length
        total_length += vector_length(coverage_vector)
        total_length += vector_length(inverted_coverage_vector)

        # Distance formula
        distance_dictionary[key] = 1 - (hit_length / total_length)
    return distance_dictionary


def d4(dict, distance_dictionary):
    '''
    Calculates the distance between two sequences using the number of identical bases between them and the HSPs length
    :param dict: Dictionary that contains the coverage vectors
    :param distance_dictionary: Dictionary that will contain the calculated distances
    :return: Returns the current distance dictionary
    '''
    for key in dict.keys():
        hit_length = 0
        sum_identities = 0

        # Sequences
        coverage_vector = dict[key]

        # Gets inverted key
        inverted_key = get_inverted_key(key)
        inverted_coverage_vector = dict[inverted_key]

        # Identities
        sum_identities += identities(coverage_vector)
        sum_identities += identities(inverted_coverage_vector)

        # Hits length
        hit_length += vector_no_zeros(coverage_vector)
        hit_length += vector_no_zeros(inverted_coverage_vector)
        if hit_length == 0:
            hit_length = 1
        total_identities = (sum_identities / hit_length)

        # Distance formula
        distance_dictionary[key] = 1 - ((2 * total_identities) / hit_length)
    return distance_dictionary


def d6(dict, distance_dictionary):
    '''
    Calculates the distance between two sequences using the number of identical bases between them and the total
    length of the two sequences
    :param dict: Dictionary that contains the coverage vectors
    :param distance_dictionary: Dictionary that will contain the calculated distances
    :return: Returns the current distance dictionary
    '''
    for key in dict.keys():
        total_length = 0
        sum_identities = 0
        dif_zero = 0

        # Sequences
        coverage_vector = dict[key]

        # Gets inverted key
        inverted_key = get_inverted_key(key)
        inverted_coverage_vector = dict[inverted_key]

        # Identities
        sum_identities += identities(coverage_vector)
        sum_identities += identities(inverted_coverage_vector)
        dif_zero += vector_no_zeros(coverage_vector)
        dif_zero += vector_no_zeros(inverted_coverage_vector)
        if dif_zero == 0:
            dif_zero = 1
        total_identities = (sum_identities / dif_zero)

        # Total length
        total_length += vector_length(coverage_vector)
        total_length += vector_length(inverted_coverage_vector)

        # Distance formula
        distance_dictionary[key] = 1 - ((2 * total_identities) / total_length)
    return distance_dictionary


def get_inverted_key(key):
    '''
    Gets the inverted key from an initial key
    :param key: key that allows us to access to an specific coverage vector
    :return: Returns the inverted key
    '''
    separated_key = key.split('-')
    inverted_key = separated_key[1] + '-' + separated_key[0]
    return inverted_key


def vector_length(coverage_vector):
    '''
    Calculates the length of an array or list
    :param coverage_vector: Vector that includes all the HSPs between two sequences
    :return: Returns the total length of the input vector
    '''
    length = len(coverage_vector)
    return length


def vector_no_zeros(coverage_vector):
    '''
    Calculate the length of the HSPs
    :param coverage_vector: Vector that includes all the HSPs between two sequences
    :return: Returns the number of positions that are not equal to zero
    '''
    dif_zero = 0
    for base in coverage_vector:
        if base != 0:
            dif_zero += 1

    return dif_zero


def identities(coverage_vector):
    '''
    Calculates the amount of characters which match exactly between two different sequences
    :param coverage_vector: Vector that includes all the HSPs between two sequences
    :return: Returns the number of identities or exact matches between two sequences
    '''
    sum_identities = 0
    for base in coverage_vector:
        sum_identities += base
    return sum_identities


def distance_matrix(dictionary):
    '''
    Creates a matrix with the distances calculated and stored in a dictionary
    :param dictionary: dictionary where the distance values are saved
    :return: Returns the final distance matrix created with all the data in the input dictionary
    '''
    key_list = []
    # Creates a list with all the possible keys
    for key in dictionary.keys():
        key_parts = key.split('-')
        if key_parts[0] not in key_list:
            key_list.append(key_parts[0])

    matrix = []
    # Create the matrix distance with the calculated distances
    for first_key in key_list:
        list = []
        for second_key in key_list:
            list.append(dictionary[first_key + '-' + second_key])
        matrix.append(list)

    return matrix, key_list


def write_csv(d_matrix, key_list, distance_file):
    '''
    Writes the content of the distance matrix in a file with extension 'csv'
    :param d_matrix: matrix that contains all the distances between sequences
    :param key_list: list with the known key of the sequences
    :param distance_file: name of the file where the content will be written
    '''
    row = 0
    aux_matrix = []

    # Creates a new matrix combining keys and distances
    k_list = [' '] + key_list
    aux_matrix.append(k_list)
    for matrix_list in d_matrix:
        matrix_list = [key_list[row]] + matrix_list
        aux_matrix.append(matrix_list)
        row += 1

    # Writes matrix in a new file
    with open(distance_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(aux_matrix)


def get_phylogenetic_tree(distance_file):
    '''
    Creates a phylogenetic tree from csv file that contained a distance matrix
    :param distance_file: name of the file where the content was written
    :return: Returns the created phylogenetic tree
    '''
    # Reads the content of a file with csv extension
    pdm1 = dendropy.PhylogeneticDistanceMatrix.from_csv(src=open(distance_file), delimiter=",")
    nj_tree = pdm1.nj_tree()
    tree_data = nj_tree.as_string("newick")  # Phylogenetic tree format
    separated_tree_data = tree_data.split(" ")

    # Create tree using Biopython
    tree = Phylo.read(StringIO(separated_tree_data[1]), "newick")
    print(tree)

    # Creates the tree
    tree = Tree(separated_tree_data[1])
    return tree


if __name__ == '__main__':

    # Configuration file
    setting_file = "settings.json"
    try:
        json_file = json.load(open(setting_file))  # Reads json file
    except IOError:
        sys.exit('Could not open settings.json')

    # Important folders
    input_folder = json_file["input_folder"]
    if not os.path.exists(input_folder):
        os.mkdir(input_folder)
    working_folder = json_file["working_folder"]
    if not os.path.exists(working_folder):
        os.mkdir(working_folder)
    output_folder = json_file["output_folder"]
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # Remaining information of the configuration file
    analysis_type = json_file["analysis_type"]
    access_ncbi_list = json_file["genome_accessions"]
    user_email = json_file["user_email"]
    e_value = json_file["e_value"]
    distance_function = json_file["distance_function"]

    # Get NCBI files
    access_ncbi(access_ncbi_list, user_email, input_folder)

    # Deletes all content from working folder
    # for file in os.listdir(working_folder):
    #    os.remove(working_folder + '/' + file)

    file_error_list = []
    # Preprocessing phase
    for file in os.listdir(input_folder):  # Navigates into the input_folder
        file_error_list = preprocessing_phase(file, file_error_list)

    # Displays a list of error detected in the preprocessing code
    display_error_messages(file_error_list)

    # Builts a database
    dictionary, distance_dictionary = make_blast_database(working_folder)

    print("This might take a few minutes...")
    # Blastp
    dictionary = blastp(working_folder, e_value, dictionary)

    # Distance
    d_matrix, distance_dictionary, key_list = get_distance(dictionary, distance_dictionary, distance_function)
    print(d_matrix)

    # Phylogenetic tree
    distance_file = "distance_matrix.csv"
    write_csv(d_matrix, key_list, distance_file)  # External file that contains the distance matrix

    # Create tree
    tree = get_phylogenetic_tree(distance_file)
    print("Tree: ", tree)


