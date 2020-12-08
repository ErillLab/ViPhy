# 2020

from Bio import SeqIO, Entrez
from Bio.Seq import Seq, reverse_complement, UnknownSeq
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastpCommandline
from Bio.Blast import NCBIXML
import numpy as np
import json
import os
import re
import sys
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
            working_file = working_folder + '/' + id_lists[0][0] + '.fasta'
            id_name = id_lists[0][0]
        else:
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


def export_fasta_as_nucleotide_folder(id_lists, protein_seq_lists, working_folder):
    print("Hello0")


def export_fasta_as_protein_folder(id_lists, protein_seq_lists, working_folder):
    '''
    Saves protein sequences from different files into a single fasta file. This file shall contain each proteins
    separated by its own identifier
    :param id_lists: list of lists. Contains the sequence identifiers
    :param protein_seq_lists: list of list of lists. Contains a group of proteins.
    :param working_folder: directory where fasta file will be saved
    '''
    working_file = working_folder + '/' + id_lists[0][0] + '.fasta'
    f = open(working_file, "w")

    for protein_double_list in protein_seq_lists:
        list_count = len(protein_double_list)
        for protein_list in protein_double_list:
            element_count = len(protein_list)
        for i in range(0, list_count):
            for j in range(0, element_count):
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
            protein_seq_lists = []
            protein_seq_lists.append(protein_seq)
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
    :return: Returns the lists with the protein names and the protein length and the current dictionary
    '''
    db_folder = "Protein_db"
    db_file = "DB.fasta"
    id_list = []
    length_list = []
    dict = {}

    if os.path.exists(db_folder):
        # Deletes old databases before creating a new one
        for file in os.listdir(db_folder):
            os.remove(db_folder + '/' + file)

    # Concatenates all the working_files
    id_list, length_list, dict = concatenate_files(working_folder, db_file, id_list, length_list, dict)

    # Create a nested dictionary
    dict = create_dictionary(id_list, length_list, dict)

    db_name = db_file.split(".")  # Separates the file name and the file extension
    # Creates a blast database
    cmd = NcbimakeblastdbCommandline(input_file=db_file, out=db_folder + '/' + db_name[0], dbtype="prot")
    os.system(str(cmd))

    return id_list, length_list, dict


def concatenate_files(directory, db_file, id_list, length_list, dict):
    '''
    Concatenates files extracted from a specific folder and saves important information that will be used to create a
    dictionary and a database.
    :param directory: folder from where we get the files we want to concatenate
    :param db_file: name of the file that will contain the information to create the database
    :param id_list: list where the protein names will be stored
    :param length_list: list where the length of the protein sequence will be stored
    :param dict: empty dictionary
    :return: Returns the lists with the protein names and length and the current dictionary with external keys
    '''
    db_content = ""
    for filen in os.listdir(directory):
        with open(directory + '/' + filen) as f:
            db_content += f.read()
        for record in SeqIO.parse(directory + '/' + filen, "fasta"):
            # Saves information that will be used to create a dictionay
            id_list.append(record.id)
            length_list.append(len(record.seq))
            # Creates the external keys of dictionary
            dict[str(id_list[-1])] = 0

    # Writes the content that will be used to create the database
    with open(db_file, 'w') as fW:
        fW.write(db_content)

    return id_list, length_list, dict


def create_dictionary(id_list, length_list, dictionary):
    '''
    Gets a normal dictionary and creates a nested dictionary. It writes a pair key:value with the information extracted
    from a several fasta files. It will have an external key and then a internal dictionary with the protein names as a
    internal key and a vector of zeros as a internal value
    :param id_list: list where the protein names are stored
    :param length_list: list where the length of the protein sequence are stored
    :param dictionary: current dictionary
    :return: Returns the current nested dictionary
    '''
    internal_dict = {}
    count = 0
    for external_loop in id_list:
        for internal_loop in id_list:
            # Creates a vector of zeros with the same size of the sequence and saves it in the internal dictionary value
            internal_dict[internal_loop] = np.zeros(length_list[count])
            count+= 1
        dictionary[external_loop] = internal_dict
        internal_dict = {}
        count = 0
    return dictionary


def blastp(working_folder, e_value, id_list, length_list, dict):
    '''
    Compares a protein sequence to a protein database
    :param working_folder: directory where fasta file were saved after the preprocessing process
    :param e_value: number of expected hits of similar quality (score) that could be found just by chance
    :param id_list: list of protein names extracted from fasta files
    :param length_list: list of the length of each protein sequence
    :param dict: dictionary where we will save a vector with the hits position and its scores
    :return:
    '''
    vector_list = []
    internal_range = []
    external_range = []
    aux_list = []
    for i in range(len(length_list)):
        aux_list.append([])
    aux_length_list, external_range = internal_dictionary_length(internal_range, external_range, length_list)
    for seq in os.listdir(working_folder):
        seq_file = working_folder + '/' + seq  # File name
        seq_file2 = seq_file.split(".")  # Removes the file extension from the file name
        seq_file_name = seq_file2[0].split("/")  # Removes the directory from the file name
        aligment_file = 'Protein_db/' + seq_file_name[1] + "_to_DB.xml"  # Output file name
        cmd = NcbiblastpCommandline(query=seq_file, db="Protein_db/DB", evalue=e_value, outfmt=5,
                                    out=aligment_file)

        os.system(str(cmd))
        vector = []
        vector = parse_file(aligment_file, vector)
        vector_list.append(vector)
        print(vector)

    update_dictionary(vector_list, id_list, length_list, dict)


def internal_dictionary_length(specific_range, range_list, length_list):
    '''
    Saves the positions of the sequences in the database in a list
    :param specific_range: list where the positions will be saved
    :param range_list: list that contains the specific_range list
    :param length_list: number of files that the database contains
    :return: Return the current list of the sequences position in the database and a list of lists with the same
    positions but empty
    '''
    piece = 0
    aux_length_list = []
    aux_length_list.append(0)  # First range start
    for length in length_list:
        aux_length_list.append(piece+length)
        piece += length
        range_list.append(specific_range)
    return aux_length_list, range_list


def parse_file(aligment_file, vector):
    '''
    Parses a blast xml file to get important information we will requiere to make a list with the hits agter blast
    :param aligment_file: file obtained after Blast
    :param vector: input list that will be modified with the hits after blast
    :return: Return a list with the position of the hits agains the database
    '''
    result_handle = open(aligment_file, "r")
    blast_records = NCBIXML.parse(result_handle)
    for rec in blast_records:
        length = rec.num_letters_in_database
        for alignment in rec.alignments:
            for hsp in alignment.hsps:
                score = hsp.score
                # evalue = hsp.expect
                hit_start = hsp.sbjct_start
                hit_end = hsp.sbjct_end
                vector = coverage(length, score, hit_start, hit_end, vector)
    result_handle.close()
    return vector


def coverage(length, score, hit_start, hit_end, hit_vector):
    '''
    Creates a list where the hits will be indicated. If two hit are mapped to the same position, we will take the
    biggest score
    :param length: length of the concatenated database sequences
    :param score: value of the hit
    :param hit_start: indicates where the hit start
    :param hit_end: indicates where the hit ends
    :param hit_vector: list with the same lenght of the database where hits are indicated
    :return: Returns list with the same lenght of the database. The 0s indicates no hit.
    '''
    for i in range(0, length):
        hit_vector.append(0)  # Creates a list full of 0s

    for i in range(hit_start, hit_end + 1, 1):
        if 0 < hit_vector[i - 1]:
            if hit_vector[i - 1] < score:
                hit_vector[i - 1] = score
        else:
            hit_vector[i - 1] = score
    return hit_vector


def update_dictionary(vector_list, id_list, length_list, dict):
    num_vector = 0
    for external_loop in id_list:
        position = 0
        for internal_loop in id_list:
            dict[external_loop][internal_loop] = vector_list[num_vector][position:(position+length_list[position])]
            position+= 1
        num_vector+= 1
    # print(dict)


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
    id_list, length_list, dict = make_blast_database(working_folder)

    print("This might take a few minutes...")
    # Blastp
    blastp(working_folder, e_value, id_list, length_list, dict)
