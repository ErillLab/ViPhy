
from Bio.Seq import Seq, reverse_complement
from src.files import read_fasta_file, read_gb_file_as_nucleotide, read_gb_file_as_protein

import os


def is_nucleotide_true(file_content):
    """
    Indicates if the input value is a nucleotide sequence by returning a boolean
    :param file_content: fragment from a file that contains a nucleotide or protein sequence
    :return: Returns a boolean, whose value is True when the input value is a nucleotide sequence
    """
    is_nucleotide = True

    for sequence in file_content:
        for element in sequence:
            if element != 'A' and element != 'C' and element != 'G' and element != 'T':
                return False

    return is_nucleotide


def translate_to_protein(file_content):
    """
    Translates a nucleotide sequence into six amino acid sequences and then concatenate them into one longer sequence
    :param file_content: list of nucleotide sequences
    :return: Returns a single amino acid sequence that come from combining the six possible amino acid sequences
    """
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
                if len(codon1) == 3:
                    codon1 = Seq(codon1)
                    protein_seq1 += codon1.translate('1')
                if len(codon2) == 3:
                    codon2 = Seq(codon2)
                    protein_seq2 += codon2.translate('1')
                if len(codon3) == 3:
                    codon3 = Seq(codon3)
                    protein_seq3 += codon3.translate('1')
            protein_seq += protein_seq1 + protein_seq2 + protein_seq3

            if j == 0:
                nucleotide_seq = reverse_complement(nucleotide_seq)  # Return the reverse complement sequence
                nucleotide_seq = str(nucleotide_seq)

        final_protein_seq.append(protein_seq)
    return final_protein_seq


def export_fasta(id_lists, protein_seq_lists, working_folder):
    """
    Saves the protein sequence into a fasta file
    :param id_lists: list of lists. Contains the sequence identifiers
    :param protein_seq_lists: list of list. Contains the combination of the six possible protein sequences.
    Obtained using "nucleotide_to_protein" function
    :param working_folder: directory where fasta file will be saved
    """
    if len(id_lists) > 0:
        str_list = str(id_lists[0])
        if str_list.count("[") > 0:
            id_name = id_lists[0][0]
            working_file = "../" + working_folder + '/' + id_name + '.fasta'  # File name
        else:
            id_name = id_lists[0]
            working_file = "../" + working_folder + '/' + id_name + '.fasta'  # File name

        file_content = ""
        for sequence_list in protein_seq_lists:
            for sequence in sequence_list:
                for element in sequence:
                    file_content += str(element)

        with open(working_file, "w") as f:
            # Writes identifiers and sequences into the fasta file
            f.write(">" + str(id_name) + '\n' + str(file_content) + "\n")


def export_fasta_as_protein_folder(id_lists, protein_seq_lists, working_folder):
    """
    Saves protein sequences from different files into a single fasta file. This file shall contain each proteins
    separated by its own identifier
    :param id_lists: list of lists. Contains the sequence identifiers
    :param protein_seq_lists: list of list of lists. Contains a group of proteins.
    :param working_folder: directory where fasta file will be saved
    """

    working_file = "../" + working_folder + '/' + id_lists[0][0] + '.fasta'  # File name
    with open(working_file, "w") as f:
        for protein_double_list in protein_seq_lists:
            list_count = len(protein_double_list)
            for i in range(0, list_count):
                for j in range(0, len(protein_double_list[i])):
                    # Writes identifiers and sequences into the fasta file
                    f.write(">" + str(id_lists[i][j]) + '\n' + str(protein_double_list[i][j]) + "\n")


class Preprocess:
    """
    Class definition for Preprocess.
    Preprocess class reads the content of fasta and genbank files and standardize their data by translating nucleotide
    sequences into amino acid sequences when necessary.
    """

    def __init__(self):
        """Constructs an object for the preprocessing phase."""
        self.error_list = []

    def preprocessing_phase(self, file_name, input_folder, sequence_type, protein_type, working_folder):
        """
        Divides the preprocessing phase depending on the analysis type.
        :param file_name: file or folder found inside input folder. It should be a fasta or a genbank file
        :param input_folder: folder where the input files are stored
        :param protein_type: wheather to obtain protein sequences from annotations or from translating DNA
        :param working_folder: folder where the standardized files will be stored
        :return: Returns the error_list with the detected errors
        """
        file_extension = os.path.splitext(file_name)[1]  # Gets the file extension

        # Checks if the file is empty and it is not a folder
        if os.stat("../" + input_folder + '/' + file_name).st_size == 0 and not os.path.isdir("../" + input_folder +
                                                                                              '/' + file_name):
            self.error_list.append([file_name, "Its an empty file"])
        else:
            # Analyse DNA
            if sequence_type.lower() == 'dna':
                self.preprocessing_as_nucleotide(file_name, file_extension, input_folder, working_folder, sequence_type)
            
            # Analyse proteins
            elif sequence_type.lower() == 'protein':
                
                # Translate proteins
                if protein_type.lower() == 'translated':
                    self.preprocessing_as_nucleotide(file_name, file_extension, input_folder, working_folder, sequence_type)
                
                # Read annotated proteins
                elif protein_type.lower() == 'annotated':
                    self.preprocessing_as_protein(file_name, file_extension, input_folder, working_folder)
        
        return self.error_list

    def preprocessing_as_nucleotide(self, file, file_extension, input_folder, working_folder, sequence_type):
        """
        Gets nucleotide sequences from fasta and genbank files. If required, it transforms
        them into protein sequences. The final sequence obtained is saved into a fasta file.
        :param file: file or folder found inside input folder. It should be a fasta or a genbank file
        :param file_extension: file extension. If it is a folder it should be fasta or genbank but it can also be a folder
        :param input_folder: folder where the input files are stored
        :param working_folder: folder where the standardized files will be stored
        :sequence_type: defines what type of sequences to use ('DNA' or 'protein')
        :return: Returns the error_list with the detected errors
        """
        if file_extension in ['.fasta', '.fas', '.fa']:
            sequence_id, sequence = read_fasta_file(input_folder, file)
            if is_nucleotide_true(sequence):  # Checks if the fasta file contains a nucleotide sequence
                if sequence_type.lower() == 'protein':
                    seq_to_export = translate_to_protein(sequence)  # Translate to proteins
                elif sequence_type.lower() == 'dna':
                    seq_to_export = [sequence]
                else:
                    raise ValueError("Please, indicate a correct sequence_type. It should be 'DNA' or 'protein'.")
                export_fasta(sequence_id, seq_to_export, working_folder)  # Creates a fasta file
            else:
                self.error_list.append([file, "Does not contain a nucleotide sequence."])
        else:
            if file_extension in ['.gb', '.gbk', '.genbank']:
                sequence_id, sequence, error_list = read_gb_file_as_nucleotide(input_folder, file, self.error_list)
                if sequence_type.lower() == 'protein':
                    seq_to_export = translate_to_protein(sequence)
                elif sequence_type.lower() == 'dna':
                    seq_to_export = [sequence]
                else:
                    raise ValueError("Please, indicate a correct sequence_type. It should be 'DNA' or 'protein'.")
                export_fasta(sequence_id, seq_to_export, working_folder)
            elif os.path.isdir('../' + input_folder + '/' + file):  # If it finds a folder rather than a file
                subdirectory = input_folder + '/' + file
                sequence_id = []
                list_of_seq = []
                sequence_list = []
                for subdirectory_file in os.listdir('../' + subdirectory):
                    extension = os.path.splitext(subdirectory_file)[1]
                    if os.path.isdir('../' + subdirectory + '/' + subdirectory_file):  # If there is a folder
                        self.error_list.append([subdirectory_file, "There should not be a folder inside a subdirectory"
                                                                   ""])
                    # Checks file extension:
                    elif extension not in ['.fasta', '.fas', '.fa', '.gb', '.gbk', '.genbank']:
                        self.error_list.append([subdirectory_file, "It is not a genbank or fasta file"])
                    else:
                        if os.stat('../' + subdirectory + '/' + subdirectory_file).st_size == 0:  # If there is a empty file
                            self.error_list.append([subdirectory_file, "It is a empty file"])
                        else:
                            identifier, seq, self.error_list = read_gb_file_as_nucleotide(subdirectory,
                                                                                          subdirectory_file,
                                                                                          self.error_list)
                            if is_nucleotide_true(seq):
                                sequence_id.append(identifier)
                                sequence_list.append(seq)  # Concatenates all sequences from the subdirectory files
                for list_element in sequence_list:
                    if sequence_type.lower() == 'protein':
                        list_of_seq.append(translate_to_protein(list_element))
                    elif sequence_type.lower() == 'dna':
                        list_of_seq.append(list_element)
                    else:
                        raise ValueError("Please, indicate a correct sequence_type. It should be 'DNA' or 'protein'.")
                list_of_seq_lists = [list_of_seq]
                # Creates a fasta file with the final translated sequence
                export_fasta_as_protein_folder(sequence_id, list_of_seq_lists,
                                               working_folder)
            else:
                self.error_list.append([file, "It is not a genbank or a fasta file."])

    def preprocessing_as_protein(self, file, file_extension, input_folder, working_folder):
        """
        Gets protein sequences from fasta and genbank files and saves them into a fasta file
        :param file: file or folder found inside input folder. It should be a fasta or a genbank file
        :param file_extension: file extension. If it is a folder it should be fasta or genbank but it can also be a
        folder
        :param input_folder: folder where the input files are stored
        :param working_folder: folder where the standardized files will be stored
        :return: Returns the error_list with the detected errors
        """
        if file_extension in ['.fasta', '.fas', '.fa']:
            sequence_id, sequence = read_fasta_file(input_folder, file)
            if not is_nucleotide_true(sequence):
                export_fasta(sequence_id, sequence,
                             working_folder)  # Creates a fasta file with the original protein sequence
            else:
                self.error_list.append([file, "Does not contain a protein sequence. "])
        else:
            if file_extension in ['.gb', '.gbk', '.genbank']:
                sequence_id, sequence, self.error_list = read_gb_file_as_protein(input_folder, file, self.error_list)
                export_fasta(sequence_id, sequence,
                             working_folder)
            else:
                if os.path.isdir('../' + input_folder + '/' + file):  # If it finds a folder rather than a file
                    subdirectory = input_folder + '/' + file
                    sequence_id = []
                    protein_seq = []
                    sequence_list = []
                    for subdirectory_file in os.listdir('../' + subdirectory):
                        extension = os.path.splitext(subdirectory_file)[1]
                        if os.path.isdir('../' + subdirectory + '/' + subdirectory_file):  # If there is a folder
                            self.error_list.append([subdirectory_file, "There should not be a folder inside a "
                                                                       "subdirectory"])
                        # Checks file extension:
                        elif extension not in ['.fasta', '.fas', '.fa', '.gb', '.gbk', '.genbank']:
                            self.error_list.append([subdirectory_file, "It is not a genbank or fasta file"])
                        else:
                            if os.stat('../' + subdirectory + '/' + subdirectory_file).st_size == 0:  # If there is a empty file
                                self.error_list.append([subdirectory_file, "It is a empty file"])
                            else:
                                identifier, seq, self.error_list = read_gb_file_as_protein(subdirectory,
                                                                                           subdirectory_file,
                                                                                           self.error_list)
                                if not is_nucleotide_true(seq):
                                    sequence_id.append(identifier)
                                    sequence_list.append(seq)
                    protein_seq.append(sequence_list)
                    export_fasta_as_protein_folder(sequence_id, protein_seq, working_folder)
                else:
                    self.error_list.append([file, "It is not a genbank or a fasta file."])
