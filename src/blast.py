
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastpCommandline, NcbiblastnCommandline
from Bio.Blast import NCBIXML

import os


class Blast:
    """
    Class definition for Blast.
    Blast class defines the second phase of the ViPhy process. During this phase a BLAST database is created using the
    available sequences. Then, the sequences will be compared with each other in order to find regions of similarity
    between each sequence pair (an input sequence and the database sequence) in a matter of seconds. More specifically,
    the the class supports BLASTP, a version of BLAST that compares two protein sequences.

    It also creates a coverage vector that includes the positions of the high score segment pairs (HSPs) between two
    aligned sequences.
    """

    def __init__(self):
        """Constructs an object for the alignment phase."""
        self.db_folder = "dbFolder"
        self.db_file = "DataBase.fasta"
        self.id_list = []  # list where the protein names are stored
        self.length_list = []  # list where the length of the protein sequence are stored
        self.dictionary = {}  # dictionary that contains the coverage vector between each sequence pair
        self.distance_dictionary = {}  # dictionary that will contain the distance between sequences

    def make_blast_database(self, sequence_type, working_folder, e_value, blast_word_size):
        """
        Concatenates all the files extracted from the working folder and creates a BLAST database
        :param working_folder: directory from where we get the fasta files with protein sequences
        :param e_value: number of expected hits of similar quality (score) that could be found just by chance
        :return: Returns two created dictionaries
        """
        if os.path.exists("../" + self.db_folder):
            # Deletes old databases before creating a new one
            for file_to_delete in os.listdir("../" + self.db_folder):
                os.remove("../" + self.db_folder + '/' + file_to_delete)

        # Concatenates all the working_files
        self.concatenate_files(working_folder)

        # Create a dictionary to store the HSPs and another dictionary that will store the distance between sequences
        self.dictionary_creation()

        db_name = self.db_file.split(".")  # Separates the file name and the file extension
        # Creates a blast database
        if sequence_type.lower() == 'protein':
            blast_db_type = "prot"
        elif sequence_type.lower() == 'dna':
            blast_db_type = "nucl"
        else:
            raise ValueError("Please, indicate a correct sequence_type. It should be 'DNA' or 'protein'.")
        
        cmd = NcbimakeblastdbCommandline(input_file="../" + self.db_file, out="../" + self.db_folder + '/' + db_name[0],
                                         dbtype=blast_db_type)
        cmd()
        print("BLAST database has been created")

        # Uses Blastp to the coverage_vector for each sequence pair
        print("The sequence alignment might take some time")
        self.coverage_vector_collection(sequence_type, working_folder, e_value, blast_word_size)

        return self.distance_dictionary, self.dictionary

    def concatenate_files(self, directory):
        """
        Concatenates files extracted from a specific folder and saves important information that will be used to create
        a dictionary.
        :param directory: folder from where we get the files we want to concatenate
        :return: Returns two lists. One with the protein names and another with the sequence length
        """
        data_base_content = ""
        for working_file in os.listdir("../" + directory):
            with open("../" + directory + '/' + working_file) as f:
                data_base_content += f.read()
                data_base_content += '\n'
            for record in SeqIO.parse("../" + directory + '/' + working_file, "fasta"):
                # Saves information that will be used to create a dictionary
                self.id_list.append(record.id)
                self.length_list.append(len(record.seq))

        # Writes the content in a separated file that will be used to create the database
        with open("../" + self.db_file, 'w') as fW:
            fW.write(data_base_content)

    def dictionary_creation(self):
        """
        Creates a dictionary. It writes a pair key:value with the information extracted from a several fasta files. It
        will have a composite key and a vector of zeros as a internal value.
        It also creates a second dictionary that will contain the distances between sequences.
        :return: Returns a dictionary with a list of zeros and another with initial distances.
        """
        count = 0

        for i in self.id_list:
            zero_list = [0] * self.length_list[count]  # List of zeros
            for j in self.id_list:
                key = i + '-' + j  # Dictionary key creation
                self.dictionary[key] = zero_list
                # The initial distance value will be 1 to indicates that two sequences are completely different
                self.distance_dictionary[key] = 1.0
            count += 1

    def coverage_vector_collection(self, sequence_type, working_folder, e_value, blast_word_size):
        """
        Compares a sequence from a fasta file to a sequence from the database
        :param working_folder: directory from where we get the fasta files
        :param e_value: number of expected hits of similar quality (score) that could be found just by chance
        :return: Returns a dictionary that contains the coverage vectors
        """
        for file in os.listdir("../" + working_folder):
            seq_file_path = "../" + working_folder + '/' + file
            seq_file2 = seq_file_path.split(".")  # Removes the file extension from the file name
            seq_file_name = seq_file2[2].split("/")  # Removes the directory from the file name
            alignment_file = "../dbFolder/" + seq_file_name[2] + "_to_DB.xml"  # Final file path

            # Alignment using BLAST of the file against the database            
            if sequence_type.lower() == 'protein':
                # blastp command
                if blast_word_size.lower() == 'default':
                    cmd = NcbiblastpCommandline(query=seq_file_path, db="../dbFolder/DataBase", evalue=e_value, outfmt=5,
                                                out="../dbFolder/" + alignment_file)
                else:
                    cmd = NcbiblastpCommandline(query=seq_file_path, db="../dbFolder/DataBase", evalue=e_value, outfmt=5,
                                                out="../dbFolder/" + alignment_file, word_size=blast_word_size)
            elif sequence_type.lower() == 'dna':
                # blastn command
                if blast_word_size.lower() == 'default':
                    cmd = NcbiblastnCommandline(query=seq_file_path, db="../dbFolder/DataBase", evalue=e_value, outfmt=5,
                                                out="../dbFolder/" + alignment_file)
                else:
                    cmd = NcbiblastnCommandline(query=seq_file_path, db="../dbFolder/DataBase", evalue=e_value, outfmt=5,
                                                out="../dbFolder/" + alignment_file, word_size=blast_word_size)
            # Run command
            cmd()

            # Gets important information from the database to create a list
            self.dictionary.update(self.parse_xml_file(alignment_file))

    def parse_xml_file(self, alignment_file):
        """
        Parses a blast xml file to get important information we will require to make a list with the hits after blast
        :param alignment_file: file obtained after Blast
        :return: Returns the current dictionary with the hits positions
        """
        result_handle = open(alignment_file, "r")
        blast_records = NCBIXML.parse(result_handle)  # Parsing blastp object

        for rec in blast_records:
            query = rec.query
            query_length = rec.query_length
            for alignment in rec.alignments:
                hit = alignment.hit_def
                # List of zeros that will be used to replace an specific value in the dictionary
                zero_list = [0] * query_length

                for hsp in alignment.hsps:
                    hit = hit.replace(' ', '')  # Deletes extra spaces in the sequence identifiers
                    identity = hsp.identities / hsp.align_length  # Fraction of identical base pairs
                    query_start = hsp.query_start
                    query_end = hsp.query_end

                    # Creates coverage vector
                    self.coverage_vector_dictionary(zero_list, query, hit, query_start, query_end, identity)

        result_handle.close()
        return self.dictionary

    def coverage_vector_dictionary(self, zero_list, query, hit, query_start, query_end, identity):
        """
        Writes in a dictionary where the hits are indicated. If two hit are mapped to the same position, we will only
        write the biggest score
        :param zero_list: list of zeros that will be modified with score values
        :param query: sequence that is compared with the database
        :param hit: sequence from the database
        :param query_start: position where the hit start
        :param query_end: position where the hit ends
        :param identity: value of the comparison of both sequences
        :return: Returns the updated dictionary and a boolean that indicates whether or not there is an overlap
        """
        # Fills the list in the positions where a hits happens with the score value
        for i in range(query_start, query_end + 1):
            if len(self.dictionary[query + '-' + hit]) >= query_end:
                if self.dictionary[query + '-' + hit][i - 1] == 0:
                    zero_list[i - 1] = identity
                elif self.dictionary[query + '-' + hit][i - 1] < identity:
                    zero_list[i - 1] = identity
                else:
                    zero_list[i - 1] = self.dictionary[query + '-' + hit][i - 1]

            # Bootstrap of the coverage vector
            self.dictionary[query + '-' + hit] = zero_list


