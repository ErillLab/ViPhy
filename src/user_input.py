import json
import sys
import os
import re


class UserInput:
    """
    Class definition for UserInput.
    UserInput class interacts with the settings.json file to access the run parameters provided by the user. It also
    checks the content of the settings.json file and sets default values for not specified or not correct parameters.
    """

    def __init__(self):
        """Constructs an object using the settings.json file."""
        setting_file = "settings.json"  # Configuration file
        self.json_file = ""
        try:
            self.json_file = json.load(open("../" + setting_file))  # Reads JSON file
        except IOError:
            sys.exit('Could not open settings.json')

    def get_input_folder(self):
        """Returns the input folder name"""
        input_folder = self.json_file["input_folder"]
        if not os.path.exists("../" + input_folder):
            os.mkdir("../" + input_folder)  # Creates a new folder
        return input_folder

    def get_working_folder(self):
        """Returns the working folder name"""
        working_folder = self.json_file["working_folder"]
        if not os.path.exists("../" + working_folder):
            os.mkdir("../" + working_folder)  # Creates a new folder
        return working_folder

    def get_output_folder(self):
        """Returns the output folder name"""
        output_folder = self.json_file["output_folder"]
        if not os.path.exists("../" + output_folder):
            os.mkdir("../" + output_folder)  # Creates a new folder
        return output_folder

    def get_sequence_type(self):
        """
        Returns the type of sequence to use for the phylogenetic analysis.
        If set to 'DNA', the whole DNA sequences are used. If set to 'protein',
        only the coding sequences are used.
        """
        sequence_type = self.json_file["sequence_type"]
        return sequence_type

    def get_protein_type(self):
        """
        Returns the protein type.
        It could be 'annotated' or 'translated'.
        """
        protein_type = self.json_file["protein_type"]
        return protein_type
    
    def get_blast_word_size(self):
        """
        Returns the word size parameter to be used when running BLAST.
        """
        blast_word_size = self.json_file["blast_word_size"]
        return blast_word_size

    def get_genome_accessions(self):
        """Returns a list of genomes that will be downloaded from ncbi"""
        access_ncbi_list = self.json_file["genome_accessions"]
        return access_ncbi_list

    def get_user_email(self):
        """
        Returns the user email address
        It also checks the format of the email address
        """
        user_email = self.json_file["user_email"]
        match = re.search(r'[\w.-]+@[\w.-]+.\w+', user_email)  # Checks if the e-mail follows the correct format
        if not match:
            print("Please, introduce a valid email address format")
            return
        else:
            return user_email

    def get_distance_function(self):
        """Returns the function that will be used to calculate the distance between sequences"""
        distance_function = self.json_file["distance_function"]
        return distance_function

    def get_replicates(self):
        """Returns the number of bootstrap replicates that will be generated"""
        replicates = self.json_file["replicates"]
        return replicates

    def get_e_value(self):
        """Returns the number of expected hits of similar quality (score) that could be found just by chance"""
        e_value = self.json_file["e_value"]
        if e_value < 0:
            print('\n', "The entered e_value is negative when it should be positive. It will be convert into a positive"
                        " value")
            e_value = e_value * -1
        return e_value

    def get_cutoff(self):
        """Returns the threshold for the values of phylogenetic branches"""
        cutoff = self.json_file["cutoff"]
        if 0 > cutoff or cutoff > 1:
            print('\n',
                  "The entered cutoff value is not correct. It should be between 0 and 1. Until you change it, cutoff"
                  " value will be 0, in other words, the default value will be used.")
            cutoff = 0
        return cutoff

    def get_phylogenetic_tree_type(self):
        """Returns the type of phylogenetic trees that will be generated"""
        phylogenetic_tree_type = self.json_file["majority_or_support_tree"]
        return phylogenetic_tree_type

    def get_original_newick_tree(self):
        """Returns a boolean depending on whether the user wants to generate a tree with the original data or not"""
        original_tree = self.json_file["get_original_newick_tree"]
        return original_tree

    def get_original_distance_matrix(self):
        """
        Returns a boolean depending on whether the user wants to generate a distance matrix with the original data or
        not.
        """
        original_distance_matrix = self.json_file["get_original_distance_matrix"]
        return original_distance_matrix

    def get_bootstrap_distance_matrix(self):
        """
        Returns a boolean depending on whether the user wants to generate a distance matrix with the data after creating
        bootstrap replicates or not.
        """
        bootstrap_distance_matrix = self.json_file["get_bootstrap_distance_matrix"]
        return bootstrap_distance_matrix
