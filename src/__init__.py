
from src.user_input import UserInput
from src.entez_utils import access_ncbi
from src.preprocess import Preprocess
from src.blast import Blast
from src.phylo import Phylogeny

import os
import warnings

#warnings.simplefilter("ignore")


def delete_folder_content(folder_name):
    """
    Deletes all files from the indicated folder
    :param folder_name: name of the folder from which the files will be deleted
    """
    for file_to_delete in os.listdir("../" + folder_name):
        os.remove("../" + folder_name + '/' + file_to_delete)


def display_error_messages(error_list):
    """
    Shows the error messages collected during the preprocessing phase
    :param error_list: list of errors that have arisen during the sequence preprocessing
    """
    print('\n', "----------------ERROR LIST------------------")
    if len(error_list) > 0:
        print("The following files have been ignored:")
        for f in error_list:
            print("-", f[0], ": ", f[1])  # File name + error
    else:
        print("No errors were found")


def check_settings(sequence_type, protein_type):
    """
    # !!! docstring here
    """
    # Check parameters correctness
    if sequence_type.lower() not in ['dna', 'protein']:
        raise ValueError("Invalid sequence_type parameter. Please chose " +
                         "between 'DNA' and 'protein'.")
    if protein_type not in ['annotated', 'translated', None]:
        raise ValueError("Invalid protein_type parameter. Please chose " +
                         "between 'annotated', 'translated', or null.")
    
    # Check parameters consistency
    if (sequence_type.lower() == 'dna' and
        protein_type != None):
        warnings.warn("protein_type was set to " + protein_type +
                      ", but it will be ignored, because sequence_type was " +
                      "set to DNA.")
    if (sequence_type.lower() == 'protein' and
        protein_type == None):
        raise ValueError("sequence_type was set to 'protein', but " +
                         "protein_type was set to null. Please chose " +
                         "'annotated' or 'translated' as protein_type.")


def go():
    """
    The entry-point for the pipeline.
    """
    u_input = UserInput()

    # Locates important folders
    input_folder = u_input.get_input_folder()
    working_folder = u_input.get_working_folder()
    output_folder = u_input.get_output_folder()

    # Remaining information of the configuration file
    sequence_type = u_input.get_sequence_type()
    protein_type = u_input.get_protein_type()
    check_settings(sequence_type, protein_type)
    accession_ncbi_list = u_input.get_genome_accessions()
    user_email = u_input.get_user_email()
    distance_function = u_input.get_distance_function()
    e_value = u_input.get_e_value()
    cutoff = u_input.get_cutoff()
    replicates = u_input.get_replicates()
    blast_word_size = u_input.get_blast_word_size()

    # Output files configuration
    majority_or_support_tree = u_input.get_phylogenetic_tree_type()
    original_newick_tree = u_input.get_original_newick_tree()
    original_distance_matrix = u_input.get_original_distance_matrix()
    bootstrap_distance_matrix = u_input.get_bootstrap_distance_matrix()

    # Deletes old content from files
    delete_folder_content(working_folder)
    # delete_folder_content(output_folder)

    # Downloads NCBI files
    access_ncbi(accession_ncbi_list, user_email, input_folder)

    # Preprocessing phase
    n_files = 0
    error_list = []
    preprocess_phase = Preprocess()
    for file in os.listdir("../" + input_folder):  # Navigates into the input_folder
        n_files += 1
        error_list = preprocess_phase.preprocessing_phase(file, input_folder, sequence_type, protein_type, working_folder)

    # Displays a list of error detected in the preprocessing code
    display_error_messages(error_list)

    if len(error_list) < n_files - 1:
        alignment = Blast()
        # Builds a database
        distance_dictionary, coverage_vector_dictionary = alignment.make_blast_database(
            sequence_type, working_folder, e_value, blast_word_size)
        print("Sequence alignment has been done")

        # Calculates distances and generates a phylogenetic tree in newick format
        phylogeny_tree = Phylogeny()
        print("Creating phylogenetic trees")
        newick_tree = phylogeny_tree.get_newick_tree(coverage_vector_dictionary, distance_dictionary, distance_function,
                                                     replicates, working_folder, output_folder,
                                                     original_distance_matrix, bootstrap_distance_matrix,
                                                     original_newick_tree)

        # Read and concatenates trees from files
        tree_list = phylogeny_tree.get_tree_list(working_folder)

        # Generates a consensus trees with or without support
        if majority_or_support_tree in ["Support", "support"]:
            phylogeny_tree.get_support_tree(newick_tree, tree_list, output_folder)
        elif majority_or_support_tree in ["Majority", "majority"]:
            phylogeny_tree.majority_consensus_tree(output_folder, tree_list, cutoff)
        else:
            if majority_or_support_tree in ["Both", "both"]:
                phylogeny_tree.get_support_tree(newick_tree, tree_list, output_folder)
                phylogeny_tree.majority_consensus_tree(output_folder, tree_list, cutoff)
            else:
                print("No majority tree consensus or support tree will be calculated")
    else:
        print('\n', "At least two correct sequences to compare are needed. Please, check the error list to solve the "
                    "detected problems and the content of the '" + input_folder + "' folder.")


if __name__ == '__main__':
    go()
