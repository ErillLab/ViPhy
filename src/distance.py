from random import randrange
from Bio import Phylo
from Bio.Phylo import TreeConstruction
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

from src.files import export_matrix

import numpy as np


def get_inverted_key(key):
    """
    Gets the inverted key from an initial key
    :param key: key that allows us to access to an specific coverage vector
    :return: Returns the inverted key
    """
    separated_key = key.split('-')
    inverted_key = separated_key[1] + '-' + separated_key[0]
    return inverted_key


def vector_length(coverage_vector):
    """
    Calculates the length of an array or list
    :param coverage_vector: Vector that includes all the HSPs between two sequences
    :return: Returns the total length of the input vector
    """
    length = len(coverage_vector)
    return length


def vector_no_zeros(coverage_vector):
    """
    Calculate the length of the HSPs
    :param coverage_vector: Vector that includes all the HSPs between two sequences
    :return: Returns the number of positions that are not equal to zero
    """
    dif_zero = 0
    for base in coverage_vector:
        if base != 0:
            dif_zero += 1

    return dif_zero


def identities(coverage_vector):
    """
    Calculates the amount of characters which match exactly between two different sequences
    :param coverage_vector: Vector that includes all the HSPs between two sequences
    :return: Returns the number of identities or exact matches between two sequences
    """
    sum_identities = 0
    for value in coverage_vector:
        sum_identities += value
    return sum_identities


def d0(coverage_vector_dictionary, distance_dictionary):
    """
    Calculates the distance between two sequences using the HSPs lenght between them and their total length
    :param coverage_vector_dictionary: Dictionary that contains the coverage vectors
    :param distance_dictionary: Dictionary that will contain the calculated distances
    :return: Returns the current distance dictionary
    """
    for key in coverage_vector_dictionary.keys():
        # Sequences
        coverage_vector = coverage_vector_dictionary[key]

        # Gets inverted key
        inverted_key = get_inverted_key(key)
        inverted_coverage_vector = coverage_vector_dictionary[inverted_key]

        # Hits length
        hit_length = vector_no_zeros(coverage_vector)
        hit_length += vector_no_zeros(inverted_coverage_vector)

        # Total length
        total_length = vector_length(coverage_vector)
        total_length += vector_length(inverted_coverage_vector)

        # Distance formula
        result = 1 - (hit_length / total_length)
        # result = "{:.8f}".format(float(result))
        distance_dictionary[key] = result

    return distance_dictionary


def d4(coverage_vector_dictionary, distance_dictionary):
    """
    Calculates the distance between two sequences using the number of identical bases between them and the HSPs length
    :param coverage_vector_dictionary: Dictionary that contains the coverage vectors
    :param distance_dictionary: Dictionary that will contain the calculated distances
    :return: Returns the current distance dictionary
    """
    for key in coverage_vector_dictionary.keys():
        # Sequences
        coverage_vector = coverage_vector_dictionary[key]

        # Gets inverted key
        inverted_key = get_inverted_key(key)
        inverted_coverage_vector = coverage_vector_dictionary[inverted_key]

        # Identities
        identities1 = identities(coverage_vector)
        identities2 = identities(inverted_coverage_vector)

        # Values of the coverage vector that are not zero
        dif_zero1 = vector_no_zeros(coverage_vector)
        dif_zero2 = vector_no_zeros(inverted_coverage_vector)

        if dif_zero1 != 0 and dif_zero2 != 0:
            # Calculates the average identity between the coverage vector and its opposite
            sum_identities = identities1 / dif_zero1
            inverted_sum_identities = identities2 / dif_zero2
            total_identities = (sum_identities + inverted_sum_identities) * 0.5

            # Hits length or sum of vales that are not zero on the coverage vector
            hit_length = dif_zero1 + dif_zero2

            # Distance formula
            result = 1 - ((2 * total_identities) / hit_length)
            if result < 0:
                result = 0.0
            # result = "{:.8f}".format(result)
            distance_dictionary[key] = result
        else:  # Both sequences are completely different
            distance_dictionary[key] = 1.0

    return distance_dictionary


def d6(coverage_vector_dictionary, distance_dictionary):
    """
    Calculates the distance between two sequences using the number of identical bases between them and the total
    length of the two sequences
    :param coverage_vector_dictionary: Dictionary that contains the coverage vectors
    :param distance_dictionary: Dictionary that will contain the calculated distances
    :return: Returns the current distance dictionary
    """
    for key in coverage_vector_dictionary.keys():
        # Sequences
        coverage_vector = coverage_vector_dictionary[key]

        # Gets inverted key
        inverted_key = get_inverted_key(key)
        inverted_coverage_vector = coverage_vector_dictionary[inverted_key]

        # Identities
        identities1 = identities(coverage_vector)
        identities2 = identities(inverted_coverage_vector)

        # Values of the coverage vector that are not zero
        dif_zero1 = vector_no_zeros(coverage_vector)
        dif_zero2 = vector_no_zeros(inverted_coverage_vector)

        if dif_zero1 != 0 and dif_zero2 != 0:
            # Calculates the average identity between the coverage vector and its opposite
            sum_identities = identities1 / dif_zero1
            inverted_sum_identities = identities2 / dif_zero2
            total_identities = (sum_identities + inverted_sum_identities) * 0.5

            # Calculate the length of the complete vector
            total_length = vector_length(coverage_vector)
            total_length += vector_length(inverted_coverage_vector)

            # Distance formula
            result = 1 - ((2 * total_identities) / total_length)
            result = "{:.8f}".format(result)
            distance_dictionary[key] = float(result)  # Limits the distance value to 5 decimals
        else:  # Both sequences are completely different
            distance_dictionary[key] = 1.0

    return distance_dictionary


def distance_matrix(dictionary, replicates, working_folder, output_folder, original_distance_matrix,
                    bootstrap_distance_matrix, original_newick_tree):
    """
    Creates a matrix with the distances calculated and stored in a dictionary
    :param dictionary: dictionary where the distance values are saved
    :param replicates: current number of the bootstrap replicates
    :param working_folder: folder from where the data is obtained
    :param output_folder: folder where the distance matrix will be stored
    :param original_distance_matrix: boolean that indicates if the user wants to generate a distance matrix with
    the original data
    :param bootstrap_distance_matrix: boolean that indicates if the user wants to generate a distance matrix with
    the data after creating bootstrap replicates
    :param original_newick_tree: boolean used to indicate if the user wants to generate a tree with original data
    :return: Returns the final distance matrix created with all the data in the input dictionary
    """
    key_list = []
    # Creates a list with all the possible keys
    for key in dictionary.keys():
        key_parts = key.split('-')
        if key_parts[0] not in key_list:
            key_list.append(key_parts[0])

    matrix = []
    # Create the matrix distance with distances calculated before
    for first_key in key_list:
        list = []
        for second_key in key_list:
            list.append(dictionary[first_key + '-' + second_key])
        matrix.append(list)
    tree = lower_triangle_matrix(matrix, key_list, replicates, working_folder, output_folder,
                                 original_newick_tree)

    if replicates == 0:
        if original_distance_matrix in ["True", "true"]:
            file_name = 'original_distance_matrix.tsv'
            export_matrix(key_list, matrix, file_name, output_folder)
    else:
        if bootstrap_distance_matrix in ["True", "true"]:
            file_name = 'all_distance_matrix.tsv'
            export_matrix(key_list, matrix, file_name, output_folder)
    return tree


def lower_triangle_matrix(d_matrix, key_list, replicates, working_folder, output_folder, original_newick_tree):
    """
    Creates a lower triangular matrix using a pre-calculated distance matrix and a list with all the sequences
    identifiers that are being compared. It also creates a tree from the created matrix.
    :param d_matrix: matrix that contains all distances between sequences pairs
    :param key_list: list where all the sequences identifiers are stored
    :param replicates: number that indicates the current bootstrap sample
    :param working_folder: folder from where the data is obtained
    :param output_folder: folder where the distance matrix and the generated trees will be stored
    :param original_newick_tree: boolean used to indicate if the user wants to generate a tree with original data
    :return: Returns a phylogenetic tree generated from a distance matrix
    """
    matrix = []
    list = []
    count = 1
    prev = 0

    d_matrix = np.array(d_matrix)
    num_sequences = len(key_list)  # Number of sequences we are comparing

    # Takes the values of the lower triangular distance matrix and saves in a list
    aux_list = []
    position = 0
    for l in d_matrix:
        internal_counter = 0
        for element in l:
            if internal_counter <= position:
                aux_list.append(element)
            internal_counter += 1
        position += 1
    array = np.array(aux_list)

    # Converts the list of distances calculated above into a matrix
    for i in range(0, num_sequences):
        for j in range(prev, count + prev):
            list.append(array[j])
        matrix.append(list)
        list = []
        prev = count + prev
        count += 1

    # Creates a lower triangular matrix object
    dm = TreeConstruction.DistanceMatrix(key_list, matrix)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # Saves the original tree into a file with newick format
    if replicates == 0:
        if original_newick_tree in ["True", " vtrue"]:
            for node in tree.get_nonterminals():
                node.name = None
            Phylo.write(tree, "../" + output_folder + "/original_tree.nwk", "newick")
    else:
        Phylo.write(tree, "../" + working_folder + "/tree_" + str(replicates) + ".nwk", "newick")

    return tree


def bootstrap(boostrap_sample_dict):
    """
    Samples a dataset with a specific number of replacements
    :param boostrap_sample_dict: dictionary that contains the coverage vectors previously calculated
    :return Returns a dictionary with the new samples created
    """
    for key in boostrap_sample_dict.keys():
        length = len(boostrap_sample_dict[key])
        aux_list = []
        for i in range(0, length):
            position = randrange(length)
            aux_list.append(boostrap_sample_dict[key][position])
        boostrap_sample_dict[key] = aux_list
    return boostrap_sample_dict
