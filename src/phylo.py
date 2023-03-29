from src.distance import d0, d4, d6, bootstrap, distance_matrix
from Bio import Phylo
from Bio.Phylo import Consensus
import sys
import os


class Phylogeny:
    """
    Class definition for Phylogeny.
    Phylogeny class defines the the code that will be used to generate phylogenetic trees following the newick format.
    User can choose between a tree made by consensus or a tree with branch support values.
    At the end, this trees will be saved in the Outputs folder.
    It is also connected with the Distance class.
    """

    def __init__(self):
        """Constructs an object for the phase of the phylogeny creation."""
        self.count = 0
        self.tree_list = []
        self.tree_format = "newick"

    def get_newick_tree(self, dictionary, distance_dictionary, distance_function, replicates, working_folder,
                        output_folder, original_distance_matrix, bootstrap_distance_matrix, original_newick_tree):
        """
        Calculates the distance between two sequences with the formula indicated by the user. It is done as many times
        as specified by the number of replicates. Once the distance values have been obtained, it is time to create a
        phylogenetic tree.
        :param dictionary: dictionary that includes the coverage vectors
        :param distance_dictionary: dictionary where the distances will be saved
        :param distance_function: indicates the formula that will be used to get the distance
        :param replicates: number of replicates of the original coverage vector using bootstrap
        :param working_folder: folder from where the data is obtained
        :param output_folder: folder where the files that contain the trees will be stored
        :param original_distance_matrix: boolean that indicates if the user wants to generate a distance matrix with
        the original data
        :param bootstrap_distance_matrix: boolean that indicates if the user wants to generate a distance matrix with
        the data after creating bootstrap replicates
        :param original_newick_tree: boolean used to indicate if the user wants to generate a tree with original data
        :return: Returns a tree following the newick format
        """
        original_dictionary = dictionary
        while self.count <= replicates:
            if self.count != 0:
                dictionary = bootstrap(original_dictionary)  # Bootstrap to all the coverage vectors
            if distance_function == 'd0':
                distance_dictionary = d0(dictionary, distance_dictionary)
            else:
                if distance_function == 'd4':
                    distance_dictionary = d4(dictionary, distance_dictionary)
                elif distance_function == 'd6':
                    distance_dictionary = d6(dictionary, distance_dictionary)
                else:
                    print("Please, introduce a valid distance function. It can be 'd0', 'd4' or 'd6'.")

            newick_tree = distance_matrix(distance_dictionary, self.count, working_folder, output_folder,
                                          original_distance_matrix, bootstrap_distance_matrix, original_newick_tree)

            if self.count == 0:
                original_newick_tree = newick_tree  # Tree before bootstrap

            self.count += 1

        return original_newick_tree

    def get_tree_list(self, working_folder):
        """
        Reads and concatenates in a single list all the bootstrap samples
        :param working_folder: directory from where we get the trees created earlier
        :return: Returns all bootstrap trees in a single list
        """
        for file in os.listdir("../" + working_folder):
            if os.path.splitext("../" + working_folder + '/' + file)[1] == '.nwk':
                # Save all bootstrap trees in a single list
                self.tree_list.append(Phylo.read("../" + working_folder + '/' + file, self.tree_format))
        return self.tree_list

    def majority_consensus_tree(self, output_folder, tree_list, cutoff):
        """
        Search majority rule consensus tree from multiple trees.
        :param output_folder: folder where the output will be stored
        :param tree_list: list of trees created with bootstrap
        :param cutoff: threshold used to compare the branches from different trees. Any clade that has <= cutoff support
        will be dropped.
        """
        majority_tree = Consensus.majority_consensus(tree_list, cutoff)
        Phylo.write(majority_tree, "../" + output_folder + "/majority_consensus_tree.nwk", self.tree_format)
        # Also print to standard output
        # print("Majority consensus tree: ")
        # Phylo.write(majority_tree, sys.stdout, self.tree_format)

    def get_support_tree(self, original_tree, trees, output_folder):
        """
        Calculate branch support for a target tree given bootstrap replicate trees
        :param original_tree: original tree created without using bootstrap samples
        :param trees: list of trees created with bootstrap
        :param output_folder: name of the file where the consensus tree will be stored
        """
        tree_with_support = Consensus.get_support(original_tree, trees)
        Phylo.write(tree_with_support, "../" + output_folder + "/support_consensus_tree.nwk", self.tree_format)
        # Also print to standard output
        # print("Support tree: ")
        # Phylo.write(tree_with_support, sys.stdout, self.tree_format)
