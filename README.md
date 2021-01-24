# ViPhy
Virus phylogeny system.


## Preprocessing process
Translates a nucleotide sequences(DNA) into a sequence of proteins.

The steps in the script are the following:

1. Opens and reads ``settings.json``

2. Downloads the genbank genomes especified in the Json file from NCBI 

3. If the analysis_type indicated in the Json file is 'nucleotide':

3.1. Navigate into input files directory until there are no more files to read.

3.1.1. When it finds a fasta file, it checks if the file contains a nucleotide sequence. If not, then returns to the step 3.1.. 

3.1.2. When it finds a genbank file, it extracts the nucleotide sequence from the file. 

3.1.3. Translates the nucleotide sequence into a protein sequence.

3.1.4. Stores the combined protein sequence in a new fasta file


4. If the analysis_type indicated in the Json file is 'protein':

4.1. Navigates into the input files directory until there are no more files to read.

4.1.1. When it finds a fasta file, it checks if the file contains a nucleotide sequence. If so, stores the sequence into a new fasta file. 

4.1.2. When it finds a genbank file, it extracts the protein sequence from the file. 

4.1.2.1. Stores the combined protein sequence in a new fasta file


## Sequence alignment
During this phase, a database formed by the combination of multiple fasta files will be created and it will be stored in the folder named ``dbFolder``.

Then, the sequences will be compared with each other in order to find regions of similarity between an input sequence and the database, in a matter of seconds. More specifically, the program will use Blastp, a version of Blast that compares two protein sequences.

Finally, the distance between them will be also calculated and saved in a distance matrix. This matrix will be saved in the ``Outputs`` folder following the phylip format.


## Phylogenetic trees
Once the distance matrix has been obtained, the code will generate a phylogenetic tree using the newick format and it will be stored in the ``Outputs`` folder.

After this, if the user wishes it, bootstrap will be used to create sub-samples from which the parameters of the model are estimated repeated times.

At the end, a consensus tree will be created with all the samples obtained and stored in the `Outputs` folder.


## Viphy_env.yml

Conda environment created using Phyton38 that contains a specific collection of conda packages installed to ensure the smooth running of the program.

To create the environment from `viphy-env.yml` file:

	$ conda env create -f viphy-env.yml

To activate the `viphy-env.yml` environment:

	$ conda activate viphy-env

To deactivate the conda environment:

	$ conda deactivate


### Dependencies

- Python 3.8.3
- Biopython 1.78


## Settings.json

This file contains important information for the correct operation of the code. Please, don't change this file path!

Structure of the JSON file:


{

	"user_email": "A.N.Other@example.com",

	"input_folder" : "genome_data",

	"genome_accessions": [["GQ919031.1", "JX182370.1"], ["NC_015464"], ["NC_042011.1"]],

	"output_folder" : "results",

	"working_folder" : "tmp",

	"analysis_type": "nucleotide",

	"e_value": 0.0001

	"distance_function": "d6",

	"replicates": 100,
	
	"cutoff": 0.7,

	"majority_or_support_tree": "both",

	"get_original_newick_tree": "True",

	"get_original_distance_matrix": "True",

	"get_bootstrap_distance_matrix": "True"
}



### Details of the JSON file: 

- "user_email": Email that will be used to identify the user in NCBI
- "analysis_type": Can be "nucleotide" or "protein"
- "input_folder": Folder where you can find the fasta or genbank files that the program will read
- "working_folder": Path for the directory in which all the other files are located
- "output_folder": Folder in which the outgoing files will be located.
- "genome_accessions": List of lists that contains the identifier of the files you want to download from Genbank database. 
- "e_value": Number of expected hits of similar score that could be found just by chance 
- "distance_function": Type of function that will be used to calculate the distance between two sequences. It should be "d0", "d4" or "d6"
- "replicates": Number of bootstrap samples 
- "majority_or_support_tree": Way to calculate a consensus tree. It can be "majority", "support" or "both"
- "cutoff": Threshold used when "majority" consensus tree is selected. It should be a number between 0 and 1.
- "get_original_newick_tree": Indicates using a boolean if you want to get the original newick tree
- "get_original_distance_matrix": Indicates using a boolean if you want to get the original distance matrix in phylip format
- "get_bootstrap_distance_matrix": Indicates using a boolean if you want to get each one of the possible distance matrix for each bootstrap sample



## Other files

- Fasta files: Inside `Inputs` folder. Documents that contain a nucleotide or amino acid sequence that the program will read. If the content is a nucleotide sequence, it will be translates into proteins.

- Genbank files: Inside `Inputs` folder. Documents where you can obtain a nucleotide or amino acid sequence, together with other relevant information. 

- DataBase files: Inside `dbFolder` folder. Database files that will appear in this folder once you have created a database.  

- DataBase.fasta: File that will formed by the combination of all the files stored in the `WorkingFolder` folder after the preprocessing process. 

- `Original_tree.nwk`: Inside `Outputs`. File that contains the original phylogenetic tree, in other words, the phylogenetic tree before using bootstrap.

- `Support_consensus_tree.nwk`: Inside `Outputs`. File that contains the final consensus phylogenetic tree in newick format with branches support.

- `Majority_consensus_tree.nwk`: Inside `Outputs`. File that contains the final consensus phylogenetic tree in newick format following the majority rule.

- `Original distance matrix.txt`: Inside `Outputs` folder. This file contains the distance matrix calculated before using bootstrap using the phylip format.

- `Bootstrap_distance_matrix.txt`: Inside `Ouputs` folder. File which contains a distance matrix for each sub-sample created using bootstrap. It follows the phylip format. All will be stored in the same file. 


## Folders

- ``Inputs``: Folder to store the set of files that the program will read.

- ``WorkingFolder``: Folder to store the resulting sequences after the translation ends. 

- ``Outputs``: Folder that stores the final results.

- ``dbFolder``: Folder where the database will be stored after its creation.

