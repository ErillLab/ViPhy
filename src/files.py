
from Bio import SeqIO
from Bio.Seq import UnknownSeq


def read_fasta_file(input_folder, fasta_file):
    """
    Reads a fasta file to obtain a nucleotide or protein sequence.
    :param input_folder: folder from where we obtain the file to read
    :param fasta_file: file that contains a nucleotide or protein sequence, should be FASTA format
    :return: Returns a list with the identifier and another list with the nucleotide sequence
    """
    record_seq = []
    record_id = []
    print('Reading ' + fasta_file)

    # This loop allows us to concatenate multiple sequences
    for record in SeqIO.parse("../" + input_folder + '/' + fasta_file, "fasta"):
        record_id.append(record.id)
        seq = record.seq.upper()  # Upper case the sequence
        record_seq.append(seq)

    return record_id, record_seq


def read_gb_file_as_nucleotide(input_folder, gb_file, error_list):
    """
    Reads and extracts a nucleotide sequence from a Genbank file
    :param input_folder: folder from where we obtain the genbank file to read
    :param gb_file: file that can store several sequences and extra information
    :param error_list: list of errors that have arisen during the sequence preprocessing
    :return: Returns a list with the identifier, another list with the nucleotide sequence and the list of errors
    """
    record_id = []
    record_seq = []
    print('Reading ' + gb_file)

    with open("../" + input_folder + '/' + gb_file, 'r') as gb_file:
        records = SeqIO.read(gb_file, "genbank")
        description = records.description

        # Modifies the file name to respect a standard format
        description = description.replace("-", "_").replace(" ", "_").replace("(", "_").replace(")", "_")
        description = description.replace("[", "_").replace("]", "_")
        description_division = description.split("/")
        description_name = description_division[0].split(",")

        sequence = records.seq
        if isinstance(sequence, UnknownSeq):  # Checks if there is a nucleotide sequence in the file
            error_list.append([gb_file, "There seems to be no sequence in this GenBank file"])
        else:
            record_id.append(records.id + '_' + str(description_name[0]))
            record_seq.append(records.seq)

    return record_id, record_seq, error_list


def read_gb_file_as_protein(input_folder, gb_file, error_list):
    """
    Reads and extracts protein sequences from a Genbank file
    :param input_folder: folder from where we obtain the genbank file to read
    :param gb_file: file that can store several sequences and extra information, including protein sequences
    :param error_list: list of errors that have arisen during the sequence preprocessing
    :return: Returns a list of identifiers and another of protein sequences obtained from the genbank file
    """
    id_list = []
    protein_list = []
    print('Reading ' + gb_file)

    gb_name = ".".join(gb_file.split('.')[:-1])
    with open("../" + input_folder + '/' + gb_file, 'r') as f:
        gb_cds = SeqIO.InsdcIO.GenBankCdsFeatureIterator(f)
        recorder = SeqIO.read("../" + input_folder + '/' + gb_file, "genbank")
        description = recorder.description

        # Modifies the file name to respect a standard format
        description_seq = description.replace("-", "_").replace(" ", "_").replace("(", "_").replace(")", "_")
        description_seq = description_seq.replace("[", "_").replace("]", "_")
        description_division = description_seq.split("/")[0]
        description_name = description_division.split(",")[0]

        for cds in gb_cds:
            if cds.seq is not None:  # Checks if there is a CDS
                if gb_name not in id_list:
                    id_list.append(gb_name + '_' + description_name)  # Final file name
                else:
                    id_list.append(str(cds.id) + '_' + description_name)
                protein_list.append(cds.seq)
            elif isinstance(cds.seq, UnknownSeq):
                error_list.append([gb_file, "There seems to be no sequence in this GenBank file"])

    return id_list, protein_list, error_list


def export_matrix(key_list, d_matrix, file_name, output_folder):
    """
    Writes in a file the values of a distance matrix following the TSV format
    :param key_list: list where all the sequences identifiers are stored
    :param d_matrix: matrix that contains all distances between sequences pairs
    :param file_name: name of the file that will contain the distance matrix in phylip format
    :param output_folder: folder where the generated file will be stored
    """
    line = ""
    count = 0
    length = len(key_list)

    # Opens a file and writes the content of the distance matrix
    with open("../" + output_folder + '/' + file_name, 'a') as f:
        # Writes the number of analyzed sequences
        f.write('{}\n'.format(str(length)))
        for i in key_list:
            for j in range(length - 1):
                line += str(d_matrix[count][j])
                line += "\t"
            line += str(d_matrix[count][length - 1])  # last element of the row
            line += '\n'
            f.write('{}\t{}'.format(i, line))  # Writes the sequences and the distances of each comparison
            count += 1
            line = ""
        f.write('\n')

