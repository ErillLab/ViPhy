import sys
import csv
import os
import json


def read_csv(file, column_number):
    """
    Reads a csv file and creates an accession list using its content
    :param file: name of the file that contains the information to create the accession list
    :param column_number: specific column that will be read
    :return: Returns a list of lists with all the accessions
    """
    l = []

    # Read the file content
    with open(file, "r", encoding='utf-8', errors='ignore') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count != 0:
                l.append(row[column_number])
            line_count += 1

        line_count = 1
        length = len(l)

        # Writes the accessions in the required format
        accessions = ('[')
        for element in l:
            accessions += '["' + str(element) + '"]'
            if line_count < length:
                accessions += ', '
            line_count += 1
        accessions += (']')

    return accessions


def result_file(accession_list):
    """
    Writes the accession list obtained before into a new file_name
    :param accession_list: List of lists with the sequence identifiers that the user wants to get
    """
    with open("../accessions_list.txt", 'w') as file:
        file.write(accession_list)


if __name__ == "__main__":
    configuration_file = "csv_configuration_file.json"
    arr = os.listdir('.')

    try:
        json_file = json.load(open(configuration_file)) # Reads json file
    except IOError:
        sys.exit('Could not open csv_configuration_file.json')

    # Reads the content of the configuration file
    file_name = json_file["file_name"]
    column_number = json_file["column_number"]

    # Checks if the file indicates has the correct format
    if os.path.splitext(file_name)[1] == '.csv':
        accessions = read_csv(file_name, column_number)

    result_file(accessions)
