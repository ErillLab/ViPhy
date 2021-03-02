
from Bio import Entrez

import os
import time


def access_ncbi(accessing_list, user_email, input_folder):
    """
    Goes through the accessing list to find the files that the user whats to download from NCBI and call
    another function to do so.
    :param accessing_list: list of lists that contains the elements we want to download
    :param user_email: email that will be used to identify the user in NCBI
    :param input_folder: folder where files will be stored
    """
    for element_list in accessing_list:
        list_length = len(element_list)
        if list_length == 1:  # One element list
            file_name = element_list[0] + ".gb"
            file_path = "../" + input_folder + '/' + file_name

            # Checks if the file was downloaded before
            if not os.path.exists(file_path):
                print('Downloading ' + file_name)
                download_file(user_email, file_path, element_list)

        if list_length > 1:  # More than one element list
            folder_path = "../" + input_folder + '/' + element_list[0]
            if not os.path.exists(folder_path):  # Checks if there is a folder with the same name
                os.mkdir(folder_path)  # Creates a folder
            for position in range(list_length):
                file_name = element_list[position] + ".gb"
                file_path = "../" + folder_path + '/' + file_name

                if not os.path.exists(file_path):  # Checks if the file was downloaded before
                    print('Downloading ' + file_name)
                    download_file(user_email, file_path, element_list)
                    time.sleep(0.25)  # Indicates the specific delay after each download


def download_file(user_email, file_path, accession_list):
    """
    Downloads files form NCBI in GenBank plain text format and saves it into the input folder
    :param user_email: email that will be used to identify the user in NCBI
    :param file_path: input folder path or subdirectory where the files will be saved
    :param accession_list: list that contains the files to download
    """
    Entrez.email = user_email  # Always tell NCBI who you are
    count = 0
    while 0 <= count < 3:  # Tries to download a file at most 3 times if it does succeed at first time
        try:
            if count >= 0:
                handle = Entrez.efetch(db="nucleotide", id=accession_list[0], rettype="gbwithparts", retmode="text")
                data = handle.read()
                f = open(file_path, "w")
                f.write(data)
                f.close()
                count = -1
        except Exception:
            count += 1
    if count == 3:
        print("Download unsuccessful! Cannot fetch " + accession_list[0])
