# generate merged reference matrix from series_matrix files obtained from GEO
# containing average beta values
# Jan 2018
# Christa Caggiano

from __future__ import print_function
import pandas as pd
import functools
import os
import argparse
import numpy as np


def parse_input(input_file):
    """
    from given file, take input files defined by user and separates them  into lists of files with and without metadata
    also create a dictionary of user-defined columns to be kept for analysis
    :param input_file: file containing paths of average beta value containing files
    :return: dict of files to remove metadata from, dict of files without metadata-- note user must define this
    """

    to_clean = {}
    cleaned = {}

    try:
        if os.path.getsize(input_file) > 0:  # checks that user-input file is not empty
            with open(input_file, "r") as f:
                for line in f:
                    split_line = line.strip("\n").split("\t")  # strip all weird newline chars and read tab-sep values

                    # if the file to be parsed needs metadata removed, add to list to be cleaned
                    if split_line[1].strip(" ").upper() == "YES":
                        to_clean[split_line[0]] = return_geo_names(split_line)  # obtain the columns to be kept

                    # otherwise, return a dictionary of already 'clean' files and relevant columns
                    else:
                        cleaned[split_line[0]] = return_geo_names(split_line)

            return to_clean, cleaned

    # checks if the file the user defines as input exists, otherwise throw error and quit
    except Exception:

        print("File input not correct. Please check input.")

        raise


def return_geo_names(split_line):
    """
    get a list of geo accessions / column values for each input file
    :param split_line: line from user file containing file name and geo ids
    :return: geo accession numbers in a list
    """
    geo = []
    if len(split_line) > 2: # checks if user gave user ID numbers
        for i in range(2, len(split_line)):  # for all data entered by the user that is not the file name or clean_matrix check
                geo.append(split_line[i])

    # if the user defines no geo accession numbers to subset on, return an empty list
    return geo


def check_file_exists(input, should_exist):
    """
    check if a given file exists in the specified path
    :param input: file path
    :param should_exist: whether the program expects the file to exist
    :return:
    """
    if should_exist:  # if it is expected, check all the inputs provided
        for file in input:
            if not os.path.exists(file) or os.path.getsize(file) == 0:  # throw an error if the file doesn't exist or is size 0
                raise FileNotFoundError
    else:  # if the file is not expected to exist, remove it if it does exist
        try:
            os.remove(input)
        except OSError:
            pass


def generate_clean_file(unclean_files):
    """
    facilitates the "cleaning" files defined by the user as needing metadata removed
    :param unclean_files: dictionary of files to be cleaned
    :return: dictionary of files with the cleaned path replaced as key
    """

    new_files = {}
    if unclean_files:  # checks to make sure dict isn't empty- if user did not define any unclean files
        check_file_exists(unclean_files.keys(), True)  # checks that files exists, given that we expect them to exist
        for file in unclean_files:  # if the files exist, iterate over and remove metadata
            clean_file_name = file.replace(".txt", "_cleaned.txt")
            remove_metadata(clean_file_name, file)
            new_files[clean_file_name] = unclean_files.get(file)  # make a copy of new dict where the clean file path is replaced as the key

        return new_files  # return this new dictionary
    else:
        return {}


def remove_metadata(clean_file, file):
    """
    opens the 'unclean' files given by generate_clean_file() and creates a file without metadata
    :param clean_file: name of new file
    :param file:  file to be cleaned
    :return:
    """

    check_file_exists(clean_file, False)  # if a cleaned file exists, delete it and start fresh

    try:
        # open file and remove the metadata. Print to a new file.
        with open(clean_file, "w") as out:
            with open(file) as f:
                start_printing = False
                for line in f:
                    if line.startswith("!Sample_geo_accession"):  # print the sample name-- standard series matrix file struct
                        print(line.rstrip(), file=out)
                    elif start_printing:  # if already begun printing, continue
                        print(line.rstrip(), file=out)
                    elif line.startswith("\"ID_REF\"") or line.startswith("ID_REF"):  # start printing where beta-values begin
                        start_printing = True
            # os.system("gzip " + file)  # zip up file with metadata because they tend to be large

    except Exception:
        raise
        # if the file is not in the proper format, throw an error
        print("ERROR: Please input a standard soft_matrix txt/tsv file.")
        print("File must include Sample_geo_accession and ID_REF lines")


def merge_file_dicts(cleaned, newly_cleaned):
    """
    takes dictionary of files not cleaned by program (assumed to be already free of metadata and merges them
    with files where metadata is removed by program
    :param cleaned: dict of files free of metadata
    :param newly_cleaned: dict of files were metadata was removed
    :return: merged dictionary
    """
    z = cleaned.copy()  # start with x's keys and values- need to make a copy because of weird python memory stuff
    z.update(newly_cleaned)  # modifies z with y's keys and values & returns None
    return z  # return reference to completed dict


def generate_merge(file_list):
    """
    controls merging of all beta value files in order to create create reference table
    :param file_list: all files to be merged
    :return: the merged data frame !
    """

    df_list = ["df_" + str(i) for i in range(len(file_list))]  # generates unique data frame names for merging on

    for i, file in zip(range(len(df_list)), file_list):  # for each file in our list load the data frame
        df_list[i] = load_df(df_list[i], file, file_list[file])  # must be df_list[i] because of annoying memory stuff - we are physically changing what datatype is at that list index

    return merge_df(df_list)  # return the merged data frame


def load_df(df, file, geo_entries):
    """
    loads the unique data frame, subsets to only columns entered by user
    :param df: data frame name
    :param file: file to be loaded into data frame
    :param geo_entries: columns to keep
    :return:
    """
    try:
        df = round_df(pd.read_table(file, delimiter="\t")) # load the tab delimited file, and round it to 4 decimal places
        return subset_df(df, geo_entries)
    except Exception:
        print("ERROR: unable to load data frame.")  # throw an error and quit if a problem is encountered
        raise


def subset_df(df, geo_entries):
    """
    subset to only columns defined by user
    :param df: data frame
    :param geo_entries: geo accession numbers/ specific columns to be restricted to
    :return: subseted dataframe
    """
    if geo_entries:
        for item in geo_entries:
            if item not in list(df):  # makes sure the column desired to be kept is actually in the data frame
                raise ValueError
        geo_entries.append("!Sample_geo_accession")
        return df[geo_entries]  # otherwise, subset
    else:
        return df  # if the user did not define any columns, do not subset and return the entire data frame


def round_df(df):
    """
    round the data frame to 4 decimals
    :param df: data frame to be rounded
    :return: rounded data frame
    """
    return df.round(decimals=4)


def merge_df(df_list):
    """
    merge all the data frames on the cpg Probe ID (sample geo accession in file)
    :param df_list: list of data frames to merge
    :return: merged final data frame
    """

    # merge the list of data frames with a left outer join on the CPG probe id, given by sample geo accession column
    df_final = functools.reduce(lambda left, right: pd.merge(left, right, how='outer', on='!Sample_geo_accession'), df_list)
    # df_final = df_final.fillna(value="NA") # fi
    return df_final


def annotate_reference(df, illumina_file):
    """
    take CpG probes and use the illumina reference to generate information about their chromosome location/position
    :param df: final merged data frame
    :param illumina_file: illumina reference file
    :return:
    """
    cpg_dict = generate_cpg_dict(illumina_file) # calls function to make dictionary of probe:chrom info
    rows_to_remove = []

    chr = []
    start = []

    for index, row in df.iterrows():  # for each row in our merged file
        if row["!Sample_geo_accession"] in cpg_dict:  # check if it is a valid probe name
            chr.append(cpg_dict[row["!Sample_geo_accession"]][0])  # append chromosome
            start.append(cpg_dict[row["!Sample_geo_accession"]][1])  # append position
        else:  # otherwise, append Na
            chr.append("NA")
            start.append("NA")
            rows_to_remove.append(index)

    df["chr"] = pd.Series(chr).values  # add this information to the data frame
    df["start"] = pd.Series(start).values
    return df.drop(df.index[rows_to_remove])


def generate_cpg_dict(illumina_file):
    """
    generate dictionary of probes to chromosome location
    :param illumina_file:
    :return: dictionary
    """
    cpg_dict = {}

    with open(illumina_file) as f:
        for line in f:
            split_line = line.split(",")
            cpg, chr, start = split_line[:3]  # known trait of file - chrom position and start are 1st 3 entries
            cpg_dict[cpg] = (chr, start)  # add info to dict

    return cpg_dict


if __name__ == "__main__":

    # parser = argparse.ArgumentParser()
    # parser.add_argument("input_file", type=str, help="txt file containing a list of soft matrix files to be merged")
    # parser.add_argument("--non_standard", type=str, help="list of files not in soft matrix format")
    # args = parser.parse_args()

    to_clean, cleaned = parse_input("file_list.txt")  # parse input given info provided by user
    all_clean_files = merge_file_dicts(cleaned, generate_clean_file(to_clean)) # generate clean files and merge into one coherent list
    merged = generate_merge(all_clean_files)  # make our merged data frame

    df_final = annotate_reference(merged, "HumanMethylationSites.txt")  # annotate with illumina reference
    df_final.to_csv('reference_matrix.txt', sep="\t")  # save as a tab-sep file 

