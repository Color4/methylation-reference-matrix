# generate merged reference matrix from soft_matrix files obtained from GEO
# containing average beta values
# Jan 2018
# Christa Caggiano

from __future__ import print_function
import pandas as pd
import functools
import os
import argparse
import numpy as np


def parse_file_list(file_list):

    files = {}
    to_clean = {}
    cleaned = {}

    with open(file_list) as f:
        for line in f:
            split_line = line.split("\t")
            if split_line[1].upper() == "YES":
                to_clean[split_line[0]] = split_line[1:]


def return_GEO_names(split_line):

    geo = []
    if len(split_line) > 2:
        for i in range(2, len(split_line)):
            geo.append(split_line[i])
    else:
        return None
    return geo



def check_file_exists(file_list):
    for file in file_list:
        if not os.path.exists(file):
            raise FileNotFoundError("ERROR " + file + "not found. Please check your input file list.")

    if os.path.exists("reference_matrix.txt"):
        os.remove("reference_matrix.txt")


def make_clean_file(file_list):

    clean_file_list = []
    for file in file_list:
        with open(file) as f:
            clean_file_name = file.replace(".txt", "_cleaned.txt")
            clean_file_list.append(clean_file_name)
            remove_metadata(clean_file_name, file)

    return clean_file_list


def remove_metadata(clean_file, file):

    os.system("rm -f " + clean_file)  # if a cleaned file exists, delete it and start fresh

    try:
        # open file and remove the metadata. Print to a new file.
        with open(clean_file, "w") as out:
            with open(file) as f:
                start_printing = False
                for line in f:
                    if line.startswith("!Sample_geo_accession"):
                        print(line.rstrip(), file=out)
                    elif start_printing:
                        print(line.rstrip(), file=out)
                    elif line.startswith("\"ID_REF\"") or line.startswith("ID_REF"):
                        start_printing = True
            # os.system("gzip " + file)  # zip up file with metadata because they tend to be large

    except Exception:
        raise
        # if the file is not in the proper format, throw an error
        print("ERROR: Please input a standard soft_matrix txt/tsv file.")
        print("File must include Sample_geo_accession and ID_REF lines")


def intiate_merge(file_list):

    df_list = ["df_" + str(i) for i in range(len(file_list))]

    for i in range(len(df_list)):
        df_list[i] = load_df(df_list[i], file_list[i])

    return merge_df(df_list)


def load_df(df, file):
    try:
        return round_df(pd.read_table(file, delimiter="\t"))

    except Exception:
        print("ERROR: unable to load data frame.")
        raise


def round_df(df):
    return df.round(decimals=4)


def merge_df(df_list):
    df_final = functools.reduce(lambda left, right: pd.merge(left, right, how='outer', on='!Sample_geo_accession'), df_list)
    df_final = df_final.fillna(value=-1)
    # df_final.to_csv('reference_matrix.txt', sep="\t")
    return df_final


def annotate_reference(df, illumina_file):
    cpg_dict = generate_cpg_dict(illumina_file)
    rows_to_remove = []

    chr = []
    start = []

    for index, row in df.iterrows():
        if row["!Sample_geo_accession"] in cpg_dict:
            chr.append(cpg_dict[row["!Sample_geo_accession"]][0])
            start.append(cpg_dict[row["!Sample_geo_accession"]][1])
        else:
            chr.append(-1)
            start.append(-1)
            rows_to_remove.append(index)

    df["chr"] = pd.Series(chr).values
    df["start"] = pd.Series(start).values
    return df.drop(df.index[rows_to_remove])


def generate_cpg_dict(illumina_file):

    cpg_dict = {}

    with open(illumina_file) as f:
        for line in f:
            split_line = line.split(",")
            cpg, chr, start = split_line[:3]
            cpg_dict[cpg] = (chr, start)

    return cpg_dict


if __name__ == "__main__":

    # parser = argparse.ArgumentParser()
    # parser.add_argument("soft_file", type=str, help="txt file containing a list of soft matrix files to be merged")
    # parser.add_argument("--non_standard", type=str, help="list of files not in soft matrix format")
    # args = parser.parse_args()

    soft_file = "file_list.txt"

    soft_file_list = [line.strip() for line in open(soft_file, 'r')]
    check_file_exists(soft_file_list)

    cleaned_file_list = make_clean_file(soft_file_list)
    merged = intiate_merge(cleaned_file_list)

    df_final = annotate_reference(merged, "HumanMethylationSites.txt")
    df_final.to_csv('reference_matrix.txt', sep="\t")

