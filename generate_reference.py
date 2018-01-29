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


def parse_input(input_file):

    to_clean = {}
    cleaned = {}

    try:
        if os.path.getsize(input_file) > 0:
            with open(input_file, "r") as f:
                for line in f:
                    split_line = line.strip("\n").split("\t")
                    if split_line[1].strip(" ").upper() == "YES":
                        to_clean[split_line[0]] = return_geo_names(split_line)
                    else:
                        cleaned[split_line[0]] = return_geo_names(split_line)

            return to_clean, cleaned

    except Exception:
        print("File input not correct. Please check input.")
        raise


def return_geo_names(split_line):

    geo = []
    if len(split_line) > 2:
        for i in range(2, len(split_line)):
            geo.append(split_line[i])
    else:
        return []
    return geo


def check_file_exists(input, should_exist):
    if should_exist:
        for file in input:
            if not os.path.exists(file) or os.path.getsize(file) == 0:
                raise FileNotFoundError
    else:
        try:
            os.remove(input)
        except OSError:
            pass


def generate_clean_file(unclean_files):

    new_files = {}
    if unclean_files:
        check_file_exists(unclean_files.keys(), True)
        for file in unclean_files:
            clean_file_name = file.replace(".txt", "_cleaned.txt")
            remove_metadata(clean_file_name, file)
            new_files[clean_file_name] = unclean_files.get(file)

        return new_files
    else:
        return {}


def remove_metadata(clean_file, file):

    check_file_exists(clean_file, False)  # if a cleaned file exists, delete it and start fresh

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


def merge_file_dicts(cleaned, newly_cleaned):
    z = cleaned.copy()  # start with x's keys and values
    z.update(newly_cleaned)  # modifies z with y's keys and values & returns None
    return z


def generate_merge(file_list):

    df_list = ["df_" + str(i) for i in range(len(file_list))]

    for i, file in zip(range(len(df_list)), file_list):
        df_list[i] = load_df(df_list[i], file, file_list[file])

    return merge_df(df_list)


def load_df(df, file, geo_entries):
    try:
        df = round_df(pd.read_table(file, delimiter="\t"))
        return subset_df(df, geo_entries)
    except Exception:
        print("ERROR: unable to load data frame.")
        raise


def subset_df(df, geo_entries):
    if geo_entries:
        for item in geo_entries:
            if item not in list(df):
                raise ValueError
        geo_entries.append("!Sample_geo_accession")
        return df[geo_entries]
    else:
        return df


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
    # parser.add_argument("input_file", type=str, help="txt file containing a list of soft matrix files to be merged")
    # parser.add_argument("--non_standard", type=str, help="list of files not in soft matrix format")
    # args = parser.parse_args()

    to_clean, cleaned = parse_input("file_list.txt")
    all_clean_files = merge_file_dicts(cleaned, generate_clean_file(to_clean))
    merged = generate_merge(all_clean_files)

    df_final = annotate_reference(merged, "HumanMethylationSites.txt")
    df_final.to_csv('reference_matrix.txt', sep="\t")

