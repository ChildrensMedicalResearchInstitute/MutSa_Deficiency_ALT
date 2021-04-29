import argparse
import json
import sys
import time

import numpy as np
import pandas as pd

from utils import *


def read_gear_vca(input_path):
    name_map = {'q_telomere': "telomere", 'inner_non_telomeric': "non_telomeric", 'p_telomere': "telomere",
                "mapq_fail_": "mapq_fail_", "other_": "other_", "qv_fail_": "qv_fail", "unmapped_": "unmapped"}
    with open(input_path) as f:
        data = json.load(f)
    data_list = list()
    for chromosome, clist in data.items():
        for item in clist:
            name = item["name"]
            for mutation, sample_data in item["mutations"].items():
                mutation_id = mutation.split(":")[0]
                mutation_list = mutation_id.split("_")
                if "INS" in mutation or "DEL" in mutation:
                    signature = "indels"
                elif len(mutation_list) == 2:
                    signature = "snp"
                else:
                    signature = "dnp"

                filter_key = mutation.split(":")[-1]
                for k, v in sample_data.items():
                    mutation_type = mutation_list[0]
                    mutation_subtype = np.nan if len(mutation_list) < 2 else mutation_list[1]
                    mutation_indel_size = np.nan if len(mutation_list) < 3 else mutation_list[2]
                    mutation_repeat_size = np.nan if len(mutation_list) < 4 else mutation_list[3]
                    data_list.append({
                        "chromosome": chromosome
                        , "name": name_map[name]
                        , "Type": mutation_type
                        , "SubType": mutation_subtype
                        , "IndelSize": mutation_indel_size
                        , "RepeatSize": mutation_repeat_size
                        , 'signature': signature
                        , "filter": filter_key
                        , "mutation_id": mutation_id
                        , "sample": k
                        , "count": v})
    return pd.DataFrame(data_list).query("count != 0")


def calculate_percentage(x):
    x["percentage"] = x["count"] / x["count"].sum()
    return x


def split_dnp(ds):
    ds["left"] = ds["Type"].split(">")[0]
    ds["right"] = ds["Type"].split(">")[1]
    return ds


def split_indels(ds):
    mutation_id = ds["mutation_id"].split("_")
    ds["Type"] = mutation_id[0]
    ds["SubType"] = mutation_id[1]
    ds["Length"] = mutation_id[2]
    ds["RepeatSize"] = mutation_id[3]
    return ds


if __name__ == "__main__":

    # store start time for benchmarking
    start_time = pd.to_datetime(time.time(), unit="s")

    # setup logger
    root_logger, log_formatter = get_cmri_logger()

    # setup arguments
    parser = argparse.ArgumentParser(description='region indexer for unknown chromosomes', epilog=epilog_text,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--silent', action='store_true', help='Starts in silent mode, no message will be output.')
    parser.add_argument('-d', '--debug', action='store_true', help='Shows debug info')
    parser.add_argument('-o', '--output', type=str, help='Output file', default="../data/")
    parser.add_argument('-i', '--input', type=str, help='Input directory', default="../data/GEAR_VCA/")

    # parse arguments and set logger
    args = parser.parse_args()
    consoleHandler = logging.StreamHandler(sys.stdout)

    if args.debug:
        root_logger.setLevel(logging.DEBUG)

    if args.silent:
        consoleHandler.setLevel(logging.ERROR)

    consoleHandler.setFormatter(log_formatter)
    root_logger.addHandler(consoleHandler)

    print_cmri_welcome("Mutation Analysis")

    directory_exists(args.input, True)

    data_list = list()
    for f in find_files(args.input, "json"):
        if os.path.basename(f)[0] == ".":
            continue
        logging.info("Processing file %s", f)
        df_tmp = read_gear_vca(f)
        df_tmp["group"] = f.split("/")[-2]
        data_list.append(df_tmp)

    if len(data_list) ==0:
        logging.error("data does not exists")
        exit(-1)

    df = pd.concat(data_list)
    df_snp = df.query("signature=='snp'").groupby(["Type", "SubType", "group"])["count"].sum().reset_index().groupby(
        "group").apply(calculate_percentage)
    file_name=os.path.join(args.output, "df_snp.csv")
    df_snp.to_csv(file_name)
    logging.info("Saved csv - %s", file_name)

    df_dnp = df.query("signature=='dnp'").groupby(["Type", "group"])["count"].sum().reset_index().groupby(
        "group").apply(calculate_percentage)
    df_dnp = df_dnp.apply(split_dnp, axis=1)
    file_name=os.path.join(args.output, "df_dnp.csv")
    df_dnp.to_csv(file_name)
    logging.info("Saved csv - %s", file_name)


    df_indels = df.query("signature=='indels'").groupby(["mutation_id", "group"])["count"].sum().reset_index().groupby(
        "group").apply(calculate_percentage)
    df_indels = df_indels.apply(split_indels, axis=1)
    file_name=os.path.join(args.output, "df_indels.csv")
    df_indels.to_csv(file_name)
    logging.info("Saved csv - %s", file_name)



    logging.info("Total computation time: %s", str(pd.to_datetime(time.time(), unit="s") - start_time))
