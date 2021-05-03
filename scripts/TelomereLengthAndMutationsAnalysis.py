import argparse
import json
import sys
import time
import gzip
import numpy as np
import pandas as pd

from utils import *


def read_gear_mutations(input_path):
    '''
    Parse a json file with GEAR file format.

    Parameters
    ----------
    input_path : str
        The json file path.

    Returns
    -------
    pandas.DataFrame
        a data frame with motif count information.
    '''

    # open the input file
    with gzip.open(input_path) as f:
        # load the data and create an empty list
        data = json.load(f)
        data_list = list()
        for i, item in enumerate(data):
            data_dict = {k: v for k, v in item.items() if not k in ["sbs", "indels"]}
            if "sbs" in item:
                for k, v in item["sbs"].items():
                    for vv in v:
                        data_list.append({**data_dict,
                                          **{"kind": "sbs", "pos": int(k), "ref_pos": vv["pos"], "ref_qv": vv["qv"],
                                             "mean_qv": vv["mean_qv"], "sbs": vv["value"],
                                             "is_trimmed": vv["is_trimmed"]}})
            if "indels" in item:
                for k, v in item["indels"].items():
                    for vv in v:
                        data_list.append({**data_dict, **{"kind": vv["kind"], "pos": int(k), "ref_pos": vv["pos"],
                                                          "mean_qv": vv["mean_qv"], "indel": vv["seq"],
                                                          "size": len(vv["seq"]), "is_trimmed": vv["is_trimmed"]}})
    return pd.DataFrame(data_list).drop_duplicates()


if __name__ == "__main__":

    # store start time for benchmarking
    start_time = pd.to_datetime(time.time(), unit="s")

    # setup arguments
    parser = argparse.ArgumentParser(description='region indexer for unknown chromosomes', epilog=epilog_text,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--silent', action='store_true', help='Starts in silent mode, no message will be output.')
    parser.add_argument('-d', '--debug', action='store_true', help='Shows debug info')
    parser.add_argument('-o', '--output', type=str, help='Output file', default="../data/")
    parser.add_argument('-i', '--input', type=str, help='Input directory', default="../data/GEAR_TELOMERE_MUTATION/")

    # parse arguments and set logger
    args = parser.parse_args()
    consoleHandler = logging.StreamHandler(sys.stdout)

    # setup logger
    root_logger, log_formatter = get_cmri_logger()

    if args.debug:
        root_logger.setLevel(logging.DEBUG)

    if args.silent:
        consoleHandler.setLevel(logging.ERROR)

    consoleHandler.setFormatter(log_formatter)
    root_logger.addHandler(consoleHandler)

    print_cmri_welcome("Mutation Analysis")

    directory_exists(args.input, True)

    data_list = list()
    for f in find_files(args.input, "json.gz",compressed=True):
        if "._" in f:
            continue
        logging.info("reading file: %s", f)
        df_tmp = read_gear_mutations(f)
        D1 = os.path.dirname(f)
        D2 = os.path.dirname(D1)
        df_tmp["sample"] = os.path.basename(D2) + os.path.basename(D1)

        data_list.append(df_tmp)
    df = pd.concat(data_list).reset_index(drop=True)

    variants = [c for c in df.columns if "G" in c]

    df_size = df[["mlen", "name", "seq", "sample"] + variants].drop_duplicates().groupby("sample")[
        variants + ["mlen"]].sum().reset_index()

    category_map = {
        "A1L1": "Con"
        , "A1L2": "Con"
        , "A2L1": "Con"
        , "A2L2": "Con"
        , "A3L1": "MSH6 KO"
        , "A3L2": "MSH6 KO"
        , "A4L1": "MSH6 KO"
        , "A4L2": "MSH6 KO"
    }
    sequence_size_map = {
        "A1L1": 178424347
        , "A1L2": 173957699
        , "A2L1": 171704752
        , "A2L2": 163338462
        , "A3L1": 171391614
        , "A3L2": 158500729
        , "A4L1": 192706327
        , "A4L2": 179876862
    }

    logging.info("Total computation time: %s", str(pd.to_datetime(time.time(), unit="s") - start_time))
