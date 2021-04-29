import argparse
import json
import sys
import time

import numpy as np
import pandas as pd
from Bio.Seq import Seq

from utils import *


def get_total_read(input_path):
    '''
    Return a map with the number of reads per sample.

    Parameters
    ----------
    input_path : str
        a path to GEAR's log output.

    Returns
    -------
    dict
        a dictionary with pairs {sample:count}.
    '''
    read_map = dict()
    # iterate all the log files
    for f in find_files(input_path, "log"):
        # iterate the log file
        with open(f) as log_file:
            # find the Total sequences analysed field
            for l in log_file:
                if "Total sequences analysed:" in l:
                    # the sample id is the GEAR directory name
                    sample_id = f.split("/")[-2]
                    read_map[sample_id] = int(l.split(" ")[-1])
    return read_map


def parse_motifs(item, field, chromosome):
    data_list = list()
    for motif, values in item[field].items():
        # temporaty dictionary, we will append the pair {quality value : count} for each qv
        data_dict = {"chromosome": chromosome, "name": item["name"], "total_reads": item["count"], "motif": motif,
                     "metric": field}
        count = 0
        tmp_dict = dict()
        # append histogram values
        for k, v in values.items():
            tmp_dict[int(k)] = v
            count += v
        # create a series (we will use it sort the values and make a cummulative sum) and append total motif count
        ds = pd.Series(tmp_dict)
        data_dict.update(ds.sort_index(ascending=False).cumsum().to_dict())
        data_dict["count"] = count
        data_list.append(data_dict)
    return data_list


def read_gear_motif_count(input_path):
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
    with open(input_path) as f:
        # load the data and create an empty list
        data = json.load(f)
        data_list = list()
        # iterate each chromosome in the file
        for chromosome, clist in data.items():
            # iterate each region in the chromosome (for example, telomere are p or q, etc).
            for item in clist:
                # total number of reads in this refion

                # motif is a histogram of motif as function of mean base quality of the motif.
                data_list += parse_motifs(item, "motifs", chromosome)
                data_list += parse_motifs(item, "regex", chromosome)

                # iterate the map of motifs.
    # return only positive count rows.
    return pd.DataFrame(data_list).query("count != 0")


def calculate_ratio(ds):
    control=ds["Con"].replace(0,np.nan).values
    n=control.shape[0]
    treatment=ds["MSH6 KO"].values
    m=treatment.shape[0]
    ratio=(treatment*np.ones((m,n))/(control*np.ones((n,m))).T).flatten()
    for i,r in enumerate(ratio):
        ds.loc[i]=r
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
    parser.add_argument('-i', '--input', type=str, help='Input directory', default="../data/GEAR_MOTIF_COUNT/")
    parser.add_argument('-q', '--quality_value_threshold', type=int, help='Set the quality value threshold for the '
                                                                          'motif filtering',
                        default=35)

    # parse arguments and set logger
    args = parser.parse_args()
    consoleHandler = logging.StreamHandler(sys.stdout)

    if args.debug:
        root_logger.setLevel(logging.DEBUG)

    if args.silent:
        consoleHandler.setLevel(logging.ERROR)

    consoleHandler.setFormatter(log_formatter)
    root_logger.addHandler(consoleHandler)

    print_cmri_welcome("Motif Count Analysis")

    qv = args.quality_value_threshold

    region_map = {'mapq_fail_': "telomeric"
        , 'qv_fail_': "qv_fail"
        , 'unmapped_': "telomeric"
        , 'other_': "other"
        , 'p_telomere': 'telomeric'
        , 'q_telomere': 'telomeric'
        , 'c_interstitial': 'interstitial'
                  }
    category_map = {"A1L1": "Con"
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

    motif_type = np.sort(
        ['CTAGGG', 'ATAGGG', 'CTAGGG', 'TTAGGC', 'GTAGGG', 'TAAGGG', 'TCAGGG', 'TTAGGA', 'TGAGGG', 'TTAAGG', 'TTACGG',
         'TTAGAG', 'TTAGCG', 'TTAGGA', 'TTAGGC', 'TTAGGG', 'TTAGGT', 'TTAGTG', 'TTATGG', 'TTCGGG', 'TTGGGG', 'TTTGGG'])

    strand_map = dict()
    for m in motif_type:
        strand_map[m] = m
        strand_map[str(Seq(m).reverse_complement())] = m

    regex_list = ["TTAGGG(G{%d})[AC]" % (x + 1) for x in np.arange(5)]
    regex_list += ["TTAGGG(G{6,})[AC]"]
    regex_list += ["TTAGGG(G{%d})TAGGG" % (x + 1) for x in np.arange(5)]
    regex_list += ["TTAGGG(G{6,})TAGGG"]

    for m in regex_list:
        strand_map[m] = m

    strand_map.update({"[TG](C{%d})CCCTAA" % (x + 1): "TTAGGG(G{%d})[AC]" % (x + 1) for x in np.arange(5)})
    strand_map.update({"[TG](C{6,})CCCTAA": "TTAGGG(G{6,})[AC]"})
    strand_map.update({"CCCTA(C{%d})CCCTAA" % (x + 1): "TTAGGG(G{%d})TAGGG" % (x + 1) for x in np.arange(5)})
    strand_map.update({"CCCTA(C{6,})CCCTAA": "TTAGGG(G{6,})TAGGG"})
    # get the read size of each sample.
    read_factor_map = {k: v / np.min(list(sequence_size_map.values())) for k, v in sequence_size_map.items()}

    motif_group_map = dict()
    for k, v in strand_map.items():
        if "(" in v:
            motif_group_map[k] = "regex"
        elif v in ["TTAGGG", "TCAGGG", "TGAGGG", "TTGGGG"]:
            motif_group_map[k] = "canonical"
        elif v[-3:] == 'GGG':
            motif_group_map[k] = "Gs"
        else:
            motif_group_map[k] = "other"

    # Iterate all the files in path
    data_list = list()
    for f in find_files(args.input, "json"):
        if os.path.basename(f)[0] == ".":
            continue
        logging.info("Processing file %s", f)
        # read the json files
        df_tmp = read_gear_motif_count(f).query("count != 0")
        # sample id is assume to be the GEAR directory name
        sample_id = f.split("/")[-3]+f.split("/")[-2]
        df_tmp["sample_id"] = sample_id
        # map each sample id to a category
        df_tmp["category"] = category_map[sample_id]

        df_tmp["motif_type"] = df_tmp["motif"].map(strand_map)
        df_tmp["region"] = df_tmp["name"].map(region_map)
        df_tmp["reads"] = sequence_size_map[sample_id]
        df_tmp["motif_group"] = df_tmp["motif"].map(motif_group_map)
        df_tmp["proportion"] = df_tmp[qv] / (df_tmp["reads"] * 150 / 6)
        df_tmp["raw_proportion"] = df_tmp[qv] / df_tmp["reads"]
        df_tmp["per_Mbp"] = 1000000 * df_tmp[qv] / (df_tmp["reads"] * 150)
        df_tmp["relative_proportion"] = df_tmp[qv] / df_tmp["total_reads"]
        df_tmp["telomere_percentage"] = df_tmp[qv] / (df_tmp.query("motif_type == 'TTAGGG'")[qv].sum())
        df_tmp["telomere_proportion"] = df_tmp[qv] / (
            df_tmp.query("motif_type == 'TTAGGG' and region=='telomeric'")[qv].sum())
        df_tmp["normalized_count"] = df_tmp[qv] * read_factor_map[sample_id]
        df_tmp["read_factor"] = read_factor_map[sample_id]
        data_list.append(df_tmp)

    df = pd.concat(data_list)


    output_file="../data/df_mutation_count.csv"
    df.to_csv(output_file)
    logging.info("File Saved: %s", output_file)

    df = pd.merge(df, df.query("motif_type == 'TTAGGG'").groupby("sample_id")[qv].sum().reset_index().rename(
        columns={qv: "TTAGGG_total"}), on="sample_id", how="left")

    df_control_treatment = df.groupby(["region", "sample_id", "category", "motif_type", "motif_group"])[
        "raw_proportion"].sum().reset_index()
    df_control_treatment = pd.pivot_table(df_control_treatment, index=["motif_type", "motif_group", "region"],
                                          columns=['category', "sample_id"], values="raw_proportion").reset_index()

    df_ratio = df_control_treatment.apply(calculate_ratio, axis=1)

    del df_ratio["Con"]
    del df_ratio["MSH6 KO"]

    df_ratio.columns = df_ratio.columns.droplevel(1)

    output_file="../data/df_ratio.csv"
    pd.melt(df_ratio, value_vars=np.arange(15), id_vars=["motif_type", "motif_group", "region"]).to_csv(output_file)
    logging.info("File Saved: %s", output_file)

    logging.info("Total computation time: %s", str(pd.to_datetime(time.time(), unit="s") - start_time))
