import argparse
import sys
import time
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from utils import *

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
    parser.add_argument('-o', '--output', type=str, help='Output file', default="../figures/fig2A.pdf")

    # parse arguments and set logger
    args = parser.parse_args()
    consoleHandler = logging.StreamHandler(sys.stdout)

    if args.debug:
        root_logger.setLevel(logging.DEBUG)

    if args.silent:
        consoleHandler.setLevel(logging.ERROR)

    consoleHandler.setFormatter(log_formatter)
    root_logger.addHandler(consoleHandler)

    print_cmri_welcome("Fig 2A")

    aetiology_map = {'Ultraviolet light exposure': np.nan,
                     'Aristolochic acid exposure': np.nan,
                     'Aflatoxin exposure': np.nan,
                     'Tobacco smoking': np.nan,
                     'Activity of APOBEC family of cytidine deaminases': "APOBEC",
                     'Damage by reactive oxygen species': "Reactive oxygen",
                     'Spontaneous deamination of 5-methylcytosine': "$^5$Me-C deamination",
                     'Defective DNA base excision repair due to NTHL1 mutations': "NTHL mutation",
                     'Tobacco smoking and other mutagens': np.nan,
                     'Defective DNA base excision repair due to MUTYH mutations': "MUTYH mutations",
                     'Platinum chemotherapy treatment': np.nan,
                     'Defective homologous recombination DNA damage repair': "HRD defects",
                     'Polymerase epsilon exonuclease domain mutations': "Polϵ mutation",
                     'Sequencing Artefacts': np.nan,
                     'Mismatch Repair': 'Mismatch repair',
                     'Unknown': 'Unknown'}

    plt.rc('axes', linewidth=3, labelsize=20)
    plt.rc('xtick', labelsize=20)  # fontsize of the x tick labels
    plt.rc('ytick', labelsize=20)  # fontsize of the y tick labels
    fig, ax = plt.subplots(figsize=(15, 8))

    df_signature_contribution = pd.read_csv("../data/df_signature_contribution.csv",index_col=0)
    df_comb_sum = df_signature_contribution.groupby(["aetiology", "metric"])["contribution"].sum().reset_index()
    df_comb_sum = pd.pivot_table(df_comb_sum, index=["aetiology"], columns="metric").fillna(0).sum(axis=1)
    df_comb_sum = df_comb_sum / df_comb_sum.sum()
    df_comb_sum = df_comb_sum.reset_index().rename(columns={0: "contribution"}).query("contribution > 0")
    df_comb_sum["aetiology"] = df_comb_sum["aetiology"].map(aetiology_map)

    sns.set_style("whitegrid")
    g = sns.barplot(
        ax=ax
        , y="aetiology"
        , x="contribution"
        , data=df_comb_sum
        , order=df_comb_sum.sort_values(by="contribution")["aetiology"].dropna()
        , color="black"
    )

    ax.set_xlim(0, 0.45)
    ax.set_xlabel("Relative contribution (%)")
    ax.set_ylabel("")
    ax.set_title("Combined Probabilistic Evidence")

    ax.tick_params(direction='out', width=3, bottom=True, left=True)
    ax.spines['bottom'].set_color('0')
    ax.spines['top'].set_color('1')
    ax.spines['right'].set_color('1')
    ax.spines['left'].set_color('0')
    ax.grid(False)

    plt.savefig(args.output, format="pdf", bbox_inches='tight', pad_inches=0.25)

    logging.info("Saved figure - %s", args.output)
    logging.info("Total computation time: %s", str(pd.to_datetime(time.time(), unit="s") - start_time))