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
    parser.add_argument('-o', '--output', type=str, help='Output file', default="../figures/fig2B.pdf")

    # parse arguments and set logger
    args = parser.parse_args()
    consoleHandler = logging.StreamHandler(sys.stdout)

    if args.debug:
        root_logger.setLevel(logging.DEBUG)

    if args.silent:
        consoleHandler.setLevel(logging.ERROR)

    consoleHandler.setFormatter(log_formatter)
    root_logger.addHandler(consoleHandler)

    print_cmri_welcome("Fig 2B")

    df_signature_contribution = pd.read_csv("../data/df_signature_contribution.csv", index_col=0)

    plt.rc('axes', linewidth=3, labelsize=20)
    plt.rc('xtick', labelsize=20)  # fontsize of the x tick labels
    plt.rc('ytick', labelsize=20)  # fontsize of the y tick labels
    fig, ax = plt.subplots(nrows=3, figsize=(15, 8))

    df_tmp = df_signature_contribution.groupby(["group", "metric", "aetiology"])["contribution"].sum().reset_index()

    sns.set_context("talk")
    sns.set(font_scale=1)
    sns.set_style("white")

    for i, (tp, title) in enumerate(zip(["SBS", "DBS", "SID"],
                                        ["Single base substitutions (SBS)", "Doublet base substitutions (DBS)",
                                         "Small insertions/deletions (SID)"])):
        ax_i = ax[i]
        df_data = df_tmp.query("metric == '%s' and contribution>0 " % (tp))
        g = sns.barplot(y="aetiology"
                        , x="contribution"
                        , data=df_data
                        # ,hue="category"
                        , ax=ax_i
                        , order=['Mismatch repair', 'Unknown']
                        , color=["b", "r", "g"][i]
                        , capsize=0.15
                        )
        ax_i.set_xlim(0, 0.8)
        ax_i.set_title(title)

        ax_i.set_ylabel("")

        if i == 2:
            ax_i.set_xlabel("Relative contribution (%)")
        else:
            ax_i.set_xlabel("")
        if i != 2:
            g.set(xticklabels=[])
            ax_i.spines['bottom'].set_color('1')
            ax_i.tick_params(direction='out', width=3, bottom=False, left=True)
        else:
            ax_i.tick_params(direction='out', width=3, bottom=True, left=True)
            ax_i.spines['bottom'].set_color('0')

        ax_i.spines['top'].set_color('1')
        ax_i.spines['right'].set_color('1')
        ax_i.spines['left'].set_color('0')
        ax_i.grid(False)

    plt.subplots_adjust(left=0.1, bottom=0.2, right=1, top=1, wspace=2, hspace=0.2)

    plt.savefig(args.output, format="pdf", bbox_inches='tight', pad_inches=0.25)

    logging.info("Saved figure - %s", args.output)
    logging.info("Total computation time: %s", str(pd.to_datetime(time.time(), unit="s") - start_time))
