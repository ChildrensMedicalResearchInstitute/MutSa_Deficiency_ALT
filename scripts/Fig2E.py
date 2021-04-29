import argparse
import sys
import time
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms

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
    parser.add_argument('-o', '--output', type=str, help='Output file', default="../figures/fig2E.pdf")

    # parse arguments and set logger
    args = parser.parse_args()
    consoleHandler = logging.StreamHandler(sys.stdout)

    if args.debug:
        root_logger.setLevel(logging.DEBUG)

    if args.silent:
        consoleHandler.setLevel(logging.ERROR)

    consoleHandler.setFormatter(log_formatter)
    root_logger.addHandler(consoleHandler)

    print_cmri_welcome("Fig 2E")

    plt.rc('axes', linewidth=3, labelsize=20)
    plt.rc('xtick', labelsize=20)  # fontsize of the x tick labels
    plt.rc('ytick', labelsize=20)  # fontsize of the y tick labels
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.set(font_scale=2)
    sns.set_style("white")

    df_sampled=pd.read_csv("../data/df_ratio.csv")

    df_tmp = df_sampled.query("region == 'telomeric' and motif_group != 'other' and motif_group != 'regex'")
    ax = sns.barplot(x="motif_type"
                     , y="value"
                     , data=df_tmp
                     , capsize=0.15
                     , palette=["#9f2830"]
                     )

    ax.set_ylabel("Motif proportion ratio\n(MSH6 KO / Con)")
    ax.set_xlabel("")
#    ax.set_ylim(0, 2.02)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontdict={"horizontalalignment": "center"})
    ax.spines['bottom'].set_color('0')
    ax.spines['top'].set_color('1')
    ax.spines['right'].set_color('1')
    ax.spines['left'].set_color('0')
    ax.tick_params(direction='out', width=3, bottom=True, left=True)

    TTAGGG_ratio = df_sampled.query("region=='telomeric' and motif_type=='TTAGGG'")["value"].mean()

    trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(1, TTAGGG_ratio, "TTAGGG", ha="left", va='bottom', transform=trans, fontsize="x-small")
    ax.axhline(TTAGGG_ratio, 0, 1, ls='--', c="#000000", lw=3)

    ax.text(1, 1, "Null effect", ha="left", va='center', transform=trans, fontsize="x-small", rotation=-90)
    ax.axhline(1, 0, 1, ls='dotted', c="#000000", lw=3)

    plt.savefig(args.output, format="pdf", bbox_inches='tight')

    logging.info("Saved figure - %s", args.output)
    logging.info("Total computation time: %s", str(pd.to_datetime(time.time(), unit="s") - start_time))
