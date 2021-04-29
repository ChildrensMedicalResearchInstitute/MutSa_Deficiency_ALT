import argparse
import sys
import time
import pandas as pd
from matplotlib.ticker import FuncFormatter
import matplotlib.transforms as transforms
import seaborn as sns
import matplotlib.pyplot as plt
from utils import *


def scientific(x, pos):
    # x:  tick value - ie. what you currently see in yticks
    # pos: a position - ie. the index of the tick (from 0 to 9 in this example)
    return '%.5f' % x


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
    parser.add_argument('-o', '--output', type=str, help='Output file', default="../figures/fig2D.pdf")

    # parse arguments and set logger
    args = parser.parse_args()
    consoleHandler = logging.StreamHandler(sys.stdout)

    if args.debug:
        root_logger.setLevel(logging.DEBUG)

    if args.silent:
        consoleHandler.setLevel(logging.ERROR)

    consoleHandler.setFormatter(log_formatter)
    root_logger.addHandler(consoleHandler)

    print_cmri_welcome("Fig 2D")


    scientific_formatter = FuncFormatter(scientific)

    groupby = ["sample_id", "motif_type"]
    groupby += ["category", "region"]

    df=pd.read_csv("../data/df_mutation_count.csv",index_col=0)
    df_tmp = df.query("motif_type == 'TTAGGG' and (region == 'telomeric' or region == 'interstitial')").groupby(groupby)[
        "per_Mbp"].sum().reset_index()

    sns.set_context("talk")
    sns.set(font_scale=2)
    sns.set_style("white")

    sns.set_context("talk")
    g = sns.catplot(x="motif_type"
                    , y="per_Mbp"
                    , data=df_tmp
                    , hue="category"
                    , col="region"
                    , palette=["#717171", "#9f2830"]
                    , kind="bar"
                    , col_order=["telomeric", "interstitial"]
                    , height=5
                    , aspect=0.5
                    , sharey=False
                    , capsize=0.15
                    )

    g.set_titles("")

    for i, ax in enumerate(g.axes.flatten()):
        # ax.set_title("Hexmer proportion (interstitial region)")
        if i == 0:
            ax.set_ylabel(" Motif count / Mbp")
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: int(x)))
        trans = transforms.blended_transform_factory(ax.transAxes, ax.transAxes)
        ax.text(0.5, -0.025, "telomeric" if i == 0 else "interstitial", ha="center", va='top', transform=trans
                # ,fontsize="x-small"
                )
        ax.set_xlabel("TTAGGG", labelpad=15)
        ax.set(xticklabels=[])
        ax.spines['bottom'].set_color('0')
        ax.spines['top'].set_color('1')
        ax.spines['right'].set_color('1')
        ax.spines['left'].set_color('0')
        ax.tick_params(direction='out', width=3, bottom=True, left=True)
        ax.grid(False)
    g._legend.set_title("")
    g._legend.set_bbox_to_anchor((1, 0.15))

    plt.subplots_adjust(bottom=0.15, wspace=0.9)

    plt.savefig(args.output, format="pdf", bbox_inches='tight')

    logging.info("Saved figure - %s", args.output)
    logging.info("Total computation time: %s", str(pd.to_datetime(time.time(), unit="s") - start_time))
