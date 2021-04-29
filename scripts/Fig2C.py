import argparse
import os.path
import sys
import time
import pandas as pd
import matplotlib.transforms as transforms
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
    parser.add_argument('-o', '--output', type=str, help='Output file', default="../figures/fig2C.pdf")

    # parse arguments and set logger
    args = parser.parse_args()
    consoleHandler = logging.StreamHandler(sys.stdout)

    if args.debug:
        root_logger.setLevel(logging.DEBUG)

    if args.silent:
        consoleHandler.setLevel(logging.ERROR)

    consoleHandler.setFormatter(log_formatter)
    root_logger.addHandler(consoleHandler)

    print_cmri_welcome("Fig 2C")

    root_file_name = os.path.splitext(args.output)[0]
    df_snp=pd.read_csv("../data/df_snp.csv",index_col=0)
    df_dnp=pd.read_csv("../data/df_dnp.csv",index_col=0)
    df_indels=pd.read_csv("../data/df_indels.csv",index_col=0)

    sns.set(font_scale=1)
    sns.set_style("whitegrid")

    plt.rc('axes', labelsize=20)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=20)  # fontsize of the x tick labels
    plt.rc('ytick', labelsize=20)  # fontsize of the y tick labels

    fig, ax = plt.subplots(ncols=6, figsize=(40, 5))

    for i, tp in enumerate(["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]):

        ax_i = ax[i]
        g = sns.barplot(y="count"
                        , x="SubType"
                        , data=df_snp.query("Type == '%s'" % (tp))
                        # ,hue="category"
                        , ax=ax_i
                        , color=sns.color_palette()[i]
                        )
        #    ax_i.set_ylim(0,0.05)
        ax_i.set_ylim(0, df_snp["count"].max())
        if i != 0:
            g.set(yticklabels=[])

        ax_i.set_xticklabels(ax_i.get_xticklabels(), rotation=90)
        ax_i.set_title(tp, fontdict={"fontweight": "bold","fontsize":20})

        ax_i.set_ylabel("")
        ax_i.set_xlabel("")
        trans = transforms.blended_transform_factory(ax_i.transAxes, ax_i.transAxes)
        rect = plt.Rectangle((0, 1), 1, 0.1, linewidth=1, edgecolor='none', facecolor=sns.color_palette()[i],
                             clip_on=False, transform=trans)
        ax_i.add_patch(rect)

    trans = transforms.blended_transform_factory(fig.transFigure, ax_i.transAxes)
    plt.text(0.1, 0.5, "Count", ha="center", va="center", transform=trans, rotation=90,fontdict={"fontsize":20})
    plt.text(0.5, -0.25, "Context", ha="center", va="top", transform=trans,fontdict={"fontsize":20})
    plt.text(0.5, 1.1, "Single base substitutions", ha="center", va="bottom", transform=trans,fontdict={"fontsize":20})


    plt.subplots_adjust(bottom=0.15, wspace=0.0)
    fig_name=root_file_name+"_1.pdf"
    plt.savefig(fig_name, format="pdf", bbox_inches='tight',pad_inches=0.25)
    logging.info("Saved figure - %s", fig_name)


    context = ["AC", "AT", "CC", "CG", "CT", "GC", "TA", "TC", "TG", "TT"]
    fig, ax = plt.subplots(ncols=len(context), figsize=(40, 5))

    for i, left in enumerate(context):
        ax_i = ax[i]
        g = sns.barplot(y="count"
                        , x="right"
                        , data=df_dnp.query("left == '%s'" % (left))
                        # ,hue="category"
                        , ax=ax_i
                        # ,palette=pal
                        , color=sns.color_palette()[i]
                        )
        ax_i.set_ylim(0, df_dnp["count"].max())
        if i != 0:
            g.set(yticklabels=[])

        ax_i.set_xticklabels(ax_i.get_xticklabels(), rotation=90)

        ax_i.set_title("%s>NN" % left, fontdict={"fontweight": "bold","fontsize":20})

        ax_i.set_ylabel("")
        ax_i.set_xlabel("")

        trans = transforms.blended_transform_factory(ax_i.transAxes, ax_i.transAxes)
        rect = plt.Rectangle((0, 1), 1, 0.1, linewidth=1, edgecolor='none', facecolor=sns.color_palette()[i],
                             clip_on=False, transform=trans)
        ax_i.add_patch(rect)

    trans = transforms.blended_transform_factory(fig.transFigure, ax_i.transAxes)
    plt.text(0.1, 0.5, "Count", ha="center", va="center", transform=trans, rotation=90,fontdict={"fontsize":20})
    plt.text(0.5, -0.25, "Context", ha="center", va="top", transform=trans,fontdict={"fontsize":20})
    plt.text(0.5, 1.1, "Doublet base substitutions", ha="center", va="bottom", transform=trans,fontdict={"fontsize":20})

    plt.subplots_adjust(bottom=0.15, wspace=0.0)
    fig_name=root_file_name+"_2.pdf"
    plt.savefig(fig_name, format="pdf", bbox_inches='tight',pad_inches=0.25)
    logging.info("Saved figure - %s", fig_name)

    colors = [
        "#ffb86a"
        , "#ff7000"
        , "#afda88"
        , "#079e20"
        , "#ffc7b8"
        , "#ff8870"
        , "#ff2223"
        , "#d40000"
        , "#cee6f3"
        , "#85c8e3"
        , "#00a0cf"
        , "#006db2"
        , "#aeb4d5"
        , "#7f89bd"
        , "#9680cb"
        , "#5e459a"
    ]

    width_ratios = [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 2, 3, 4]
    fig, ax = plt.subplots(ncols=len(width_ratios), figsize=(40, 5),
                           gridspec_kw={'width_ratios': width_ratios})
    for i, (Type, SubType, Length) in enumerate([
        ("DEL", "C", "1")
        , ("DEL", "T", "1")
        , ("INS", "C", "1")
        , ("INS", "T", "1")
        , ("DEL", "repeats", "2")
        , ("DEL", "repeats", "3")
        , ("DEL", "repeats", "4")
        , ("DEL", "repeats", "5+")
        , ("INS", "repeats", "2")
        , ("INS", "repeats", "3")
        , ("INS", "repeats", "4")
        , ("INS", "repeats", "5+")
        , ("DEL", "MH", "2")
        , ("DEL", "MH", "3")
        , ("DEL", "MH", "4")
        , ("DEL", "MH", "5+")
    ]):
        ax_i = ax[i]
        g = sns.barplot(y="count"
                        , x="RepeatSize"
                        , data=df_indels.query(
                "Type == '%s' and SubType == '%s' and Length == '%s'" % (Type, SubType, Length))
                        , ax=ax_i
                        # ,hue="category"
                        # ,palette=pal
                        , color=colors[i]
                        )
        ymax = df_indels["count"].max()
        ax_i.set_ylim(0, 1.1 * ymax)
        ax_i.set_title("%s" % (SubType if SubType in ["C", "T"] else Length),
                       fontdict={"fontweight": "bold", "fontsize": 20})

        if i == 0:
            ax_i.set_ylabel("")
        else:
            ax_i.set_ylabel("")
            g.set(yticklabels=[])

        ax_i.set_xlabel("")
        ax_i.set_xticklabels(ax_i.get_xticklabels(), rotation=90)

        trans = transforms.blended_transform_factory(ax_i.transAxes, ax_i.transAxes)
        rect = plt.Rectangle((0, 1), 1, 0.1, linewidth=1, edgecolor='none', facecolor=colors[i], clip_on=False,
                             transform=trans)
        ax_i.add_patch(rect)

        trans = transforms.blended_transform_factory(ax_i.transAxes, ax_i.transAxes)
        label = "insertion" if Type == "INS" else "deletion"
        if i in [0, 2, 5, 9, 14]:

            if i in [0, 2]:
                ax_i.text(1, 1.1, "1bp " + label, ha="center", va='bottom', transform=trans, fontdict={"fontsize": 20})
            if i in [5, 9]:
                ax_i.text(1, 1.1, ">1bp " + label, ha="center", va='bottom', transform=trans, fontdict={"fontsize": 20})
            if i == 14:
                ax_i.text(1, 1.1, label, ha="center", va='bottom', transform=trans, fontdict={"fontsize": 20})
                ax_i.text(1, 1.2, "micro-homologies", ha="center", va='bottom', transform=trans, fontdict={"fontsize": 20})

    trans = transforms.blended_transform_factory(fig.transFigure, ax_i.transAxes)
    plt.text(0.1, 0.5, "Count", ha="right", va="center", transform=trans, rotation=90, fontdict={"fontsize": 20})
    plt.text(0.5, -0.25, "Number of repeat units", ha="center", va="top", transform=trans, fontdict={"fontsize": 20})
    plt.text(0.5, 1.2, "Small insertions/deletions", ha="center", va="bottom", transform=trans,
             fontdict={"fontsize": 20})

    plt.subplots_adjust(bottom=0.15, wspace=0.0)
    fig_name=root_file_name+"_3.pdf"
    plt.savefig(fig_name, format="pdf", bbox_inches='tight',pad_inches=0.25)
    logging.info("Saved figure - %s", fig_name)


    logging.info("Total computation time: %s", str(pd.to_datetime(time.time(), unit="s") - start_time))
