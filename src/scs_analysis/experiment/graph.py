import warnings

warnings.filterwarnings("ignore")
import os
from typing import List
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.ticker import (
    FuncFormatter,
    LogLocator,
    MultipleLocator,
    ScalarFormatter,
    StrMethodFormatter,
)

sns.set_theme()
sns.set_context("paper")
sns.set_style("whitegrid")


from scs_analysis.experiment.distance_calculator import ORDERING

from scs_analysis.experiment.experiment import SCS_FAST, BCDG, BCDN, MCS


METHOD_MAP = {SCS_FAST: "SCS", BCDG: "BCD (GSCM)", BCDN: "BCD (No GSCM)", MCS: "MCS"}


def format_title(header, root):
    if "dcm_source_trees" in root:
        try:
            taxa, max_subprob = map(int, root.split("/")[-3::2])
            suffix = f" SCS DCM Exact (n={taxa}, m={max_subprob})"
        except:
            (taxa,) = map(int, root.split("/")[-2::2])
            suffix = f" SCS DCM Exact (n={taxa})"
    elif "iq_source_trees" in root:
        try:
            taxa, max_subprob = map(int, root.split("/")[-3::2])
            suffix = f" SCS DCM IQ-TREE (n={taxa}, m={max_subprob})"
        except:
            (taxa,) = map(int, root.split("/")[-2::2])
            suffix = f" SCS DCM IQ-TREE (n={taxa})"
    elif "SuperTriplets" in root:
        try:
            d, k = map(lambda x: int(x[1:]), root.split("/")[-2:])
            suffix = f" SuperTriplets (d={d}, k={k})"
        except:
            (d,) = map(lambda x: int(x[1:]), root.split("/")[-1:])
            suffix = f" SuperTriplets (d={d})"
    elif "SMIDGenOutgrouped" in root:
        if "10000" in root:
            suffix = f" SMIDGenOG-5500"
        else:
            try:
                taxa, scaffold = map(int, root.split("/")[-2:])
                suffix = f" SMIDGenOG (n={taxa}, s={scaffold}%)"
            except:
                (taxa,) = map(int, root.split("/")[-1:])
                suffix = f" SMIDGenOG (n={taxa})"
    else:
        raise ValueError
    return header + " Results for" + suffix


def graph_experiment(image_folder: str, root: str, distance_files: List[str]):
    print(root)
    distance_files = sorted(
        distance_files, key=lambda x: (ORDERING.get(x[:-12], float("inf")), x)
    )
    methods = [file[:-27] for file in distance_files]
    distance_file_paths = [root + "/" + file for file in distance_files]

    header = [
        "Model Tree",
        "Source Tree",
        "Wall Time",
        "CPU Time",
        "RF Distance",
        "Matching Cluster Distance",
        "F1 Score",
        "Bifurcating RF Distance",
        "Bifurcating Matching Cluster Distance",
        "Bifurcating F1 Score",
        "Tree",
    ]
    drops = ["Model Tree", "Source Tree", "Tree"]
    dfs = []
    for method, distance_file in zip(methods, distance_file_paths):
        if method not in METHOD_MAP:
            continue

        df = pd.read_csv(distance_file, delimiter="\t", names=header)
        df = df.drop(columns=drops)
        df["Method"] = METHOD_MAP[method]
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)

    image_directory = image_folder + "/".join(root.split("/")[1:]) + "/"

    if not os.path.exists(image_directory):
        os.makedirs(image_directory)

    for col in set(header).difference(drops):
        plt.figure()
        if "Time" in col:
            scaled, unit = make_time_scale_reasonable(df, col)
            min_log_y_tick = scaled[col].min()
            g = sns.boxplot(scaled, x="Method", y=col)
            plt.ylabel(col + f" ({unit})")
        else:
            min_log_y_tick = df[col].min()
            g = sns.boxplot(df, x="Method", y=col)
        if "F1" not in col:
            plt.ylim(0, None)
        plt.title(format_title(col, root))
        # sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.savefig(image_directory + f"{col}_all.pdf")
        try:
            g.set_yscale("log")
            plt.ylim(min_log_y_tick / 2, None)
            plt.yticks()
            g.yaxis.set_major_formatter(FuncFormatter(format_tick))
            g.yaxis.set_minor_formatter(ScalarFormatter())
            g.yaxis.set_minor_formatter(FuncFormatter(format_tick))
            g.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
            g.yaxis.set_minor_locator(LogLocator(base=10, subs=(1, 2, 5)))
            # sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.savefig(image_directory + f"{col}_all_log.pdf")
        except:
            pass
        plt.close()

        mcs_mask = df["Method"] == METHOD_MAP[MCS]

        mcs_data = df[mcs_mask]
        if len(mcs_data) > 0:
            plt.figure()
            if "Time" in col:
                scaled, unit = make_time_scale_reasonable(mcs_data, col)
                min_log_y_tick = scaled[col].min()
                g = sns.boxplot(scaled, x="Method", y=col)
                plt.ylabel(col + f" ({unit})")
            else:
                min_log_y_tick = df[col].min()
                g = sns.boxplot(mcs_data, x="Method", y=col)
            if "F1" not in col:
                plt.ylim(0, None)
            plt.title(format_title(col, root))
            # sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.savefig(image_directory + f"{col}_mcs.pdf")
            try:
                g.set_yscale("log")
                plt.ylim(min_log_y_tick / 2, None)
                plt.yticks()
                g.yaxis.set_major_formatter(FuncFormatter(format_tick))
                g.yaxis.set_minor_formatter(ScalarFormatter())
                g.yaxis.set_minor_formatter(FuncFormatter(format_tick))
                g.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
                g.yaxis.set_minor_locator(LogLocator(base=10, subs=(1, 2, 5)))
                # sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
                plt.tight_layout()
                plt.savefig(image_directory + f"{col}_mcs_log.pdf")
            except:
                pass
            plt.close()

        no_mcs_data = df[~mcs_mask & (df["Method"] != METHOD_MAP[BCDN])]
        plt.figure()
        if "Time" in col:
            scaled, unit = make_time_scale_reasonable(no_mcs_data, col)
            min_log_y_tick = scaled[col].min()
            g = sns.boxplot(scaled, x="Method", y=col)
            plt.ylabel(col + f" ({unit})")
        else:
            min_log_y_tick = df[col].min()
            g = sns.boxplot(no_mcs_data, x="Method", y=col)
        if "F1" not in col:
            plt.ylim(0, None)
        plt.title(format_title(col, root))
        # sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.savefig(image_directory + f"{col}_others.pdf")
        try:
            g.set_yscale("log")
            plt.ylim(min_log_y_tick / 2, None)
            plt.yticks()
            g.yaxis.set_major_formatter(FuncFormatter(format_tick))
            g.yaxis.set_minor_formatter(ScalarFormatter())
            g.yaxis.set_minor_formatter(FuncFormatter(format_tick))
            g.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
            g.yaxis.set_minor_locator(LogLocator(base=10, subs=(1, 2, 5)))
            # sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.savefig(image_directory + f"{col}_others_log.pdf")
        except:
            pass
        plt.close()

        scs_mask = df["Method"] == METHOD_MAP[SCS_FAST]

        scs_data = df[scs_mask]
        plt.figure()
        if "Time" in col:
            scaled, unit = make_time_scale_reasonable(scs_data, col)
            min_log_y_tick = scaled[col].min()
            g = sns.boxplot(scaled, x="Method", y=col)
            plt.ylabel(col + f" ({unit})")
        else:
            min_log_y_tick = df[col].min()
            g = sns.boxplot(scs_data, x="Method", y=col)
        if "F1" not in col:
            plt.ylim(0, None)
        plt.title(format_title(col, root))
        # sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.savefig(image_directory + f"{col}_scs.pdf")
        try:
            g.set_yscale("log")
            plt.ylim(min_log_y_tick / 2, None)
            plt.yticks()
            g.yaxis.set_major_formatter(FuncFormatter(format_tick))
            g.yaxis.set_minor_formatter(ScalarFormatter())
            g.yaxis.set_minor_formatter(FuncFormatter(format_tick))
            g.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
            g.yaxis.set_minor_locator(LogLocator(base=10, subs=(1, 2, 5)))
            # sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.savefig(image_directory + f"{col}_scs_log.pdf")
        except:
            pass
        plt.close()


def make_time_scale_reasonable(df, column):
    df = df.copy()
    if df[column].max() < 600:
        return df, "seconds"
    df[column] /= 60
    if df[column].max() < 600:
        return df, "minutes"
    df[column] /= 60
    return df, "hours"


def load_data(folder):
    distance_files = list(
        filter(lambda x: x.endswith("_with_distances.tsv"), os.listdir(folder))
    )
    distance_files = sorted(
        distance_files, key=lambda x: (ORDERING.get(x[:-12], float("inf")), x)
    )
    methods = [file[:-27] for file in distance_files]
    distance_file_paths = [folder + "/" + file for file in distance_files]

    header = [
        "Model Tree",
        "Source Tree",
        "Wall Time",
        "CPU Time",
        "RF Distance",
        "Matching Cluster Distance",
        "F1 Score",
        "Bifurcating RF Distance",
        "Bifurcating Matching Cluster Distance",
        "Bifurcating F1 Score",
        "Tree",
    ]
    drops = ["Model Tree", "Source Tree", "Tree"]
    dfs = []
    for method, distance_file in zip(methods, distance_file_paths):
        if method not in METHOD_MAP:
            continue

        df = pd.read_csv(distance_file, delimiter="\t", names=header)
        df = df.drop(columns=drops)
        df["Method"] = METHOD_MAP[method]
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def graph_smidgenog_combined(image_folder: str):
    densities = (20, 50, 75, 100)
    for taxa in (100, 500, 1000):
        dfs = []
        for density in densities:
            df = load_data(f"results/SMIDGenOutgrouped/{taxa}/{density}")
            df["Scaffold Factor"] = density
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)

        graph_combined(
            image_folder + f"combined/SMIDGenOutgrouped/{taxa}/",
            df,
            f"results/SMIDGenOutgrouped/{taxa}",
            "Scaffold Factor",
        )


def graph_supertriplets_combined(image_folder: str):
    ks = (10, 20, 30, 40, 50)
    for d in (25, 50, 75):
        dfs = []
        for k in ks:
            df = load_data(f"results/SuperTripletsBenchmark/d{d}/k{k}")
            df["k"] = k
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)

        graph_combined(
            image_folder + f"combined/SuperTripletsBenchmark/d{d}/",
            df,
            f"results/SuperTripletsBenchmark/d{d}",
            "k",
        )


def graph_dcm_combined(image_folder: str):
    ms = (50, 100)
    for taxa in (500, 1000, 2000, 5000, 10000):
        dfs = []
        for m in ms:
            df = load_data(f"results/birth_death/{taxa}/dcm_source_trees/{m}/")
            df["m"] = m
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)

        graph_combined(
            image_folder + f"combined/birth_death/{taxa}/dcm_source_trees/",
            df,
            f"results/birth_death/{taxa}/dcm_source_trees",
            "m",
        )


def graph_iq_combined(image_folder: str):
    ms = (50, 100)
    for taxa in (500, 1000, 2000, 5000, 10000):
        dfs = []
        for m in ms:
            df = load_data(f"results/birth_death/{taxa}/iq_source_trees/{m}/")
            df["m"] = m
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)

        graph_combined(
            image_folder + f"combined/birth_death/{taxa}/iq_source_trees/",
            df,
            f"results/birth_death/{taxa}/iq_source_trees",
            "m",
        )


def graph_combined(image_directory, df, root, x_col):
    if not os.path.exists(image_directory):
        os.makedirs(image_directory)

    header = [
        "Model Tree",
        "Source Tree",
        "Wall Time",
        "CPU Time",
        "RF Distance",
        "Matching Cluster Distance",
        "F1 Score",
        "Bifurcating RF Distance",
        "Bifurcating Matching Cluster Distance",
        "Bifurcating F1 Score",
        "Tree",
    ]
    drops = ["Model Tree", "Source Tree", "Tree"]

    for col in set(header).difference(drops):
        plt.figure()
        if "Time" in col:
            scaled, unit = make_time_scale_reasonable(df, col)
            min_log_y_tick = scaled[col].min()
            g = sns.boxplot(scaled, x=x_col, hue="Method", y=col)
            plt.ylabel(col + f" ({unit})")
        else:
            min_log_y_tick = df[col].min()
            g = sns.boxplot(df, x=x_col, hue="Method", y=col)
        if "F1" not in col:
            plt.ylim(0, None)
        plt.title(format_title(col, root))
        sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.savefig(image_directory + f"{col}_all.pdf")
        try:
            g.set_yscale("log")
            plt.ylim(min_log_y_tick / 2, None)
            plt.yticks()
            g.yaxis.set_major_formatter(FuncFormatter(format_tick))
            g.yaxis.set_minor_formatter(ScalarFormatter())
            g.yaxis.set_minor_formatter(FuncFormatter(format_tick))
            g.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
            g.yaxis.set_minor_locator(LogLocator(base=10, subs=(1, 2, 5)))
            sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.savefig(image_directory + f"{col}_all_log.pdf")
        except:
            pass
        plt.close()

        mcs_mask = df["Method"] == METHOD_MAP[MCS]

        mcs_data = df[mcs_mask]
        if len(mcs_data) > 0:
            plt.figure()
            if "Time" in col:
                scaled, unit = make_time_scale_reasonable(mcs_data, col)
                min_log_y_tick = scaled[col].min()
                g = sns.boxplot(scaled, x=x_col, hue="Method", y=col)
                plt.ylabel(col + f" ({unit})")
            else:
                min_log_y_tick = mcs_data[col].min()
                g = sns.boxplot(mcs_data, x=x_col, hue="Method", y=col)
            if "F1" not in col:
                plt.ylim(0, None)
            plt.title(format_title(col, root))
            sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.savefig(image_directory + f"{col}_mcs.pdf")
            try:
                g.set_yscale("log")
                plt.ylim(min_log_y_tick / 2, None)
                plt.yticks()
                g.yaxis.set_major_formatter(FuncFormatter(format_tick))
                g.yaxis.set_minor_formatter(ScalarFormatter())
                g.yaxis.set_minor_formatter(FuncFormatter(format_tick))
                g.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
                g.yaxis.set_minor_locator(LogLocator(base=10, subs=(1, 2, 5)))
                sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
                plt.tight_layout()
                plt.savefig(image_directory + f"{col}_mcs_log.pdf")
            except:
                pass
            plt.close()

        no_mcs_data = df[~mcs_mask & (df["Method"] != METHOD_MAP[BCDN])]
        plt.figure()
        if "Time" in col:
            scaled, unit = make_time_scale_reasonable(no_mcs_data, col)
            min_log_y_tick = scaled[col].min()
            g = sns.boxplot(scaled, x=x_col, hue="Method", y=col)
            plt.ylabel(col + f" ({unit})")
        else:
            min_log_y_tick = no_mcs_data[col].min()
            g = sns.boxplot(no_mcs_data, x=x_col, hue="Method", y=col)
        if "F1" not in col:
            plt.ylim(0, None)
        plt.title(format_title(col, root))
        sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.savefig(image_directory + f"{col}_others.pdf")
        try:
            g.set_yscale("log")
            plt.ylim(min_log_y_tick / 2, None)
            plt.yticks()
            g.yaxis.set_major_formatter(FuncFormatter(format_tick))
            g.yaxis.set_minor_formatter(ScalarFormatter())
            g.yaxis.set_minor_formatter(FuncFormatter(format_tick))
            g.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
            g.yaxis.set_minor_locator(LogLocator(base=10, subs=(1, 2, 5)))
            sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.savefig(image_directory + f"{col}_others_log.pdf")
        except:
            pass
        plt.close()

        scs_mask = df["Method"] == METHOD_MAP[SCS_FAST]

        scs_data = df[scs_mask]
        plt.figure()
        if "Time" in col:
            scaled, unit = make_time_scale_reasonable(scs_data, col)
            min_log_y_tick = scaled[col].min()
            g = sns.boxplot(scaled, x=x_col, hue="Method", y=col)
            plt.ylabel(col + f" ({unit})")
        else:
            min_log_y_tick = df[col].min()
            g = sns.boxplot(scs_data, x=x_col, hue="Method", y=col)
        if "F1" not in col:
            plt.ylim(0, None)
        plt.title(format_title(col, root))
        sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.savefig(image_directory + f"{col}_scs.pdf")
        try:
            g.set_yscale("log")
            plt.ylim(min_log_y_tick / 2, None)
            plt.yticks()
            g.yaxis.set_major_formatter(FuncFormatter(format_tick))
            g.yaxis.set_minor_formatter(ScalarFormatter())
            g.yaxis.set_minor_formatter(FuncFormatter(format_tick))
            g.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
            g.yaxis.set_minor_locator(LogLocator(base=10, subs=(1, 2, 5)))
            sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
            plt.tight_layout()
            plt.savefig(image_directory + f"{col}_scs_log.pdf")
        except:
            pass
        plt.close()


def format_tick(value, pos):
    if float.is_integer(value):
        return str(int(value))
    return f"{value:.2f}"


def graph_results(image_folder: str, results_folder: str):
    print("GRAPHING", "SMIDGenOG")
    graph_smidgenog_combined(image_folder)
    print("GRAPHING", "SuperTriplets")
    graph_supertriplets_combined(image_folder)
    print("GRAPHING", "DCM")
    graph_dcm_combined(image_folder)
    print("GRAPHING", "IQ")
    graph_iq_combined(image_folder)

    for root, subdirs, files in os.walk(results_folder):
        distance_files = list(
            filter(lambda x: x.endswith("_with_distances.tsv"), files)
        )
        if len(distance_files) > 0:
            print("GRAPHING", root)
            graph_experiment(image_folder, root, distance_files)
