import os

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

from scs_analysis.experiment.distance_calculator import ORDERING
from scs_analysis.experiment.experiment import BCDG, BCDN, MCS, SCS_FAST


sns.set_theme()
sns.set_context("paper")
sns.set_style("whitegrid")

METHOD_MAP = {SCS_FAST: "SCS", BCDG: "BCD (GSCM)", BCDN: "BCD (No GSCM)", MCS: "MCS"}


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


def graph_smidgenog_combined(image_folder: str, methods):
    # Original SMIDGenOG
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
            "Scaffold Factor",
            f"SMIDGenOG Dataset ({taxa} taxa)",
            methods,
        )

    # SMIDGenOG-5500
    df = load_data(f"results/SMIDGenOutgrouped/10000/0")
    graph_combined(
        image_folder + f"combined/SMIDGenOutgrouped/10000/0/",
        df,
        None,
        "SMIDGenOG-5500 Dataset",
        methods,
    )


def graph_supertriplets_combined(image_folder: str, methods):
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
            "k",
            f"Supertriplets Dataset ({d=}%)",
            methods,
        )


def graph_dcm_combined(image_folder: str, methods):
    ms = (50, 100)
    for taxa in (500, 1000, 2000, 5000, 10000):
        dfs = []
        for m in ms:
            df = load_data(f"results/birth_death/{taxa}/dcm_source_trees/{m}/")
            df["Max Subproblem Size"] = m
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)

        graph_combined(
            image_folder + f"combined/birth_death/{taxa}/dcm_source_trees/",
            df,
            "Max Subproblem Size",
            f"SCS DCM Exact Dataset ({taxa} taxa)",
            methods,
        )


def graph_iq_combined(image_folder: str, methods):
    ms = (50, 100)
    for taxa in (500, 1000, 2000, 5000, 10000):
        dfs = []
        for m in ms:
            df = load_data(f"results/birth_death/{taxa}/iq_source_trees/{m}/")
            df["Max Subproblem Size"] = m
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)

        graph_combined(
            image_folder + f"combined/birth_death/{taxa}/iq_source_trees/",
            df,
            "Max Subproblem Size",
            f"SCS DCM IQ-TREE Dataset ({taxa} taxa)",
            methods,
        )


def generate_hue_order(hues):
    order = ["SCS", "BCD (GSCM)", "BCD (No GSCM)", "MCS"]
    assert all(h in order for h in hues)
    return [o for o in order if o in hues]


def include_code(methods):
    codes = {
        METHOD_MAP[SCS_FAST]: "S",
        METHOD_MAP[BCDG]: "G",
        METHOD_MAP[BCDN]: "N",
        METHOD_MAP[MCS]: "M",
    }
    code = ""
    for method in methods:
        code += codes[method]
    return code


def graph_combined(image_directory, df, x_col, suptitle, methods):
    if not os.path.exists(image_directory):
        os.makedirs(image_directory)

    df = df[df["Method"].isin(methods)]

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(11, 7))

    fig.suptitle(suptitle)

    hue_order = generate_hue_order(df["Method"].unique())
    palette = {
        METHOD_MAP[SCS_FAST]: "C0",
        METHOD_MAP[BCDG]: "C1",
        METHOD_MAP[BCDN]: "C2",
        METHOD_MAP[MCS]: "C3",
    }

    # RF Distance
    rf_graph = sns.boxplot(
        df,
        x=x_col,
        y="RF Distance",
        hue="Method",
        hue_order=hue_order,
        palette=palette,  # type: ignore
        ax=ax[0, 0],
    )
    rf_graph.set_ylim(0, None)
    rf_graph.set_title("Robinson Foulds Distance (Lower is Better)")

    # F1 Score
    f1_graph = sns.boxplot(
        df,
        x=x_col,
        y="F1 Score",
        hue="Method",
        hue_order=hue_order,
        palette=palette,  # type: ignore
        ax=ax[0, 1],
    )
    f1_graph.set_ylim(None, 1)
    f1_graph.set_title("F1 Score (Higher is Better)")

    # Matching Cluster Distance
    mc_graph = sns.boxplot(
        df,
        x=x_col,
        y="Matching Cluster Distance",
        hue="Method",
        hue_order=hue_order,
        palette=palette,  # type: ignore
        ax=ax[1, 0],
    )
    mc_graph.set_ylim(0, None)
    mc_graph.set_title("Matching Cluster Distance (Lower is Better)")

    # CPU Time
    time_graph = sns.boxplot(
        df,
        x=x_col,
        y="CPU Time",
        hue="Method",
        hue_order=hue_order,
        palette=palette,  # type: ignore
        ax=ax[1, 1],
    )
    time_graph.set_yscale("log")
    time_graph.set_title("CPU Time (Lower is Better)")

    min_time = df["CPU Time"].min()
    max_time = df["CPU Time"].max()

    time_ticks = log_time_ticks(min_time, max_time)
    time_graph.minorticks_off()
    time_graph.set_yticks(time_ticks)

    time_graph.set_ylim(time_ticks[0], time_ticks[-1])
    time_graph.yaxis.set_major_formatter(FuncFormatter(format_time_tick))

    # Legend
    handles, labels = time_graph.get_legend_handles_labels()
    rf_graph.legend().remove()
    f1_graph.legend().remove()
    mc_graph.legend().remove()
    time_graph.get_legend().remove()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=len(df["Method"].unique()),
        bbox_to_anchor=(0.5, -0.03),
    )

    fig.tight_layout()

    fig.savefig(
        image_directory
        + f"_{include_code(methods)}_{suptitle.lower().replace(' ','_').replace('(','').replace(')','')}.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)


def log_time_ticks(min_value, max_value, padding=0):
    times = [
        0.01,  # ms
        0.02,
        0.05,
        0.1,
        0.2,
        0.5,
        1,  # s
        2,
        4,
        8,
        15,
        30,
        60,
        2 * 60,  # m
        4 * 60,
        8 * 60,
        15 * 60,
        30 * 60,
        60 * 60,
        2 * 60 * 60,  # h
        4 * 60 * 60,
        8 * 60 * 60,
        16 * 60 * 60,
        24 * 60 * 60,  # d
        36 * 60 * 60,
        48 * 60 * 60,
    ]

    min_index, max_index = None, None
    for i, val in enumerate(times):
        if min_index is None and val > min_value:
            min_index = i - 1
        if max_index is None and val > max_value:
            max_index = i + 1
            break
    assert min_index is not None
    assert max_index is not None

    return times[max(min_index - padding, 0) : max_index + padding]


def format_time_tick(value, pos):
    value = float(value)
    if value < 60:
        return (str(int(value)) if float.is_integer(value) else str(value)) + "s"
    assert value % 60 == 0
    value //= 60
    if value < 60:
        return str(int(value)) + "m"
    assert value % 60 == 0
    value //= 60
    if value < 24:
        return str(int(value)) + "h"
    assert value % 24 == 0 or value % 24 == 12
    value /= 24
    return (str(int(value)) if float.is_integer(value) else str(value)) + "d"


def method_combinations():
    methods = []
    methods.append(
        [METHOD_MAP[SCS_FAST], METHOD_MAP[BCDG], METHOD_MAP[BCDN], METHOD_MAP[MCS]]
    )
    methods.append([METHOD_MAP[SCS_FAST], METHOD_MAP[BCDG], METHOD_MAP[BCDN]])
    methods.append([METHOD_MAP[SCS_FAST], METHOD_MAP[BCDG]])
    methods.append([METHOD_MAP[SCS_FAST], METHOD_MAP[MCS]])
    return methods


def graph_results(image_folder: str, verbosity: int = 0):
    method_combs = method_combinations()

    if verbosity:
        print("GRAPHING", "SMIDGenOG")

    for methods in method_combs:
        graph_smidgenog_combined(image_folder, methods)

    if verbosity:
        print("GRAPHING", "SuperTriplets")
    for methods in method_combs:
        graph_supertriplets_combined(image_folder, methods)

    if verbosity:
        print("GRAPHING", "DCM")
    for methods in method_combs:
        graph_dcm_combined(image_folder, methods)

    if verbosity:
        print("GRAPHING", "IQ")
    for methods in method_combs:
        graph_iq_combined(image_folder, methods)
