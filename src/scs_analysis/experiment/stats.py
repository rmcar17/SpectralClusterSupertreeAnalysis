from .graph import load_data

ALLOWED = ("BCD (GSCM)", "SCS")


def iq_stats(col, taxa, m):
    print("\t\t" + str(taxa).ljust(5), "&", str(m).ljust(3), "& ", end="")

    file = f"results/birth_death/{taxa}/iq_source_trees/{m}/"

    df = load_data(file)

    values = []
    index = 0
    for method, sub_df in df.groupby("Method"):
        if method not in ALLOWED:
            continue
        assert method == ALLOWED[index]
        index += 1

        values.append(sub_df[col].median())
        # print(method.ljust(10), round(sub_df[col].median(), 2))
    for i, value in enumerate(values):
        if values[i] < values[1 - i]:
            print(f"\\textbf{{{round(value, 2)}}}", "& ", end="")
        else:
            print(round(value, 2), "& ", end="")
    print(str(round(values[0] / values[1], 2)) + "$\\times$ \\\\")


def smidgenog_stats(col, taxa, density):
    print(
        "\t\t" + str(taxa).ljust(5),
        "&",
        str(density if density != 0 else "-").ljust(3),
        "& ",
        end="",
    )

    file = f"results/SMIDGenOutgrouped/{taxa}/{density}/"

    df = load_data(file)

    values = []
    index = 0
    for method, sub_df in df.groupby("Method"):
        if method not in ALLOWED:
            continue
        assert method == ALLOWED[index]
        index += 1

        values.append(sub_df[col].median())
        # print(method.ljust(10), round(sub_df[col].median(), 2))
    for i, value in enumerate(values):
        if values[i] < values[1 - i]:
            print(f"\\textbf{{{round(value, 2)}}}", "& ", end="")
        else:
            print(round(value, 2), "& ", end="")
    print(str(round(values[0] / values[1], 2)) + "$\\times$ \\\\")


def supertriplets_stats(col, d, k):
    print(
        "\t\t" + str(d).ljust(2),
        "&",
        str(k).ljust(2),
        "& ",
        end="",
    )

    file = f"results/SuperTripletsBenchmark/d{d}/k{k}/"

    df = load_data(file)

    values = []
    index = 0
    for method, sub_df in df.groupby("Method"):
        if method not in ALLOWED:
            continue
        assert method == ALLOWED[index]
        index += 1

        values.append(sub_df[col].median())
        # print(method.ljust(10), round(sub_df[col].median(), 2))
    for i, value in enumerate(values):
        if values[i] < values[1 - i]:
            print(f"\\textbf{{{round(value, 2)}}}", "& ", end="")
        else:
            print(round(value, 2), "& ", end="")
    print(str(round(values[0] / values[1], 2)) + "$\\times$ \\\\")


def table_prelude():
    print(
        """\\begin{table}[]
\t\\centering
\t\\begin{tabular}{c|c|c|c|c}"""
    )


def table_prologue(caption, label):
    print(
        """\t\\end{tabular}
\t\\caption{Median """
        + caption
        + """}
\t\\label{tab:"""
        + label
        + """}
\\end{table}
"""
    )


LABELS = {
    "iq": "scs-dcm-iq",
    "smidgenog": "smidgenog",
    "supertriplets": "super_triplets",
}


def table(func):
    def table_wrapper(*args, **kwargs):
        table_prelude()
        func(*args, **kwargs)
        table_prologue(args[0], LABELS[func.__name__.split("_")[0]])

    return table_wrapper


@table
def iq_stats_from_results(col):
    print(
        """\t\tTaxa  & $m$ & BCD (GSCM) (seconds) & SCS (seconds) & Improvement \\\\
\t\t\\hline"""
    )
    for taxa in (500, 1000, 2000, 5000, 10000):
        for m in (50, 100):
            iq_stats(col, taxa, m)


@table
def smidgenog_from_results(col):
    print(
        """\t\tTaxa  & Density & BCD (GSCM) (seconds) & SCS (seconds) & Improvement \\\\
\t\t\\hline"""
    )
    for taxa in (100, 500, 1000, 10000):
        for density in (20, 50, 75, 100) if taxa != 10000 else (0,):
            smidgenog_stats(col, taxa, density)


@table
def supertriplets_from_results(col):
    print(
        """\t\t$d$ & $k$ & BCD (GSCM) (seconds) & SCS (seconds) & Improvement \\\\
\t\t\\hline"""
    )
    for taxa in (25, 50, 75):
        for density in (10, 20, 30, 40, 50):
            supertriplets_stats(col, taxa, density)


def stats_from_results():
    iq_stats_from_results("CPU Time")
    iq_stats_from_results("Matching Cluster Distance")

    smidgenog_from_results("CPU Time")
    smidgenog_from_results("Matching Cluster Distance")

    supertriplets_from_results("CPU Time")
    supertriplets_from_results("Matching Cluster Distance")
