import sys
from typing import List
from cogent3.core.tree import PhyloNode
from cogent3 import make_tree
import numpy as np
from min_cut_supertree import min_cut_supertree


def parse_trees(file_path: str) -> List[PhyloNode]:
    trees = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                trees.append(make_tree(line.strip()))
    return trees


if __name__ == "__main__":
    input_trees = parse_trees(sys.argv[-1])

    if "-s" in sys.argv:
        index = sys.argv.index("-s")
        seed = int(sys.argv[index + 1])
    else:
        seed = None

    pcg_weighting = "branch"
    if sys.argv[1].startswith("data/SuperTripletsBenchmark/"):
        pcg_weighting = "depth"
    supertree = min_cut_supertree(
        input_trees,
        pcg_weighting=pcg_weighting,
        contract_edges=False,
        normalise_pcg_weights=False,
        depth_normalisation=False,
    )
    print(supertree)
