import sys
from typing import List
from cogent3 import make_tree
from cogent3.core.tree import TreeNode, PhyloNode
import numpy as np

from sc_supertree import construct_supertree


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

    contract_edges = "-c" in sys.argv

    pcg_weighting = "branch"
    if sys.argv[1].startswith("data/SuperTripletsBenchmark/"):
        pcg_weighting = "depth"
    supertree = construct_supertree(
        input_trees,
        pcg_weighting=pcg_weighting,
        contract_edges=contract_edges,
        random_state=np.random.RandomState(seed),
    )
    print(supertree)
