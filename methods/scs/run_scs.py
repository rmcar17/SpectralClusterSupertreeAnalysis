import sys
from typing import List
from cogent3 import make_tree
from cogent3.core.tree import TreeNode, PhyloNode

from spectral_cluster_supertree import spectral_cluster_supertree


def parse_trees(file_path: str) -> List[PhyloNode]:
    trees = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                trees.append(make_tree(line.strip()))
    return trees


if __name__ == "__main__":
    input_trees = parse_trees(sys.argv[1])
    supertree = spectral_cluster_supertree(
        input_trees,
        pcg_weighting="branch",
        contract_edges=False,
        normalise_pcg_weights=sys.argv[2] == "False",
        depth_normalisation=sys.argv[2] == "True",
    )
    print(supertree)