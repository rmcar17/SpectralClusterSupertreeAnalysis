import sys
from typing import List
from cogent3.core.tree import PhyloNode
from cogent3 import make_tree
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
    input_trees = parse_trees(sys.argv[1])
    supertree = min_cut_supertree(input_trees)
    print(supertree)
