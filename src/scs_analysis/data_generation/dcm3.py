from typing import List, Set
from cogent3.core.tree import PhyloNode
import heapq


def ucs_descending(node: PhyloNode) -> Set:
    frontier = []
    return set()


def compute_short_subtree(internal_node: PhyloNode) -> Set:
    # Compute the short subtree for a specific internal node
    taxa = set()

    return taxa


def compute_short_subtrees(tree: PhyloNode) -> List[Set]:
    # For a cogent3 tree, any internal node excluding the root
    # corresponds to an edge (between it and its parent) that
    # can be broken to generate 4 trees (need to be careful at root)
    internal_nodes = tree.nontips(include_self=False)

    short_subtrees = []
    for internal_node in internal_nodes:
        short_subtrees.append(compute_short_subtree(internal_node))

    return short_subtrees


def compute_short_subtree_graph(short_subtrees: List[Set]) -> None:
    pass


def partition_short_subtree_graph(
    guide_tree, short_subtrees, short_subtree_graph
) -> List[Set[str]]:
    return []


def split_tree(guide_tree: PhyloNode) -> List[PhyloNode]:
    short_subtrees = compute_short_subtrees(guide_tree)
    short_subtree_graph = compute_short_subtree_graph(short_subtrees)

    partition = partition_short_subtree_graph(
        guide_tree, short_subtrees, short_subtree_graph
    )

    subtrees = []
    for taxa_group in partition:
        subtrees.append(guide_tree.get_sub_tree(taxa_group))

    return subtrees


def dcm3(
    guide_tree: PhyloNode,
    max_problem_size: int,
) -> List[PhyloNode]:
    subtrees = split_tree(guide_tree)

    if len(subtrees) == 1:
        # DCM3 cannot split the guide tree anymore
        if len(subtrees[0].tips()) > max_problem_size:
            print(
                "Warning... DCM3 cannot split a guide tree to less than the max subproblem size."
            )
        return subtrees

    result = []

    for subtree in subtrees:
        if len(subtree.tips()) > max_problem_size:
            result.extend(dcm3(subtree, max_problem_size))
        else:
            result.append(subtree)

    return result
