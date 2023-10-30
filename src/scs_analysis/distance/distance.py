"""
The Generalized Robinson-Foulds Distance for Phylogenetic Trees

https://www.liebertpub.com/doi/10.1089/cmb.2021.0342
"""

from typing import Any, FrozenSet, Set

import networkx as nx
from cogent3.core.tree import TreeNode

from .day_distance import (
    ClusterTable,
    make_psw,
    rename_trees,
)


def rooted_rf_distance(tree_1: TreeNode, tree_2: TreeNode) -> float:
    tree_1_tips = set(tree_1.get_tip_names())
    tree_2_tips = set(tree_2.get_tip_names())

    assert tree_1_tips == tree_2_tips, (
        str(tree_1_tips.difference(tree_2_tips))
        + " "
        + str(tree_2_tips.difference(tree_1_tips))
    )

    tree_1 = tree_1.deepcopy()
    tree_2 = tree_2.deepcopy()
    inverse = rename_trees([tree_1, tree_2])

    psws = list(map(make_psw, [tree_1, tree_2]))

    cluster_tables = list(map(ClusterTable, psws))

    num_clusters = list(map(lambda x: x.number_of_clusters(), cluster_tables))
    intersection = 0

    cluster_table = cluster_tables[0]
    psw = psws[1]
    S = []
    psw.treset()
    v, w = psw.nvertex()
    while v != -1:
        if w == 0:
            S.append((cluster_table.encode(v), cluster_table.encode(v), 1, 1))
        else:
            L, R, N, W = float("inf"), 0, 0, 1
            while w != 0:
                Ls, Rs, Ns, Ws = S.pop()
                L, R, N, W = min(L, Ls), max(R, Rs), N + Ns, W + Ws
                w = w - Ws
            S.append((L, R, N, W))
            if N == R - L + 1 and cluster_table.is_clust(L, R):
                intersection += 1
        v, w = psw.nvertex()
    return sum(num_clusters) - 2 * intersection


def get_clusters(tree: TreeNode) -> Set[FrozenSet[Any]]:
    clusters = set()
    for node in tree.postorder(include_self=True):
        clusters.add(frozenset(node.get_tip_names()))
    return clusters


def cluster_matching_distance(tree_1: TreeNode, tree_2: TreeNode) -> float:
    tree_1_tips = set(tree_1.get_tip_names())
    tree_2_tips = set(tree_2.get_tip_names())
    assert tree_1_tips == tree_2_tips, (
        str(tree_1_tips.difference(tree_2_tips))
        + " "
        + str(tree_2_tips.difference(tree_1_tips))
    )

    clusters_1 = get_clusters(tree_1)
    clusters_2 = get_clusters(tree_2)

    difference_1 = clusters_1.difference(clusters_2)
    difference_2 = clusters_2.difference(clusters_1)

    assert len(difference_1) == len(difference_2)

    graph = nx.Graph()

    graph.add_nodes_from(difference_1)
    graph.add_nodes_from(difference_2)

    for node_1 in difference_1:
        for node_2 in difference_2:
            graph.add_edge(
                node_1, node_2, weight=len(node_1.symmetric_difference(node_2))
            )

    distance = 0
    matching_edges = nx.min_weight_matching(graph)
    for edge in matching_edges:
        distance += graph.edges[edge]["weight"]

    return distance
