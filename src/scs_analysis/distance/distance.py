from typing import Any, List, Set
import networkx as nx
from cogent3.core.tree import TreeNode
from .day_distance import ClusterTable, make_psw, rename_trees


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


def get_clusters(tree: TreeNode) -> List[Set[Any]]:
    clusters = []
    for node in tree.postorder(include_self=True):
        clusters.append(set(node.get_tip_names()))
    return clusters


def matching_cluster_distance(tree_1: TreeNode, tree_2: TreeNode) -> int:
    """Matching cluster distance between rooted trees.

    Source:
    https://sciendo.com/article/10.2478/amcs-2013-0050

    Args:
        tree_1 (TreeNode): A rooted tree
        tree_2 (TreeNode): A rooted tree

    Returns:
        int: The matching cluster distance between the trees.
    """
    tree_1_tips = set(tree_1.get_tip_names())
    tree_2_tips = set(tree_2.get_tip_names())
    assert tree_1_tips == tree_2_tips, (
        str(tree_1_tips.difference(tree_2_tips))
        + " "
        + str(tree_2_tips.difference(tree_1_tips))
    )

    clusters_1 = get_clusters(tree_1)
    clusters_2 = get_clusters(tree_2)

    while len(clusters_1) < len(clusters_2):
        clusters_1.append(set())
    while len(clusters_2) < len(clusters_1):
        clusters_2.append(set())

    cluster_num = 0
    cluster_ids_1 = {}
    for cluster in clusters_1:
        cluster_ids_1[cluster_num] = cluster
        cluster_num += 1
    cluster_ids_2 = {}
    for cluster in clusters_2:
        cluster_ids_2[cluster_num] = cluster
        cluster_num += 1

    graph = nx.Graph()

    graph.add_nodes_from(cluster_ids_1.keys())
    graph.add_nodes_from(cluster_ids_2.keys())

    for node_1, cluster_1 in cluster_ids_1.items():
        for node_2, cluster_2 in cluster_ids_2.items():
            graph.add_edge(
                node_1, node_2, weight=len(cluster_1.symmetric_difference(cluster_2))
            )

    distance = 0
    matching_edges = nx.min_weight_matching(graph)
    for edge in matching_edges:
        distance += graph.edges[edge]["weight"]

    return distance
