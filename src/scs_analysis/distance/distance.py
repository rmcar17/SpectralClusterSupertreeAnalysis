import time

from typing import Any, List, Set

import numpy as np

from cogent3.core.tree import TreeNode
from scipy.optimize import linear_sum_assignment

from .day_distance import ClusterTable, make_psw, rename_trees


def rooted_rf_distance(tree_1: TreeNode, tree_2: TreeNode) -> int:
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


def matching_cluster_distance(
    tree_1: TreeNode, tree_2: TreeNode, remove_identical=True
) -> int:
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
    if remove_identical:
        intersection = set(map(frozenset, clusters_1)).intersection(
            map(frozenset, clusters_2)
        )
        clusters_1 = list(
            map(set, set(map(frozenset, clusters_1)).difference(intersection))
        )
        clusters_2 = list(
            map(set, set(map(frozenset, clusters_2)).difference(intersection))
        )

    while len(clusters_1) < len(clusters_2):
        clusters_1.append(set())
    while len(clusters_2) < len(clusters_1):
        clusters_2.append(set())

    adjacency = np.zeros(shape=(len(clusters_1), len(clusters_2)))

    for i, cluster_1 in enumerate(clusters_1):
        for j, cluster_2 in enumerate(clusters_2):
            adjacency[i, j] = len(cluster_1.symmetric_difference(cluster_2))

    row_ind, col_ind = linear_sum_assignment(adjacency)
    distance = int(adjacency[row_ind, col_ind].sum())

    return distance


def rooted_f1_distance(tree_1: TreeNode, tree_2: TreeNode) -> float:
    """
    A variation of the F1 accuracy defined in https://watermark.silverchair.com/msx191.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAA3UwggNxBgkqhkiG9w0BBwagggNiMIIDXgIBADCCA1cGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMQUCxmQl5ruQu2990AgEQgIIDKCAJoYQG5_dta6NgjGrJr4l3V1c8cTGkro6X9OUGkxsHS1opCg9NZ-Qx-NSHBr10AZUrnq1p2CnoKp4L1fXsBZXS4_rvv_-UW_xKXPJx2PepzfosFMnf2QIyFauc5MKCD9D8SwoMb6ZBTeX1KgXiHODHjE2L-8VYOzdmgYJALKdDc8xd6Y8xhb2n9gx-Lj_1AFvawAYe_uMktQfA4w5WSQXvOgXSO7g_21uqvTApHNLQk04m31bygsJrj0Po2pR4mFiEcWkyMMHhbkCiSx6bnVWyVddjTNQDKJ7-g-KfHWbQukdsWmDIcEi62_bhgJ3BYp8lLDmDB1lgb59LvH2Cmd2pnG_De-lY2diojqFcWcJ_Mxs_L2zpSUbd6YPaV3Loo4F15-3kPYyVxFMz842orzdvsblGTCuZmGmNQmGFjdtVYYxRd97Y509fw701gy8hHau5W5p2Wg8aCGum2GWVxoaA1uk55qAoNr-M6SEOslKb9-G0OUtLZkLLTTG8fiybf25txG0GWEKWx-ITY-f01SDoRKeiCTE4LIqoLTpDJLg_X_7WpkseWcoqwUzL3ihUjcRQ5Ht2SrRyqAfYtKqUGVi25Hkn5NZUN21mlaByThUVbi0AGy1u43JdYv9LkiT5XGWBJDzT6ZKSUf68VvRhwkdX3DA8iEwM_0rHwsIIkhBQ9_FPrCnHiVXYVVJSFs_NB-v-F956g1qJ1kiGVGQjmAlutvQr2QnHsml7rxnx8rA00xNxOuIxUT_xpccHAlaIZQF9EULCg48u7NvkF6mPuW95hcBJC074t_8a0AaF3zQKra96UYzGGain4A3GmRXWwuwDrZkcr7V81mw89bayk9Rlwl2HgefSsPafURgxNMX8p4itZFMG3pGQsy9M1IvFWzyMh_xIAqil2zZ6Y_jvAhbWP9nRUGdvsNIkNwbDhV2xTYYzZwPthNuFyREI5mgKOecMTv9GShN8OUqBhMO1NasOWrQ-7X1R7V7SKhzF6HMm5dBX_SRmQBTNCSvdzRntRari3Tws_EZB0pbAlMFs6Xo35saNuBe9pb0CmM0Uapk9mo2DT3a3VCY
    taking into account information about the location of the root.

    This is done my measuring f1 accuracy of clades

    Args:
        tree_1 (TreeNode): A tree
        tree_2 (TreeNode): A tree

    Returns:
        float: The f1 score between the clusters of the two trees.
    """
    tree_1_clusters = set(map(frozenset, get_clusters(tree_1)))
    tree_2_clusters = set(map(frozenset, get_clusters(tree_2)))

    # Note that the distance definition is symmetric so
    # which one is the model tree strictly doesn't matter

    # WLOG, assume tree_1 is the model tree

    # True positives are clusters that appear in both the model and estimated tree
    tp = len(tree_1_clusters.intersection(tree_2_clusters))

    # False positives are clusters that appear in the estimated but not the model
    fp = len(tree_2_clusters.difference(tree_1_clusters))

    # False negatives are clusters that appear in the model but not the estimated
    fn = len(tree_1_clusters.difference(tree_2_clusters))

    f1 = 2 * tp / (2 * tp + fp + fn)
    return f1
