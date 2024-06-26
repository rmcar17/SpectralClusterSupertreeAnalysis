from typing import Dict, Literal, Optional, Sequence, Set, Tuple

from cogent3.core.tree import TreeNode

from sc_supertree.scs import (
    _component_to_names_set,
    _connect_trees,
    _contract_proper_cluster_graph,
    _denamify,
    _generate_induced_trees_with_weights,
    _get_all_tip_names,
    _get_graph_components,
    _proper_cluster_graph_edges,
    _tip_names_to_tree,
)


def min_cut_supertree(
    trees: Sequence[TreeNode],
    weights: Optional[Sequence[float]] = None,
    pcg_weighting: Literal["one", "branch", "depth"] = "one",
    contract_edges: bool = True,
) -> TreeNode:
    """
    Min-Cut Supertree (MCS).

    Constructs a supertree from a collection of input trees. The supertree
    method is inspired by Min-Cut Supertree (Semple & Steel, 2000), using
    spectral clustering instead of min-cut to improve efficiency.

    The set of input trees must overlap, the optional weights parameter
    allows the biasing of some trees over others.

    Args:
        trees (Sequence[TreeNode]): Overlapping subtrees.
        weights (Optional[Sequence[float]]): Optional weights for the trees.

    Returns:
        TreeNode: The supertree containing all taxa in the input trees.
    """

    assert len(trees) >= 1, "there must be at least one tree"

    assert pcg_weighting in ["one", "branch", "depth"]

    # Input trees are of equal weight if none is specified
    if weights is None:
        weights = [1.0 for _ in range(len(trees))]

    assert len(trees) == len(weights), "trees and weights must be of same length"

    if len(trees) == 1:  # If there is only one tree left, we can simply graft it on
        _denamify(trees[0])
        return trees[0]

    # The vertices of the proper cluster graph
    # are the names of the tips of all trees
    all_names = _get_all_tip_names(trees)

    # If there are less than or only two names, can instantly return a tree
    if len(all_names) <= 2:
        tree = _tip_names_to_tree(all_names)
        return tree

    pcg_vertices = set((name,) for name in all_names)

    (
        pcg_edges,
        pcg_weights,
        taxa_ocurrences,
        taxa_co_occurrences,
    ) = _proper_cluster_graph_edges(
        pcg_vertices,
        trees,
        weights,
        pcg_weighting,
    )

    components = _get_graph_components(pcg_vertices, pcg_edges)

    if len(components) == 1:
        if contract_edges:
            # Modifies the proper cluster graph inplace
            _contract_proper_cluster_graph(
                pcg_vertices,
                pcg_edges,
                pcg_weights,
                taxa_ocurrences,
                taxa_co_occurrences,
            )
        components = min_cut_partition(pcg_vertices, pcg_weights)

    # The child trees corresponding to the components of the graph
    child_trees = []

    for component in components:
        component = _component_to_names_set(component)
        # Trivial case for if the size of the component is <=2
        # Simply add a tree expressing that
        if len(component) <= 2:
            child_trees.append(_tip_names_to_tree(component))
            continue

        # Otherwise, need to induce the trees on each compoment
        # and recursively call SCS

        # Note, inducing could possible remove trees.
        new_induced_trees, new_weights = _generate_induced_trees_with_weights(
            component, trees, weights
        )

        # Find the supertree for the induced trees
        child_trees.append(
            min_cut_supertree(
                new_induced_trees,
                new_weights,
                pcg_weighting,
                contract_edges,
            )
        )

        missing_tips = component.difference(_get_all_tip_names(new_induced_trees))

        # In this case, treat these tips as individual subtrees
        child_trees.extend(map(lambda x: _tip_names_to_tree((x,)), missing_tips))

    # Connect the child trees by making adjacent to a new root.
    supertree = _connect_trees(child_trees)
    return supertree


def min_cut_partition(vertices: Set, edge_weights: Dict):
    import networkx as nx

    graph = nx.Graph()
    graph.add_nodes_from(vertices)

    weighted_edges = []
    for edge in edge_weights:
        weighted_edges.append((*edge, edge_weights[edge]))
    graph.add_weighted_edges_from(weighted_edges)

    cut_value, partition = nx.stoer_wagner(graph)

    return list(map(set, partition))
