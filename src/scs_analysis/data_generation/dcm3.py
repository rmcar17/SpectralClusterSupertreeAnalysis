from typing import List, Set
from cogent3.core.tree import PhyloNode
import heapq
from dataclasses import dataclass, field


@dataclass(order=True)
class DistanceNode:
    distance: int
    node: PhyloNode = field(compare=False)


def ucs_descending(initial_node: PhyloNode) -> Set:
    frontier = []
    closest_leaves = set()
    shortest_distance_to_tip = float("inf")

    heapq.heappush(frontier, DistanceNode(0, initial_node))
    while len(frontier) > 0:
        distance_node: DistanceNode = heapq.heappop(frontier)

        # Used to include ties for closest tips
        if distance_node.distance > shortest_distance_to_tip:
            break

        if distance_node.node.is_tip():
            closest_leaves.add(distance_node.node.name)
            shortest_distance_to_tip = distance_node.distance
        else:
            for child in distance_node.node.children:
                heapq.heappush(
                    frontier, DistanceNode(distance_node.distance + child.length, child)
                )
    return closest_leaves


@dataclass(order=True)
class DistanceNodeBanned:
    distance: int
    node: PhyloNode = field(compare=False)
    banned: PhyloNode = field(compare=False)


def ucs_ascending(initial_node: PhyloNode) -> Set:
    frontier = []
    closest_leaves = set()
    shortest_distance_to_tip = float("inf")

    heapq.heappush(
        frontier, DistanceNodeBanned(0, initial_node.parent.parent, initial_node.parent)
    )
    while len(frontier) > 0:
        distance_node_banned: DistanceNodeBanned = heapq.heappop(frontier)

        # Used to include ties for closest tips
        if distance_node_banned.distance > shortest_distance_to_tip:
            break

        if distance_node_banned.node.is_tip():
            closest_leaves.add(distance_node_banned.node.name)
            shortest_distance_to_tip = distance_node_banned.distance
        else:
            for neighbour in distance_node_banned.node._getNeighboursExcept(
                distance_node_banned.banned
            ):
                if neighbour is distance_node_banned.node.parent:
                    # The nieghbour is the node's parent need to add the node's length
                    heapq.heappush(
                        frontier,
                        DistanceNodeBanned(
                            distance_node_banned.distance
                            + distance_node_banned.node.length,
                            neighbour,
                            distance_node_banned.node,
                        ),
                    )
                else:
                    # The neighbour is a child of the ndoe. Add the child's length, banned is the parent node
                    heapq.heappush(
                        frontier,
                        DistanceNodeBanned(
                            distance_node_banned.distance + neighbour.length,
                            neighbour,
                            distance_node_banned.node,
                        ),
                    )

    return closest_leaves


def compute_short_subtree(internal_node: PhyloNode) -> Set:
    # Compute the short subtree for a specific internal node
    assert internal_node.parent is not None  # The internal node is not the root.
    # TODO I can enforce bifurcating, though can also extend it to handle more general.
    short_subtree = set()

    for child in internal_node.children:
        short_subtree.update(ucs_descending(child))

    if internal_node.parent.is_root():
        # The two edges that go across the root are treated as if they were a single edge (treat the tree as unrooted)
        # Need to go through the sibling's children
        for sibling in internal_node.siblings():
            for sibling_child in sibling.children:
                short_subtree.update(ucs_descending(sibling_child))
    else:
        # We must handle
        # 1. The internal node's sibling
        # 2. UCS going up the tree starting from the internal node's parent

        # Case 1
        for sibling in internal_node.siblings():
            short_subtree.update(ucs_descending(sibling))
        # Case 2
        short_subtree.update(ucs_ascending(internal_node))

    return short_subtree


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