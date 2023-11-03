from typing import Any, Dict, Iterable, List, Set
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


class ShortSubtreeGraph:
    def __init__(self, short_subtrees: Iterable[Set]) -> None:
        self.vertices = set()
        self.edges: Dict[Any, Set] = {}

        for short_subtree in short_subtrees:
            self.add_vertices(short_subtree)
            self.add_edges(short_subtree)

    def add_vertices(self, short_subtree: Set) -> None:
        self.vertices.update(short_subtree)
        for vertex in short_subtree:
            if vertex not in self.edges:
                self.edges[vertex] = set()

    def add_edges(self, short_subtree: Set) -> None:
        for vertex in short_subtree:
            self.edges[vertex].update(short_subtree.difference((vertex,)))

    def get_maximal_cliques(self) -> List[List]:
        all_cliques = []

        chromatic_number = 1
        vertex_list = list(self.vertices)
        vertex_indices = {vertex_list[i]: i for i in range(len(vertex_list))}
        s = {vertex: 0 for vertex in vertex_list}

        for i in range(len(vertex_list)):
            v = vertex_list[i]
            potential_clique = [
                x for x in self.edges[v] if vertex_indices[v] < vertex_indices[x]
            ]

            if len(self.edges[v]) == 0:
                all_cliques.append([v])
            if len(potential_clique) == 0:
                continue

            u = vertex_list[min([vertex_indices[x] for x in potential_clique])]
            s[u] = max(s[u], len(potential_clique) - 1)

            if s[v] < len(potential_clique):
                all_cliques.append(potential_clique + [v])
                chromatic_number = max(chromatic_number, 1 + len(potential_clique))

        return all_cliques

    def compute_components_with_separator(self, separator: Iterable) -> List[Set]:
        vertices_backup = self.vertices.copy()
        self.vertices.difference_update(separator)

        edges_backup = {}
        for vertex, connections in self.edges.items():
            edges_backup[vertex] = connections.copy()

        for vertex in separator:
            del self.edges[vertex]

        for edge in self.edges:
            self.edges[edge].difference_update(separator)

        components = []
        while self.vertices:
            vertex = self.vertices.pop()
            new_component = set([vertex])

            frontier = set(self.edges[vertex])
            while frontier:
                vertex = frontier.pop()
                new_component.add(vertex)
                self.vertices.remove(vertex)

                frontier.update(self.edges[vertex].difference(new_component))
            components.append(new_component)
        self.vertices = vertices_backup
        self.edges = edges_backup
        return components


def centroid_heuristic_separator(tree: PhyloNode) -> Set:
    """
    Find the edge which most separates the taxa, and return
    the short subtree around that edge.

    Args:
        tree (PhyloNode): A tree.

    Returns:
        Set: The short subtree for the edge which most separates the taxa.
    """
    # Cache the number of descendants through a postorder traversal
    for node in tree.postorder():
        if node.is_tip():
            setattr(node, "_tmp_num_descendants", 1)
        else:
            setattr(
                node,
                "_tmp_num_descendants",
                sum([getattr(n, "_tmp_num_descendants") for n in node.children]),
            )
    # Possible to go faster by starting from the root, don't need to iterate through all nodes, just move in
    # direction of most descendants. But this is the formula it is minimising anyway.
    most_balanced = min(
        tree.nontips(include_self=False),
        key=lambda node: abs(
            getattr(
                node, "_tmp_num_descendants"
            )  # The number of desecndants for a node
            - (  # The number of other taxa in the tree
                getattr(tree, "_tmp_num_descendants")
                - getattr(node, "_tmp_num_descendants")
            )
        ),
    )

    for node in tree.postorder():
        delattr(node, "_tmp_num_descendants")

    return compute_short_subtree(most_balanced)


def find_optimal_partition(short_subtree_graph: ShortSubtreeGraph) -> List[Set]:
    best_partition = [short_subtree_graph.vertices.copy()]
    best_partition_score = len(best_partition[0])

    for potential_separator in short_subtree_graph.get_maximal_cliques():
        components = short_subtree_graph.compute_components_with_separator(
            potential_separator
        )
        if len(components) == 1:
            continue

        for component in components:
            component.update(potential_separator)

        partition_score = max(map(len, components))
        if partition_score < best_partition_score:
            best_partition = components
            best_partition_score = partition_score

    return best_partition


def partition_short_subtree_graph(
    guide_tree: PhyloNode,
    short_subtree_graph: ShortSubtreeGraph,
) -> List[Set]:
    centroid_separator = centroid_heuristic_separator(guide_tree)

    ssg_components = short_subtree_graph.compute_components_with_separator(
        centroid_separator
    )
    if len(ssg_components) <= 1:
        return find_optimal_partition(short_subtree_graph)

    for component in ssg_components:
        component.update(centroid_separator)  # Add the separator to all components
    return ssg_components


def split_tree(guide_tree: PhyloNode) -> List[PhyloNode]:
    short_subtrees = compute_short_subtrees(guide_tree)
    short_subtree_graph = ShortSubtreeGraph(short_subtrees)

    partition = partition_short_subtree_graph(guide_tree, short_subtree_graph)

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
