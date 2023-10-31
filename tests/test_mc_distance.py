import pytest
from scs_analysis.distance.distance import matching_cluster_distance
from cogent3 import make_tree


def test_mc_distance():
    a = make_tree("((a,b),(c,d))")
    b = make_tree("((a,c,b),d)")

    assert matching_cluster_distance(a, b) == 3
    assert matching_cluster_distance(b, a) == 3

    a = make_tree("(((a,b,x),c,d),e)")
    b = make_tree("(((a,b),x,c,d),e)")  # should be 1

    assert matching_cluster_distance(a, b) == 1
    assert matching_cluster_distance(b, a) == 1

    # should be n-1
    a = make_tree("((a,b),c)")
    b = make_tree("((a,c),b)")  # should be 2
    assert matching_cluster_distance(a, b) == 2
    assert matching_cluster_distance(b, a) == 2

    a = make_tree("((a,(x,y)),c)")
    b = make_tree("((a,c),(x,y))")  # should be 3
    assert matching_cluster_distance(a, b) == 3
    assert matching_cluster_distance(b, a) == 3

    a = make_tree("((a,(x,y)),(m,(n,o)))")
    b = make_tree("((a,(m,(n,o))),(x,y))")  # should be 5
    assert matching_cluster_distance(a, b) == 5
    assert matching_cluster_distance(b, a) == 5

    a = make_tree("((a,(x,y)),(m,n,o))")
    b = make_tree("((a,(m,n,o)),(x,y))")  # should be 5
    assert matching_cluster_distance(a, b) == 5
    assert matching_cluster_distance(b, a) == 5

    a = make_tree("((((a,d),e),(b,c)),f)")
    b = make_tree("((((e,c),a),b),(d,f))")
    assert matching_cluster_distance(a, b) == 7
    assert matching_cluster_distance(b, a) == 7
