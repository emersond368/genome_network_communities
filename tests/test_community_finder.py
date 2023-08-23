

from numpy import array

from Bin import Bin
from Community import Community
from community_finder import CommunityFinder

# Inputs used for all tests
counts_array = array([
    [3, 3, 0, 0, 0, 0],
    [3, 3, 0, 0, 0, 0],
    [0, 0, 4, 4, 4, 0],
    [0, 0, 4, 4, 4, 0],
    [0, 0, 4, 4, 4, 0],
    [0, 0, 0, 0, 0, 1]
])

regional_bin_map = [
    Bin('chr1', 'BIN0', 0, 100),
    Bin('chr1', 'BIN1', 100, 200),
    Bin('chr1', 'BIN2', 200, 300),
    Bin('chr1', 'BIN3', 300, 400),
    Bin('chr1', 'BIN4', 400, 500),
    Bin('chr1', 'BIN5', 500, 600),
]

community_finder = CommunityFinder(
    counts_array, regional_bin_map, number_partitions=4)


def test_partition_network():

    gamma = 1.0
    partitions = community_finder._partition_network(gamma)

    expected_partitions = array([
        [1, 1, 2, 2, 2, 3],
        [1, 1, 2, 2, 2, 3],
        [1, 1, 2, 2, 2, 3],
        [1, 1, 2, 2, 2, 3],
    ])

    comparison = partitions == expected_partitions
    assert comparison.all()


def test_get_similarity_consensus():

    partitions = array([
        [1, 1, 1, 2, 2, 3],
        [1, 1, 2, 2, 2, 3],
        [1, 1, 2, 2, 2, 2],
        [1, 1, 2, 2, 2, 3],
    ])

    expected_consensus = array(
        [1, 1, 2, 2, 2, 3]
    )

    consensus = community_finder._get_similarity_consensus(partitions)
    comparison = consensus == expected_consensus
    assert comparison.all()


def test_get_community_genomic_coordinates():

    consensus = array(
        [1, 1, 2, 2, 2, 3]
    )

    expected_communities = [
        Community('chr1', 0, 200),
        Community('chr1', 200, 500),
        Community('chr1', 500, 600)
    ]

    communities = community_finder._get_community_genomic_coordinates(
        consensus)

    assert len(communities) == len(expected_communities)

    for i in range(len(communities)):
        assert communities[i] == expected_communities[i]
