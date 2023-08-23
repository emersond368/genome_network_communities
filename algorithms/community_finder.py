from typing import List

import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score

from .Community import Community
from .greedy_louvain import greedy_louvain
#from .leiden import leiden

class CommunityFinder:

    def __init__(self, A: np.ndarray, region_start: int, chrom:str,number_partitions: int = 10,resolution: int = 10000):
        self.A = A
        self.region_start = region_start
        self.number_partitions = number_partitions
        self.resolution = resolution
        self.chrom = chrom
        self.network_size = len(A)

        # Newman Girvan null model
        # null_model_ij = ki*kj/(2*m) where ki is the degree of node i, kj is the degree of node j, and m is the total sum of edge weights for the network
        self.null_model = np.outer(np.sum(A, 1), np.sum(A, 1))/np.sum(A)

    def _partition_network(self, gamma: float,search_func) -> np.ndarray:
        modularity_matrix = (self.A - gamma * self.null_model)/np.sum(self.A)
        # each row represents a community partition
        partitions = np.zeros((self.number_partitions, self.network_size))
        for i in range(self.number_partitions):
            #community_partition, _modularity_score = greedy_louvain(
            community_partition, _modularity_score = search_func(
                modularity_matrix)
            # sort community partition such that community numbers are increasing
            value, loc, counts = np.unique(
                community_partition, return_index=True, return_counts=True)
            counts_ordered = counts[np.argsort(loc)]
            sorted_community_partition = np.repeat(
                range(1, len(value)+1), counts_ordered)
            partitions[i, :] = sorted_community_partition

        return partitions

    def _get_similarity_consensus(self, partitions) -> np.ndarray:
        pairwise_similarity = np.zeros(
            (self.number_partitions, self.number_partitions))

        # Calculate pairwise similarities between each partition
        for i in range(self.number_partitions):
            for j in range(self.number_partitions):
                pairwise_similarity[i, j] = adjusted_rand_score(
                    partitions[i, :], partitions[j, :])

        # make matrix symmetric
        pairwise_similarity = (pairwise_similarity +
                               np.transpose(pairwise_similarity))/2

        # compute average pairwise similarity
        avg_pairwise_similarity = np.average(pairwise_similarity, axis=0)

        # find the partition that has the maximal average similarity
        best_index = np.argmax(avg_pairwise_similarity)

        consensus = partitions[best_index, :]

        return consensus

    def _get_community_genomic_coordinates(self, consensus, include_region_edges=True) -> List[Community]:

        coordinates = []
        boundaries = np.where(consensus[1:] != consensus[:-1])[0] + 1
        #print("consensus: ",consensus)
        #print("boundaries: ",boundaries)

        # Include the first and last coordinate of the region as boundaries
        if include_region_edges:
            boundaries = np.insert(boundaries, 0, 0)
            boundaries = np.insert(boundaries, len(boundaries), len(consensus))

        for i in range(len(boundaries) - 1):
            # community_start to community_stop includes the bin numbers contained within the community
            start_bin = boundaries[i]
            stop_bin = boundaries[i+1] - 1
            start = (self.region_start + start_bin)*self.resolution
            end = (self.region_start + stop_bin)*self.resolution
            chrom = self.chrom
            if end >  start:
                community = Community(chrom, start, end)
                coordinates.append(community)

        return coordinates

    def find_communities(self, gamma: float,search_func = greedy_louvain) -> List[Community]:
        partitions = self._partition_network(gamma,search_func)
        consensus = self._get_similarity_consensus(partitions)
        community_coordinates = self._get_community_genomic_coordinates(
            consensus)

        return community_coordinates
