import argparse
import os
from parsers.parser import Parser
from typing import List

import numpy as np
import scipy.sparse

from algorithms.community_finder import CommunityFinder
from util.write_results_utils import write_community_nested_dict_coordinates


def handle_options():

    parser = argparse.ArgumentParser()
    parser.add_argument('--npz', type=str, nargs='?',
                         help='Path to Hi-C sparse counts matrix .npz.', default='input_files/H1hESC2.5_balanced_chr21_10000.npz')
    parser.add_argument('--chrom', type=str, nargs='?',
                         help='chromosome of matrix')
    parser.add_argument('--res', type=int, nargs='?',
                         help='resolution of heatmap in basepair')


    OPTIONS = parser.parse_args()

    return OPTIONS


def main():

    OPTIONS = handle_options()
    counts_parser = Parser(OPTIONS.npz,OPTIONS.chrom,OPTIONS.res)
    counts_arrays = counts_parser.parse_matrix()

    output_label = str(OPTIONS.npz).split('/')[1][:-4]
    
    gammas = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,
              1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]

    #plot_dir = 'output/plots/'
    #if not os.path.exists(plot_dir):
    #    os.makedirs(plot_dir)

    boundaries_dir = 'output/boundary_files/'
    if not os.path.exists(boundaries_dir):
        os.makedirs(boundaries_dir)

    # Iterate over Hi-C genomic regions
    community_coordinates = {}
    for region in counts_arrays:
        community_coordinates[region] = {}
        community_finder = CommunityFinder(counts_arrays[region],int(region),OPTIONS.chrom,resolution=OPTIONS.res)
        # Iterate over resolutions for partitioning the network
        for gamma in gammas:
            print(
                f'finding communities for region {region} with resolution {gamma}')
            community_coordinates[region][gamma] = community_finder.find_communities(gamma)
            if len(community_coordinates[region][gamma]) == 0:  #empty region, skip to next region
                break

    community_file = f"{boundaries_dir}{output_label}.bed"
    write_community_nested_dict_coordinates(community_coordinates, community_file)

    

if __name__ == "__main__":
    main()
