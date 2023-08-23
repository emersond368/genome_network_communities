import matplotlib as mpl
mpl.use('Agg')

from lib5c.plotters.extendable import ExtendableHeatmap
from lib5c.parsers.bed import load_features
import lib5c.plotters
import scipy.sparse
import numpy as np
import argparse
import os

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('heatmap_name',type=str,help="pathway + name of heatmap")
    parser.add_argument('domains',type=str,help="pathway + name of domain file")
    parser.add_argument('chr', type=str, help = "chromosome")
    parser.add_argument('start', type=int, help = "start")
    parser.add_argument('end', type=int, help = "end")
    parser.add_argument('resolution',type=int, help = "resolution")
    parser.add_argument('max_color',default = 100,type=float,help="maximum color scale")

    args = parser.parse_args()

    matrix_pre = scipy.sparse.load_npz(args.heatmap_name)
    matrix_pre_csr = matrix_pre.tocsr()
    size = matrix_pre_csr.shape[0]

    domains = lib5c.parsers.bed.load_features(args.domains)
    domains_chrom_list = domains[args.chr]

    #slice section of matrix based on resolution 
    matrix = matrix_pre_csr[int(args.start/args.resolution):int(args.end/args.resolution),int(args.start/args.resolution):int(args.end/args.resolution)].todense()
    some_square_matrix = np.triu(matrix) + np.triu(matrix).T - np.diag(np.diag(matrix))

    #set up heatmap object h
    h = ExtendableHeatmap(array=some_square_matrix,grange_x={'chrom': args.chr, 'start': args.start,'end': args.end},colorscale=(0, args.max_color),colormap='Reds')
 
    #add domain
    h.outline_domains(domains_chrom_list,color='green')

    output_label = args.domains.split('/')[2][:-4]

    plot_dir = 'output/plots/'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    h.save("output/plots/plots_" + output_label + '_' + args.chr  + '_' + str(args.start)  + '_' + str(args.end) + '_' +  str(args.max_color) + '.png')


if __name__ == "__main__":
    main()
     
