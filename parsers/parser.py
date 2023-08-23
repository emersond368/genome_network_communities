
from .Bin import Bin
from os import path
import numpy as np
import scipy.sparse

class Parser:
    """ Used for parsing input counts and bed files

    A bed file defines the genomic coordinates of a list of evenly-spaced genomic bins.
    The bed file should be in the following format:
        chrom\tstart\tend\tgenomic_bin
    where genomic_bin is in the following format:
        region_BIN_number
    For example: HiCchr1800037_BIN_0003


    A counts file stores a list of interaction scores between two different genomic loci as defined by their genomic bins. 
    The counts file should be in the following format:
        forward_bin\treverse_bin\tinteraction_score
    where forward_bin and reverse_bin are in the same format as the bed file genomic_bin

    
    """
    def __init__(self, matrix_fname: str,counts_fname = "None", bed_fname =  "None"):
        if not path.exists(matrix_fname):
            raise ValueError(f"matrix file {bed_fname} does not exist")
        self.matrix = matrix_fname
        self.counts = counts_fname
        self.bed = bed_fname


    def parse_bed(self):
        """Parses bed file into a bin_map, which defines the genomic coordinates of each bin per region.
        
        Returns:
            bin_map : dictionary where key is region and value is a list of genomic Bins sorted on their start coordinate
        """

        bin_map = {}
        with open(self.bed, 'r') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                bed_bin, region = self._parse_bed_line(line)
                bin_map.setdefault(region, [])
                bin_map.get(region).append(bed_bin)

        for region in bin_map:
            bin_map.get(region).sort(key=lambda x: x.start)

        return bin_map


    def parse_counts(self, bin_map=None):
        """Parses counts file into an adjacency interaction matrix per genomic region
        
        Keyword Arguments:
            bin_map {dict} -- dictionary where key is region and value is a list of genomic Bins sorted on their start coordinate (default: {None})
        
        Raises:
            ValueError: Region names and region sizes must match between counts and bed fiels
        
        Returns:
            counts_as_arrays -- dictionary where key is region and value is a 2D numpy array adjaceny matrix of interaction scores
        """

        if not bin_map:
            bin_map = self.parse_bed()
        counts_as_arrays = {}
        for region in bin_map:
            region_size = len(bin_map.get(region))
            counts_as_arrays.setdefault(region, np.zeros((region_size, region_size), dtype=np.float64))

        with open(self.counts, 'r') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                region, forward_bin_number, \
                    reverse_bin_number, interaction_score = self._parse_counts_line(line)
                if region not in bin_map:
                    raise ValueError("Regions do not match between counts and bed files")
                region_size = len(bin_map.get(region))
                if forward_bin_number > region_size - 1 or reverse_bin_number > region_size -1:
                    raise ValueError("Region sizes do not match between counts and bed files")
                
                counts_as_arrays[region][forward_bin_number, reverse_bin_number] = interaction_score
                counts_as_arrays[region][reverse_bin_number, forward_bin_number] = interaction_score

        return counts_as_arrays

    def parse_matrix(self,resolution=10000):
        """Parses the chromosome  2D numpy array adjaceny matrix of interaction scores into overlapping chunked regions

        Keyword Argument:
            resolution -- integer for basepair size of matrix

        Returns:
            counts_as_arrays  -- dictionary where key is region and value is a 2D numpy array adjaceny matrix of interaction scores
        """

        counts_matrix = scipy.sparse.load_npz(self.matrix).tocsr()
        counts_matrix = counts_matrix.todense() 
       
        counts_matrix_regions = {}   

        scale = int(10000/resolution)       

        start = 0
        size = 600*scale
        overlap = 400*scale
        end = start + size
        while end < counts_matrix.shape[0]:
            counts_matrix_regions[start] = counts_matrix[start:end,start:end]
            start = start + overlap
            end = start + size

        return counts_matrix_regions



    def _parse_bed_line(self, bed_line):
        pieces = bed_line.strip().split('\t')
        if len(pieces) != 4:
            raise ValueError("bed file not formatted as expected")
        chrom = pieces[0]
        try:
            start = int(pieces[1])
        except ValueError:
            raise ValueError("start bin cannot be parsed to an integer")
        try:
            end = int(pieces[2])
        except ValueError:
            raise ValueError("end bin cannot be parsed to an integer")
        name = pieces[3]
        bin_name_pieces = name.split('_')
        if len(bin_name_pieces) != 3:
            raise ValueError("bed file not formatted as expected")
        region = bin_name_pieces[0]
        if bin_name_pieces[1]!= 'BIN':
            raise ValueError('bed file region name not formatted as expected')
        bed_bin = Bin(chrom, name, start, end)

        return bed_bin, region


    def _parse_counts_line(self, counts_line):

        pieces = counts_line.strip().split('\t')
        if len(pieces) != 3:
            raise ValueError("counts file not formatted as expected")
        region = pieces[0].split('_')[0]
        try:
            forward_bin = pieces[0].split('_')[-1]
            forward_bin_number = int(forward_bin)
        except ValueError:
            raise ValueError(f"Could not parse bin number from {forward_bin}")
        try:
            reverse_bin = pieces[1].split('_')[-1]
            reverse_bin_number = int(reverse_bin)
        except ValueError:
            raise ValueError(f"Could not parse bin number from {reverse_bin}")
        try:
            interaction_score = float(pieces[2])
        except ValueError:
            raise ValueError(f"Could not cast interaction score {pieces[2]} to float")

        return region, forward_bin_number, reverse_bin_number, interaction_score



        

        




        



        

    

