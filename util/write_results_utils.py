from typing import List
from typing import Dict

from algorithms.Community import Community


def write_community_coordinates(coordinates: List[Community], fname: str):
    out_file = open(fname, 'w')
    header = "#Chrom, Start, End\n"
    out_file.write(header)
    for coord in coordinates:
        line = f"{coord.chrom}\t{coord.start}\t{coord.end}\n"
        out_file.write(line)

    out_file.close()

def write_community_nested_dict_coordinates(coordinates: Dict[str,Dict[str,List[Community]]], fname: str):
    """coordinates are a nested dictionary of list of communities per chunked region and gamma    
    """

    out_file = open(fname, 'w')
    header = "#Chrom, Start, End, Gamma\n"
    out_file.write(header)
    for region in coordinates:
        for gamma in coordinates[region]:
            for coord in coordinates[region][gamma]:
                line = f"{coord.chrom}\t{coord.start}\t{coord.end}\t{gamma}\n"
                out_file.write(line)

    out_file.close()


