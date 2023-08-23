from ..parser import Parser


def test_parse_correctly_formatted():
    bed_fname = './example_files/example.bed'
    counts_fname = './example_files/example.counts'
    parser = Parser(counts_fname, bed_fname)
    counts_array = parser.parse_counts()
    expected_region_size = 151
    for region in counts_array:
        region_size = len(counts_array.get(region))
        assert region_size == expected_region_size


#TODO: parse counts incorrect format






