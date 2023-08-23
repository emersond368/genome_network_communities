from ..Bin import Bin
import pytest


def test_bin_init():
    chrom = 'chr19'
    start = 1010
    end = 2010
    name = 'test'
    _test_bin = Bin(chrom, name, start, end)


def test_test_bin_end_less_than_start():
    chrom = 'chr19'
    start = 10010
    end = 2010
    name = 'test'
    with pytest.raises(ValueError, match="Bin end must be greater than bin start"):
        _test_bin = Bin(chrom, name, start, end)
