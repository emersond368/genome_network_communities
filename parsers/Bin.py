from dataclasses import dataclass

@dataclass
class Bin:
    chrom: str
    name: str
    start: int
    end: int

    def __init__(self, chrom, name, start, end):

        if end <= start:
            raise ValueError("Bin end must be greater than bin start")
        self.chrom = chrom
        self.name = name
        self.start = start
        self.end = end
