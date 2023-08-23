from dataclasses import dataclass

@dataclass
class Community:
    chrom: str
    start: int
    end: int

    def __init__(self, chrom, start, end):

        if end <= start:
            raise ValueError("Bin end must be greater than bin start")
        self.chrom = chrom
        self.start = start
        self.end = end
