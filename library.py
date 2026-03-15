from typing import NewType

# Static-only aliases (mypy catches misuse, zero runtime cost)
CAIScore = NewType("CAIScore", float)
GenbankAccession = NewType("GenbankAccession", str)

# Runtime-validated types (bad data raises immediately)
class DNASequence(str):
    """A validated, uppercase ACGT-only DNA string."""
    VALID = frozenset("ACGT")
    def __new__(cls, value: str) -> "DNASequence":
        v = value.strip().upper()
        if invalid := frozenset(v) - cls.VALID:
            raise ValueError(f"Invalid nucleotides: {invalid!r}")
        return super().__new__(cls, v)

    @property
    def gc_percentage(self) -> float:
        return (self.count("G") + self.count("C")) / len(self) * 100
    
    @property
    def at_percentage(self) -> float:
        return (self.count("A") + self.count("T")) / len(self) * 100

    def reverse_complement(self) -> "DNASequence":
        table = str.maketrans("ACGT", "TGCA")
        return DNASequence(self.translate(table)[::-1])


class AminoAcidSequence(str):
    """A validated single-letter IUPAC amino acid string."""
    VALID = frozenset("ARNDCEQGHILKMFPSTWYV")

    def __new__(cls, value: str) -> "AminoAcidSequence":
        v = value.strip().upper()
        if invalid := frozenset(v) - cls.VALID:
            raise ValueError(f"Invalid residues: {invalid!r}")
        return super().__new__(cls, v)