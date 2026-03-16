# https://github.com/jreusser/RosalindSolutions

from library import DNASequence, Fasta
from typing import Optional
import os
if __name__ == "__main__":
    # read a file
    fastas: list[Fasta] = []
    current_fasta_name: str = ""
    current_fasta_lines: list[str] = []
    lines = None
    with open(os.path.join(os.path.dirname(__file__), 'rosalind_004.txt')) as f:
        lines = f.readlines()

    for line in lines:
        l = line.strip()
        if l.startswith(">"):
            new_fasta_name = l[1:]
            is_first_loop = current_fasta_name == ""
            if is_first_loop:
                current_fasta_name = new_fasta_name

            if len(current_fasta_lines):
                fastas.append(Fasta(current_fasta_name, DNASequence("".join(current_fasta_lines))))
        
            current_fasta_name = l[1:]
            current_fasta_lines = []
        else:
            current_fasta_lines.append(l)
    
    if len(current_fasta_lines) :
        fastas.append(Fasta(current_fasta_name, DNASequence("".join(current_fasta_lines))))

    # find the largest fasta
    max_fasta = max(fastas, key=lambda f: f.sequence.gc_percentage)
    print(max_fasta.id)
    gc_pct: float = max_fasta.sequence.gc_percentage
    print(f"{gc_pct:.6f}")