import sys
from Bio.Seq import Seq

seq = sys.argv[1]

long = len(seq)

for nt in seq:
    A = (seq.count("A")/long)*100
    T = (seq.count("T")/long)*100
    C = (seq.count("C")/long)*100
    G = (seq.count("G")/long)*100

print(f"El porcentage de cada A, T , C, G {A:.2f}% {T:.2f}% {C:.2f}% {G:.2f} %")