from Bio.Seq import Seq
import sys
fasta_file = sys.argv[1]

name_seq = []
sequence = []
seq = ""

# Abrir y leer el archivo
with open(fasta_file) as file:
    for line in file:
        line = line.strip("\n")  # Eliminar saltos de línea
        if line.startswith(">"):  # Encabezado de la secuencia
            name = line[1:].split(" ")[0]
            name_seq.append(name)
        else:  # Parte de la secuencia
            seq += line  # Acumulamos la secuencia

# Al final, procesamos la última secuencia
if seq:
    sequence.append(seq)

    nc = seq.count("C")
    per = (nc / len(seq)) * 100
    print(f"Nombre: {name_seq}")
    print(f"Secuencia: {seq}")
    print(f"Longitud: {len(seq)}")
    print(f"Porcentaje de C: {per:.4f}%")

    # traduzco a proteina

    my_dna = Seq(seq)
    my_prot= my_dna.translate()

    print(f"Poteina teórica que se traduciría: {my_prot}")



            