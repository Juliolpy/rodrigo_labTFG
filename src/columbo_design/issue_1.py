from Bio import SeqIO
from Bio.Seq import Seq
import sys

fasta_file = sys.argv[1]

# se lee la secuencia usando SeqIO
record = next(SeqIO.parse(fasta_file, "fasta")) 
# definimos el nombre y la secuencia
name = record.id
seq = record.seq
# para calcular el porcetaje de citosina
nc = seq.count("C")
per = (nc / len(seq))*100

# paso opcional que hago para ver las funciones de Biopython
proteina = seq.translate() # asi de facil xd

# sacar por pantalla los resultados
print(f"Nombre de la secuencia: {name} y longitud: {len(seq)}")
print(f"Secuencia: {seq}")
print(f"Porcentaje de C: {per:.4f}%")
print(f"Proteína teórica traducida: {proteina}")