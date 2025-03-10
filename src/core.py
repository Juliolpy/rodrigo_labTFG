# Lean los contenidos de un fichero fasta.
# Localicen las posiciones dentro de la secuencia donde tengamos NGG (recuerda que Nsignifica cualquier base).
# Devuelvan esas posiciones.

# importamos el modulo sys y Seq de Biopython
import sys
from Bio.Seq import Seq
from Bio import SeqIO

fasta_file = sys.argv[1]

def read_fasta(fastafile: str) -> dict:
   # creamos un diccionario donde guardemos las secuencias
   seq_name = {}
   # nombre y secuencia
   for record in SeqIO.parse(fasta_file, "fasta"):
      seq_name[record.id] = str(record.seq)
   return seq_name

def find_NGG_motivs(seq_name: dict) -> dict:
   # diccionario donde guardamos las posiciones de los motivos NGG
   positions = {}

   for seq_id, seq in seq_name.items():
      # esta línea: i + 1 for i in range(len(seq)-2) --> hace que no nos salgamos del indice 
      position_list = [i + 1 for i in range(len(seq)-2) if seq[i + 1 : i + 3] == "GG"] # Si es un motivo NGG en la sencuenia
      positions[seq_id] = position_list
   return positions

# llamamos a las funciones
seq_name = read_fasta(fasta_file)
NGG_positions = find_NGG_motivs(seq_name)
# vemos los resultados

for seq_id, pos in NGG_positions.items():
   print(f"La Secuencia correspondiente con nombre: {seq_id}, tiene {len(pos)} motivos NGG en las posiciones {pos}")
