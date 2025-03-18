# Lean los contenidos de un fichero fasta.
# Localicen las posiciones dentro de la secuencia donde tengamos NGG (recuerda que Nsignifica cualquier base).
# Devuelvan esas posiciones.

# importamos el modulo sys y Seq de Biopython
import json
import pickle
from Bio.Seq import Seq
from Bio import SeqIO

def read_fasta(fastafile: str) -> dict:
   # creamos un diccionario donde guardemos las secuencias
   seq_name = {}
   # nombre y secuencia
   for record in SeqIO.parse(fastafile, "fasta"):
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

# funcion para serializar los ouput de las funciones 
def output_NGG_json(positions: dict) -> dict:
# JSON
   # guardamos los datos en json en file
      with open("output_NGG.json", "w") as file:
         json.dump(positions, file, indent=4)
      # los leemos y cargamos
      with open("output_NGG.json", "r") as file:
         json_positions = json.load(file)


def output_NGG_pickle(positions: dict) -> dict:
# pickle
    # guardamos los datos en pickle en file
   with open("output_NGG.pkl", "wb") as file:
      pickle.dump(positions, file)
   # los leemos y cargamos
   with open("output_NGG.pkl", "rb") as file:
      pickle_positions = pickle.load(file)


