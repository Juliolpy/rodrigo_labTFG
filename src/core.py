# Lean los contenidos de un fichero fasta.
# Localicen las posiciones dentro de la secuencia donde tengamos NGG (recuerda que Nsignifica cualquier base).
# Devuelvan esas posiciones.

# para calcular la melting temperature
from Bio.SeqUtils import MeltingTemp as mt # mt es la función de la herramienta meltingtemp dentro del modulo sequtils de BIO
# importamos el modulo sys y Seq de Biopython
import json
import pickle
from Bio.Seq import Seq
from Bio import SeqIO

# vamos a definir la clase ColumboParts   
class ColumboParts:
   # definimos constructor --> debe aceptar un string o un objeto Bio.Seq.Seq
   def __init__(self, sequence: Seq, position : int): # sequence es un objeto seq y la positicion un numero entero donde se encuentra cada nt
      self.PAM = sequence[position : position + 3] if position + 3 <= len(sequence) else ""
      self.protospacer = sequence[position - 20 : position + 3] if position >= 20 else ""
      self.position = position
      self.scores = self.calcular_scores ()
      self.score_medio = sum(self.scores) / len(self.scores) if self.scores else 0 # basicamente la suma de los números que da el return de la funcion scores entre su lenght que siempre será 4
      self.tm = self.calcular_tm()

      # definimos la función del score
   def calcular_scores(self):
      # ponemos la seed_region
      seed_region = self.protospacer[-12:-4] # cogemos los nucleotidos y nos aseguramos que no coja el PAM
      # En sistemas CRISPR-Cas9, la región semilla es la parte más crítica para la hibridación entre la guía de ARN y el ADN diana
      # Esta región está situada cerca del PAM (NGG) y juega un papel clave en la especificidad y estabilidad de la unión -> es donde se calculan los scores
      # los scores funcionan 1 para verdad y 0 para falso y hacemos una lista con ellos para tenerlos encapsulados
      scores = [
         1 if (seed_region.count("G") and seed_region.count("C")) / len(seed_region) < 0.8 else 0, # el contenido de la seed_region de GC ha de ser meno de 80%
         1 if self.protospacer[-5] == "G" else 0, # que no haya una G antes de la PAM NGG
         1 if self.protospacer[-7] != "T" else 0, # que ni haya una T 2 posiciones upstream de la PAM 5' -----TXNGG---> 3'
         1 if not any(x * 3 in seed_region for x in "ATCG") else 0 # que no haya ningun triplete de la misma base, AAA, CCC , TTT , GGG   
      ]
      return scores
      
   def calcular_tm(self):
      melting = mt.Tm_NN(self.protospacer) if self.protospacer else 0 # para calcular la tm del protospcaer
      return melting
      # transformar el self en diccionairo para qeu pueda ser transfomrado en json
   def to_dict(self):
     return {
        "PAM" : self.PAM,
        "Protospacer" : self.protospacer,
        "Position" : self.position,
        "Scores" : self.scores,
        "Score_medio" : self.score_medio,
        "Tm" : self.tm
     }

# estas funciones behind the musgo, fuera de la clase
def read_fasta(fastafile: str) -> dict:
   # creamos un diccionario donde guardemos las secuencias
   seq_name = {}
   # nombre y secuencia
   for record in SeqIO.parse(fastafile, "fasta"):
      seq_name[record.id] = str(record.seq)
   return seq_name

def find_NGG_motivs(seq_name: dict) -> dict:
   # diccionario donde guardamos las posiciones de los motivos NGG
   motifs = {}

   for seq_id, seq in seq_name.items():
      # esta línea: i + 1 for i in range(len(seq)-2) --> hace que no nos salgamos del indice 
      motifs_list = [i + 1 for i in range(len(seq)-2) if seq[i + 1 : i + 3] == "GG"] # Si es un motivo NGG en la sencuenia
      motifs[seq_id] = motifs_list
   return motifs  

# funcion para procesar el genoma y convertirlo en columboparts
#  procesa el genoma identificando los sitios donde aparece el motivo NGG y
#  crea objetos ColumboParts para cada uno, los evalúa y devuelve los 20 mejores candidatos ordenados por su score_medio.
def process_genome(seq_name: dict, motifs: dict) -> list:
   columbo_list = []
   for seq_id, positions in motifs.items():
      for pos in positions:
         if pos >= 20: # ya que el protospacer no puede ser menor de 20, nuestro codigo no funcionaria: index out of range
            columbo_list.append(ColumboParts(seq_name[seq_id], pos)) # objeto seq y un int
            # ahora para quedarnos con los 20 mejores solo
            columbo_list.sort(key= lambda x : x.score_medio, reverse= True)
   return columbo_list[:20]

# funcion para serializar los ouput de las funciones 
def output_NGG_json(motifs: dict) -> dict:
# JSON
   # guardamos los datos en json en file
      with open("output_NGG.json", "w") as file:
         json.dump([obj.to_dict() for obj in motifs], file, indent=4)
      # los leemos y cargamos
      with open("output_NGG.json", "r") as file: # los cargamos, no hace falta return 
         json_positions = json.load(file)


def output_NGG_pickle(motifs: dict) -> dict:
# pickle
    # guardamos los datos en pickle en file
   with open("output_NGG.pkl", "wb") as file:
      pickle.dump([obj.to_dict() for obj in motifs], file)
   # los leemos y cargamos
   with open("output_NGG.pkl", "rb") as file:
      pickle_positions = pickle.load(file)