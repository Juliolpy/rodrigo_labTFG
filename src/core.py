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
      self.protospacer = sequence[position - 20 : position + 3] if position >= 20 else "" # es la seed_region también
      self.position = position
      self.scores = self.calcular_scores ()
      self.score_medio = sum(self.scores) / len(self.scores) if self.scores else 0 # basicamente la suma de los números que da el return de la funcion scores entre su lenght que siempre será 4
      self.tm = self.calcular_tm()
   """
    >--- Funciones separadas / Consejo del Luksgrin ---<
    seed_region = self.protospacer-->  el propio protospacer no coge el PAM
    En sistemas CRISPR-Cas9, la región semilla es la parte más crítica para la hibridación entre la guía de ARN y el ADN diana
    Esta región está situada cerca del PAM (NGG) y juega un papel clave en la especificidad y estabilidad de la unión -> es donde se calculan los scores
    los scores funcionan 1 para verdad y 0 para falso
   """
   def calc_1(self): # 1
      """
      Función 1
      -->El contenido de la seed_region de GC ha de ser meno de 80%
      
      """
      score_1 = 1 if (self.protospacer.count("G") + self.protospacer.count("C")) / len(self.protospacer) < 0.8 else 0
      return score_1 # primer score
   
   def calc_2(self): # 2
      """
      Función 2
      --> Que no haya una G antes de la PAM NGG
      
      """
      score_2 = 0 if self.protospacer[-5] == "G" else 1
      return score_2 # segundo score
   
   def calc_3(self): # 3
      """
      Función 3
      --> Que no haya una T 2 posiciones upstream de la PAM, ej: 5' -----TXNGG---> 3'
      
      """
      score_3 = 1 if self.protospacer[-7] != "T" else 0
      return score_3 # tercer score
   
   def calc_4(self): # 4
      """
      Función 4
      --> Que no haya ningun triplete de la misma base, AAA, CCC , TTT , GGG
      
      """
      score_4 = 1 if not any(x * 3 in self.protospacer for x in "ATCG") else 0
      return score_4 # cuarto score
   

   # definimos la función del score_global
   def calcular_scores(self): # GLOBAL
      """
      Función global que suma todas las funciones definidas anteriores
      Args : seed_region, el propio protospacer, no hace falta que definamos un valor de seed_region
      return : lista de valores, score
      
      """
      score_list = [ self.calc_1() , self.calc_2(), self.calc_3(), self.calc_4()] # recuerda usar () para ejecutar las funciones
      return score_list
    
   # Función para calcular la temperatura de melting del protospacer
   def calcular_tm(self):
      """
      Args : objeto seq del constructor, el protospacer
      Return : int correspondiente a la Tm_NN
      """
      melting = mt.Tm_NN(self.protospacer) if self.protospacer else 0 # para calcular la tm del protospcaer
      return melting
   
   # transformar el self en diccionairo para qeu pueda ser transfomrado en json
   def to_dict(self):
     """
     Args : objetos seq del constructor
     Return : diccionario con [keys] títulos de cada cosa y [values] los objetos de self definidos en __init__
     
     """
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