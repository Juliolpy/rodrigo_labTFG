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
import re
# colores para que uno se pueda aclarar
RED = "\033[91m"
CIAN = "\033[96m"
GREEN = "\033[92m"
GR = "\033[92m"
RESET = "\033[0m"
MAG = "\033[35m"
YELL = "\033[33m"

# crea una CLASE para hacer un RAISE para crear excepciones
class Out_of_frame_ERROR(Exception):
   """
   Exception raised when a PAM or protospacer falls outside the bounds of the sequence.

   :param message: Error message to display
   :type message: str

   """
   def __init__(self, mensaje = "El PAM o Protospacer no tienen el tamaño mínimo requerido"):
      super().__init__(mensaje) # El método super() en Python se usa para llamar a métodos de la superclase (clase padre) desde una subclase.
      # si no pusieras super, Exception no sabe que existe mensaje, porque nunca lo pasamos a Exception

# vamos a definir la clase ColumboParts   
class ColumboParts:
   """
   Represents a CRISPR-Cas9 protospacer candidate (Columbo part). Stores sequence, scores, and melting temperature.

   :param sequence: Full DNA sequence as a Biopython Seq object
   :type sequence: Seq
   :param position: Zero-based index of the PAM start (NGG) within `sequence`
   :type position: int

   """
   def __init__(self, sequence: Seq, position : int): 
      """
      Initialize a ColumboParts instance. Validates PAM and protospacer frame and computes scores.

      :raises OutOfFrameError: If PAM or protospacer exceed sequence boundaries or invalid PAM
      :raises TypeError: If incorrect types are passed
      
      """
      if position + 3 > len(sequence):
         # saca por pantalla el error de la clase definida arriba
         raise Out_of_frame_ERROR("El PAM no puede calcularse porque se excede la longitud de la secuencia")
      pam = sequence[position : position + 3]
      if not pam[1] == "G" and pam[2] == "G":
         raise Out_of_frame_ERROR(f"PAM inválido en la posicion {position}:{pam} no tiene estructura NGG ") # descartar PAM inválidos
      self._pam = pam
      # si la posicion es menor de 20
      if position < 25: 
         # saca por pantalla el error definido en la clase de arriba
         raise Out_of_frame_ERROR("El protospacer no tiene el tamaño requerido --> 20 nt upstream PAM")
      self._protospacer = sequence[position - 25 : position + 3]
      self._beacon_site = sequence[position - 50: position + 3]
      self._region = sequence[position - 100: position + 100]
      self._position = position
      self._scores = self.calcular_scores ()
      self._score_medio = sum(self._scores) / len(self._scores) if self._scores else 0 # basicamente la suma de los números que da el return de la funcion scores entre su lenght que siempre será 4
      self._tm = self.calcular_tm()
   
   # redefinimos todos los atributos usando @property
   @property
   def beacon_site(self) -> Seq: # SITIO DE UNIÓN BEACON
      """
      Sequence region used for beacon design (50 nt upstream + PAM).

      :return: Beacon site sequence
      :rtype: Seq
      
      """
      return self._beacon_site
   
   @property
   def region(self) -> Seq: # REGION
      """
      Flanking region around protospacer (200 nt upstream to 300 nt downstream).

      :return: Sequence region
      :rtype: Seq
      
      """
      return self._region
   
   @property
   def pam(self) -> Seq: # PAM
      """
      PAM sequence (3 nt NGG).

      :return: PAM sequence
      :rtype: Seq
      
      """
      return self._pam
   
   @property
   def position(self) -> int: # POSITION
      """
      Position of PAM within the original sequence.

      :return: Zero-based index of PAM
      :rtype: int
      
      """
      return self._position
   @property
   def protospacer(self) -> Seq: # PROTOSPACER
      """
      Protospacer sequence (25 nt + PAM).

      :return: Protospacer
      :rtype: Seq
      
      """
      return self._protospacer
   
   @property
   def scores(self) -> float: #SCORE
      """
      Mean score across all features.

      :return: Average score
      :rtype: float
      
      """
      return self._scores
   
   @property
   def score_medio(self) -> list[int]:  #SCORE_MEDIO
      """
      List of individual feature scores [gc_content, no_G_before_pam, no_T_2_upstream, no_homopolymer].

      :return: Individual scores
      :rtype: list[int]
      
      """
      return self._score_medio
   
   @property
   def tm(self) -> float:  # TM
      """
      Melting temperature of the protospacer (in °C).

      :return: Melting temperature
      :rtype: float

      """
      return self._tm
      
   # Funciones separadas
   def calc_1(self) -> int: # 1
      """
      Score 1: 1 if GC content of protospacer < 80%, else 0.

      :return: 1 or 0
      :rtype: int

      """
      score_1 = 1 if (self._protospacer.count("G") + self._protospacer.count("C")) / len(self._protospacer) < 0.8 else 0
      return score_1 # primer score
   
   def calc_2(self) -> int: # 2
      """
      Score 2: 1 if the base immediately upstream of PAM is not 'G', else 0.

      :return: 1 or 0
      :rtype: int
      
      """
      score_2 = 0 if self._protospacer[-5] == "G" else 1
      return score_2 # segundo score
   
   def calc_3(self) -> int: # 3
      """
      Score 3: 1 if the base two positions upstream of PAM is not 'T', else 0.

      :return: 1 or 0
      :rtype: int
      
      """
      score_3 = 1 if self._protospacer[-7] != "T" else 0
      return score_3 # tercer score
   
   def calc_4(self) -> int: # 4
      """
      Score 4: 1 if no homopolymer triplet (AAA, CCC, GGG, TTT), else 0.

      :return: 1 or 0
      :rtype: int
      
      """
      score_4 = 1 if not any(x * 3 in self._protospacer for x in "ATCG") else 0
      return score_4 # cuarto score
   

   # definimos la función del score_global
   def calcular_scores(self) ->list[int]: # GLOBAL
      """
      Compute all four individual scores.

      :return: List of four scores
      :rtype: list[int]
      
      """
      score_list = [ self.calc_1() , self.calc_2(), self.calc_3(), self.calc_4()] # recuerda usar () para ejecutar las funciones
      return score_list
    
   # Función para calcular la temperatura de melting del protospacer
   def calcular_tm(self) -> float:
      """
      Compute melting temperature using nearest-neighbor method.

      :return: Melting temperature in °C
      :rtype: float

      """
      melting = mt.Tm_NN(self._protospacer) if self._protospacer else 0 # para calcular la tm del protospcaer
      return melting
   
   # transformar el self en diccionairo para qeu pueda ser transfomrado en json
   def to_dict(self) -> dict:
     """
     Serialize instance to a dictionary.

     :return: Dict with keys ['PAM','Protospacer','Position','Scores','Score_medio','Tm']
     :rtype: dict

     """
     return {
        "PAM" : self.pam,
        "Protospacer" : self.protospacer,     # añadido el underscore
        "Position" : self.position,
        "Scores" : self.scores,
        "Score_medio" : self.score_medio,
        "Tm" : self.tm
     }
   def to_json(self):
      """
      Serialize instance to JSON string.

      :return: JSON representation
      :rtype: str
 
      """
      return json.dumps(self.to_dict()) # json lee los datos y los almacena
   
   @classmethod
   def from_parts(cls, sequence, data):
      """
      Reconstruct ColumboParts from a dict of properties.

      :param sequence: Original full sequence
      :type sequence: Seq
      :param data: Dict with keys matching to_dict()
      :type data: dict
      :return: New ColumboParts instance
      :rtype: ColumboParts

      """
      return cls(sequence, data["Position"])
   @classmethod
   def from_json(cls, sequence, json_str): # cls en vez de self
      """
      Reconstruct ColumboParts from a JSON string.

      :param sequence: Original full sequence
      :type sequence: Seq
      :param json_str: JSON string as produced by to_json()
      :type json_str: str
      :return: New ColumboParts instance
      :rtype: ColumboParts

      """
      data = json.loads(json_str)
      return cls.from_parts(sequence, data)

# estas funciones behind the musgo, fuera de la clase
def read_fasta(fastafile: str) -> dict:
   """
   Read a FASTA file into a dict of clean DNA sequences.

   :param fastafile: Path to FASTA file
   :type fastafile: str
   :return: Dict mapping record IDs to uppercase ATCG-only strings
   :rtype: dict[str, str]
   
   """
   
   seq_name = {}
   # nombre y secuencia
   for record in SeqIO.parse(fastafile, "fasta"):
      clean_seq = re.sub(r'[^ATCG]', '', str(record.seq).upper()) # limpiamos secuencia
      seq_name[record.id] = clean_seq
   return seq_name

def find_NGG_motivs(seq_name: dict) -> dict:
   """
   Find all NGG PAM motifs in each sequence.

   :param sequences: Dict of ID -> sequence
   :type sequences: dict
   :return: Dict of ID -> list of zero-based PAM start positions
   :rtype: dict[str, list[int]]
   
   """
   # diccionario donde guardamos las posiciones de los motivos NGG
   motifs = {} # seq_nombre = número de motivos

   for seq_id, seq in seq_name.items():
      # esta línea: i + 1 for i in range(len(seq)-2) --> hace que no nos salgamos del indice 
      motifs_list = [i + 1 for i in range(len(seq)-2) if seq[i + 1 : i + 3] == "GG"] # Si es un motivo NGG en la sencuenia
      motifs[seq_id] = motifs_list
   return motifs  

def process_genome(seq_name: dict, motifs: dict) -> list:
   """
   Build ColumboParts for each valid PAM and return top 20 by average score.

   :param sequences: Dict of ID -> sequence
   :param motifs: Dict of ID -> list of PAM positions
   :return: Top 20 ColumboParts instances, sorted by score_medio descending
   :rtype: list[ColumboParts]
   
   """
   columbo_list = []
   for seq_id, positions in motifs.items():
      sequence = seq_name[seq_id]
      for pos in positions:
         # ya que el protospacer no puede ser menor de 20, nuestro codigo no funcionaria: index out of range
         if pos >= 25 and pos + 3 <= len(sequence):
            columbo_list.append(ColumboParts(sequence, pos)) # objeto seq y un int
         else:
            print(f"{RED} ⚠️ [WARNING]⚠️ {RESET} posición {YELL}{pos}{RESET} descartada en {CIAN}{seq_id}{RESET}:{RED}{Out_of_frame_ERROR}{RESET}")
   # ahora para quedarnos con los 20 mejores solo
   columbo_list.sort(key= lambda x : x.score_medio, reverse= True)
   return columbo_list[:20]
 
def output_NGG_json(motifs: dict) -> dict:  #JSON
      """
      Write top ColumboParts to output_NGG.json and return list of dicts.

      :param parts: List of ColumboParts
      :type parts: list
      :return: List of serialized dicts
      :rtype: list[dict]
      
      """
   # guardamos los datos en json en file
      with open("output_NGG.json", "w") as file:
         json.dump([obj.to_dict() for obj in motifs], file, indent=4)
      # los leemos y cargamos
      with open("output_NGG.json", "r") as file: # los cargamos, no hace falta return 
         json_positions = json.load(file)


def output_NGG_pickle(motifs: dict) -> dict:    # PICKLE
   """
   Write top ColumboParts to output_NGG.pkl and return list of dicts.

   :param parts: List of ColumboParts
   :type parts: list
   :return: List of serialized dicts
   :rtype: list[dict]
      
   """

    # guardamos los datos en pickle en file
   with open("output_NGG.pkl", "wb") as file:
      pickle.dump([obj.to_dict() for obj in motifs], file)
   # los leemos y cargamos
   with open("output_NGG.pkl", "rb") as file:
      pickle_positions = pickle.load(file)
