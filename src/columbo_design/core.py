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
   >--Excepción personalizada dentro de una clase para indicar que nuestro archivo tipo Seq no tiene el tamaño requerido--<

   :param mensaje: Mensaje de error por tamaño
   :type mensaje: str

   """
   def __init__(self, mensaje = "El PAM o Protospacer no tienen el tamaño mínimo requerido"):
      super().__init__(mensaje) # El método super() en Python se usa para llamar a métodos de la superclase (clase padre) desde una subclase.
      # si no pusieras super, Exception no sabe que existe mensaje, porque nunca lo pasamos a Exception

# vamos a definir la clase ColumboParts   
class ColumboParts:
   """
   Clase donde se almacenan las instancias correspondientes al score y tm

   :param sequence: Objeto Seq a analizar
   :type sequence: Seq
   :param position: Posición de la repetición NGG (sitio PAM)
   :type position: int

   """
   def __init__(self, sequence: Seq, position : int): 
      """
      Función Constructora donde sequence es un objeto seq y la posición un número entero donde se encuentra cada nt
      si la posición + 3 es mayor que la longitud de la secuencia

      :raises TypeError: Error obj Seq &/o int no reconocido

      ..note::
         Este programa solo admite objetos Seq: .fasta

      ..warning::
         No usar floats
      
      """
      if position + 3 > len(sequence):
         # saca por pantalla el error de la clase definida arriba
         raise Out_of_frame_ERROR("El PAM no puede calcularse porque se excede la longitud de la secuencia")
      self._pam = sequence[position : position + 3]
      # si la posicion es menor de 20
      if position < 20: 
         # saca por pantalla el error definido en la clase de arriba
         raise Out_of_frame_ERROR("El protospacer no tiene el tamaño requerido --> 20 nt upstream PAM")
      self._protospacer = sequence[position - 20 : position + 3] 
      self._position = position
      self._scores = self.calcular_scores ()
      self._score_medio = sum(self._scores) / len(self._scores) if self._scores else 0 # basicamente la suma de los números que da el return de la funcion scores entre su lenght que siempre será 4
      self._tm = self.calcular_tm()
   
   # redefinimos todos los atributos usando @property
   
   @property
   def pam(self): # PAM
      """
      Getter para el atributo privado _pam

      :return: Valor del atributo pam
      :rtype: int
      
      """
      return self._pam
   @property
   def position(self): # POSITION
      """
      Getter para el atributo privado _position

      :return: Valor del atributo position
      :rtype: int
      
      """
      return self._position
   @property
   def protospacer(self): # PROTOSPACER
      """
      Getter para el atributo privado _protospacer

      :return: Valor del atributo protospacer
      :rtype: str
      
      """
      return self._protospacer
   
   @property
   def scores(self): #SCORE
      """
      Getter para el atributo privado _scores

      :return: Valor del atributo _scores
      :rtype: int
      
      """
      return self._scores
   
   @property
   def score_medio(self):  #SCORE_MEDIO
      """
      Getter para el atributo privado _score_medio

      :return: Valor del atributo score_medio
      :rtype: int
      
      """
      return self._score_medio
   
   @property
   def tm(self):  # TM
      """
      Getter para el atributo privado _tm

      :return: Valor del atributo tm
      :rtype: str
      
      """
      return self._tm
      
   # Funciones separadas
   def calc_1(self): # 1
      """
      Función 1 -->El contenido de la seed_region de GC ha de ser meno de 80%
      Suma 1 si el contendo de G + C dividido entre la longitud del protospacer es menor a 0.8

      :return: primer score
      :rtype: int

      
      """
      score_1 = 1 if (self._protospacer.count("G") + self._protospacer.count("C")) / len(self._protospacer) < 0.8 else 0
      return score_1 # primer score
   
   def calc_2(self): # 2
      """
      Función 2--> Que no haya una G antes de la PAM NGG
      Suma 1 en el caso de que no haya un G antes de la PAM upstream, si no es 0


      :return: segundo score
      :rtype: int
      
      """
      score_2 = 0 if self._protospacer[-5] == "G" else 1
      return score_2 # segundo score
   
   def calc_3(self): # 3
      """
      Función 3--> Que no haya una T 2 posiciones upstream de la PAM, ej: 5' -----TXNGG---> 3'
      Suma 1 en el caso de que no haya una T 2 posiciones upstream, si no es 0

      :return: tercer score
      :rtype: int
      
      """
      score_3 = 1 if self._protospacer[-7] != "T" else 0
      return score_3 # tercer score
   
   def calc_4(self): # 4
      """
      Función 4--> Que no haya ningun triplete de la misma base, AAA, CCC , TTT , GGG
      Suma 1 en el caso de que no hayan repeticiones de trinucleótidos de la misma base, si no es 0

      :return: cuarto score
      :rtype: int
      
      """
      score_4 = 1 if not any(x * 3 in self._protospacer for x in "ATCG") else 0
      return score_4 # cuarto score
   

   # definimos la función del score_global
   def calcular_scores(self): # GLOBAL
      """
      Función global que suma todas las funciones definidas anteriores
      :return: lista de scores
      :rtype: list
      
      """
      score_list = [ self.calc_1() , self.calc_2(), self.calc_3(), self.calc_4()] # recuerda usar () para ejecutar las funciones
      return score_list
    
   # Función para calcular la temperatura de melting del protospacer
   def calcular_tm(self):
      """
      Devuelve la temperatura de melting del protospacer

      :return: Temperatura de melting del prospacer
      :rtype: str

      """
      melting = mt.Tm_NN(self._protospacer) if self._protospacer else 0 # para calcular la tm del protospcaer
      return melting
   
   # transformar el self en diccionairo para qeu pueda ser transfomrado en json
   def to_dict(self):
     """
    Convierte la instancia de la clase en un diccionario
     
     :return: diccionario con [keys] títulos de cada cosa y [values] los objetos de self definidos en __init__
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
      Serializa la instancia en formato JSON

      :return: JSONI con objetos ColumboParts
      :rtype: str

      
      """
      return json.dumps(self.to_dict()) # json lee los datos y los almacena
   
   @classmethod
   def from_parts(cls, sequence, data):
      """
      Creamos una instancia en la clase ColumboParts a partir del diccionario
      :param sequence: Secuencia completa donde se encuentran el PAM y protospacer.
      :type sequence: Seq
      :param data: Diccionario con los datos necesarios para reconstruir la instancia.
      :type data: dict
      :return: Instancia de ColumboParts
      :rtype: ColumboParts

      """
      return cls(sequence, data["Position"])
   @classmethod
   def from_json(cls, sequence, json_str): # cls en vez de self
      """
      Creamos una instancia en la clase ColumboParts a partir un str JSON

      :param sequence: Secuencia completa donde se encuentran el PAM y protospacer.
      :type sequence: Seq
      :param json_str: JSON con los datos de la instancia.
      :type json_str: str
      :return: Instancia de ColumboParts
      :rtype: ColumboParts

      """
      data = json.loads(json_str)
      return cls.from_parts(sequence, data)

# estas funciones behind the musgo, fuera de la clase
def read_fasta(fastafile: str) -> dict:
   """
   Función que se encarga de leer el archivo con extensión .fasta que se introduce

   :param fastafile: Archivo fasta a analizar
   :type fastafile: str

   :return: String con los nombres y su secuencia correspondiente
   :rtype: str
   
   """
   
   seq_name = {}
   # nombre y secuencia
   for record in SeqIO.parse(fastafile, "fasta"):
      clean_seq = re.sub(r'[^ATCG]', '', str(record.seq).upper()) # limpiamos secuencia
      seq_name[record.id] = clean_seq
   return seq_name

def find_NGG_motivs(seq_name: dict) -> dict:
   """
   Función que se encarga de encontrar las PAM (motivos NGG) dentro de la sencuencia del fasta

   :param seq_name: Diccionario que contiene nombres y sus secuencias correspondientes
   :type seq_name: dict

   :return: Diccionario con las posiciones de los motivos NGG para cada secuencia
   :rtype: dict 
   
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
   Función que procesa un fasta, identifica los sitios NGG y crea  objetos ColumboParts 
   para cada uno, los evalúa y devuelve los 20 mejores candidatos ordenados por su score_medio
   
   :param seq_name: Diccionario con el nombre [keys] y su correspondiente secuencia de fasta [values]
   :type seq_name: dict
   :param motifs: Diccionario con las posiciones de cada motivo NGG

   :return: Lista creada tras evaluar los diccionarios con las instancias de ColumboParts
   :rtype: list
   
   """
   columbo_list = []
   for seq_id, positions in motifs.items():
      sequence = seq_name[seq_id]
      for pos in positions:
         # ya que el protospacer no puede ser menor de 20, nuestro codigo no funcionaria: index out of range
         if pos >= 20 and pos + 3 <= len(sequence):
            columbo_list.append(ColumboParts(sequence, pos)) # objeto seq y un int
         else:
            print(f"{RED} ⚠️ [WARNING]⚠️ {RESET} posición {YELL}{pos}{RESET} descartada en {CIAN}{seq_id}{RESET}:{RED}{Out_of_frame_ERROR}{RESET}")
   # ahora para quedarnos con los 20 mejores solo
   columbo_list.sort(key= lambda x : x.score_medio, reverse= True)
   return columbo_list[:20]
 
def output_NGG_json(motifs: dict) -> dict:  #JSON
      """
      Devuelve un archivo json con los objetos creados por ColumboParts

      :param motifs: Diccionario con las posiciones de los motivos NGG
      :type: dict

      :return: Diccionario en formato .json
      :rtype: dict
      
      """
   # guardamos los datos en json en file
      with open("output_NGG.json", "w") as file:
         json.dump([obj.to_dict() for obj in motifs], file, indent=4)
      # los leemos y cargamos
      with open("output_NGG.json", "r") as file: # los cargamos, no hace falta return 
         json_positions = json.load(file)


def output_NGG_pickle(motifs: dict) -> dict:    # PICKLE
   """
      Devuelve un archivo pickle con los objetos creados por ColumboParts
      
      :param motifs: Diccionario con las posiciones de los motivos NGG
      :type: dict

      :return: Diccionario en formato .pickle
      :rtype: dict
      
      """

    # guardamos los datos en pickle en file
   with open("output_NGG.pkl", "wb") as file:
      pickle.dump([obj.to_dict() for obj in motifs], file)
   # los leemos y cargamos
   with open("output_NGG.pkl", "rb") as file:
      pickle_positions = pickle.load(file)
