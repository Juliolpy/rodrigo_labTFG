# PARTE 2 del core.py para una mayor organización, no tener todo amontonado en un solo script
# Contiene las FUNCIONES para el diseño de primers, ... 

# importamos el módulo subprocess para lanzar programas externos dentro de nuestro script
import subprocess
import re
import json
import pickle
# Para que primer3 pueda trabajar necesita un tipo de archivo especial
def primer3_design_columbo(sequence: str, pam_position: int) -> str:
   """
   Función que utiliza programa primer3  invocado con subprocess para diseñar primers, toma como input: una secuencia y la posicion del pam

   :param sequence: secuencia a pasar a primer3
   :type sequence: str

   :return: diseño de primers Forward y Reverse para las regiones upstream y downstream del protospacer
   :rtype: str
   
   """
   # todo en mayúsculas
   sequence = sequence.upper()
   # limpiamos secuencia, eliminamos cualquier cosa que no sea A, T, C o G
   sequence = re.sub(r'[^ATCG]', '', sequence) 
   # Creamos el contenido con formato válido
   primer3_input = f"""SEQUENCE_ID=test
SEQUENCE_TEMPLATE={sequence} 
SEQUENCE_TARGET={pam_position - 20}, {23}
PRIMER_PRODUCT_SIZE_RANGE=65-165
PRIMER_TASK=generic
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=22
PRIMER_NUM_RETURN=5
=
"""
   
   with open("entry.txt", "w") as f:
      f.write(primer3_input)

   process = subprocess.Popen(
      ["primer3_core", "entry.txt"], # le damos como input el fasta (str)
      stdout= subprocess.PIPE,
      stderr= subprocess.PIPE
   )
   stdout, stderr= process.communicate()
   # si no devuelve nada
   if process.returncode != 0: # si devuelve algo que no sea 0
      raise RuntimeError(f"❌ Error al ejecutar primer3_core:\n{stderr.decode()}")
   # si primer3_core nos devuelve texto plano, la devolvemos como str
   output = stdout.decode()
   # procesamos dicho str y lo convertimos en un diccionario

   return output

# ------ df∈[27,42], dr∈[15,100] -------
   # min = distancia mínima de forward al PAM + tamaño del protospacer + 3 (longitud PAM) + distancia mínima de reverse al PAM
   # max = distancia máxima de forward al PAM + tamaño del protospacer + 3 + distancia máxima de reverse al PAM
   # protospacer+PAM ocupa 23 nt, y quieres df∈[27,42], dr∈[15,100], entonces:
   # PRIMER_PRODUCT_SIZE_RANGE= (27 + 23 + 15, 42 + 23 + 100) = (65, 165)

def score_primers(output: str, pam_pos: int, protospacer_len: int) -> float:
   """
   Función que evalua la bondad de los primers y las posiciones en las que tienen que unirse
   Dado el dict parseado de Primer3 y la posición del PAM, devuelve un score 0 - 1 combinando:

      - Cumplimiento de distancias

      - Contenido GC (~50 %)

      - Tm (proporcionado por Primer3)

   :param output: dict con los primers que saca primer3
   :type output: dict
   :param pam_pos: posicion de las PAM evaluadas por Columbo
   :type pam_pos: int
   :param protospacer_len: longitud del protospacer + PAM (23 nt)
   :type protospacer_len: int

   :return: score de primers
   :rtype: float
   
   output:
   - 'PRIMER_LEFT_0': (start, length)
   - 'PRIMER_RIGHT_0': (start, length)
   - 'PRIMER_LEFT_0_TM', 'PRIMER_RIGHT_0_TM'
   - 'PRIMER_LEFT_0_GC_PERCENT', 'PRIMER_RIGHT_0_GC_PERCENT'
   pam_pos = índice 0-based donde empieza el PAM en la secuencia.
   protospacer_len = longitud de protospacer+PAM (por defecto 23).
   
   """
   # posiciones
   left_start, left_len = output["PRIMER_LEFT_0"]
   right_start, right_len = output["PRIMER_RIGHT_0"]
   # (1) calculamos distacias 
   df = pam_pos - left_start
   dr = (right_start + right_len - 1) - (pam_pos + protospacer_len - 1)
   # (2) puntuación por distancia (1 si esta en el rango y 0 si no)
   s_df = 1.0 if 27 <= df <= 42 else 0.0
   s_dr = 1.0 if 15 <= dr <= 100 else 0.0
   # (3) puntuación por contenido GC (máx ->50 %), la da el programa primer3
   def gc_score(gc_percent):
      return 1.0 - abs(gc_percent - 50)/50 # 1.0 si gc = 50, decae linealmente
   s_gc = ((gc_score(output['PRIMER_LEFT_0_GC_PERCENT'])) + (gc_score(output["PRIMER_RIGHT_0_GC_PERCENT"]))) / 2
   # (4) puntuación por la tm, la da el programa -> ideal 58 - 62 ºC
   def tm_score(tm):
      if 58 <= tm <= 62:
         return 1.0
      return max(0, 1 - min(abs(tm-60)/10, 1)) # 60 al ser la media entre 58 y 62, entre 10 para escalarlo entre 0 y 1 y que no pase de 1. Funciones anidadas
   s_tm = ( tm_score(output['PRIMER_LEFT_0_TM']) + tm_score(output['PRIMER_RIGHT_0_TM']) ) / 2

   # score global, dandole más importancia a la distancia
   score_global = 0.5 * (s_df * s_dr) + 0.25 * s_gc + 0.25 * s_tm
   return score_global # a la funcion score_primers



# funcion parser para gestionar las salidas del programa primer3
def parse_primers_output(output: str) -> dict:
   """
   Función que recoge los datos del archivo generado por primer3 y genera un diccionario rápidamente legible

   :param output: output del programa primer3
   :type output: str

   :return: dict con los primers forward y reverse
   :rtype: dict
   
   """
   results = {}

   for line in output.strip().split("\n"):
        if line.startswith("PRIMER_LEFT_0="):
            coords = line.split("=")[1].strip()
            start, length = map(int, coords.split(","))
            results["PRIMER_LEFT_0"] = (start, length)
        elif line.startswith("PRIMER_RIGHT_0="):
            coords = line.split("=")[1].strip()
            start, length = map(int, coords.split(","))
            results["PRIMER_RIGHT_0"] = (start, length)
        elif line.startswith("PRIMER_LEFT_0_TM="):
            results["PRIMER_LEFT_0_TM"] = float(line.split("=")[1].strip())
        elif line.startswith("PRIMER_RIGHT_0_TM="):
            results["PRIMER_RIGHT_0_TM"] = float(line.split("=")[1].strip())
        elif line.startswith("PRIMER_LEFT_0_GC_PERCENT="):
            results["PRIMER_LEFT_0_GC_PERCENT"] = float(line.split("=")[1].strip())
        elif line.startswith("PRIMER_RIGHT_0_GC_PERCENT="):
            results["PRIMER_RIGHT_0_GC_PERCENT"] = float(line.split("=")[1].strip())

   return results

# pasar los archivos generados del programa a json
def output_pimers_json(output: str) -> dict:  #JSON
      """
      Devuelve un archivo json con los objetos creados por primer3_design

      :param output: output del programa primer3
      :type: str

      :return: Diccionario en formato .json
      :rtype: dict
      
      """
   # guardamos los datos en json en file
      with open("output_PRIMERS.json", "w") as file:
         json.dump(output, file, indent=4)
      # los leemos y cargamos
      with open("output_PRIMERS.json", "r") as file: # los cargamos, no hace falta return 
         json_primers = json.load(file)
      return json_primers

# pasar los archivos generados del programa a pickle
def output_primers_pickle(output: str) -> dict:    # PICKLE
   """
      Devuelve un archivo pickle con los objetos creados por primer3_design
      
      :param output: output del programa primer3
      :type: str

      :return: Diccionario en formato .pickle
      :rtype: dict
      
      """

    # guardamos los datos en pickle en file
   with open("output_PRIMERS.pkl", "wb") as file:
      pickle.dump(output, file)
   # los leemos y cargamos
   with open("output_PRIMERS.pkl", "rb") as file:
      pickle_primers = pickle.load(file)
   return pickle_primers
   