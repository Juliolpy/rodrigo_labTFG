# PARTE 2 del core.py para una mayor organización, no tener todo amontonado en un solo script
# Contiene las FUNCIONES para el diseño de primers, ... 

# importamos el módulo subprocess para lanzar programas externos dentro de nuestro script
import subprocess
import json
import pickle
from core import read_fasta
# Para que primer3 pueda trabajar necesita un tipo de archivo especial
def primer3_design(fastafile: str) -> dict:
   """
   Función que utiliza programa primer3 invocado con subprocess, toma como input un fasta (str) y saca los diseños de primers (stdout)

   :param fastafile: Archivo fasta como input para primer3
   :type fastafile: str

   :return: diseño de primers Forward y Reverse para PCR
   :rtype: dict
   
   """
   seq_name = read_fasta(fastafile)
   for seq_id, seq in seq_name.items():
      seq_name[seq_id] = seq
   # Creamos el contenido con formato válido
   primer3_input = f"""SEQUENCE_ID=test
SEQUENCE_TEMPLATE={seq}
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
      if line.startswith("PRIMER_LEFT_0_SEQUENCE"):
         results["Primer_Forward"] = line.split("=")[1].strip()
      if line.startswith("PRIMER_RIGHT_0_SEQUENCE"):
         results["Primer_Reverse"] = line.split("=")[1].strip()
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
   