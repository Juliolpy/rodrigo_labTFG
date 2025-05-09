# PARTE 2 del core.py para una mayor organización, no tener todo amontonado en un solo script
# Contiene las FUNCIONES para el diseño de primers, ... 

import re
import json
import pickle
from Bio.Seq import Seq
import primer3

def primer3_design_columbo(sequence: str, pam_position: int) -> dict:
   """
   Función que utiliza programa primer3  invocado con subprocess para diseñar primers, toma como input: una secuencia y la posicion del pam

   :param sequence: secuencia a pasar a primer3
   :type sequence: str

   :return: diseño de primers Forward y Reverse para las regiones upstream y downstream del protospacer
   :rtype: dict
   
   """

   # Clean and prepare the sequence
   sequence = sequence.upper()
   sequence = re.sub(r'[^ATCG]', '', sequence)

   # Prepare Primer3 parameters
   primer3_params = {
      'SEQUENCE_ID': 'test',
      'SEQUENCE_TEMPLATE': sequence,
      'SEQUENCE_TARGET': [pam_position - 20, 23],
      'PRIMER_PRODUCT_SIZE_RANGE': [[65, 165]],
      'PRIMER_TASK': 'generic',
      'PRIMER_MIN_SIZE': 18,
      'PRIMER_MAX_SIZE': 22,
      'PRIMER_NUM_RETURN': 5
   }
   global_args = {}

   # Call Primer3 and return the result
   result = primer3.bindings.design_primers(primer3_params, global_args)
   return result

def score_primers(output, pam_relative_pos, protospacer_len=23):
    # Check if primers were found
    if "PRIMER_LEFT_0" not in output or "PRIMER_RIGHT_0" not in output:
        raise ValueError("No primers found by Primer3 for scoring.")

    left_start = output["PRIMER_LEFT_0"][0]
    right_start = output["PRIMER_RIGHT_0"][0]

    # (1) calculamos distancias
    df = pam_relative_pos - left_start
    dr = (right_start + protospacer_len - 1) - (pam_relative_pos + protospacer_len - 1)
    
    # (2) puntuación por distancia (1 si esta en el rango y 0 si no)
    s_df = 1.0 if 27 <= df <= 42 else 0.0
    s_dr = 1.0 if 15 <= dr <= 100 else 0.0
    
    # (3) puntuación por contenido GC (máx ->50 %), la da primer3
    def gc_score(gc_percent):
        return 1.0 - abs(gc_percent - 50)/50 # 1.0 si gc = 50, decae linealmente
    s_gc = (gc_score(output['PRIMER_LEFT_0_GC_PERCENT']) + gc_score(output["PRIMER_RIGHT_0_GC_PERCENT"])) / 2
    
    # (4) puntuación por la tm, la da primer3 -> ideal 58 - 62 ºC
    def tm_score(tm):
        if 58 <= tm <= 62:
            return 1.0
        return max(0, 1 - min(abs(tm-60)/10, 1))
    s_tm = (tm_score(output['PRIMER_LEFT_0_TM']) + tm_score(output['PRIMER_RIGHT_0_TM'])) / 2

    # score global, dandole más importancia a la distancia
    score_global = 0.5 * (s_df * s_dr) + 0.25 * s_gc + 0.25 * s_tm
    return score_global

def parse_primers_output(output, region):
    # Check if primers were found
    if "PRIMER_LEFT_0" not in output or "PRIMER_RIGHT_0" not in output:
        raise ValueError("No primers found by Primer3 for the given region.")

    start, length = output["PRIMER_LEFT_0"]
    right_start, right_length = output["PRIMER_RIGHT_0"]

    left_primer = region[start:start+length]
    right_primer = region[right_start:right_start+right_length]

    return {
        "left_primer": left_primer,
        "right_primer": right_primer
    }

def output_pimers_json(output: dict) -> dict:
   """
   Devuelve un archivo json con los objetos creados por primer3_design

   :param output: output del programa primer3
   :type output: dict

   :return: Diccionario en formato .json
   :rtype: dict
   """
   with open("output_PRIMERS.json", "w") as file:
      json.dump(output, file, indent=4)
   with open("output_PRIMERS.json", "r") as file:
      json_primers = json.load(file)
   return json_primers

def output_primers_pickle(output: dict) -> dict:
   """
   Devuelve un archivo pickle con los objetos creados por primer3_design
   
   :param output: output del programa primer3
   :type output: dict

   :return: Diccionario en formato .pickle
   :rtype: dict
   """
   with open("output_PRIMERS.pkl", "wb") as file:
      pickle.dump(output, file)
   with open("output_PRIMERS.pkl", "rb") as file:
      pickle_primers = pickle.load(file)
   return pickle_primers