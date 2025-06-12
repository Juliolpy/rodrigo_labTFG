# PARTE 2 del core.py para una mayor organización, no tener todo amontonado en un solo script
# Contiene las FUNCIONES para el diseño de primers, ... 

import re
import json
import pickle
from Bio.Seq import Seq
import primer3

# le damos ya los parametros definidos de serie y solo necesitaria la seq y la pam_pos ( un int)
def primer3_design_columbo(seq: str, pam_pos: int) -> dict:
   """
    Design forward and reverse primers around a CRISPR PAM site using Primer3.

    This function takes a DNA sequence and the zero-based index of the PAM
    site, then calls Primer3 to select primers flanking the protospacer region
    (25 nt upstream + PAM).

    :param seq: Full DNA sequence to search for primers.
    :type seq: str
    :param pam_pos: Zero-based index of the first base of the PAM motif in `seq`.
    :type pam_pos: int

    :return: Primer3 raw output dictionary containing primer positions and metrics.
    :rtype: dict

    :raises ValueError: If `seq` contains no valid bases after cleaning.
    :raises RuntimeError: When Primer3 fails to find any primer pairs.
   
   """
   # Clean and prepare the sequence
   seq = seq.upper()
   seq = re.sub(r'[^ATCG]', '', seq)
   # parametros
   protospacer_len = 25   # tu longitud upstream

   # Prepare Primer3 parameters
   primer3_params = {
      'SEQUENCE_ID': 'amp',
      'SEQUENCE_TEMPLATE': seq,
      'PRIMER_PRODUCT_SIZE_RANGE': [[65, 92]],
      'SEQUENCE_INCLUDED_REGION': [pam_pos - 45,  45*2 + 3],
      'PRIMER_PICK_LEFT_PRIMER':   1,
      'PRIMER_PICK_RIGHT_PRIMER':  1,
      'PRIMER_MIN_SIZE':           16,            
      'PRIMER_MAX_SIZE':           25,
      'PRIMER_MIN_TM':             55,
      'PRIMER_OPT_TM':             60,
      'PRIMER_MAX_TM':             65,
      'PRIMER_MIN_GC':             30,
      'PRIMER_MAX_GC':             70,
      'PRIMER_NUM_RETURN':         1,
      'PRIMER_TASK': 'generic',
   }
   global_args = {}

   # Call Primer3 and return the result
   out = primer3.bindings.design_primers(primer3_params, global_args)
   return out

def score_primers(output, pam_relative_pos, protospacer_len=25):
    """
    Compute a composite score for primer pairs based on distance, GC content, and Tm.

    :param output: Primer3 output dictionary from `primer3_design_columbo`.
    :type output: dict
    :param pam_relative_pos: Position of PAM relative to the start of the provided sequence.
    :type pam_relative_pos: int
    :param protospacer_len: Length of the protospacer (default 25 nt).
    :type protospacer_len: int

    :return: Composite score in [0.0, 1.0].
    :rtype: float

    :raises ValueError: If expected Primer3 keys are missing in `output`.

    """
    # Check if primers were found
    if "PRIMER_LEFT_0" not in output or "PRIMER_RIGHT_0" not in output:
        raise ValueError("No primers found by Primer3 for scoring.")

    left_start = output["PRIMER_LEFT_0"][0]
    right_start = output["PRIMER_RIGHT_0"][0]

    # (1) calculamos distancias
    df = pam_relative_pos - left_start
    dr = (right_start + protospacer_len - 1) - (pam_relative_pos + protospacer_len - 1)
    
    # (2a) puntuación por distancia (1 si esta en el rango y 0 si no)
    s_df = 1.0 if 20 <= df <= 45 else 0.0
    s_dr = 1.0 if 14 <= dr <= 100 else 0.0
    # (2b) solapamiento entre el cebador 
    # El protospacer empezará en la región PAM - 25 --> pam_relative_pos - protospacer_len
    prot_start = pam_relative_pos - protospacer_len
    prot_end = pam_relative_pos + 3
    primer_end = left_start + output["PRIMER_LEFT_0"][1]
    # calculo de solapamiento
    overlap_start = max(left_start, prot_start)
    overlap_end = min(primer_end, prot_end)
    overlap_len = max(0, overlap_end - overlap_start)
    # score de solapamiento: 0 si = 0 o > 7, lineal 1..7 -> 1.0
    if overlap_len == 0 or overlap_len > 10:
        s_ov = 0.0
    else:
        s_ov = overlap_len / 10.0

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
    score_global = 0.4 * (s_df * s_dr) + 0.2 * s_gc + 0.2 * s_tm + 0.2 * s_ov
    return score_global

def parse_primers_output(output, region):
    """
    Extract primer sequences and their complements from Primer3 output.

    :param output: Primer3 output dictionary.
    :type output: dict
    :param region: DNA substring used for primer design.
    :type region: str

    :return: Dictionary with left/right primers and complement/reverse-complement.
    :rtype: dict

    :raises ValueError: If no primers are found in `output`.
    """
    # Check if primers were found
    if "PRIMER_LEFT_0" not in output or "PRIMER_RIGHT_0" not in output:
        raise ValueError("No primers found by Primer3 for the given region.")

    start, length = output["PRIMER_LEFT_0"]
    right_start, right_length = output["PRIMER_RIGHT_0"]

    left_primer = region[start:start+length]
    right_primer = region[right_start:right_start+right_length]

    return {
        "left_primer": left_primer,
        "left_primer_comp": str(Seq(left_primer).complement()),
        "right_primer": right_primer,
        "right_primer_revcomp" : str(Seq(right_primer).reverse_complement())
    }

def output_pimers_json(output: dict) -> dict:
   """
    Write Primer3 output to JSON and load it back.

    :param output: Primer3 output dictionary.
    :type output: dict

    :return: Parsed JSON as Python dict.
    :rtype: dict
   """
   with open("output_PRIMERS.json", "w") as file:
      json.dump(output, file, indent=4)
   with open("output_PRIMERS.json", "r") as file:
      json_primers = json.load(file)
   return json_primers

def output_primers_pickle(output: dict) -> dict:
   """
    Serialize Primer3 output to a pickle file and reload it.

    :param output: Primer3 output dictionary.
    :type output: dict

    :return: Unpickled Python dict.
    :rtype: dict
    """
   with open("output_PRIMERS.pkl", "wb") as file:
      pickle.dump(output, file)
   with open("output_PRIMERS.pkl", "rb") as file:
      pickle_primers = pickle.load(file)
   return pickle_primers