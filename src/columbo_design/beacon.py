# TERCERA PARTE DEL SCRIPT CORE donde se diseña el Beacon con sus scores correspondientes

# Importamos los módulos necesarios
import RNA
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO
from Bio.Seq import Seq
import re

# Función de folding
def fold_beacon(seq: str) -> str:
    """
    Fold a DNA beacon sequence into its RNA secondary structure and compute its minimum free energy (MFE).

    This function converts the input DNA sequence to RNA by replacing Ts with Us, then uses ViennaRNA
    to calculate the secondary structure in dot-bracket notation and the associated MFE.

    :param seq: DNA sequence of the beacon
    :type seq: str
    :return: tuple containing the secondary structure (dot-bracket) and its MFE in kcal/mol
    :rtype: (str, float)
    """
    # Convert DNA to RNA (replace T with U)
    rna = seq.replace("T", "U")
    # Use ViennaRNA to fold the RNA sequence
    struct, mfe = RNA.fold(rna)
    return struct, mfe

# Función de energía libre en la hibridación
def hybridation_energy(beacon: str, target:str) -> float:
    """
    Compute the hybridization free energy between a beacon and its target sequence.

    Both input sequences are first converted to RNA (T -> U) before calling RNA.duplexfold.

    :param beacon: DNA sequence of the beacon
    :type beacon: str
    :param target: DNA sequence of the target (complementary displaced strand)
    :type target: str
    :return: hybridization energy ΔG in kcal/mol (negative value indicates favorable binding)
    :rtype: float
    """
    # hacer RNA la secuencia del beacon
    r_beacon = beacon.replace("T", "U")
    r_target = target.replace("T", "U")

    return RNA.duplexfold(r_beacon, r_target).energy

# reescalado de energía, score de energía libre
def score_energy(dG: float, thr: float = -15.0) -> float:
    """
    Rescale a free energy value into a [0,1] score.

    For ΔG >= 0, returns 0.0. For ΔG <= 2*thr, returns 1.0. Linearly interpolates in between.

    :param dG: free energy change (normally negative)
    :type dG: float
    :param thr: threshold free energy (negative) such that 2*thr maps to score 1.0
    :type thr: float
    :return: normalized energy score in [0.0, 1.0]
    :rtype: float
    """
    # 15 kcal/mol 
    strenght = -dG
    if strenght <= 0:
        return 0.0
    if strenght >= -thr*2:
        return 1.0
    return strenght / (-thr*2)

# vemos que el hairpin tiene efecticamente 11 nt en stem y 9 en el loop
def hairpin_correct(struct: str, ideal_stem=11 , ideal_loop:int=9) -> bool: # true or false
    r"""
    Evaluate how close a folded beacon hairpin matches the ideal stem-loop-stem pattern.

    The function checks the lengths of the initial run of '(' (stem1),
    the central run of '.' (loop), and the final run of ')' (stem2). Each is normalized
    by its ideal length and averaged to produce a score in [0,1].

    :param struct: RNA secondary structure in dot-bracket notation
    :type struct: str
    :param ideal_stem: ideal number of base pairs in each stem
    :type ideal_stem: int
    :param ideal_loop: ideal number of unpaired bases in the loop
    :type ideal_loop: int
    :return: average normalized score of (stem1, loop, stem2)
    :rtype: float
    """
    # patron = r"^\({%d,}\.{%d,}\){%d,}" % (stem_len, loop_len, stem_len)
    # redex desglosado:
    # parentesis ( desde el incio del primer stem_len
    left = len(re.match(r'^\(+', struct).group(0)) if re.match(r'^\(+', struct) else 0
    # puntos del loop_len
    middle = len(re.match(r'^\.+', struct[left:-1:]).group(0)) if re.match(r'^\.+', struct[left:-1:]) else 0
    # parentesis de la derecha, stem_len final
    right = len(re.match(r'\)+', struct[::-1]).group(0)) if re.match(r'\)+', struct[::-1]) else 0
    # ahora normalizamos cada parte
    s_left  = min(left,  ideal_stem) / ideal_stem
    s_mid   = min(middle,   ideal_loop) / ideal_loop
    s_right = min(right, ideal_stem) / ideal_stem

    # y devolvemos la media de los tres (valor en [0,1])
    return (s_left + s_mid + s_right) / 3.0
    
def beacon_tm(tm: float, t_min:float=50.0, t_floor:float=45.0) -> float:
    """
    Activation function for melting temperature: maps Tm to [0,1].

    0 if tm < t_floor, linearly increases to 1 at tm >= t_min.

    :param tm: melting temperature in °C
    :type tm: float
    :param t_min: temperature at which score == 1.0
    :type t_min: float
    :param t_floor: temperature below which score == 0.0
    :type t_floor: float
    :return: normalized Tm score between 0 and 1
    :rtype: float
    """
    if tm < t_floor:
        return 0.0
    if tm >= t_min:
        return 1.0
    return (tm - t_floor)/(t_min - t_floor)

# para ver unicamente la tm
def melting_temperature(seq: str) -> float:
    """
    Calculate the DNA melting temperature (Tm) of a sequence using nearest-neighbor thermodynamics.

    :param seq: nucleotide sequence (DNA or RNA)
    :type seq: str
    :return: melting temperature in °C
    :rtype: float
    """
    seq = seq.upper().replace("U", "T")  # Asegurar secuencia RNA válida
    return float(mt.Tm_NN(seq))

# para encapsular la funcion de la complemetariedad de bases
def ratio_complement(seq1: str, seq2: str) -> float:
    """
    Compute the fraction of complementary base pairs (A-T, G-C) between two aligned sequences.

    :param seq1: first nucleotide sequence
    :type seq1: str
    :param seq2: second nucleotide sequence
    :type seq2: str
    :return: fraction of positions with canonical complementarity
    :rtype: float
    """
    base_pair = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}
    min_len = min(len(seq1), len(seq2))
    matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if base_pair.get(a.upper()) == b.upper())
    return matches / min_len if min_len else 0.0

def score_beacon(beacon_seq: str, beacon_site: str, gRNA_seq: str ,tm_floor:float=50.0) -> float:
    """
    Compute the overall molecular beacon performance score and its components.

    Combines:
      1. hairpin structural correctness
      2. melting temperature activation F_tm
      3. complementarity to displaced strand R_bn
      4. lack of complementarity to gRNA R_gr
      5. hybridization free energy score F_e

    :param beacon_seq: beacon DNA sequence
    :type beacon_seq: str
    :param beacon_site: displaced target site sequence
    :type beacon_site: str
    :param gRNA_seq: guide RNA sequence (protospacer + tracr)
    :type gRNA_seq: str
    :param tm_floor: minimum acceptable Tm for activation
    :type tm_floor: float
    :return: tuple (global_score, R_bn, R_gr, F_tm, F_e, hairpin_score)
    :rtype: (float, float, float, float, float, float)

    """
    # 1. Estructura
    beacon_region = beacon_site[:len(beacon_seq)]
    struct, beacon_mfe = fold_beacon(beacon_seq)
    hp = hairpin_correct(struct)

    # 2. Tm
    tm_val = mt.Tm_NN(beacon_seq)
    F = beacon_tm(tm_val, t_min=tm_floor)

    # 3. Ratio beacon:nicada (simulación simple: fraction of complementarity)
    
    R_bn = ratio_complement(beacon_seq, str(Seq(beacon_region)))

    # 4. Ratio gRNA libre = 1 - fraction complementary to beacon

    R_gr = 1.0 - ratio_complement(beacon_seq, gRNA_seq[:len(beacon_seq)])

    # 5 Energía mínima para la unión beacon con hebra desplazada
    dG = hybridation_energy(beacon_seq, str(Seq(beacon_region).complement()))
    F_e = score_energy(dG, thr=-15.0)
    # componentes de todo el score en una lista
    score_components = [hp, F, R_bn, R_gr, F_e]
    # Score global
    score = 1.0
    for elements in score_components:
        score *= elements if elements > 0 else 1e-9
        score **= (1.0/len(score_components))

    return score, R_bn, R_gr, F, F_e, hp


def design_beacon(beacon_site: str, stem_len:int=11, loop_len:int=9) -> str:
    """
    Design a DNA molecular beacon sequence from a target site.

    Beacon = [stem1 | loop | stem2], where:
      - stem1 = first stem_len nt of the displaced strand reverse
      - loop   = complementary of central region
      - stem2 = inverse complement of stem1

    :param beacon_site: target region sequence (displaced strand)
    :type beacon_site: str
    :param stem_len: length of each stem arm
    :type stem_len: int
    :param loop_len: length of the loop (unused, inferred from site)
    :type loop_len: int
    :return: designed beacon DNA sequence
    :rtype: str
    """
    # 1) hebra desplazada (la que se emparejará con beacon)
    beacon_region = beacon_site[:32] # aquí es donde ajustamos donde queremos que vaya el beacon de nuestra región más grande

    target = str(Seq(beacon_region[::-1]))

    # 2) definimos stem₁
    stem1 = target[:stem_len]

    # 3) stem₂ = complemento inverso de stem₁
    stem2 = str(Seq(stem1).complement()[::-1])

    # 4) todo el trozo intermedio (justo lo que queda entre stems) es el loop
    loop = (target[len(stem1) : len(target)-len(stem1)])
    loop_real = str(Seq(loop).complement()[::-1])
    # 4) stem₂ = complemento inverso de stem₁
    
    return stem1 + loop_real + stem2

# funcion auxiliar para contar --> CHATGPT 
def count_hairpin(struct: str):
    """
    Count stem1, loop, and stem2 lengths in a dot-bracket structure.

    :param struct: RNA secondary structure in dot-bracket
    :type struct: str
    :return: tuple (stem1_length, loop_length, stem2_length)
    :rtype: (int, int, int)
    
    """
    # 1) stem₁ = run inicial de '('
    m1 = re.match(r'^\(+', struct)
    stem1 = len(m1.group(0)) if m1 else 0

    # 2) loop = run de '.' justo después de ese stem₁ (y antes del último ')')
    middle_region = struct[stem1:len(struct)-stem1]
    m2 = re.match(r'^\.+', middle_region)
    loop = len(m2.group(0)) if m2 else 0

    # 3) stem₂ = run final de ')'
    m3 = re.search(r'\)+$', struct)
    stem2 = len(m3.group(0)) if m3 else 0

    return stem1, loop, stem2