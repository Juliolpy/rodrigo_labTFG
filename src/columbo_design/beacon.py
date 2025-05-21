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
    Función que recoge un str de ADN y lo convierte a ARN para poder usar ViennaRNA sobre él.

    :param seq: secuencia de ADN
    :type seq: str

    :return: devuelve la estructura secundaria en notación de paréntesis (dot-bracket)
    :rtype: str
    """
    # Convert DNA to RNA (replace T with U)
    rna = seq.replace("T", "U")
    # Use ViennaRNA to fold the RNA sequence
    struct, mfe = RNA.fold(rna)
    return struct, mfe

# Función de energía libre en la hibridación
def hybridation_energy(beacon: str, target:str) -> float:
    """
    Utilizamos RNA.duplexfold
    Tendremos que convertir primero el beacon a RNA para calcular dicha energía de interacción

    duplexfold devuelve un objeto .energy !! -> el target tiene que ser la non - target strand, la hebra desplazada, la reverse complement de mi protospacer

    """
    # hacer RNA la secuencia del beacon
    r_beacon = beacon.replace("T", "U")
    r_target = target.replace("T", "U")

    return RNA.duplexfold(r_beacon, r_target).energy

# reescalado de energía, score de energía libre
def score_energy(dG: float, thr: float = -15.0) -> float:
    """
    Convierte un ΔG (normalmente negativo) en un valor entre 0 y 1,
    de modo que para ΔG >= 0 → 0, y para ΔG <= 2*thr → 1, lineal en medio.
    """
    # 15 kcal/mol 
    strenght = -dG
    if strenght <= 0:
        return 0.0
    if strenght >= -thr*2:
        return 1.0
    return strenght / (-thr*2)

# vemos que el hairpin tiene efecticamente 11 nt en stem y 9 en el loop
def hairpin_correct(struct: str, stem_len=11 , loop_len:int=9) -> bool: # true or false
    r"""
    Función que comprueba en dot-bracket que hay un stem de al menos stem_len pares,
    un loop de al menos loop_len puntos, y luego el stem inverso.
    Buscamos el patrón: ^\( {stem,}\.{loop,}\){stem,}

    :param struct: estructura generada por la función fold_beacon
    :type struct: str
    :param stem_len: longitud que tiene que tener la secuencia stem, 8 nt
    :type stem_len: int
    :param loop_len: longitud que debe tener la secuencia loop, 6 nt
    :type loop_len: int

    :return: un boleano
    :rtype: bool
    """
    # patron = r"^\({%d,}\.{%d,}\){%d,}" % (stem_len, loop_len, stem_len)
    # redex desglosado:
    # parentesis ( desde el incio del primer stem_len
    left = len(re.match(r'^\(+', struct).group(0)) if re.match(r'^\(+', struct) else 0
    # puntos del loop_len
    middle = len(re.match(r'^\.+', struct[left:-1:]).group(0)) if re.match(r'^\.+', struct[left:-1:]) else 0
    # parentesis de la derecha, stem_len final
    right = len(re.match(r'\)+', struct[::-1]).group(0)) if re.match(r'\)+', struct[::-1]) else 0
    return left >= stem_len and right >= stem_len and middle >= loop_len # True si al menos len de parentesis y puntos (stem_len, loop_len)
    
def beacon_tm(tm: float, t_min:float=50.0, t_floor:float=45.0) -> float:
    """
    Función de activación para Tm:
      0 si tm < t_floor,
      lineal de 0→1 entre t_floor→t_min,
      1 si tm ≥ t_min
    """
    if tm < t_floor:
        return 0.0
    if tm >= t_min:
        return 1.0
    return (tm - t_floor)/(t_min - t_floor)

# para ver unicamente la tm
def melting_temperature(seq: str) -> float:
    """
    Calcula la temperatura de melting (Tm) de una secuencia de DNA.
    """
    seq = seq.upper().replace("U", "T")  # Asegurar secuencia RNA válida
    return float(mt.Tm_NN(seq))

# para encapsular la funcion de la complemetariedad de bases
def ratio_complement(seq1: str, seq2: str) -> float:
    """
    Calcula el porcentaje de pares complementarios (A-T, G-C) entre dos secuencias alineadas por posición.

    :param seq1: Primera secuencia (ej. beacon)
    :param seq2: Segunda secuencia (ej. amplicón o gRNA)
    :return: Fracción [0.0 - 1.0] de posiciones con pares complementarios
    """
    base_pair = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}
    min_len = min(len(seq1), len(seq2))
    matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if base_pair.get(a.upper()) == b.upper())
    return matches / min_len if min_len else 0.0

def score_beacon(beacon_seq: str, protospacer: str, gRNA_seq: str ,tm_floor:float=50.0) -> float:
    """
  Calcula el score del molecular beacon.

    Combina:
    - R_beacon:nicada → complementariedad con la hebra objetivo.
    - R_gRNA_libre    → falta de complementariedad con el gRNA.
    - F(Tm)           → función de activación según Tm mínima requerida.

    :param beacon_seq: secuencia del beacon diseñado
    :param protospacer: la reversa complementaria para hacer el target del beacon
    :param gRNA_seq: guía CRISPR (protospacer + tracrRNA)
    :param tm_floor: Tm mínima aceptada
    :return: score global, R_beacon:nicada, R_gRNA_libre, F_tm

    """
    seq_amplicon = str(Seq(protospacer).reverse_complement())
    # 1. Estructura
    struct, beacon_mfe = fold_beacon(beacon_seq)
    if not hairpin_correct(struct):
        return 0.0, 0.0, 0.0, 0.0, 0.0

    # 2. Tm
    tm_val = mt.Tm_NN(beacon_seq)
    F = beacon_tm(tm_val, t_min=tm_floor)

    # 3. Ratio beacon:nicada (simulación simple: fraction of complementarity)
    
    R_bn = ratio_complement(beacon_seq, seq_amplicon)

    # 4. Ratio gRNA libre = 1 - fraction complementary to beacon

    R_gr = 1.0 - ratio_complement(beacon_seq, gRNA_seq[:len(beacon_seq)])

    # 5 Energía mínima para la unión beacon con hebra desplazada
    dG = hybridation_energy(beacon_seq, str(Seq(protospacer).reverse_complement()))
    F_e = score_energy(dG, thr=-15.0)

    # Score global
    score = (R_bn * R_gr * F * F_e )**(1/4)
    return score, R_bn, R_gr, F, F_e


def design_beacon(protospacer: str, stem_len:int=11, loop_len:int=9) -> str:
    """
    Concatena dominios:  s (stem), r (loop), t* (loop), s* (stem)
    protospacer ≡ n+x+p (n = nicado seed; x=3nt; p=PAM)
    Deja el fluoróforo en 3' en primera posición del stem.
    """
    # definimos dominios arbitrarios: aquí s = complement(stem of protospacer[0:stem_len]) y el target
    target = str(Seq(protospacer).reverse_complement())
    s_1 = complement(protospacer[:stem_len])[::-1]
    s_2 = s_1[::-1].translate(str.maketrans("ATCG","TAGC"))
    # loop interno r* y t* elegimos secuencias neutrales (A/T rico)
    loop_seq = target[stem_len:stem_len + loop_len]
    # beacon = 5'– s + l + s* –3'
    beacon = s_1 + loop_seq + s_2
    return beacon

def complement(seq:str)->str:
    return seq.translate(str.maketrans("ATCG","TAGC"))