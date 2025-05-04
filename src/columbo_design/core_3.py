# TERCERA PARTE DEL SCRIPT CORE donde se diseña el Beacon con sus scores correspondientes

# Importamos los módulos necesarios
import subprocess
from tempfile import NamedTemporaryFile
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO
from Bio.Seq import Seq
import re
import RNA  # Make sure this import is at the top of your file

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
    rna = seq.replace('T', 'U')
    # Use ViennaRNA to fold the RNA sequence
    structure, _ = RNA.fold(rna)
    return structure

# vemos que el hairpin tiene efecticamente 8 nt en stem y 6 en el loop
def hairpin_correct(struct: str, stem_len:int=8, loop_len:int=6) -> bool: # true or false
    """
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
    patron = r"^\({%d,}\.{%d,}\){%d,}" % (stem_len, loop_len, stem_len)
    return re.match(patron, struct) is not None

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
    seq = seq.upper().replace("U", "T")  # Asegurar secuencia DNA válida
    return float(mt.Tm_NN(seq))

def score_beacon(beacon_seq: str,amplicon_seq: str,gRNA_seq: str,tm_floor:float=50.0) -> float:
    """
    Calcula R_beacon:nicada y R_gRNA_libre de forma simplificada
    (a falta de NUPACK real, aquí suponemos ideales si no hay
    motivos de homología extensiva), y F_tm.
    """
    seq_amplicon = Seq(amplicon_seq).reverse_complement()
    # 1. Estructura
    struct = fold_beacon(beacon_seq)
    if not hairpin_correct(struct):
        return 0.0, 0.0, 0.0, 0.0

    # 2. Tm
    tm_val = mt.Tm_NN(beacon_seq)
    F = beacon_tm(tm_val, t_min=tm_floor)

    # 3. Ratio beacon:nicada (simulación simple: fraction of complementarity)
    #    = número de nt complementarios / len(beacon)
    comp = sum(1 for a,b in zip(beacon_seq, seq_amplicon) if (a=="A" and b=="T") or (a=="T" and b=="A") or (a=="C" and b=="G") or (a=="G" and b=="C"))
    R_bn = (comp/len(beacon_seq))

    # 4. Ratio gRNA libre = 1 - fraction complementary to beacon
    comp2 = sum(1 for a,b in zip(beacon_seq, gRNA_seq) if (a=="A" and b=="T") or (a=="T" and b=="A") or (a=="C" and b=="G") or (a=="G" and b=="C"))
    R_gr = 1.0 - comp2/len(beacon_seq)

    # Score global
    score = R_bn * R_gr * F
    return score, R_bn, R_gr, F

def design_beacon(protospacer: str,pam_seq: str = "NGG",stem_len:int=8,loop_len:int=6) -> str:
    """
    Concatena dominios:  s (stem), r (loop), t* (loop), s* (stem)
    protospacer ≡ n+x+p (n = nicado seed; x=3nt; p=PAM)
    Deja el fluoróforo en 3' en primera posición del stem.
    """
    # definimos dominios arbitrarios: aquí s = complement(stem of protospacer[0:stem_len])
    s_1 = complement(protospacer[:stem_len])[::-1]
    s_2 = s_1[::-1].translate(str.maketrans("ATCG","TAGC"))
    # loop interno r* y t* elegimos secuencias neutrales (A/T rico)
    r_1 = "A"*loop_len
    t_2 = "T"*loop_len
    # beacon = 5'– s + r* + t* + s* –3'
    beacon = s_1 + r_1 + t_2 + s_2
    return beacon

def complement(seq:str)->str:
    return seq.translate(str.maketrans("ATCG","TAGC"))