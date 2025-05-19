# Pytest para ver que fallos estamos acumulando en el diseño del Beacon beacon.py
import pytest
import RNA
# importamos todas las funciones de script
from columbo_design.beacon import fold_beacon, hairpin_correct, beacon_tm, melting_temperature, score_beacon, design_beacon, complement

# 1 Función que comprueba que se haga el fold correctamente
def test_fold_beacon():
    """
    Coja una secuencia ramdon y vea como se produce el fold con Vienna
    
    
    """

    DNA_seq = "GCGCGCATAAAAAAATATGCGCGC"
    fold_structure = fold_beacon(DNA_seq)
    assert isinstance(fold_structure, str)
    assert set(fold_structure).issubset({"(",")", "."})

# 2 hacemos 2 funciones diferentes para testear la funcion hairping_correct
# esta estrucutra si forma haipir
@pytest.mark.parametrize("dot_bracket",["((((((((......))))))))",
    "((((((((......))))))))"]) # lista de folds
def test_hairpin_correct(dot_bracket):

    r"""
    Función que comprueba en dot-bracket que hay un stem de al menos stem_len pares,
    un loop de al menos loop_len puntos, y luego el stem inverso.
    Buscamos el patrón: ^\( {stem,}\.{loop,}\){stem,}
    
    """
    assert hairpin_correct(dot_bracket, stem_len=8, loop_len=6)

# esta estructura no tiene que forma hairpin
@pytest.mark.parametrize("dot_bracket", ["...(((...)))...", "....((..))...."])
def test_hairpin_incorrect(dot_bracket):
    assert not hairpin_correct(dot_bracket, stem_len=8, loop_len=6)

# 3. Test para comprobar los scores de beacon_tm en diferentes rangos
def test_beacon_tm():
    assert beacon_tm(44) == 0.0
    assert beacon_tm(50) == 1.0
    assert beacon_tm(47.5) == 0.5

# 4. Test para la melting temperature
def test_melting_temperature():
    DNA_seq = "ATCGATCGATCGATCGATCG"
    tm = melting_temperature(DNA_seq)
    assert isinstance(tm, float)
    assert tm > 0

# 5. Test para el score del beacon
def test_score_beacon():
    protospacer = "ATGCGTACGTAGCTAGCTAAGG" # secuencia ramdon modelo protospacer 20 nt + GG final = 22 nt (PAM)
    beacon = design_beacon(protospacer)
    amplicon = "CCTTAGCTAGCTACGTACGCAT"  # reversa complementaria del protospacer
    tracrRNA = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"
    gRNA = protospacer + tracrRNA

    score, R_bn, R_gr, F = score_beacon(beacon, gRNA, protospacer)
    assert score > 0, "score global que tiene que ser mayor que 0"
    assert 0 <= R_bn <= 1, "R_beacon debe ser entre 0 y 1"
    assert 0 <= R_gr <= 1, "R_gRNA dene estar entre 0 y 1"
    assert 0 <= F <= 1, "F_tm debe estar enre 0 y 1"

# 6. Test para el diseño del beacon
def test_beacon_design():
    protospacer = "ATGCGTACGTAGCTAGCTAAGG"
    beacon = design_beacon(protospacer)
    assert isinstance(beacon, str)
    assert len(beacon) == 8 + 6 + 8 # stem + loop + loop + stem


# 7. Test de la funcion complement
def test_complement():
    assert complement("ATCG") == "TAGC"
    assert complement("GCTA") == "CGAT"

