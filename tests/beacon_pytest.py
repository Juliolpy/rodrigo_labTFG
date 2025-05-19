# Pytest para ver que fallos estamos acumulando en el diseño del Beacon beacon.py
import pytest
import subprocess
# importamos todas las funciones de script
from columbo_design.beacon import fold_beacon, hairpin_correct, beacon_tm, melting_temperature, score_beacon, design_beacon, complement

def fold_beacon_test():
    """
    Coja una secuencia ramdon y vea como se produce el fold con Vienna
    
    
    """

    DNA_seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTA"
    fold_structure = fold_beacon(DNA_seq)
    assert isinstance(fold_structure, str)
    assert set(fold_structure).issubset({"(",")", "."})

# hacemos 2 funciones diferentes para testear la funcion hairping_correct
# esta estrucutra si forma haipir
@pytest.mark.parametrize("dot_bracket",["(((((((......)))))))",
    "((((((((......))))))))"]) # lista de folds
def hairpin_correct_test(dot_bracket):

    """
    Función que comprueba en dot-bracket que hay un stem de al menos stem_len pares,
    un loop de al menos loop_len puntos, y luego el stem inverso.
    Buscamos el patrón: ^\( {stem,}\.{loop,}\){stem,}
    
    """
    assert hairpin_correct(dot_bracket, stem_len=8, loop_len=6)

# esta estructura no tiene que forma hairpin
@pytest.mark.parametrize("dot_bracket", ["...(((...)))...", "....((..))...."])
def hairpin_incorrect_test(dot_bracket):
    assert not hairpin_correct(dot_bracket, stem_len=8, loop_len=6)

# 4. Test para comprobar los scores de beacon_tm en diferentes rangos
def beacon_tm_test():
    assert beacon_tm(44) == 0.0
    assert beacon_tm(50) == 1.0
    assert beacon_tm(47.5) == 0.5

# 5. Test para la melting temperature
def melting_temperature_test():
    DNA_seq = "ATCGATCGATCGATCGATCG"
    tm = melting_temperature(DNA_seq)
    assert isinstance(tm, float)
    assert tm > 0