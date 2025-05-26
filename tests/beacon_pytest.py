# tests/beacon_pytest.py
import pytest
import RNA
from Bio.Seq import Seq

from columbo_design.beacon import (
    fold_beacon,
    hybridation_energy,
    score_energy,
    hairpin_correct,
    beacon_tm,
    melting_temperature,
    ratio_complement,
    count_hairpin,
    design_beacon,
    score_beacon,
)

# 1) fold_beacon → devuelve (struct, mfe)
def test_fold_beacon():
    dna = "GCGCGCATAAAAAAATATGCGCGC"
    struct, mfe = fold_beacon(dna)
    assert isinstance(struct, str)
    assert set(struct).issubset({"(", ")", "."})
    assert isinstance(mfe, float)

# 2) count_hairpin → cuenta stems y loop
def test_count_hairpin():
    struct = "(((((((((((((......)))))))))))))"
    stem1, loop, stem2 = count_hairpin(struct)
    assert stem1 == 13
    assert loop == 6
    assert stem2 == 13

# 3) hairpin_correct → devuelve media normalizada entre 0 y 1
@pytest.mark.parametrize("struct, exp", [
    ("(((((((((((.........)))))))))))", 1.0),   # ideal 13/11,6/9,13/11 → normalized = 1
    ("(((((....)))))", pytest.approx((5/11 + 4/9 + 5/11)/3)),
])
def test_hairpin_correct(struct, exp):
    score = hairpin_correct(struct, ideal_stem=11, ideal_loop=9)
    assert pytest.approx(score) == exp

# 4) beacon_tm en rangos
def test_beacon_tm():
    assert beacon_tm(44) == 0.0
    assert beacon_tm(50) == 1.0
    assert beacon_tm(47.5) == 0.5

# 5) melting_temperature
def test_melting_temperature():
    dna = "ATCGATCGATCGATCGATCG"
    tm = melting_temperature(dna)
    assert isinstance(tm, float)
    assert tm > 0

# 7) score_energy
def test_score_energy():
    assert score_energy(0) == 0.0
    assert score_energy(-30, thr=-15) == 1.0
    # intermedio: strenght = 10, max=30 → 10/30 = 0.3333
    assert score_energy(-10, thr=-15) == pytest.approx(10/30)

# 8) hybridation_energy debe ser bastante negativo para un emparejamiento perfecto
def test_hybridation_energy():
    beacon = "ACGCCATTCGGTATGGTCACCGAATGGCGT"
    target = str(Seq(beacon).reverse_complement())
    dG = hybridation_energy(beacon, target)
    assert isinstance(dG, float)
    assert dG < -15.0

# 9) design_beacon y su longitud
def test_design_beacon():
    prot = "ATGCGTACGTAGCTAGCTAAGG"
    b = design_beacon(prot, stem_len=8, loop_len=6)
    assert isinstance(b, str)
    assert len(b) == 8 + 6 + 8

# 10) score_beacon firma actual (6 valores) y rangos
def test_score_beacon():
    prot = "ATGCGTACGTAGCTAGCTAAGG"
    beacon = design_beacon(prot, stem_len=8, loop_len=6)
    # construimos target y gRNA
    beacon_region = prot[:len(beacon)]
    target = str(Seq(beacon_region).reverse_complement())
    tracr = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"
    gRNA = prot + tracr

    score, R_bn, R_gr, F_tm, F_e, hp = score_beacon(beacon, target, gRNA, tm_floor=50.0)
    assert 0.0 <= score <= 1.0
    assert 0.0 <= R_bn <= 1.0
    assert 0.0 <= R_gr <= 1.0
    assert 0.0 <= F_tm <= 1.0
    assert 0.0 <= F_e <= 1.0
    assert 0.0 <= hp   <= 1.0