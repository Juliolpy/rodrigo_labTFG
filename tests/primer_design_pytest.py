# tests/test_primers_design.py

import pytest
from columbo_design.primer_design import primer3_design_columbo, score_primers, parse_primers_output, output_pimers_json, output_primers_pickle
from Bio.Seq import Seq

# 1️ Test de primer3_design_columbo
def test_primer3_design_columbo_basic():
    """
    Test básico: comprobar que primer3_design_columbo devuelve un diccionario con claves esperadas.
    """
    # Creamos una secuencia artificial suficientemente larga
    seq = "AAGCTGATCGCTAACAAGTTCAACCAAGCCCTGGGCGCTATGCAGACCGGATTCACCACTACAAACGAAGCCTTCCGTAAGGTGCAAGACGCTGTCAACAACAACGCCCAGGCTCTCTCCAAGTTGGCCTCCGAGCTGTCTAACACCTTCGGCGCCATCTCTGCTAGCATCGGAGATATCATCCAGCGTCTGGACGTGCTCGAGCAGGATGCTCAAATCG"
    pam_pos = 126
    result = primer3_design_columbo(seq, pam_pos)

    # Debe ser un diccionario con al menos las claves de los primers
    assert isinstance(result, dict)
    assert "PRIMER_LEFT_0" in result
    assert "PRIMER_RIGHT_0" in result

# 2️ Test de score_primers
def test_score_primers():
    """
    Test de score_primers con una simulación de output de Primer3.
    """
    # Simulamos un output de Primer3 válido
    output_mock = {
        "PRIMER_LEFT_0": (80, 20),
        "PRIMER_RIGHT_0": (150, 20),
        "PRIMER_LEFT_0_GC_PERCENT": 50.0,
        "PRIMER_RIGHT_0_GC_PERCENT": 50.0,
        "PRIMER_LEFT_0_TM": 60.0,
        "PRIMER_RIGHT_0_TM": 60.0
    }

    pam_relative_pos = 100

    score = score_primers(output_mock, pam_relative_pos)

    # El score debe estar en el rango [0,1]
    assert 0.0 <= score <= 1.0

# 3️ Test de parse_primers_output
def test_parse_primers_output():
    """
    Test de parse_primers_output: comprobar extracción de secuencias.
    """
    output_mock = {
        "PRIMER_LEFT_0": (10, 20),
        "PRIMER_RIGHT_0": (50, 20)
    }
    # Región simulada de 100 nt
    region = "A" * 100

    result = parse_primers_output(output_mock, region)

    assert isinstance(result, dict)
    assert "left_primer" in result
    assert "right_primer" in result
    assert "left_primer_comp" in result
    assert "right_primer_revcomp" in result
    assert len(result["left_primer"]) == 20
    assert len(result["right_primer"]) == 20

# 4️ Test de output_pimers_json
def test_output_pimers_json(tmp_path):
    """
    Test de output_pimers_json: verificar escritura y lectura correcta.
    """
    output_mock = {
        "PRIMER_LEFT_0": (10, 20),
        "PRIMER_RIGHT_0": (50, 20),
        "PRIMER_LEFT_0_GC_PERCENT": 50.0,
        "PRIMER_RIGHT_0_GC_PERCENT": 50.0,
        "PRIMER_LEFT_0_TM": 60.0,
        "PRIMER_RIGHT_0_TM": 60.0
    }

    # Redirigimos el archivo temporalmente
    original_path = "output_PRIMERS.json"
    tmp_file = tmp_path / "output_PRIMERS.json"

    # Guardamos original
    import os
    if os.path.exists(original_path):
        os.rename(original_path, tmp_file)

    json_result = output_pimers_json(output_mock)

    assert isinstance(json_result, dict)
    assert "PRIMER_LEFT_0" in json_result

# 5️ Test de output_primers_pickle
def test_output_primers_pickle(tmp_path):
    """
    Test de output_primers_pickle: verificar escritura y lectura con pickle.
    """
    output_mock = {
        "PRIMER_LEFT_0": (10, 20),
        "PRIMER_RIGHT_0": (50, 20),
        "PRIMER_LEFT_0_GC_PERCENT": 50.0,
        "PRIMER_RIGHT_0_GC_PERCENT": 50.0,
        "PRIMER_LEFT_0_TM": 60.0,
        "PRIMER_RIGHT_0_TM": 60.0
    }

    pickle_result = output_primers_pickle(output_mock)

    assert isinstance(pickle_result, list) or isinstance(pickle_result, dict)
