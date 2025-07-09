# tests/test_core.py

import pytest
from Bio.Seq import Seq
from columbo_design.core import read_fasta, find_NGG_motivs, process_genome, ColumboParts, Out_of_frame_ERROR

# 1️ Test de read_fasta
def test_read_fasta(tmp_path):
    """
    Test for `read_fasta`.

    Verifies that the function correctly parses a FASTA file and returns a 
    dictionary mapping sequence IDs to their nucleotide sequences

    """
    fasta_content = ">seq1\nATGCCGTAGCCGTTAGC\n"
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    result = read_fasta(str(fasta_file))

    assert "seq1" in result
    assert result["seq1"] == "ATGCCGTAGCCGTTAGC"

# 2️ Test de find_NGG_motivs
def test_find_NGG_motivs():
    """
    Test for `find_NGG_motivs`.

    Checks that the function identifies NGG motifs (PAM sites) in a given 
    sequence dictionary and returns their positions.

    """
    # ejemplo con una secuencia sencilla que tenga un NGG
    seq_dict = {"seq1": "ATCGTACCGGCATCGGCTAG"}
    motifs = find_NGG_motivs(seq_dict)

    # Debe encontrar la posición del primer NGG
    assert "seq1" in motifs
    # En esta secuencia hay GG en posiciones 8 y 14 (index+1 porque usas i+1)
    assert motifs["seq1"] == [8, 14]

# 3️ Test de process_genome
def test_process_genome():
    """
    Test for `process_genome`.

    Validates that the function generates a list of ColumboParts objects
    from sequences containing NGG motifs and checks key attributes such as
    PAM, protospacer length, and position.
    """
    # Creamos un FASTA ficticio que tenga al menos un NGG válido
    seq_dict = {"seq1": "A" * 30 + "AGG" + "A" * 50}  # Posición 31
    motifs = find_NGG_motivs(seq_dict)

    # Ejecutamos process_genome
    columbo_parts = process_genome(seq_dict, motifs)

    # Debe devolver una lista con al menos un ColumboParts
    assert isinstance(columbo_parts, list)
    assert len(columbo_parts) > 0
    assert isinstance(columbo_parts[0], ColumboParts)

    # Verificamos que los atributos sean coherentes
    part = columbo_parts[0]
    assert str(part.pam) == "GGA"
    assert part.position == 31
    assert len(part.protospacer) == 28  # 25 + PAM

# 4️ Test de ColumboParts: excepciones
def test_out_of_frame_error():
    """
    Test for `ColumboParts` exception.

    Ensures that the `Out_of_frame_ERROR` is raised when a PAM position is
    too close to the end of a sequence to extract a full protospacer.
    """
    # Caso en el que la posición del PAM excede la longitud de la secuencia
    seq = Seq("A" * 10)  # Demasiado corta

    with pytest.raises(Out_of_frame_ERROR):
        ColumboParts(seq, 8)

# 5️ Test de ColumboParts: cálculo de scores
def test_columbo_scores():
    """
    Test for `ColumboParts.scores` and related properties.

    Confirms that the `ColumboParts` object calculates a list of binary 
    scores, the mean score, and a valid melting temperature (Tm).
    """
    # Creamos una secuencia que cumpla los criterios
    seq = Seq("A" * 50 + "AGG" + "C" * 50)
    position = 50
    part = ColumboParts(seq, position)

    # Debe devolver una lista de 4 scores
    scores = part.scores
    assert isinstance(scores, list)
    assert len(scores) == 4
    for s in scores:
        assert s in [0, 1]

    # Verificamos que el score medio sea correcto
    expected_score_medio = sum(scores) / len(scores)
    assert part.score_medio == expected_score_medio

    # Temperatura de melting debe ser un número positivo
    assert part.tm >= 0
