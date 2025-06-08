# tests/test_cli.py

import subprocess
import pytest
import os
from pathlib import Path

# 1️ Test básico: ejecución CLI con output JSON
def test_cli_json_output(tmp_path):
    """
    Test de integración básico: comprobar que cli.py se ejecuta con output JSON.
    """

    # Creamos un FASTA de prueba (valida y con un PAM NGG)
    fasta_content = ">seq1\n" + "A" * 100 + "AGG" + "T" * 100
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    # Lanzamos el CLI con subprocess
    result = subprocess.run(
        ["python3", "src/columbo_design/cli.py", str(fasta_file), "--output", "json"],
        capture_output=True,
        text=True
    )

    # Verificamos que el script se ejecutó sin error
    assert result.returncode == 0

    # Verificamos que se imprimió el mensaje esperado de guardado en JSON
    assert "output_NGG.json" in result.stdout
    assert "output_PRIMERS.json" in result.stdout

    # Verificamos que los archivos se generaron en el directorio de trabajo
    assert os.path.exists("output_NGG.json")
    assert os.path.exists("output_PRIMERS.json")

# 2️ Test CLI con output PICKLE
def test_cli_pickle_output(tmp_path):
    """
    Test de integración: ejecución de cli.py con output en formato pickle.
    """

    fasta_content = ">seq1\n" + "A" * 100 + "AGG" + "T" * 100
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    result = subprocess.run(
        ["python3", "src/columbo_design/cli.py", str(fasta_file), "--output", "pickle"],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0
    assert "output_NGG.pkl" in result.stdout
    assert "output_PRIMERS.pkl" in result.stdout

    assert os.path.exists("output_NGG.pkl")
    assert os.path.exists("output_PRIMERS.pkl")

# 3️ Test de error: Fichero FASTA no existe
def test_cli_file_not_found():
    """
    Test de control: lanzar cli.py con un archivo FASTA que no existe.
    El CLI debería fallar con error controlado.
    """

    result = subprocess.run(
        ["python3", "src/columbo_design/cli.py", "nonexistent.fasta", "--output", "json"],
        capture_output=True,
        text=True
    )

    # Debe fallar (código de error diferente de 0)
    assert result.returncode != 0
    assert "No such file or directory" in result.stderr or "No such file" in result.stderr or "No existe" in result.stderr