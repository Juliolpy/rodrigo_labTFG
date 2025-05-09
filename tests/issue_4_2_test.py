import pytest
from src.issue_4_2 import issue_4_2
import RNA

# definimos la función ya que vamos a probar si funciona en una secuencia que creemos
# las demás funciones de vienna_RNA están definidas en el otro script que importamos

def test_vienna_RNA():

    test_seq = "AUGCGCGGCGCACACGAGCAU" # secuencia de prueba

    resultado = vienna_RNA(test_seq) # llamamos a la funcion del otro script

    assert resultado["Secuencia"] == "AUGCGCGGCGCACACGAGCAU"
    assert resultado["secondary_structure"] == "((((.((.......)).))))"
    assert resultado["min_energy"] == -5.10 # tiene que se run float, que antes lo habia puesto como string
    assert len(resultado["secondary_structure"]) == len(test_seq)
