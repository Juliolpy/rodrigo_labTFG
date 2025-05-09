import pytest
import subprocess
from src.core import read_fasta
from pathlib import Path
from src.core import ColumboParts, find_NGG_motivs, process_genome

# vamos a crear un test para ver que nuestro script funciona
# para ello creamos un fasta ficticio
def test_read_fasta(tmp_path):
    """
    Vamos a ver que la función read_fasta en el script core.py de verdad coge y nos hace
    un diccionario donde los guarda los nombres de las secuencias en keys y las secuencias en 
    los values 
    tal que así:
    
    leer = read_fasta(str(fasta_file))
    --> leer["seq1"] == "ATGCCGTAGCCGTTAGC"


    """
    # contenido del fasta
    fasta_content = ">seq1\nATGCCGTAGCCGTTAGC\n"
    fasta_file = tmp_path / "test.fasta"
    # añadimos el contenido de fasta content dentro de fasta file
    fasta_file.write_text(fasta_content)

    # ejecutamos la función aqui de core.py para que nos lea el fasta

    results = read_fasta(str(fasta_file))

    # comprobamos todo

    # comprobamos que el diccionario tenga la key correcta y la secuencia correspondiente
    assert "seq1" in results
    assert results["seq1"] == "ATGCCGTAGCCGTTAGC"
    # como la función read_fasta nos devuelve un dict con las claves como nombre y secuencias como value
    # con assert y este test nos aseguramos que nos devuelve seq1 como key y la secuencia ATGCCGTAGCCGTTAGC como value
    
    # otro test
    def test_parse_analysis(tmp_path):
        """
        Vemos que parse_analysis.py se ejecuta como debe y devuelve salida
        Para ello creamos el mismo fasta ficticio que en la funcion test_read_fasta
        con el módulo subprocess podemos importar todo el script como si estuvieramos ejecutandolo desde la terminal

        !! La secuencia que se muestra en fasta_content tiene todos los sitios de restriccion
        de todas las enzimas

        """
        fasta_content = ">seq1\nATGCCGTAGCCGTTAGCGCGGCCGCGGATCCGAATTCGTCGACAAGCTTGAATGC\n"
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)

        # ejecutamos el script entero de parse_analysis.py con el archivo de prueba
        result = subprocess.run(["python3", "src/parse_analysis.py", str(fasta_file), "--enzyme", "BamHI"], capture_output=True, text= True)

        # vemos que el script no haya dado ningun error

        assert result.returncode == 0
        assert "Analizando secuencia" in result.stdout