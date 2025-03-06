import pytest
from src.issue_4_1 import issue_4_1

# definimos la función y creamos un archivo fasta para ver como funciona el script del Issue_4.1

def test_analyze_fasta(tmp_path):
    # y aqui creamos el archivo fasta temporal:
    fasta_content = ">Nombre\nATGCCGTAGCCGTTAGC\n"
    fasta_file = tmp_path / "test_fasta"
    fasta_file.write_text(fasta_content) # añadimos el string de fasta_content a nuestro archivo temporal

    # ejecutamos la función importada del otro script
    results = analyze_fasta(str(fasta_file))


    # comprobamos

    assert results["name"] == "Nombre"
    assert results["longitud"] == 17
    assert results["sequence"] == "ATGCCGTAGCCGTTAGC"
    assert results["percentage_C"] == round((5/17)*100,4)
    assert results["translated_protein"] == "MP*P*"

