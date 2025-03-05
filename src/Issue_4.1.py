from Bio import SeqIO
from Bio.Seq import Seq

def analyze_fasta(fasta_file: str) -> dict:
    """
    definimos la funcion para analizar un archivo fasta y devolver
    la información que nos piden

    --> Parámetros: ruta del archivo a analizar 
    --> Output : return, crea un dict con la siguiente información:
        -"name" : ID de la secuencia
        -"length" : Longitud de la secuencia
        -"sequence" : La secuencia
        -"percentage_C": porcentaje de citosina en la secuencia
        -"transalated_protein": la proteina teórica producida
    
    """
    # se lee la secuencia usando SeqIO
    record = next(SeqIO.parse(fasta_file, "fasta")) 
    # definimos el nombre y la secuencia
    name = record.id
    seq = record.seq
    # para calcular el porcetaje de citosina
    nc = seq.count("C")
    per = (nc / len(seq))*100
    # paso opcional que hago para ver las funciones de Biopython
    proteina = seq.translate() # asi de facil xd


    return {
        "name": name,
        "longitud": len(seq),
        "sequence": str(seq),
        "percentage_C": round(per, 4),
        "translated_protein": str(proteina)
    }
if __name__ == "__main__":
    import sys
    fasta_file = sys.argv[1]
    results = analyze_fasta(fasta_file)

    print(f"Nombre de la secuencia: {results['name']} y longitud: {results['longitud']}")
    print(f"Secuencia: {results['sequence']}")
    print(f"Porcentaje de C: {results['percentage_C']}%")
    print(f"Proteína teórica traducida: {results['translated_protein']}")