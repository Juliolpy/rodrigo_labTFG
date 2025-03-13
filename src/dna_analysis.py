import re
import sys

sequence = input(str("introduce la secuecia de interés para encontrar sitios de restriccion BamHI: "))

def search_restriction_sites(sequence: str) -> dict:
    """
    Busca los sitios de corte de varias enzimas de restricción en una secuencia de ADN.
    
    Args:
        sequence (str): Secuencia de ADN a analizar.
    
    Returns:
        dict: Diccionario con las posiciones donde se encuentra el patrón de cada enzima.
    """
    # Definimos los patrones para cada enzima
    pattern = 'GGATCC' # Sitio de corte de BamHI
    sitios = re.findall(pattern,sequence) # he usado findall, creo que es más sencillo

    return sitios
    
    # Buscamos las posiciones de cada patrón en la secuencia
     # Has visto esta manera de construir una lista? Echale un vistazo a https://ellibrodepython.com/list-comprehension-python

# Aquí tendrías que implementar (o importar!) una función que lea un fichero fasta
# ...

def main(file_path: str) -> None:
    # Leemos la secuencia de ADN desde el archivo, usando tu función
    # ... sequence = ....
    
    # Buscamos las posiciones de los sitios de corte de las enzimas
    results = search_restriction_sites(sequence)
    
    # Imprimimos las posiciones encontradas
    print(f"Se encontraron {len(results)} posiciones con el sitio de corte de BamHI: {results}")

if __name__ == "__main__":
    main(sys.argv[1])