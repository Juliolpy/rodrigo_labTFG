import re

sequence = input(str("introduce la secuecia de interés para encontrar sitios de restriccion BamHI, Recuerda que es  < 'GGATCC'> : "))

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
    numero = re.findall(pattern,sequence) # he usado findall, creo que es más sencillo

    lugares = [match.start() for match in re.finditer(pattern, sequence)] # lista con los lugares donde esta


    return numero, lugares
    
    # Buscamos las posiciones de cada patrón en la secuencia
     # Has visto esta manera de construir una lista? Echale un vistazo a https://ellibrodepython.com/list-comprehension-python

# Aquí tendrías que implementar (o importar!) una función que lea un fichero fasta
# ...

def main():
    # Leemos la secuencia de ADN desde el archivo, usando tu función
    # ... sequence = ....
    
    # Buscamos las posiciones de los sitios de corte de las enzimas
    cuantos, sitio = search_restriction_sites(sequence)

    encontro = "Se encontró" if len(cuantos) == 1 else "Se encontraron"
    if len(cuantos) == 1:
        coincidencias = "coincidencia"
    else:
        coincidencias = "coincidencias"

    posicion = "la posición" if len(sitio) == 1 else "las posiciones"

    
    # Imprimimos las posiciones encontradas

    print(f"{encontro} {len(cuantos)} {coincidencias} en {posicion} {sitio} con el sitio de corte de BamHI: GGATCC")

if __name__ == "__main__":
    main()