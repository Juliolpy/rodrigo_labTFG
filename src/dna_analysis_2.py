import re

def search_restriction_sites(sequence: str, enzyme_name: str) -> tuple:
    """
    Busca los sitios de corte de una enzima de restricción específica en una secuencia de ADN.
    
    Args:
        sequence (str): Secuencia de ADN a analizar.
        enzyme_name (str): Nombre de la enzima a buscar.

    Returns:
        tuple: Número de coincidencias y lista de posiciones donde se encuentra el patrón.
    """
    # Definimos los sitios de corte de algunas enzimas de restricción
    enzymes = {
        "BamHI": "GGATCC",
        "HindIII": "AAGCTT",
        "SalI": "GTCGAC",
        "EcoRI": "GAATTC",
        "NotI": "GC.NGC",  # Se usa "." en regex para "N" (cualquier base)
        "BsmI": "GA.STC"   # Se usa "." en regex para "S" (G o C)
    }

    # Verificamos si la enzima está en el diccionario
    if enzyme_name not in enzymes:
        print("ERROR: Enzyme not found or incorrect name provided")
        return None, None

    # Obtener el patrón de la enzima
    pattern = enzymes[enzyme_name]

    # Buscar coincidencias en la secuencia
    matches = list(re.finditer(pattern, sequence))

    # Extraer las posiciones donde se encuentran las coincidencias
    positions = [match.start() for match in matches]

    return len(matches), positions

def main():
    # Pedimos la secuencia al usuario
    sequence = input("Introduce la secuencia de ADN de interés: ").upper()

    # Pedimos la enzima que quiere buscar
    enzyme_name = input("¿Qué enzima deseas buscar? (BamHI, HindIII, SalI, EcoRI, NotI, BsmI): ")

    # Buscamos los sitios de corte
    num_matches, positions = search_restriction_sites(sequence, enzyme_name)

    # Si hubo error en la búsqueda, terminamos el programa
    if num_matches is None:
        return

    # Mostramos el resultado
    if num_matches == 0:
        print(f"No se encontraron sitios de corte para la enzima {enzyme_name}.")
    else:
        coincidencias = "coincidencia" if num_matches == 1 else "coincidencias"
        posicion_texto = "posición" if num_matches == 1 else "posiciones"
        print(f"Se encontraron {num_matches} {coincidencias} en las {posicion_texto} {positions} para la enzima {enzyme_name}.")

if __name__ == "__main__":
    main()