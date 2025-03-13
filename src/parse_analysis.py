import argparse
from core import read_fasta 
from dna_analysis_2 import search_restriction_sites

# para analizar un fasta y decir si tiene secuencias de restricción

def get_args() -> argparse.Namespace:

    # añadimos la variable
    parser = argparse.ArgumentParser(description= "Analizador para encontrar coincidencias con secuencias de enzimas en archivos fasta")
    # añadimos los argumentos
    parser.add_argument("fasta", type= str, help= "introduzca el archivo fasta para encontrar coincidencias")
    parser.add_argument("--enzyme", choices=["BamHI", "HindII","SalI", "EcoRI", "NotI", "BsmI"], default= "BamHI", help= "Elija una de las enzimas disponibles en nuestra biblioteca porfavor")

    # guardamos los valores en la función
    return parser.parse_args()

# ahora definimos la main function que funciona como argparse y se ejecuta directamente

def main() -> None:
    # llamamos a la función de arriba
    args = get_args()
    # llamamos a la funcion de core.py para leer el archivo fasta y obtener su secuencia, PERO DA ERROR POR SER UN DICT, expected string or bytes - like object
    sequences =  read_fasta(args.fasta)
    # por si el archivo fasta tiene mas de una secuencia
    for seq_id, sequence in sequences.items():
        print(f"\nAnalizando secuencia: {seq_id}")

    if not sequences:
        return # abortamos si hubo algun error
    
    # llamamos a la funcion de dna_analysis_2.py para obtener cuantos y donde estan los sitios de restriccion
    matches, position = search_restriction_sites(sequence, args.enzyme) # buscamos sitios de la secuencia con la enzyma elegida
    if matches is None:
        return

    # Mostramos el resultado
    if matches == 0:
        print(f"No se encontraron sitios de corte para la enzima {args.enzyme}.")
    else:
        coincidencias = "coincidencia" if matches == 1 else "coincidencias"
        posicion_texto = "posición" if matches == 1 else "posiciones"
        print(f"Se encontraron {matches} {coincidencias} en las {posicion_texto} {position} para la enzima {args.enzyme}.")




if __name__ == "__main__":
    main()

    
