# entrypoint file
# añadimos parser que le de al usuario la facilidad de pasar un path al fichero fasta

# El fichero cli.py debe importar las funciones de core.py y darle los argumentos adecuados a las mismas
# para que estas realicen los cálculos, y luego desde cli.py hay que sacar por pantalla los resultados para que los vea el usuario.
import argparse
# importamos las funciones del codigo
from core import read_fasta 
from core import find_NGG_motivs

# definimos nuestra funcion parser
def get_parse() -> argparse.Namespace:

    # variable
    parser = argparse.ArgumentParser(description= "esto es un programa para encontrar motivos NGG en secuencias de ADN en virus")

    # añadimos los argumentos que necesitará nuestro comando de argumentos
    parser.add_argument("fasta", type = str, help="Especique el archivo .fasta a analizar")
    
    # guardamos los valores en la funcion
    return parser.parse_args()

# para ejecutar la función de manera manual

def main() -> None:

    # reclutamos la función get_parse()
    args = get_parse()

    # asignamos la función de leer el fasta y encontrar el ngg
    sequences =  read_fasta(args.fasta)
    NGG_positions = find_NGG_motivs(sequences)
    # imprimir resultados
    for seq_id, pos in NGG_positions.items():
        print(f"La Secuencia correspondiente con nombre: {seq_id}, tiene {len(pos)} motivos NGG en las posiciones {pos}")

# para ejecutar la función como principal y que no de error

if __name__ == "__main__":
    main()