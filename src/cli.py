# entrypoint file
# añadimos parser que le de al usuario la facilidad de pasar un path al fichero fasta

# El fichero cli.py debe importar las funciones de core.py y darle los argumentos adecuados a las mismas
# para que estas realicen los cálculos, y luego desde cli.py hay que sacar por pantalla los resultados para que los vea el usuario.
import argparse
# importamos las funciones del codigo
from core import read_fasta, find_NGG_motivs, process_genome, output_NGG_json, output_NGG_pickle

# definimos nuestra funcion parser
def get_parse() -> argparse.Namespace:

    # variable
    parser = argparse.ArgumentParser(description= "esto es un programa para encontrar motivos NGG en secuencias de ADN en virus")

    # añadimos los argumentos que necesitará nuestro comando de argumentos
    parser.add_argument("fasta", type = str, help="Especique el archivo .fasta a analizar")
    parser.add_argument("--output", choices=["json", "pickle"], default= "json", help= "Especifica el tipo de archivo de salida que deseas obtener")
    # guardamos los valores en la funcion
    return parser.parse_args()

# para ejecutar la función de manera manual

def main() -> None:

    # reclutamos la función get_parse()
    args = get_parse()

    # asignamos la función de leer el fasta y encontrar el ngg
    sequences =  read_fasta(args.fasta)
    NGG_positions = find_NGG_motivs(sequences)
    top_candidates = process_genome(sequences, NGG_positions) # vemos cuales son los resultados de la clase y los almacenamos en la variable top_candidates
    # especificamos el archivo que queremos
    if args.output == "pickle":
        pickle_file = output_NGG_pickle(top_candidates)
    else:
        json_file = output_NGG_json(top_candidates)
    # imprimir resultados
    for seq_id, pos in NGG_positions.items():
        print(f"La Secuencia correspondiente con nombre: {seq_id}, tiene {len(pos)} motivos NGG en las posiciones:")
        print("  ")
        print(">" * 170)
        print(f"{pos}")
        print(("<" * 170))
        # en este caso no ponemos ningun return porque la funcion main() no la vamos a ejecutar en ningun lado y "and" me iba a devolver solo NGG_positions
        # El operador and en Python devuelve el primer valor falsy que encuentre o, si no hay ninguno, devuelve el último valor evaluado.
        print("\n🔬 Resultados de análisis CRISPR-Cas9:\n")
        counter = 0
        for obj in top_candidates:
            counter += 1
            print(f"-----ColumboPart nº {counter}-----")
            print(f" PAM: {obj.PAM}")
            print(f" Protospacer: {obj.protospacer}")
            print(f" Localización en el genoma: {obj.position}")
            print(f" Temperatura de melting (Tm): {obj.tm:.2f}°C")
            print(f" Scores individuales: {obj.scores}")
            print(f" Score global: {obj.score_medio:.2f}")
            print("-" * 50)  # Separador entre cada resultado

# para ejecutar la función como principal y que no de error

if __name__ == "__main__":
    main()