# entrypoint file
# añadimos parser que le de al usuario la facilidad de pasar un path al fichero fasta

# El fichero cli.py debe importar las funciones de core.py y darle los argumentos adecuados a las mismas
# para que estas realicen los cálculos, y luego desde cli.py hay que sacar por pantalla los resultados para que los vea el usuario.
import argparse
# importamos las funciones del codigo
from core import read_fasta, find_NGG_motivs, process_genome, output_NGG_json, output_NGG_pickle
from core_2 import primer3_design_columbo, parse_primers_output, output_pimers_json, output_primers_pickle, score_primers

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

    primers_for_obj = {}
    for seq_id, seq in sequences.items():
        for obj in top_candidates:
            # definir region alrrederor del PAM para diseñar los primers ( x nt upstream y downstream)
            flank = 100
            pam_pos = obj._position
            full_seq = sequences[seq_id] # solo la secuencia del dict
            # controlar los bordes del genoma
            start = max(0, pam_pos - flank)
            end = min(len(full_seq), pam_pos + 23 + flank)  # 23 = protospacer+PAM
            region = full_seq[start:end]
            try:
                # diseñar dichos primers
                raw_output = primer3_design_columbo(region, pam_pos - start)
                parsed_primers = parse_primers_output(raw_output, region)
                score = score_primers(parsed_primers, pam_pos-start, protospacer_len=23)
                primers_for_obj[pam_pos] = {
                    "primers": parsed_primers,
                    "score": round(score, 2)
                }
            except Exception as exc:
                primers_for_obj[pam_pos] = {"error": str(exc)}
    # especificamos el archivo que queremos
    if args.output == "pickle":
        pickle_file = output_NGG_pickle(top_candidates)
        pickle_primers = output_primers_pickle(raw_output)
        print("📦 Archivos guardados: output_NGG.pkl, output_PRIMERS.pkl")
    else:
        json_file = output_NGG_json(top_candidates)
        json_primers = output_pimers_json(raw_output)
        print("📁 Archivos guardados: output_NGG.json, output_PRIMERS.json")
    # imprimir resultados
    for seq_id, pos in NGG_positions.items():
        print(f"La Secuencia correspondiente con nombre: {seq_id}, tiene {len(pos)} motivos NGG en las posiciones:")
        print("  ")
        print(">" * 170)
        print(f"{pos}")
        print(("<" * 170))
        # en este caso no ponemos ningun return porque la funcion main() no la vamos a ejecutar en ningun lado y "and" me iba a devolver solo NGG_positions
        # El operador and en Python devuelve el primer valor falsy que encuentre o, si no hay ninguno, devuelve el último valor evaluado.
        print("\n🔬 Resultados de análisis COLUMBO:\n")
        for i, obj in enumerate(top_candidates, 1):
            print(f"-----ColumboPart nº {i}-----")
            print(f"creada con el archivo: 📁 route -> {args.fasta}")
            print(f" PAM: {obj._pam}")
            print(f" Protospacer: {obj._protospacer}")
            print(f" Localización del protospacer en el genoma: {obj._position}")
            print(f" Temperatura de melting (Tm): {obj._tm:.2f}°C")
            print(f" Scores individuales: {obj._scores}")
            print(f" Score global: {obj._score_medio:.2f}")

            primer_data = primers_for_obj.get(obj._position, {})
            if "error" in primer_data:
                print(f"❌ Error al diseñar los primers: {primer_data['error']}")
            else:
                print(f"Primers diseñados: {primer_data['primers']}")
                print(f"Score de primer: {primer_data['score']:.2f}")
            print("-" * 50)

# para ejecutar la función como principal y que no de error

if __name__ == "__main__":
    main()