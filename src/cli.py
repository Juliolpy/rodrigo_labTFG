# entrypoint file
# añadimos parser que le de al usuario la facilidad de pasar un path al fichero fasta

# El fichero cli.py debe importar las funciones de core.py y darle los argumentos adecuados a las mismas
# para que estas realicen los cálculos, y luego desde cli.py hay que sacar por pantalla los resultados para que los vea el usuario.
import argparse
# importamos las funciones del codigo
from core import read_fasta, find_NGG_motivs, process_genome, output_NGG_json, output_NGG_pickle
from core_2 import primer3_design_columbo, parse_primers_output, output_pimers_json, output_primers_pickle, score_primers
from core_3 import design_beacon, score_beacon, fold_beacon, beacon_tm, melting_temperature

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

    # colores y el tracRNA que es siempre el mismo
    RED = "\033[91m"
    CIAN = "\033[96m"
    GREEN = "\033[92m"
    GR = "\033[92m"
    RESET = "\033[0m"
    MAG = "\033[35m"
    YELL = "\033[33m"
    tracrRNA = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"

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
            # funcion beacon
            prot = obj._protospacer
            gRNA = prot + tracrRNA
            beacon = design_beacon(str(prot))
            beacon_score = score_beacon(beacon, str(prot), gRNA)
            beacon_struct = fold_beacon(beacon)
            beacon_melting = melting_temperature(beacon)
            beacon_tm_score = beacon_tm(beacon_melting)
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
        print(f"💾 Archivos guardados: {GREEN}output_NGG.pkl, output_PRIMERS.pkl{RESET}")
    else:
        json_file = output_NGG_json(top_candidates)
        json_primers = output_pimers_json(raw_output)
        print(f"💾 Archivos guardados: {RED}output_NGG.json, output_PRIMERS.json{RESET}")
    # imprimir resultados
    for seq_id, pos in NGG_positions.items():
        print(f"La Secuencia correspondiente con nombre: {YELL}{seq_id}{RESET}, tiene {YELL}{len(pos)}{RESET} motivos NGG en las posiciones:")
        print("  ")
        print(">" * 170)
        print(f"{pos}")
        print(("<" * 170))
        # en este caso no ponemos ningun return porque la funcion main() no la vamos a ejecutar en ningun lado y "and" me iba a devolver solo NGG_positions
        # El operador and en Python devuelve el primer valor falsy que encuentre o, si no hay ninguno, devuelve el último valor evaluado.
        print("\n🔬 Resultados de análisis COLUMBO:\n")
        for i, obj in enumerate(top_candidates, 1):
            print(f">>>>>> ColumboPart nº {i} <<<<<<")
            print(f"Creada con el archivo: 📁 route -> {args.fasta}")
            print(f" PAM: {GREEN}{obj._pam}{RESET}")
            print(f" Protospacer: {RED}{obj._protospacer}{RESET}")
            print(f" Localización del protospacer en el genoma: {YELL}{obj._position}{RESET}")
            print(f" Temperatura de melting (Tm): {YELL}{obj._tm:.2f}°C{RESET}")
            print(f" Scores individuales: {YELL}{obj._scores}{RESET}")

            print(f" TracrRNA Streptococcus pyogenes Cas9 {CIAN}{tracrRNA}{RESET} :  ")
            print(f" Guide RNA (gRNA) resultante {GR}5'--{RESET}{RED}{obj._protospacer}{RESET}{CIAN}{tracrRNA}{RESET}{GR}--3'{RESET}")
            print(f" Score global: {YELL}{obj._score_medio:.2f}{RESET}")

            print(f" Diseño de Beacon: {GREEN}{beacon}{RESET} con Score {YELL}{beacon_score:.3f}{RESET}")
            print(f" RNAfold hairping generado {GREEN}{beacon_struct}{RESET} y una Temperatura de melting {YELL}{beacon_melting:.2f}ºC{RESET} Score tm:{YELL}{beacon_tm_score}{RESET} ")
            primer_data = primers_for_obj.get(obj._position, {})
            if "error" in primer_data:
                print(f"{RED}❌ ERROR: diseño de primers interrumpido{RESET} {primer_data['error']}")
            else:
                print(f" Primers diseñados: {MAG}{primer_data['primers']}{RESET}")
                print(f" Score de primer: {YELL}{primer_data['score']:.2f}{RESET}")
            print("-" * 50)

# para ejecutar la función como principal y que no de error

if __name__ == "__main__":
    main()