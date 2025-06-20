# Entrypoint file
# añadimos parser que le de al usuario la facilidad de pasar un path al fichero fasta

# El fichero cli.py debe importar las funciones de core.py y darle los argumentos adecuados a las mismas
# para que estas realicen los cálculos, y luego desde cli.py hay que sacar por pantalla los resultados para que los vea el usuario.
import argparse
from Bio.Seq import Seq
# importamos las funciones del codigo
from .core import read_fasta, find_NGG_motivs, process_genome, output_NGG_json, output_NGG_pickle
from .primer_design import primer3_design_columbo, parse_primers_output, output_pimers_json, output_primers_pickle, score_primers
from .beacon import design_beacon, score_beacon, fold_beacon, beacon_tm, melting_temperature

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

    # definimos distancias
    max_protospacer_lenght = 30
    initial_protospacer_lenght = 23 # Protospacer + PAM

    beacons_for_obj = {}
    primers_for_obj = {}
    for seq_id, seq in sequences.items():
        for obj in top_candidates:
            # posicion del pam y la secuencia completa
            pam_pos = obj._position
            full_seq = sequences[seq_id] # solo la secuencia del dict
            # para diferentes longitudes de protospacer si el diseño falla
            prot_len = initial_protospacer_lenght
            primers_ok = False
            primer_data = {"error" : "No se pudieron encontrar primers para estas regiones"}

            while not primers_ok and prot_len <= max_protospacer_lenght:
                start_ps = pam_pos - prot_len
                if start_ps < 0:
                    break
                new_protospacer = full_seq[start_ps: pam_pos + 3]
                 # controlar los bordes del genoma
                flank = 100
                start = max(0, start_ps - flank)
                end = min(len(full_seq), pam_pos + 3 + flank)  # 23 = protospacer+PAM
                region = full_seq[start:end]

                try:
                    raw_output = primer3_design_columbo(region, pam_pos - start)
                    parsed_primers = parse_primers_output(raw_output, region)
                    score = score_primers(raw_output, pam_pos-start, protospacer_len=len(new_protospacer))
                    primer_data = {
                        "primers": parsed_primers,
                        "score": round(score, 2),
                        "used_prot_len": prot_len
                    }
                    primers_ok = True
                # except ValueError as exc:
                    # Handle the case where no primers are found for this region
                    # primers_for_obj[pam_pos] = {"error": str(exc)}
                except Exception:
                # Handle other errors in primer design
                # primers_for_obj[pam_pos] = {"error": str(exc)}
                    prot_len += 1

            primers_for_obj[pam_pos] = primer_data

            # funcion beacon
            gRNA = str(obj._protospacer) + tracrRNA
            beacon = design_beacon(obj._protospacer)
            beacon_struct = fold_beacon(beacon)
            beacon_melting = melting_temperature(beacon)
            beacon_tm_score = beacon_tm(beacon_melting)
            beacon_score, R_bn, R_gr, F_tm = score_beacon(beacon, obj._protospacer, gRNA)
        
            beacons_for_obj[obj._position] = {
            "beacon":          beacon,
            "beacon_struct":   beacon_struct,
            "tm":              beacon_melting,
            "tm_score":        beacon_tm_score,
            "score":           beacon_score,
            "R_bn":            R_bn,
            "R_gr":            R_gr,
            "F_tm":            F_tm,
            }

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
            print(f" PAM: {YELL}3'-- {RESET}{GREEN}{obj._pam}{RESET} {YELL}--5'{RESET}")
            print(f" Amplicon:                {YELL}3'-- {RESET}{CIAN}{str(Seq(obj._protospacer).reverse_complement())}{RESET}{YELL} --5'{RESET} 🎯 Target del Beacon (reversa complementaria)")
            print(f"                                         ⇄               ")
            print(f" Hebra complementaria:    {YELL}5'--{RESET} {GREEN}{str(Seq(obj._protospacer).complement())}{RESET} {YELL}--3'{RESET} --> Hebra complementaria del Protospacer" )
            print(f"                               |||||||||||||||||||||||      ")
            print(f" Protospacer:             {YELL}3'--{RESET} {RED}{obj._protospacer}{RESET} {YELL}--5'{RESET}  de longitud = {YELL}{len(obj._protospacer)} nt {RESET} ")
            print("                                                          ")
            print(f" Localización del protospacer en el genoma: {YELL}{obj._position}{RESET}")
            print(f" Temperatura de melting (Tm): {YELL}{obj._tm:.2f}°C{RESET}")
            print(f" Scores individuales = {YELL}{obj._scores}{RESET}")
            print(f" Score global protospacer = {YELL}{obj._score_medio:.2f}{RESET}")

            print(f" TracrRNA Streptococcus pyogenes Cas9 {CIAN}{tracrRNA}{RESET} de longitud = {YELL}{len(tracrRNA)} nt{RESET}")
            print(f" Guide RNA (gRNA) resultante {GR}5'--{RESET}{RED}{obj._protospacer}{RESET}{CIAN}{tracrRNA}{RESET}{GR}--3'{RESET} de longitud = {YELL}{len(tracrRNA)+len(obj._protospacer)} nt{RESET}")
            print(f" gRNA (corte): {MAG}{gRNA[:len(beacon)]}{RESET}")
            print("                                                          ")
            print(f" Diseño de Beacon: {GREEN}{beacons_for_obj[obj._position]['beacon']}{RESET} con Score GLobal = {YELL}{beacons_for_obj[obj._position]['score']:.2f}{RESET}")
            print(f" Score R_beacon = {YELL}{beacons_for_obj[obj._position]['R_bn']:.2f}{RESET}, Score R_guide = {YELL}{beacons_for_obj[obj._position]['R_gr']:.2f}{RESET}")
            print(f" RNAfold hairping generado {GREEN}{beacons_for_obj[obj._position]['beacon_struct']}{RESET} y una Temperatura de melting {YELL}{beacons_for_obj[obj._position]['tm']:.2f}ºC{RESET} Score tm = {YELL}{beacons_for_obj[obj._position]['tm_score']:.4f}{RESET}")
            primer_data = primers_for_obj.get(obj._position, {})
            if "error" in primer_data:
                print(f"{RED}❌ ERROR: diseño de primers interrumpido{RESET} {primer_data['error']}")
                print(f" Score de primer: {YELL}{0}{RESET}")
            else:
                print(f" Primers diseñados: {MAG}{primer_data['primers']}{RESET}")
                print(f" Score de primer = {YELL}{primer_data['score']:.2f}{RESET}")
            print("-" * 50)

# para ejecutar la función como principal y que no de error

if __name__ == "__main__":
    main()