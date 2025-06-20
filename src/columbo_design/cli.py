# Entrypoint file
# añadimos parser que le de al usuario la facilidad de pasar un path al fichero fasta

# El fichero cli.py debe importar las funciones de core.py y darle los argumentos adecuados a las mismas
# para que estas realicen los cálculos, y luego desde cli.py hay que sacar por pantalla los resultados para que los vea el usuario.
import argparse
import time
from Bio.Seq import Seq
# importamos las funciones del codigo
from .core import read_fasta, find_NGG_motivs, process_genome, output_NGG_json, output_NGG_pickle
from .primer_design import primer3_design_columbo, output_pimers_json, output_primers_pickle, score_primers, parse_primers_output
from .beacon import design_beacon, score_beacon, fold_beacon, beacon_tm, melting_temperature, score_energy, hybridation_energy, count_hairpin

# definimos nuestra funcion parser
def get_parse() -> argparse.Namespace:
    """
    Build and return the command-line argument parser.

    Prepares two arguments:
      - `fasta` (positional): path to the input FASTA file.
      - `--output` (optional): one of "json" or "pickle" (default: "json"), 
        specifying desired output format.

    :return: Parsed arguments namespace.
    :rtype: argparse.Namespace
    """

    # variable
    parser = argparse.ArgumentParser(description= "esto es un programa para encontrar motivos NGG en secuencias de ADN en virus")

    # añadimos los argumentos que necesitará nuestro comando de argumentos
    parser.add_argument("fasta", type = str, help="Especique el archivo .fasta a analizar")
    parser.add_argument("--output", choices=["json", "pickle"], default= "json", help= "Especifica el tipo de archivo de salida que deseas obtener")
    # guardamos los valores en la funcion
    return parser.parse_args()

# para ejecutar la función de manera manual

def main() -> None:
    """
    Main entry point for the Columbo CLI tool.

    1. Parses command-line arguments.
    2. Loads sequences from the given FASTA.
    3. Finds NGG PAM sites in each sequence.
    4. Builds ColumboParts objects and ranks top candidates.
    5. Designs molecular beacons and primers around each PAM.
    6. Writes results to JSON or pickle files.
    7. Prints a summary report to stdout.

    :raises ValueError: If flanking window around a PAM goes out of sequence bounds.
    """

    # ANSI colores y el tracRNA que es siempre el mismo
    RED = "\033[91m"
    AZU = "\033[34m"
    CIAN = "\033[96m"
    GREEN = "\033[92m"
    GR = "\033[92m"
    RESET = "\033[0m"
    MAG = "\033[35m"
    YELL = "\033[33m"
    tracrRNA = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU"
    tracrDNA = tracrRNA.upper().replace("U", "T") # Transformamos a DNA para ver complementariedad DNA // DNA

    time_start_total = time.time()

    # reclutamos la función get_parse()
    args = get_parse()

    # asignamos la función de leer el fasta y encontrar el ngg
    sequences =  read_fasta(args.fasta)
    NGG_positions = find_NGG_motivs(sequences)
    top_candidates = process_genome(sequences, NGG_positions) # vemos cuales son los resultados de la clase y los almacenamos en la variable top_candidates

    beacons_for_obj = {}
    primers_for_obj = {}
    #  -- Datos para el benchmarking -- 
    columbo_times = []      # lista de tiempos de cada ColumboPart
    primer_scores = []      # lista de scores obtenidos para cada ColumboPart
    primer_failure = 0     # cuantas ColumboParts no pudieron diseñar primers
    total_parts = 0
    for seq_id, seq in sequences.items():
        for obj in top_candidates:
            part_start = time.time()
            total_parts += 1
            # definir region alrrederor del PAM para diseñar los primers ( x nt upstream y downstream)
            pam_pos = obj._position
            full_seq = sequences[seq_id] # solo la secuencia del dict
            flank = 45
            start = max(0,pam_pos - flank)
            end = min(len(full_seq), pam_pos + 25 + flank)
            if start < 0 or end > len(full_seq):
                raise ValueError("El flanqueo sale de los límites de la secuencia.")
            # objetos para la funcion de primers
            subseq = full_seq[start:end]
            #region_for_primers = str(Seq(obj._primer_region))
            # funcion beacon
            gRNA = str(obj._protospacer[:-3]) + tracrRNA
            gDNA = str(Seq(obj._protospacer[:8]).replace("C","T")) + str(Seq(obj._protospacer[8:-3])) + tracrDNA
            beacon = design_beacon(obj._beacon_site.replace(" ", ""))
            beacon_region = obj._beacon_site[:len(beacon)]
            beacon_struct, beacon_mfe = fold_beacon(beacon)
            stem1_len, loop_len, stem2_len = count_hairpin(beacon_struct)
            hybridation_e = hybridation_energy(beacon, str(Seq(beacon_region).complement()))
            score_hybridation = score_energy(hybridation_e)
            beacon_melting = melting_temperature(beacon)
            beacon_tm_score = beacon_tm(beacon_melting)
            beacon_score, R_bn, R_gr, F_tm, F_e, hp = score_beacon(beacon, str(Seq(obj._beacon_site)), gDNA) # score hecho con DNA
            beacons_for_obj[obj._position] = {
            "beacon":          beacon,
            "beacon_struct":   beacon_struct,
            "beacon_mfe"   :   beacon_mfe,
            "hybridation_e" :  hybridation_e,
            "score_hybridation": score_hybridation,
            "tm":              beacon_melting,
            "tm_score":        beacon_tm_score,
            "score":           beacon_score,
            "R_bn":            R_bn,
            "R_gr":            R_gr,
            "F_tm":            F_tm,
            "F_e" :            F_e,
            "hp":              hp,
            "stem1_len":       stem1_len,
            "stem2_len":       stem2_len,
            "loop_len":        loop_len
            }
            try:
                # diseñar dichos primers
                raw_output = primer3_design_columbo(subseq, pam_pos - start)
                try:
                    parsed_primers = parse_primers_output(raw_output, subseq)
                    score = score_primers(raw_output, pam_pos - start, protospacer_len = 25)
                    primers_for_obj[pam_pos] = {
                        "primers": parsed_primers,
                        "score": round(score, 2),
                        "all_output": raw_output
                    }
                    primer_scores.append(score)  # guardamos el score en la lista
                except ValueError as exc:
                    # Handle the case where no primers are found for this region
                    primers_for_obj[pam_pos] = {"error": str(exc)}
                    primer_failure += 1 # contador le añadimos + 1 cada vez que no pueda diseñar primers
            except Exception as exc:
                # Handle other errors in primer design
                primers_for_obj[pam_pos] = {"error": str(exc)}
                part_end = time.time()
            finally:
                # Aquí garantizamos que SIN EXCEPCIÓN o NO, midamos el tiempo
                part_end = time.time()
                columbo_times.append(part_end - part_start)

    # Al finalizar todos los bucles, medimos tiempo total
    time_end_total = time.time()
    total_time = time_end_total - time_start_total  # <-- Métrica (A)

    # especificamos el archivo que queremos
    if args.output == "pickle":
        pickle_file = output_NGG_pickle(top_candidates)
        pickle_primers = output_primers_pickle(primers_for_obj)
        print(f"💾 Archivos guardados: {GREEN}output_NGG.pkl, output_PRIMERS.pkl{RESET}")
    else:
        json_file = output_NGG_json(top_candidates)
        json_primers = output_pimers_json(primers_for_obj)
        print(f"💾 Archivos guardados: {RED}output_NGG.json, output_PRIMERS.json{RESET}")

    # --- benchmarking ---
    numero_parts = total_parts
    numero_éxitos = numero_parts - primer_failure
    porcentaje_éxitos = (numero_éxitos / numero_parts) * 100 if numero_parts else 0.0
    porcentaje_fallas = 100 - porcentaje_éxitos

    # imprimir resultados del output ColumboParts
    for seq_id, pos in NGG_positions.items():
        print(f"La Secuencia correspondiente con nombre: {YELL}{seq_id}{RESET}, tiene {YELL}{len(pos)}{RESET} motivos NGG en las posiciones:")
        print("  ")
        print(">" * 170)
        print(f"{pos}")
        print(("<" * 170))
        # en este caso no ponemos ningun return porque la funcion main() no la vamos a ejecutar en ningun lado y "and" me iba a devolver solo NGG_positions
        # El operador and en Python devuelve el primer valor falsy que encuentre o, si no hay ninguno, devuelve el último valor evaluado.
        print("\n🔬 Resultados de análisis COLUMBO:\n")
        right_paren = beacon_struct.count(")")
        loop_dot = beacon_struct.count(".")
        left_paren = beacon_struct.count("(")
        for i, obj in enumerate(top_candidates, 1):
            print(f">>>>>> ColumboPart nº {i} <<<<<<")
            print(f"Creada con el archivo: 📁 route -> {args.fasta}")
            print("                                           ")
            print(f" Region visual: ")
            print(f" {RED}{obj._region}{RESET}")
            print(f" {GREEN}{str(Seq(obj._region).complement())}{RESET}")
            print("                                           ")
            print(f" PAM: {YELL}3'-- {RESET}{GREEN}{obj._pam}{RESET} {YELL}--5'{RESET}")
            print(f"                                                      ")
            print(f"Región donde hibrida el Beacon: {YELL}5'--{RESET}{RED}{obj._beacon_site}{RESET}{YELL}--3'{RESET}")
            print(f" Protospacer:             {YELL}5'--{RESET} {RED}{obj._protospacer}{RESET} {YELL}--3'{RESET}  de longitud {YELL}{len(obj._protospacer)} nt {RESET} ")
            print(f" Hebra complementaria:    {YELL}3'--{RESET} {GREEN}{str(Seq(obj._protospacer).complement())}{RESET} {YELL}--5'{RESET} --> Hebra complementaria del Protospacer" )
            print("                                                          ")
            print(f" Localización del protospacer en el genoma: {YELL}{obj._position}{RESET}")
            print("                                           ")
            print(f" Temperatura de melting (Tm): {YELL}{obj._tm:.2f}°C{RESET}")
            print(f" Scores individuales = {YELL}{obj._scores}{RESET}")
            print(f" Score global protospacer = {YELL}{obj._score_medio:.2f}{RESET}")
            print(f" TracrRNA Streptococcus pyogenes Cas9 {CIAN}{tracrRNA}{RESET} de longitud = {YELL}{len(tracrRNA)} nt{RESET}")
            print(f" Guide RNA (gRNA) resultante {GR}5'--{RESET}{RED}{str(Seq(obj._protospacer[:8]).replace("T", "U").replace("C","U"))}{RESET}{RED}{str(Seq(obj._protospacer[8:-3]).replace("T", "U"))}{RESET}{CIAN}{tracrRNA}{RESET}{GR}--3'{RESET} de longitud = {YELL}{len(tracrRNA)+len(obj._protospacer)} nt{RESET}")
            print("                                                          ")
            print(f" Diseño de Beacon: {GREEN}{beacons_for_obj[obj._position]['beacon']}{RESET} con Score Global = {YELL}{beacons_for_obj[obj._position]['score']:.2f}{RESET}")
            print(f" RNAfold hairping generado {GREEN}{beacons_for_obj[obj._position]['beacon_struct']} --energía libre estrucutural--> {RESET}{GREEN}{beacons_for_obj[obj._position]['beacon_mfe']:.4f} kcal/mol{RESET} y una Temperatura de melting {YELL}{beacons_for_obj[obj._position]['tm']:.2f}ºC{RESET}")
            print(f" Score structural Beacon: {YELL}{beacons_for_obj[obj._position]['hp']:.2f}{RESET}")
            print(f" Score tm = {YELL}{beacons_for_obj[obj._position]['tm_score']:.4f}{RESET}")
            print(f" Unión Beacon con Target mfe score = {YELL}{beacons_for_obj[obj._position]['F_e']:.4f}{RESET}")
            print(f" Score de hibridation (mfe) --> comprobante = {YELL}{beacons_for_obj[obj._position]['score_hybridation']:.4f}{RESET}")
            print(f" Score R_beacon = {YELL}{beacons_for_obj[obj._position]['R_bn']:.2f}{RESET}, Score R_guide = {YELL}{beacons_for_obj[obj._position]['R_gr']:.2f}{RESET}")
            print("                                                          ")
            print(f" Estructura stem -> {GREEN}( <-- {RESET} {YELL}{beacons_for_obj[obj._position]['stem1_len']}{RESET} se junta con {GREEN} --> ){RESET} {YELL}{beacons_for_obj[obj._position]['stem2_len']}{RESET} con un loop de puntos {GREEN} --> .{RESET} {YELL}{beacons_for_obj[obj._position]['loop_len']}{RESET}                        Beacon --> {GREEN}{beacons_for_obj[obj._position]['beacon']}{RESET}")
            print(f" Energía libre de hibridación = {GREEN}{beacons_for_obj[obj._position]['hybridation_e']:.4f} kcal/mol{RESET}, con la hebra complementaria al protospacer(HEBRA DESPLAZADA): {YELL}5'--{RESET} {AZU}{str(Seq(obj._beacon_site))}{RESET} {YELL}--3'{RESET}")
            print(f"                                                                                                                      {AZU}{obj._beacon_site[:32 ]}{RESET}")
            primer_data = primers_for_obj.get(obj._position, {})
            if "error" in primer_data:
                print(f"{RED}❌ ERROR: diseño de primers interrumpido{RESET} {primer_data['error']}")
                print(f" Score de primer: {YELL}{0}{RESET}")
            else:
                print(f" Primers diseñados: {MAG}{primer_data['primers']}{RESET}")
                print(f" Score de primers = {YELL}{primer_data['score']:.2f}{RESET}")
            print("-" * 50)

    # imprimir resultados
    print("\n\n===== RESUMEN DE BENCHMARKING =====")
    print(f" --> creado con : {GREEN}{args.fasta}{RESET} <--")
    print(f" Correspondiente a la secuencia : {YELL}{seq_id}{RESET}")
    print(f"• Tiempo total de ejecución: \t{total_time:.2f} segundos")
    print(f"• Total ColumboParts analizadas: \t{numero_parts}")
    print(f"• ColumboParts con primers válidos: {GREEN}\t{numero_éxitos} ({porcentaje_éxitos:.1f} %){RESET}")
    print(f"• ColumboParts que fallaron: {RED}\t{primer_failure} ({porcentaje_fallas:.1f} %){RESET}")
    if columbo_times:
        print(f"• Tiempo medio por ColumboPart:\t{sum(columbo_times)/numero_parts:.3f} s")
        print(f"• Tiempo mínimo:\t{min(columbo_times):.3f} s")
        print(f"• Tiempo máximo:\t{max(columbo_times):.3f} s")

    # --- (Opcional) Volcar lista de scores y tiempos en disco para análisis posterior ---
    # Podrías hacer algo como:
    with open(f"benchmark_{seq_id}_scores.txt", "w") as f:
        for s in primer_scores:
            f.write(f"{s:.4f}\n")
    with open(f"benchmark_{seq_id}_times.txt", "w") as f:
        for t in columbo_times:
            f.write(f"{t:.4f}\n")

    # --- Generamos un gráfico desde la CLI ---
    # Usando matplotlib
    try:
        import matplotlib.pyplot as plt

        # 1) Histograma de scores de primers
        plt.figure(figsize=(6,4))
        plt.hist(primer_scores, bins=20, edgecolor='black')
        plt.title("Distribución de Scores de Primers")
        plt.xlabel("Score")
        plt.ylabel("Número de ColumboParts")
        plt.tight_layout()
        plt.savefig(f"histograma_{seq_id}_scores.png")
        plt.close()

        # 2) Histograma de tiempos por ColumboPart
        plt.figure(figsize=(6,4))
        plt.hist(columbo_times, bins=20, edgecolor='black')
        plt.title("Distribución de Tiempos por ColumboPart")
        plt.xlabel("Tiempo (s)")
        plt.ylabel("Cantidad de ColumboParts")
        plt.tight_layout()
        plt.savefig(f"histograma_{seq_id}_tiempos.png")
        plt.close()

        print(f"✅ Gráficos guardados: {YELL}histograma_{seq_id}_scores.png, histograma_{seq_id}_tiempos.png{RESET}")
    except ImportError:
        print(f"⚠️ {RED}No se pudo generar gráficos {RESET}: {YELL}falta instalar matplotlib.{RESET}")

    print("=================================\n\n")

# para ejecutar la función como principal y que no de error

if __name__ == "__main__":
    main()