import argparse
# defino la fockin funcion
def get_parse() -> argparse.Namespace:
    # ahora vamos a agregar parsers
    parser = argparse.ArgumentParser(description= "Esto es una calculadora mejorada")

    parser.add_argument("num1", type= int, help= "Primer número de la operación")
    parser.add_argument("num2", type= int, help= "Segundo número de la operación")
    parser.add_argument("--Operation", choices=["suma", "resta"], help= "elige que operación realizar: suma o resta")

    # ahora almacenamos los valores en la función
    return parser.parse_args()


def main() -> None:
    args = get_parse()

    if args.Operation == "suma":
        result = args.num1 + args.num2
    else:
        result = args.num1 - args.num2

# sacamos por pantalla
    print(f"EL resultado de la {args.Operation} es: {result}")

if __name__ == "__main__":
    main()