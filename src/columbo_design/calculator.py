import argparse
# defino la fockin funcion
def get_parse() -> argparse.Namespace:
    # ahora vamos a agregar parsers
    parser = argparse.ArgumentParser(description= "Esto es una calculadora mejorada")

    parser.add_argument("num1", type= int, help= "Primer número de la operación")
    parser.add_argument("num2", type= int, help= "Segundo número de la operación")
    parser.add_argument("--operation", choices=["suma", "resta","multiplicacion", "division"], default= "suma", help= "elige que operación realizar: suma o resta & multiplcación y división")

    # ahora almacenamos los valores en la función
    return parser.parse_args()

# esta función la definimos para que solo se ejecute si corremos el codigo manualmente, como principal, sin importarlo
def main() -> None: # no devuelve nada
    
    # reclutamos la funcion get_parse y definimos una variable que nos recopile todos los objetos definidos en ella
    args = get_parse()

    # en vez de usar if statement voy a crear un diccionario con las operaciones:
    operaciones = {
        "suma" : lambda x, y: x + y, # lambda es una forma de crear funciones rápidas y anónimas en una sola línea, sin necesidad de usar def
        "resta" : lambda x, y: x - y,
        "multiplicacion" : lambda x, y: x * y,
        "division" : lambda x, y: x / y if y != 0 else "Error: No puedes dividir por 0 meloncio, da infinito y cosas raras"
    }
    # agrupamos los valores que obtenemos segun lo que pongamos por la linea de comandos, con el default de suma
    result = operaciones[args.operation](args.num1, args.num2) # de esta manera argparse no tomará el valor de None si no le especificamos argumento

    print(f"El resultado de tu {args.operation} es : {result}")

if __name__ == "__main__":
    main()