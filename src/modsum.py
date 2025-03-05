# modsum con el script proporcionado por luksgrin
import sys

def modsum(a: int, b: int, c: int = 1) -> int:
    """
    Realiza la suma en módulo ``c`` de dos números ``a`` y ``b``.

    En caso de que el usuario no proporcione la base ``c``, devuelve la suma habitual.

    :param a: Primer número a sumar.
    :type a: int
    :param b: Segundo número a sumar.
    :type b: int
    :param c: Base para calcular el módulo. Default: 1.
    """
    return (a + b)%c

if __name__ == "__main__":
    a, b, c = map(int, sys.argv[1:4])
    print(
        f"La suma en módulo {c} de {a} y {b} es {modsum(a, b, c)}"
    )