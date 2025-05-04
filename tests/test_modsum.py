from columbo_design.modsum import modsum

import sys
print(sys.path)

def test_mod2() -> None:
    """
    Comprueba que la suma en módulo 2 funciona.

    Para ello voy a tomar como argumentos los números 24 y 67.
    """
    assert modsum(24, 67, 2) == 1

def test_mod3() -> None:
    """
    Comprueba que la suma en módulo 3 funciona.

    Para ello voy a tomar como argumentos los números 24 y 67.
    """
    assert modsum(24, 67, 3) == 1
# aqui esta el error
def test_normal_sum() -> None:
    # entiendo que aquí querías decir en módulo normal, cuando no especificas valor del módulo y asume que es 1, siendo == 0 el resultado
    """
    Comprueba que la suma en módulo 3 funciona.

    Para ello voy a tomar como argumentos los números 24 y 67.
    """
    # para corregirlo unicamente tienes que igualarlo a 0
    assert modsum(24, 67) == 0