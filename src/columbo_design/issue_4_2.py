import RNA
import sys

def vienna_RNA(RNA_seq: str) -> dict:
    """
    Función que analiza una secuencia de ARN y retorna su estructura secundaria y energía mínima.
    
    Parámetros:
        RNA_seq (str): Secuencia de ARN en formato string.
    
    Retorno:
        dict: Diccionario con la estructura secundaria y la energía mínima de la secuencia.
    """
    
    # Creamos el objeto fold compound
    fc = RNA.fold_compound(RNA_seq)
    
    # Obtener estructura secundaria (ss) y energía mínima (mfe)
    ss, mfe = fc.mfe()
    
    return {
        "secondary_structure": ss,
        "min_energy": mfe
    }

seq = sys.argv[1]  # Obtenemos la secuencia de ARN desde la terminal
resultado = vienna_RNA(seq)

# Imprimimos el resultado
print(f"Secuencia: {seq}")
print(f"Estructura secundaria: {resultado['secondary_structure']}")
print(f"Energía mínima: {resultado['min_energy']:.2f} kcal/mol")