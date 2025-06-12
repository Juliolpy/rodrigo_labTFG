if __name__ == "__main__":
    import RNA
    import sys
    if len(sys.argv) != 2:
        print("Inserte dos argumentos, de tal manera que:")
        print(": python3 name_script.py [secuencia de ARN]")

    seq = sys.argv[1]

    # creamos el objeto fold compound

    fc = RNA.fold_compound(seq)

    # obtener su structura (ss) y su energía mínima (mfe)
    (ss, mfe) = fc.mfe()

    print(f"{seq}\n{ss} [ {mfe:6.2f} ] kcal/mol")