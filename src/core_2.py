import sys
# Aquí van tus otros imports

fasta_file = sys.argv[1]

# Aquí van tus funciones

def main() -> None:
    positions = find_NGG_motivs(read_fasta(sys.argv[1]))
    print(f"Estas son las posiciones de los NGG encontrados en {sys.argv[1]}:")
    for k,v in positions.items():
        print(f"{k}: {v}")

if __name__ == "__main__":
    main()