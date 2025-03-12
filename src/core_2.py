import sys
from Bio import SeqIO
# Aquí van tus otros imports

fastafile = sys.argv[1]

# Aquí van tus funciones
def read_fasta(fastafile: str) -> dict:
   # creamos un diccionario donde guardemos las secuencias
   seq_name = {}
   # nombre y secuencia
   for record in SeqIO.parse(fastafile, "fasta"):
      seq_name[record.id] = str(record.seq)
   return seq_name

def find_NGG_motivs(seq_name: dict) -> dict:
   # diccionario donde guardamos las posiciones de los motivos NGG
   positions = {}

   for seq_id, seq in seq_name.items():
      # esta línea: i + 1 for i in range(len(seq)-2) --> hace que no nos salgamos del indice 
      position_list = [i + 1 for i in range(len(seq)-2) if seq[i + 1 : i + 3] == "GG"] # Si es un motivo NGG en la sencuenia
      positions[seq_id] = position_list
   return positions

def main() -> None:
    # obtener dic con secuencias llamando a la funcion"
    seq_dict = read_fasta(fastafile)
    positions = find_NGG_motivs(seq_dict)
    print(f"Estas son las posiciones de los NGG encontrados en {fastafile}:")
    for k,v in positions.items():
        print(f"{k}: {v}")

if __name__ == "__main__":
    main()