import pickle
import json

def klkmanin():
    return "Un manguito clásico"

if __name__ == "__main__":
    with open("pickle.pkl", "wb") as file:
        pickle.dump(klkmanin, file)

    with open("jason.json", "w") as file:
        json.dump(klkmanin, file, indent=4)



# no se pueden guardar funciones en JSON
    

