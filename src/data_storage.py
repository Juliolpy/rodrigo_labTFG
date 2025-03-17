import json
import pickle
import time

# datos que vamos a guardar en archivos pickle y json
data = {
    "enzima": "BamHI",
    "sitio": 30,
    "meta": "BamHI es una enzima de restricción de tipo II producida por el microorganismo Bacillus amyloliquefaciens que posee una diana de restricción en el ADN de cadena doble dependiente de una secuencia no metilada, palindrómica y asimétrica, sobre la cual su actividad catalítica hidrolasa genera extremos cohesivos"
}

# hacemos la función que guarde los datos en el archivo json
def save_json(data: dict, file_path : str) -> None: # esta funcion no devuelve nada, solo guarda el archivo

    """
    Args:
    esta funcion va a guardar los datos de el diccionario data en un archivo json
    la funcion coge los datos del diccioanrio y una ruta de archivo donde guardara los datos de json, 
    osea un STRING donde guardara la ruta del archivo

    como no queremos usar datos ni leer datos en este caso, 
    la función solo guarda los datos en un archivo JSON. No es necesario devolver nada, 
    ya que su única tarea es escribir el archivo.

    """

    with open(file_path, "w") as file: # abrimos el archivo file_path para escribir y le asignamos la variable file
        json.dump(data, file, indent= 4) # con indent = 4 hacemos que sea mas legible

        # aqui no ponemos return porque su unica función es escribir un archivo, no calcular ni devolver datos.

    
# hacemos la función que lea el archivo json
def read_json(file_path: str) -> dict:
    """
    Args:
    esta función en cambio lee el archivo json, file_path, la ruta del archvivo JSON
    Return:
    un diccionario con los datos cargados
    
    """

    with open(file_path, "r") as file:
        return json.load(file)





### ahora lo mismo pero con pickle 

def save_pickle(data: dict, file_path : str) -> None: # esta funcion no devuelve nada, solo guarda el archivo

    """
    Args:
    esta funcion va a guardar los datos de el diccionario data en un archivo pickle
    la funcion coge los datos del diccioanrio y una ruta de archivo donde guardara los datos de json, 
    osea un STRING donde guardara la ruta del archivo

    como no queremos usar datos ni leer datos en este caso, 
    la función solo guarda los datos en un archivo JSON. No es necesario devolver nada, 
    ya que su única tarea es escribir el archivo.

    """

    with open(file_path, "wb") as file: # abrimos el archivo file_path para escribir y le asignamos la variable file
        pickle.dump(data, file) # con indent = 4 hacemos que sea mas legible

        # aqui no ponemos return porque su unica función es escribir un archivo, no calcular ni devolver datos.

    
# hacemos la función que lea el archivo json
def read_pickle(file_path: str) -> dict:
    """
    Args:
    esta función en cambio lee el archivo json, file_path, la ruta del archvivo JSON
    Return:
    un diccionario con los datos cargados
    
    """

    with open(file_path, "rb") as file:
        return pickle.load(file)
    

    # la funcion que solo se ejecutara si ejecutamos el script como main

def main() -> None:

    # llamamos a la función de guardar datos de json y pickle
    # JSON
    save_json(data, "data.json")
    print(f"Se estan cargando los datos en json")
    delay = time.sleep(1)
    print(f". . .")
    time.sleep(2)
    print(f"Datos cargados correctamente")

    # pickle
    save_pickle(data, "data.pkl")
    print(f"Se estan cargando los datos en pickle")
    delay = time.sleep(1)
    print(f".  .  .")
    time.sleep(2)
    print(f"Datos cargados correctamente")

    # los cargamos y leemos ambos ya que ambas funciones reciben un str, no un diccionario
    loaded_json = read_json("data.json")
    print(f"Datos de archivo json: {loaded_json}")

    loaded_pickle = read_pickle("data.pkl")
    print(f"Datos del archivo pickle: {loaded_pickle}")

if __name__ == "__main__":
    main()


        

        


