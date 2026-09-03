# =============================================================================
# CONSTRUIR PATH DE TODOS LOS TAX DE eDNA USANDO LA API DE WoRMS
# Andrea Moriñigo — TFM Bioinformática
# =============================================================================


import pandas as pd
import requests
import time
import os
import re

# =============================================================================
# BLOQUE 1 — RUTAS
# =============================================================================
ruta_in  = "/mnt/c/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_r"   # ruta donde están los CSVs exportados desde R
ruta_out = "/mnt/c/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_python"   # ruta para guardar la tabla TAX final



# =============================================================================
# BLOQUE 2 — CARGAR LOS NOMBRES DE TAXONES DE CADA NIVEL
# =============================================================================
# Exportados con write.csv desde R

class_names   = pd.read_csv(os.path.join(ruta_in, "nombres_cox1_class.csv"))
order_names   = pd.read_csv(os.path.join(ruta_in, "nombres_cox1_order.csv"))
family_names  = pd.read_csv(os.path.join(ruta_in, "nombres_cox1_family.csv"))
genus_names   = pd.read_csv(os.path.join(ruta_in, "nombres_cox1_genus.csv"))
species_names = pd.read_csv(os.path.join(ruta_in, "nombres_cox1_species.csv"))

# Crear una lista única con todos los taxones de todos los niveles
# (sin duplicados porque algunos nombres pueden aparecer en varios niveles)
all_taxa = pd.concat([class_names, order_names, family_names, genus_names, species_names])
all_taxa = all_taxa["taxon"].drop_duplicates().tolist()  

# se ha colado un taxon mal (si sigue estando en la lista se quita
taxon_malo = "Nemertea.NemBar0082..Nemertea.NemBar0082..Nemertea.NemBar0082..order..."
if taxon_malo in all_taxa:
    all_taxa.remove(taxon_malo)

# =============================================================================
# BLOQUE 3 — FUNCIÓN PARA CONSULTAR WORMS
# =============================================================================
# La API de WoRMS tiene dos pasos:
#   Paso 1: buscar el taxón por nombre → devuelve un AphiaID (identificador único)
#   Paso 2: con ese AphiaID, pedir la clasificación completa (Phylum, Class, Order...)
#
# URL del paso 1: https://www.marinespecies.org/rest/AphiaRecordsByName/{nombre}?like=false&marine_only=false
# URL del paso 2: https://www.marinespecies.org/rest/AphiaClassificationByAphiaID/{aphia_id}
#
# IMPORTANTE: hay que poner un time.sleep(0.3) entre llamadas para no saturar la API

def get_taxonomy(taxon_name):
    """
    Consulta WoRMS y devuelve la clasificación taxonómica completa de un taxón.
    Devuelve un diccionario con: Phylum, Class, Order, Family, Genus, Species
    Si no encuentra el taxón, devuelve un diccionario con todos los valores como None
    """
    
    # Resultado vacío por defecto (por si falla la búsqueda)
    resultado = {                               # {}: es para diccionarios
        "taxon":   taxon_name,
        "Phylum":  None,
        "Class":   None,
        "Order":   None,
        "Family":  None,
        "Genus":   None,
        "Species": None
    }
    
    try:                                      # si algo falla en vez de romper el script salta a la siguiente
# --- PASO 1: buscar el AphiaID por nombre ---
        # Los nombres llegan desde R con formato sintactico (puntos en vez de espacios):
        # "Echinocardium.flavescens" -> "Echinocardium flavescens". Sin esta limpieza
        # WoRMS no encuentra ninguna especie y la columna Species sale vacia.
        # Mismo criterio que en midori_comparison.py e indicator_comparison.py.
        nombre_limpio = re.sub(r"\s+", " ", taxon_name.replace(".", " ")).strip()
        nombre_url = nombre_limpio.replace(" ", "%20")
        url_busqueda = f"https://www.marinespecies.org/rest/AphiaRecordsByName/{nombre_url}?like=false&marine_only=false"
        
        respuesta = requests.get(url_busqueda)    # requests.get(url) es la función que hace la consulta a WoRMS (abre el URL)
        
        # Si la respuesta está vacía, devuelve el resultado vacío
        if respuesta.status_code != 200 or not respuesta.json():
            print(f"     NO ENCONTRADO en WoRMS: {taxon_name} (buscado como '{nombre_limpio}')")
            return resultado
        
        # Coger el primer resultado (el más relevante)
        # La respuesta es una lista de diccionarios, coge el primero
        primer_resultado = respuesta.json()[0]
        aphia_id = primer_resultado["AphiaID"]   # extrae el identificador único de WoRMS
        
# --- PASO 2: pedir la clasificación completa con el AphiaID ---
        url_clasificacion = f"https://www.marinespecies.org/rest/AphiaClassificationByAphiaID/{aphia_id}"
        
        respuesta_tax = requests.get(url_clasificacion)
        
        if respuesta_tax.status_code != 200:
            return resultado
        
        clasificacion = respuesta_tax.json()
        
        # La respuesta de clasificación es un árbol anidado con esta estructura:
        # {"rank": "Kingdom", "scientificname": "Animalia", "child": 
        #     {"rank": "Phylum", "scientificname": "Annelida", "child": 
        #         {"rank": "Class", ...}}}
        # Hay que recorrerlo hasta llegar al nivel que nos interesa
        # recorriendo el árbol y extrayendo Phylum, Class, Order, Family, Genus, Species
        # se usa un while loop que vaya bajando por "child" mientras exista
        
        nodo = clasificacion
        while nodo is not None:
            rango = nodo.get("rank", "")
            nombre = nodo.get("scientificname", "")
            
            if rango == "Phylum":
                resultado["Phylum"] = nombre
            elif rango == "Class":
                resultado["Class"] = nombre
            elif rango == "Order":
                resultado["Order"] = nombre
            elif rango == "Family":
                resultado["Family"] = nombre
            elif rango == "Genus":
                resultado["Genus"] = nombre
            elif rango == "Species":
                resultado["Species"] = nombre
            
            # Bajar al siguiente nivel del árbol
            nodo = nodo.get("child", None)
        
    except Exception as e:
        # Si hay cualquier error (conexión, formato...) imprime el error y devuelve vacío
        print(f"  Error con {taxon_name}: {e}")
    
    return resultado



# =============================================================================
# BLOQUE 4 — CONSULTAR LA API PARA TODOS LOS TAXONES
# =============================================================================
# Aquí recorres la lista all_taxa y llamas a get_taxonomy() para cada uno
# Hay que guardar los resultados en una lista y luego convertirla a DataFrame

resultados = []
total = len(all_taxa)

for i, taxon in enumerate(all_taxa):
    print(f"  [{i+1}/{total}] Consultando: {taxon}")
    
    # Llamar a la función que acabas de definir
    tax_info = get_taxonomy(taxon)
    resultados.append(tax_info)
    
    # Pausa entre llamadas para no saturar la API
    # (sin esto WoRMS puede bloquearte temporalmente)
    time.sleep(0.5)   



# =============================================================================
# BLOQUE 5 — GUARDAR LA TABLA TAX
# =============================================================================
# Convertir la lista de diccionarios en un DataFrame
tax_cox1 = pd.DataFrame(resultados)

print(f"\nTaxones consultados: {len(tax_cox1)}")
print(f"Con Class asignada: {tax_cox1['Class'].notna().sum()}")
print(f"Sin Class (no encontrados): {tax_cox1['Class'].isna().sum()}")
print(f"Con Species asignada: {tax_cox1['Species'].notna().sum()}")   # antes del arreglo era 0

# Aviso explicito si algo sigue sin resolverse
sin_nada = tax_cox1[tax_cox1[["Phylum", "Class", "Order", "Family", "Genus", "Species"]].isna().all(axis=1)]
if len(sin_nada) > 0:
    print(f"\nAVISO: {len(sin_nada)} taxones no se resolvieron en WoRMS:")
    print(sin_nada["taxon"].tolist())
print("\nPrimeras filas:")
print(tax_cox1.head(10))

# Guardar como CSV
tax_cox1.to_csv(os.path.join(ruta_out, "tax_cox1.csv"), index=False)
print(f"\nGuardado en: {ruta_out}tax_cox1.csv")
