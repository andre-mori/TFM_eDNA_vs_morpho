#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comparar especies ausentes en eDNA (exclusivas de morfología) con la base de
datos de referencia MIDORI253darn, para ver si su ausencia en eDNA tiene
sentido (no hay secuencia de referencia) o no (debería haberse detectado).
TFM Andrea Moriñigo — MetaBridge

Versión completa con nombres de columnas en_MIDORI_especie / en_MIDORI_genero
(las que has estado usando en tu notebook). El midori_comparison.py original
se ha dejado intacto en esta misma carpeta.
"""

import pandas as pd
import gzip
import os

# =============================================================================
# 1. Cargar lista de especies ausentes en eDNA (exclusivas de morfología)
# =============================================================================
# Este CSV lo genera wilcoxon_comparison.R (bloque "Exportar especies exclusivas de
# morfología para midori_comparison.py"). Regenéralo desde R antes de correr
# este script si has cambiado algo en el pipeline de R.
os.chdir('/mnt/c/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_r')

especies = pd.read_csv("sp_morpho_compare.csv", sep=";", decimal=",")
lista_especies = especies["taxon"].tolist()

print(f"Especies exclusivas de morfología a comprobar: {len(lista_especies)}")

# =============================================================================
# 2. Cargar MIDORI
# =============================================================================
os.chdir('/mnt/c/Users/andre/Desktop/TFM_AZTI/comparison/archivos/Abundance_tables_Andrea_Metabridge')

# Aspecto del archivo MIDORI (separado por tabuladores):
#   MH535936  \t  Cellular organisms;Eukaryota;...;Amoebozoa sp.  \t  Amoebozoa sp.
#   [0]                          [1]                                    [2]
midori_especies = set()
midori_generos  = set()

with gzip.open("MIDORI253DARN_nds.tsv.gz", "rt", encoding="utf-8") as f:
    for linea in f:
        partes = linea.strip().split("\t")          # separado por tabuladores
        if len(partes) >= 3:                          # necesita al menos 3 columnas
            nombre = partes[2].strip()                 # columna [2] = nombre de la especie
            midori_especies.add(nombre)
            genero = nombre.split(" ")[0]               # primera palabra = género
            midori_generos.add(genero)

print(f"Especies en MIDORI: {len(midori_especies)} | Géneros en MIDORI: {len(midori_generos)}")

# =============================================================================
# 3. Buscar cada especie en MIDORI
# =============================================================================
resultados = []

for sp in lista_especies:
    sp_limpio = sp.replace(".", " ")          # "Adontorhina.similis" -> "Adontorhina similis"
    genero = sp_limpio.split(" ")[0]

    en_midori_especie = sp_limpio in midori_especies
    en_midori_genero  = genero in midori_generos
    especies_del_genero = [e for e in midori_especies if e.startswith(genero + " ")]

    resultados.append({
        "especie"                  : sp,
        "genero"                   : genero,
        "en_MIDORI_especie"        : en_midori_especie,
        "en_MIDORI_genero"         : en_midori_genero,
        "n_especies_genero_MIDORI" : len(especies_del_genero),
        "especies_genero_MIDORI"   : ", ".join(especies_del_genero)
    })

df_resultados = pd.DataFrame(resultados)
print(df_resultados.head())

os.chdir("/mnt/c/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_python")
df_resultados.to_csv("compare_MIDORI_genus_sp.csv", index=False, sep=";")
print("Guardado: compare_MIDORI_genus_sp.csv")

# =============================================================================
# 4. Añadir Class / Order / Family (desde el shadow report de morfología)
#    y comprobar la familia en MIDORI
# =============================================================================
os.chdir('/mnt/c/Users/andre/Desktop/TFM_AZTI/comparison/archivos/Abundance_tables_Andrea_Metabridge')
tax = pd.read_csv("Metabridge_shadow_report_morpho_TAX.tsv", sep="\t")

tax["Species_full"] = tax["Genus"] + " " + tax["Species"]
tax_cols = tax[["Species_full", "Family", "Order", "Class"]].drop_duplicates(subset="Species_full")

df_resultados["especie_limpio"] = df_resultados["especie"].str.replace(".", " ", regex=False)  # mismo arreglo de puntos, para que el merge sí encuentre coincidencias
df_all_tax = df_resultados.merge(tax_cols, left_on="especie_limpio", right_on="Species_full", how="left")
df_all_tax = df_all_tax.drop(columns=["Species_full", "especie_limpio"])
df_all_tax = df_all_tax.rename(columns={"especie": "Specie", "genero": "Genus"})

cols_finales = [
    "Specie", "Genus",
    "Class", "Order", "Family",
    "en_MIDORI_especie", "en_MIDORI_genero",
    "n_especies_genero_MIDORI", "especies_genero_MIDORI"
]
df_all_tax = df_all_tax[cols_finales]

# --- Comparación de family con MIDORI ---
# Todas las familias animales terminan en "idae", así que buscamos ese patrón
# en la ruta taxonómica completa (columna [1] de MIDORI).
midori_family = set()

with gzip.open("MIDORI253DARN_nds.tsv.gz", "rt", encoding="utf-8") as f:
    for linea in f:
        partes = linea.strip().split("\t")
        if len(partes) >= 3:
            taxonomia = partes[1].split(";")
            for nivel in taxonomia:
                if nivel.strip().endswith("idae"):
                    midori_family.add(nivel.strip())
                    break

family_en_MIDORI  = []
n_familia_MIDORI  = []
sp_familia_MIDORI = []

for _, fila in df_all_tax.iterrows():
    familia = fila["Family"]
    en_MIDORI = familia in midori_family if pd.notna(familia) else False

    # Especies de esa misma familia (según morfología) que sí están en MIDORI
    generos_familia_tax = tax_cols[tax_cols["Family"] == familia]["Species_full"].tolist()
    generos_en_midori = [sp for sp in generos_familia_tax if sp in midori_especies]

    family_en_MIDORI.append(en_MIDORI)
    n_familia_MIDORI.append(len(generos_en_midori))
    sp_familia_MIDORI.append(", ".join(generos_en_midori))

df_all_tax["family_en_MIDORI"]  = family_en_MIDORI
df_all_tax["n_familia_MIDORI"]  = n_familia_MIDORI
df_all_tax["sp_familia_MIDORI"] = sp_familia_MIDORI

cols_finales = [
    "Specie", "Genus",
    "Class", "Family",
    "en_MIDORI_especie", "en_MIDORI_genero", "family_en_MIDORI",
    "n_especies_genero_MIDORI", "n_familia_MIDORI",
    "especies_genero_MIDORI", "sp_familia_MIDORI"
]
df_all_tax = df_all_tax[cols_finales]

os.chdir("/mnt/c/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_python")
df_all_tax.to_csv("compare_MIDORI_con_family.csv", index=False, sep=";")
print("Guardado: compare_MIDORI_con_family.csv")

# =============================================================================
# 5. Añadir Order y comprobar el orden en MIDORI
# =============================================================================
# NOTA: en el script original este bloque solo imprimía comprobaciones sueltas
# (algunos ejemplos concretos) sin llegar a guardar el resultado. Aquí lo he
# completado siguiendo el mismo patrón que Class/Family, para que quede un
# CSV final también a nivel de orden.
orden_map = tax_cols.set_index("Species_full")["Order"]
df_all_tax["Order"] = df_all_tax["Specie"].str.replace(".", " ", regex=False).map(orden_map)  # mismo arreglo de puntos que en el merge anterior

os.chdir('/mnt/c/Users/andre/Desktop/TFM_AZTI/comparison/archivos/Abundance_tables_Andrea_Metabridge')
midori_order = set()

with gzip.open("MIDORI253DARN_nds.tsv.gz", "rt", encoding="utf-8") as f:
    for linea in f:
        partes = linea.strip().split("\t")
        if len(partes) >= 3:
            taxonomia = partes[1].split(";")
            for j, nivel in enumerate(taxonomia):
                if nivel.strip().endswith("idae") and j > 0:
                    # el orden es el nivel justo antes de la familia
                    orden = taxonomia[j - 1].strip()
                    midori_order.add(orden)
                    break

print(f"Órdenes distintos encontrados en MIDORI: {len(midori_order)}")

df_all_tax["order_en_MIDORI"] = df_all_tax["Order"].apply(
    lambda o: o in midori_order if pd.notna(o) else False
)

cols_finales = [
    "Specie", "Genus",
    "Class", "Family", "Order",
    "en_MIDORI_especie", "en_MIDORI_genero", "family_en_MIDORI", "order_en_MIDORI",
    "n_especies_genero_MIDORI", "n_familia_MIDORI",
    "especies_genero_MIDORI", "sp_familia_MIDORI"
]
df_all_tax = df_all_tax[cols_finales]

os.chdir("/mnt/c/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_python")
df_all_tax.to_csv("compare_MIDORI_con_order.csv", index=False, sep=";")
print("Guardado: compare_MIDORI_con_order.csv")

print(f"\nListo. Especies procesadas: {len(df_resultados)}")
