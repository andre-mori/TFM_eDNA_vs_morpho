# =============================================================
# Comparar especies ausentes en eDNA con base de datos MIDORI
# =============================================================
import os
import pandas as pd
import gzip
base        = "/mnt/c/Users/andre/Desktop/TFM_AZTI/comparison/archivos"
ruta_r      = os.path.join(base, "download_from_r")
ruta_midori = os.path.join(base, "Abundance_tables_Andrea_Metabridge")
ruta_d      = os.path.join(base, "download_from_python")
ruta_ind    = os.path.join(base, "indicaciones_anders")

def resumen_midori(nombre, serie_bool):
    """Imprime cuántos True (en MIDORI) y cuántos False (no) hay en la columna."""
    si = int(serie_bool.sum())
    no = int((~serie_bool).sum())
    print(f"  [{nombre}]  En MIDORI: {si}   |   NO en MIDORI: {no}   (de {len(serie_bool)})")

# --- 1. Cargar lista de especies ausentes en eDNA ---
especies = pd.read_csv(os.path.join(ruta_r, "sp_morpho_compare.csv"), sep=";", decimal=",")
# Extraer solo los nombres de taxón
lista_especies = especies["taxon"].tolist()
print(f"Especies exclusivas de morfología a comprobar: {len(lista_especies)}")
# --- 2. Cargar MIDORI ---
midori_especies = set()
midori_generos  = set()
midori_family     = set()
midori_order    = set()
with gzip.open(os.path.join(ruta_midori, "MIDORI253DARN_nds.tsv.gz"), "rt", encoding="utf-8") as f:
    for linea in f:
        partes = linea.strip().split("\t")
        if len(partes) >= 3:
            nombre = partes[2].strip()
            midori_especies.add(nombre)
            # Extraer el género (primera palabra)
            genero = nombre.split(" ")[0]
            midori_generos.add(genero)
print(f"Especies únicas en MIDORI: {len(midori_especies)}")
print(f"Géneros únicos en MIDORI:  {len(midori_generos)}")
# --- 3. Buscar cada especie en MIDORI ---
resultados = []
for sp in lista_especies:
    sp_limpio = str(sp).replace(".", " ").strip()
    genero = sp_limpio.split(" ")[0]

    en_midori_especie = sp_limpio in midori_especies
    en_midori_genero  = genero in midori_generos

    # Cuántas especies del mismo género hay en MIDORI
    especies_del_genero = [e for e in midori_especies if e.startswith(genero + " ")]

    resultados.append({
        "especie"               : sp_limpio,
        "genero"                : genero,
        "en_MIDORI_especie"     : en_midori_especie,
        "en_MIDORI_genero"      : en_midori_genero,
        "n_especies_genero_MIDORI" : len(especies_del_genero),
        "especies_genero_MIDORI": ", ".join(especies_del_genero)
    })
# --- 4. Guardar resultados ---
df_resultados = pd.DataFrame(resultados)
print(df_resultados)
df_resultados.to_csv(os.path.join(ruta_d, "midori_sp.csv"), index=False, sep=";")
# --- 5. Resumen de los datos ---
# 1. Cuántas especies están en MIDORI
df_resultados["en_MIDORI_especie"].sum()
# 2. Cuántos géneros están en MIDORI
df_resultados["en_MIDORI_genero"].sum()
# 3. De las especies que NO están, cuántas sí tienen el género en MIDORI
((~df_resultados["en_MIDORI_especie"]) & (df_resultados["en_MIDORI_genero"])).sum()
# 4. Especies sin especie NI género en MIDORI
df_resultados.loc[                                                                # loc = seleccionar filas y columnas por etiqueta
    (~df_resultados["en_MIDORI_especie"]) & (~df_resultados["en_MIDORI_genero"]),
    "especie"
]

# ===================== AÑADIDO (recuento apartado ESPECIE / GÉNERO) =====================
print("\n--- Apartado ESPECIE / GÉNERO ---")
resumen_midori("especie", df_resultados["en_MIDORI_especie"])
resumen_midori("genero",  df_resultados["en_MIDORI_genero"])
# =========================================== fin añadido ================================

# --- 5. Comprobar si family está en MIDORI
# Cargar la ruta completa de las especies de morfo
tax_path = pd.read_csv(os.path.join(ruta_midori, "Metabridge_shadow_report_morpho_TAX.tsv"), sep="\t")
tax_path["Species_full"] = tax_path["Genus"] + " " + tax_path["Species"]
tax_cols = tax_path[["Species_full", "Family", "Order", "Class"]].drop_duplicates(subset="Species_full")   # subset = sobre que columna quere,os trabaajr; así drop_duplicates funciona solo si esa columna esta duplicada, no el resto
# df_resultados["especie_limpio"] = df_resultados["especie"].str.replace(".", " ", regex=False)  # mismo arreglo de puntos, para que el merge sí encuentre coincidencias
# Juntamos los resultados de comparar sp+genero con la ruta tax completa de morfo
df_all_tax = df_resultados.merge(tax_cols, left_on="especie", right_on="Species_full", how="left")
df_all_tax = df_all_tax.rename(columns={"especie": "Specie", "genero": "Genus"})
cols_finales = [
    "Specie", "Genus",
    "Class", "Order", "Family",
    "en_MIDORI_especie", "en_MIDORI_genero",
    "n_especies_genero_MIDORI", "especies_genero_MIDORI"
]
df_all_tax = df_all_tax[cols_finales]  # Guarda en df_all_taxa solo las columnas que estan en cols_finales
# Comprobar las familias en MIDORI
midori_taxa = set()
with gzip.open(os.path.join(ruta_midori, "MIDORI253DARN_nds.tsv.gz"), "rt", encoding="utf-8") as f:
    for linea in f:
        partes = linea.strip().split("\t")
        if len(partes) >= 3:
            for nivel in partes[1].split(";"):        # phylum, class, order, family, genus...
                    midori_taxa.add(nivel.strip())

# ===================== AÑADIDO (columna family_en_MIDORI + recuento apartado FAMILIA) ===
df_all_tax["family_en_MIDORI"] = df_all_tax["Family"].apply(
    lambda f: f in midori_taxa if pd.notna(f) else False)
print("\n--- Apartado FAMILIA ---")
resumen_midori("familia", df_all_tax["family_en_MIDORI"])
# =========================================== fin añadido ================================

def estado_familia(fam):
    if pd.isna(fam) or str(fam).strip() == "":
        return "SIN FAMILIA EN LA RUTA"           # la especie no tiene familia asignada
    return "family SÍ en MIDORI" if fam in midori_taxa else "family NO en MIDORI"
sin_gen_sp = df_all_tax[(~df_all_tax["en_MIDORI_especie"]) &
                        (~df_all_tax["en_MIDORI_genero"])].copy()
sin_gen_sp["family_estado"] = sin_gen_sp["Family"].apply(estado_familia)
print(sin_gen_sp[["Specie", "Family", "family_estado"]].to_string(index=False))
print(sin_gen_sp["family_estado"].value_counts())
# --- 7. Guardar resultados ---
df_resultados = pd.DataFrame(resultados)
print(df_resultados)
df_resultados.to_csv(os.path.join(ruta_d, "midori_with_family.csv"), index=False, sep=";")
os.chdir("/mnt/c/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_python")
df_all_tax.to_csv("compare_MIDORI_con_family.csv", index=False, sep=";")
print("Guardado: compare_MIDORI_con_family.csv")
# --- 8.  Comprobar orden ---
# ===================== AÑADIDO (apartado ORDEN: esto es lo que faltaba) =================
# Set de ÓRDENES de MIDORI = nivel justo antes de la familia "...idae" en la ruta
midori_order = set()
with gzip.open(os.path.join(ruta_midori, "MIDORI253DARN_nds.tsv.gz"), "rt", encoding="utf-8") as f:
    for linea in f:
        partes = linea.strip().split("\t")
        if len(partes) >= 3:
            taxonomia = partes[1].split(";")
            for j, nivel in enumerate(taxonomia):
                if nivel.strip().endswith("idae") and j > 0:
                    midori_order.add(taxonomia[j - 1].strip())   # el orden es el nivel anterior a la familia
                    break

# Columna booleana: True si el orden está en MIDORI (False si no, o si no hay orden)
df_all_tax["order_en_MIDORI"] = df_all_tax["Order"].apply(
    lambda o: o in midori_order if pd.notna(o) else False)

print("\n--- Apartado ORDEN ---")
resumen_midori("orden", df_all_tax["order_en_MIDORI"])
print("  (de los NO, sin Orden en la ruta taxonómica:", int(df_all_tax["Order"].isna().sum()), ")")

# El CSV de orden necesita family_en_MIDORI. La recalculamos aquí para que este
# bloque funcione aunque no hayas añadido el recuento de FAMILIA (o si re-ejecutaste
# la línea df_all_tax = df_all_tax[cols_finales] y la perdiste). Usa 'midori_taxa'.
df_all_tax["family_en_MIDORI"] = df_all_tax["Family"].apply(
    lambda f: f in midori_taxa if pd.notna(f) else False)

# Guardar con las columnas ordenadas (ya incluye family_en_MIDORI y order_en_MIDORI)
cols_orden = [
    "Specie", "Genus",
    "Class", "Order", "Family",
    "en_MIDORI_especie", "en_MIDORI_genero", "family_en_MIDORI", "order_en_MIDORI",
    "n_especies_genero_MIDORI", "especies_genero_MIDORI"
]
df_all_tax = df_all_tax[cols_orden]
df_all_tax.to_csv(os.path.join(ruta_d, "compare_MIDORI_con_order.csv"), index=False, sep=";")
print("Guardado: compare_MIDORI_con_order.csv")
# =========================================== fin añadido ================================

# --- COMPROBAR ÍNDICES BIÓTICOS ---


# ¿Cuántos de los taxones ausentes en eDNA son INDICADORES de índices bióticos?
#
# Reunión del 24-07 (Miriam / Anders):
#   - Tabla 5 = de las 632 especies ausentes en eDNA, cuántas son indicadoras;
#               de los 399 géneros ausentes en eDNA, cuántos son indicadores.
#   - Solo NIVEL ESPECIE y NIVEL GÉNERO. La fila de FAMILIA se elimina.
#   - Match EXACTO: "Abra" NO puede casar con "Abra alba".
#
# ENTRADAS (las tres cargas que hacían falta):
#   A) download_from_r/sp_morpho_compare.csv    -> las 632 especies ausentes
#      download_from_r/gen_morpho_compare.csv   -> los 399 géneros ausentes
#   B) indicaciones_anders/actual_indicator_species_NSI_AMBI_only.csv
#      indicaciones_anders/actual_indicator_genus_list_NSI_AMBI_only.csv
#   C) download_from_python/compare_MIDORI_con_order.csv  -> cobertura MIDORI
#
# SALIDAS (en download_from_python/):
#   indicadores_especies_ausentes.csv   listado COMPLETO de las 632
#   indicadores_generos_ausentes.csv    listado COMPLETO de los 399
#   tabla5_indicadores_resumen.csv      la Tabla 5 del TFM
# =============================================================================

# Cargar file especies
sp = pd.read_csv(os.path.join(ruta_r, "sp_morpho_compare.csv"), sep=";", decimal=",")
print("sp_morpho_compare:", sp.shape)


# Quitar el "." que separa las especies
def norm(x):
    return str(x).replace(".", " ").strip().lower()


# Cargar file géneros
gen = pd.read_csv(os.path.join(ruta_r, "gen_morpho_compare.csv"), sep=";", decimal=",")
print("gen_morpho_compare:", gen.shape)


# Cargar listas de indicadores bioticos
ind_sp = pd.read_csv(os.path.join(ruta_ind, "actual_indicator_species_NSI_AMBI_only.csv"), 
                     sep="\t")


ind_gen = pd.read_csv(os.path.join(ruta_ind, "actual_indicator_genus_list_NSI_AMBI_only.csv"), 
                      sep="\t")


# =============================================================================
# EL MATCH — ¿es indicador? (match exacto contra el set)
# =============================================================================
sp_aus["es_indicador"]  = sp_aus["taxon_norm"].isin(ind_sp)
gen_aus["es_indicador"] = gen_aus["taxon_norm"].isin(ind_gen)

# Valores AMBI / NSI. drop_duplicates evita que un taxón duplicado en la lista
# (Callianassa, Chironomus) duplique filas al hacer el merge.
sp_aus = sp_aus.merge(ind_sp.drop_duplicates("key")[["key", "AMBI", "NSI"]],
                      left_on="taxon_norm", right_on="key", how="left").drop(columns="key")
gen_aus = gen_aus.merge(
    ind_gen.drop_duplicates("key")[["key", "AMBI", "NSI.value"]].rename(columns={"NSI.value": "NSI"}),
    left_on="taxon_norm", right_on="key", how="left").drop(columns="key")


# =============================================================================
# CARGA C — COBERTURA DE MIDORI (salida de midori_comparison.py)
# =============================================================================
print("\n[C] COBERTURA MIDORI")
midori = cargar(os.path.join(ruta_d, "compare_MIDORI_con_order.csv"), sep=";")

# R escribe los booleanos como texto "True"/"False". Sin convertirlos, Python
# los trata como cadenas no vacías = siempre True, y los recuentos salen mal.
for c in ["en_MIDORI_especie", "en_MIDORI_genero", "family_en_MIDORI", "order_en_MIDORI"]:
    midori[c] = midori[c].astype(str).str.strip().str.lower().map({"true": True, "false": False})

midori["key_sp"]  = midori["Specie"].map(norm)
midori["key_gen"] = midori["Genus"].map(norm)

sp_aus = sp_aus.merge(
    midori[["key_sp", "Class", "Order", "Family",
            "en_MIDORI_especie", "en_MIDORI_genero", "family_en_MIDORI", "order_en_MIDORI"]],
    left_on="taxon_norm", right_on="key_sp", how="left").drop(columns="key_sp")

# Para géneros hay que colapsar: en compare_MIDORI cada género aparece repetido
# una vez por cada especie suya. groupby + max sobre el booleano = "algún hit".
gen_mid = (midori.groupby("key_gen")
           .agg(Class=("Class", "first"), Order=("Order", "first"), Family=("Family", "first"),
                en_MIDORI_genero=("en_MIDORI_genero", "max"),
                family_en_MIDORI=("family_en_MIDORI", "max"))
           .reset_index())
gen_aus = gen_aus.merge(gen_mid, left_on="taxon_norm", right_on="key_gen",
                        how="left").drop(columns="key_gen")

print(f"  -> especies con info MIDORI: {sp_aus['en_MIDORI_especie'].notna().sum()}/{len(sp_aus)}"
      f" | géneros: {gen_aus['en_MIDORI_genero'].notna().sum()}/{len(gen_aus)}")
print("     (los géneros sin info son los que no tienen ninguna especie entre las 632)")


# =============================================================================
# RESULTADO — TABLA 5
# =============================================================================
def resumen(df, etiqueta):
    n_tot = len(df)
    n_ind = int(df["es_indicador"].sum())
    fila = {"level": etiqueta, "n_ausentes": n_tot, "n_indicador": n_ind,
            "n_no_indic": n_tot - n_ind, "pct_indicador": round(n_ind / n_tot * 100, 1),
            "n_con_AMBI": int(df.loc[df["es_indicador"], "AMBI"].notna().sum()),
            "n_con_NSI":  int(df.loc[df["es_indicador"], "NSI"].notna().sum())}
    print(f"  {etiqueta:8s} ausentes {n_tot:4d} | indicadores {n_ind:4d} "
          f"({fila['pct_indicador']:5.1f}%) | con AMBI {fila['n_con_AMBI']:4d} | con NSI {fila['n_con_NSI']:4d}")
    return fila

print("\n" + "=" * 74)
print("TABLA 5 — Taxones ausentes en eDNA que son indicadores (AMBI / NSI)")
print("=" * 74)
tabla5 = pd.DataFrame([resumen(sp_aus, "Species"), resumen(gen_aus, "Genus")])


# --- Cruce indicador x MIDORI: el argumento del esfuerzo de barcoding --------
print("\n" + "=" * 74)
print("CRUCE INDICADOR x MIDORI (especies)")
print("=" * 74)
print(pd.crosstab(sp_aus["es_indicador"], sp_aus["en_MIDORI_especie"],
                  rownames=["es_indicador"], colnames=["en MIDORI"], dropna=False))
sin_ref = sp_aus[sp_aus["en_MIDORI_especie"] == False]
print(f"\n  De las {len(sin_ref)} especies SIN referencia en MIDORI, "
      f"{int(sin_ref['es_indicador'].sum())} son indicadoras "
      f"({round(sin_ref['es_indicador'].mean()*100, 1)}%)  <- prioridad de barcoding")


# =============================================================================
# GUARDAR — listados COMPLETOS (632 y 399) para el apéndice
# =============================================================================
cols_sp = ["level", "taxon", "Class", "Order", "Family", "mean_morpho",
           "es_indicador", "AMBI", "NSI",
           "en_MIDORI_especie", "en_MIDORI_genero", "family_en_MIDORI", "order_en_MIDORI"]
cols_gen = ["level", "taxon", "Class", "Order", "Family", "mean_morpho",
            "es_indicador", "AMBI", "NSI", "en_MIDORI_genero", "family_en_MIDORI"]

sp_out  = sp_aus[[c for c in cols_sp if c in sp_aus.columns]].sort_values("taxon")
gen_out = gen_aus[[c for c in cols_gen if c in gen_aus.columns]].sort_values("taxon")

sp_out.to_csv(os.path.join(ruta_d, "indicadores_especies_ausentes.csv"), index=False, sep=";")
gen_out.to_csv(os.path.join(ruta_d, "indicadores_generos_ausentes.csv"), index=False, sep=";")
tabla5.to_csv(os.path.join(ruta_d, "tabla5_indicadores_resumen.csv"), index=False, sep=";")

print(f"\nGuardado en download_from_python/:")
print(f"  indicadores_especies_ausentes.csv  ({len(sp_out)} especies, listado completo)")
print(f"  indicadores_generos_ausentes.csv   ({len(gen_out)} géneros, listado completo)")
print(f"  tabla5_indicadores_resumen.csv     (Tabla 5)")


# =============================================================================
# VERIFICACIÓN MANUAL (lo que pidió Miriam: "make absolutely sure of this")
# =============================================================================
# Descomenta para volcar las 632 completas a pantalla y revisarlas una a una.
# pd.set_option("display.max_rows", None)
# print(sp_out[["taxon", "es_indicador", "AMBI", "NSI"]].to_string(index=False))

print("\n--- 15 indicadoras detectadas (comprobar a mano que son correctas) ---")
print(sp_out.loc[sp_out["es_indicador"], ["taxon", "AMBI", "NSI"]].head(15).to_string(index=False))
print("\n--- 15 NO indicadoras (comprobar que efectivamente no están en AMBI/NSI) ---")
print(sp_out.loc[~sp_out["es_indicador"], ["taxon"]].head(15).to_string(index=False))
