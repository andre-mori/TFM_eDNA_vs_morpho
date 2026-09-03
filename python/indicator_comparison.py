# =============================================================================
# indicator_comparison.py
# -----------------------------------------------------------------------------
# ¿Cuántos de los taxones ausentes en eDNA son INDICADORES de índices bióticos?
#
# Reunión del 24-07 (Miriam / Anders):
#   - Tabla 5 = de las 632 especies ausentes en eDNA, cuántas son indicadoras;
#               de los 399 géneros ausentes en eDNA, cuántos son indicadores.
#   - Solo NIVEL ESPECIE y NIVEL GÉNERO. La fila de FAMILIA se elimina.
#   - Match EXACTO: "Abra" NO puede casar con "Abra alba".
#
# ENTRADAS:
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

import os
import pandas as pd

# ---- Rutas (detecta sola si estás en WSL o en Windows) ----------------------
CANDIDATAS = [
    "/mnt/c/Users/andre/Desktop/TFM_AZTI/comparison/archivos",   # WSL
    r"C:\Users\andre\Desktop\TFM_AZTI\comparison\archivos",      # Windows
]
base = next((p for p in CANDIDATAS if os.path.isdir(p)), None)
if base is None:
    raise FileNotFoundError(
        "No encuentro la carpeta 'archivos'. Rutas probadas:\n  " + "\n  ".join(CANDIDATAS))
print(f"Base: {base}")

ruta_r   = os.path.join(base, "download_from_r")
ruta_ind = os.path.join(base, "indicaciones_anders")
ruta_d   = os.path.join(base, "download_from_python")


def cargar(ruta, **kw):
    """Lee un CSV y AVISA si el fichero no existe, en vez de fallar en silencio."""
    if not os.path.exists(ruta):
        raise FileNotFoundError(f"No existe el fichero:\n  {ruta}")
    df = pd.read_csv(ruta, **kw)
    print(f"  OK  {os.path.basename(ruta):52s} {df.shape[0]:6d} filas x {df.shape[1]} col.")
    return df


def limpiar(x):
    """Nombre LEGIBLE: quita los puntos que mete R, respetando mayúsculas.

    R sustituye los espacios por puntos al usar los nombres como columnas
    (Acanthocardia.echinata). Esta versión es la que va a las tablas y al TFM.
    """
    return str(x).replace(".", " ").strip()


def norm(x):
    """Clave de COMPARACIÓN: como limpiar() pero en minúsculas.

    Las minúsculas evitan perder un indicador por una mayúscula distinta.
    El match sigue siendo EXACTO sobre el nombre completo: 'Abra' no casa
    con 'Abra alba'. Esta columna es interna, nunca se exporta.
    """
    return limpiar(x).lower()


# =============================================================================
# CARGA A — LAS 632 ESPECIES Y LOS 399 GÉNEROS AUSENTES EN eDNA
# =============================================================================
print("\n[A] TAXONES AUSENTES EN eDNA (salida de comparison_v2.R)")

def cargar_taxones(fichero, nivel, present_in_val, col_abund):
    # R escribe con write.csv2 -> separador ";" y decimal ","
    # present_in_val: "only_morphology" (ausentes en eDNA) u "only_eDNA" (detectados solo por eDNA)
    # col_abund: qué columna de abundancia conservar ("mean_morpho" o "mean_cox1")
    df = cargar(os.path.join(ruta_r, fichero), sep=";", decimal=",")
    df = df[df["present_in"] == present_in_val].copy()
    df["taxon"]      = df["taxon"].map(limpiar)   # nombre legible, sin puntos
    df["taxon_norm"] = df["taxon"].map(norm)      # clave interna de comparación
    df["level"] = nivel
    return df[["level", "taxon", "taxon_norm", col_abund]].reset_index(drop=True)

def cargar_ausentes(fichero, nivel):
    return cargar_taxones(fichero, nivel, "only_morphology", "mean_morpho")

sp_aus  = cargar_ausentes("sp_morpho_compare.csv",  "Species")
gen_aus = cargar_ausentes("gen_morpho_compare.csv", "Genus")

print(f"  -> {len(sp_aus)} especies ausentes | {len(gen_aus)} géneros ausentes")
if len(sp_aus) != 632:
    print(f"  AVISO: esperaba 632 especies y hay {len(sp_aus)}. ¿Has re-ejecutado comparison_v2.R?")
if len(gen_aus) != 399:
    print(f"  AVISO: esperaba 399 géneros y hay {len(gen_aus)}. ¿Has añadido el export de géneros?")

# --- RECÍPROCO: taxones detectados SOLO por eDNA (ausentes en morfología) ----
# [NUEVO] Para la Tabla 4: de las especies/géneros que el eDNA aporta y la
# morfología no ve, cuántos son también indicadores AMBI/NSI. Mismas listas
# de indicadores (actual_indicator_*_NSI_AMBI_only.csv), mismo match exacto.
sp_edna  = cargar_taxones("sp_edna_compare.csv",  "Species", "only_eDNA", "mean_cox1")
gen_edna = cargar_taxones("gen_edna_compare.csv", "Genus",   "only_eDNA", "mean_cox1")

print(f"  -> {len(sp_edna)} especies solo eDNA | {len(gen_edna)} géneros solo eDNA")
if len(sp_edna) != 105:
    print(f"  AVISO: esperaba 105 especies solo-eDNA y hay {len(sp_edna)}. ¿Has re-ejecutado comparison_v2.R?")
if len(gen_edna) != 97:
    print(f"  AVISO: esperaba 97 géneros solo-eDNA y hay {len(gen_edna)}. ¿Has re-ejecutado comparison_v2.R?")


# =============================================================================
# CARGA B — LAS LISTAS DE ÍNDICES BIÓTICOS (AMBI + NSI, Anders 24-07)
# =============================================================================
# OJO: la extensión es .csv pero el separador real es TABULADOR.
# Si usas sep=";" pandas devuelve UNA sola columna y el match da 0.
print("\n[B] LISTAS DE ÍNDICES BIÓTICOS (indicaciones_anders)")

ind_sp  = cargar(os.path.join(ruta_ind, "actual_indicator_species_NSI_AMBI_only.csv"), sep="\t")
ind_gen = cargar(os.path.join(ruta_ind, "actual_indicator_genus_list_NSI_AMBI_only.csv"), sep="\t")

ind_sp["key"]  = ind_sp["Species"].map(norm)
ind_gen["key"] = ind_gen["Genus"].map(norm)

# set() en vez de lista: la búsqueda es instantánea (si no, 632 x 9903 comparaciones)
set_sp_ind  = set(ind_sp["key"])
set_gen_ind = set(ind_gen["key"])
print(f"  -> {len(set_sp_ind)} especies indicadoras únicas | {len(set_gen_ind)} géneros indicadores únicos")

# --- Control de calidad que pidió Anders ("double check that this is complete")
print("\n  Control de calidad de las listas:")
n_pal_sp  = ind_sp["Species"].str.split().str.len()
n_pal_gen = ind_gen["Genus"].str.split().str.len()
print(f"    especies: 2 palabras {(n_pal_sp==2).sum()} | 1 palabra {(n_pal_sp==1).sum()} "
      f"| 3+ {(n_pal_sp>=3).sum()} | duplicadas {ind_sp['key'].duplicated().sum()}")
print(f"    géneros : 1 palabra {(n_pal_gen==1).sum()} | 2+ palabras {(n_pal_gen>=2).sum()} "
      f"| duplicados {ind_gen['key'].duplicated().sum()}")

# Entradas de rango alto (phylum/clase/orden) coladas en la lista de géneros.
# Pasan el filtro de "una sola palabra", así que hay que detectarlas aparte.
rango_alto = {"animalia","polychaeta","oligochaeta","bivalvia","gastropoda","nemertea",
              "sipuncula","nematoda","priapulida","echiura","phoronida","bryozoa",
              "platyhelminthes","turbellaria","anthozoa","porifera","ascidiacea",
              "hemichordata","chaetognatha","ophiuroidea","asteroidea","echinoidea",
              "holothuroidea","amphipoda","isopoda","cumacea","tanaidacea",
              "pycnogonida","harpacticoida"}
colados_gen = sorted(set_gen_ind & rango_alto)
colados_sp  = sorted(set_sp_ind  & rango_alto)
print(f"    entradas de rango alto en la lista de GÉNEROS: {len(colados_gen)}")
if colados_gen:
    print(f"      {', '.join(colados_gen)}")
if colados_sp:
    print(f"    AVISO: rango alto en la lista de ESPECIES: {', '.join(colados_sp)}")

# ¿Alguna de esas entradas sucias produce realmente un falso positivo?
fp_gen = sorted(set(gen_aus["taxon_norm"]) & set(colados_gen))
print(f"    falsos positivos reales por esas entradas: {len(fp_gen)}"
      + (f" -> {', '.join(fp_gen)}" if fp_gen else "  (ninguno: el recuento es fiable)"))


# =============================================================================
# EL MATCH — ¿es indicador? (match exacto contra el set)
# =============================================================================
def marcar_indicador(df, set_ind, ind_df, nsi_col):
    """Añade es_indicador (bool) + columnas AMBI/NSI a df, por match exacto
    sobre taxon_norm contra el set de indicadores. Se reutiliza para los
    taxones ausentes en eDNA y para los detectados solo por eDNA."""
    df = df.copy()
    df["es_indicador"] = df["taxon_norm"].isin(set_ind)
    df = df.merge(
        ind_df.drop_duplicates("key")[["key", "AMBI", nsi_col]].rename(columns={nsi_col: "NSI"}),
        left_on="taxon_norm", right_on="key", how="left").drop(columns="key")
    return df

sp_aus  = marcar_indicador(sp_aus,  set_sp_ind,  ind_sp,  "NSI")
gen_aus = marcar_indicador(gen_aus, set_gen_ind, ind_gen, "NSI.value")

# [NUEVO] mismo cruce para los taxones detectados solo por eDNA
sp_edna  = marcar_indicador(sp_edna,  set_sp_ind,  ind_sp,  "NSI")
gen_edna = marcar_indicador(gen_edna, set_gen_ind, ind_gen, "NSI.value")

print(f"\n  [eDNA] especies indicadoras: {int(sp_edna['es_indicador'].sum())}/{len(sp_edna)}"
      f" | géneros indicadores: {int(gen_edna['es_indicador'].sum())}/{len(gen_edna)}")


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
def resumen(df, etiqueta, direccion):
    # direccion: "only_morphology" (perdidos/ausentes en eDNA) u
    # "only_eDNA" (detectados solo por eDNA, ausentes en morfología) [NUEVO]
    n_tot = len(df)
    n_ind = int(df["es_indicador"].sum())
    fila = {"level":         etiqueta,
            "direccion":     direccion,
            "n_taxones":     n_tot,
            "n_indicador":   n_ind,
            "n_no_indic":    n_tot - n_ind,
            "pct_indicador": round(n_ind / n_tot * 100, 1),
            "n_con_AMBI":    int(df.loc[df["es_indicador"], "AMBI"].notna().sum()),
            "n_con_NSI":     int(df.loc[df["es_indicador"], "NSI"].notna().sum())}
    print(f"  {etiqueta:8s} [{direccion:16s}] n {n_tot:4d} | indicadores {n_ind:4d} "
          f"({fila['pct_indicador']:5.1f}%) | con AMBI {fila['n_con_AMBI']:4d} "
          f"| con NSI {fila['n_con_NSI']:4d}")
    return fila

print("\n" + "=" * 74)
print("TABLA 5 — Taxones indicadores (AMBI / NSI), taxones perdidos y taxones")
print("aportados en exclusiva por cada método")
print("=" * 74)
tabla5 = pd.DataFrame([
    resumen(sp_aus,  "Species", "only_morphology"),
    resumen(gen_aus, "Genus",   "only_morphology"),
    resumen(sp_edna,  "Species", "only_eDNA"),   # [NUEVO] 105 especies solo eDNA
    resumen(gen_edna, "Genus",   "only_eDNA"),   # [NUEVO] 97 géneros solo eDNA
])


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
# taxon_norm NO se exporta: es la clave interna en minúsculas.
cols_sp = ["level", "taxon", "Class", "Order", "Family", "mean_morpho",
           "es_indicador", "AMBI", "NSI",
           "en_MIDORI_especie", "en_MIDORI_genero", "family_en_MIDORI", "order_en_MIDORI"]
cols_gen = ["level", "taxon", "Class", "Order", "Family", "mean_morpho",
            "es_indicador", "AMBI", "NSI", "en_MIDORI_genero", "family_en_MIDORI"]

# El filtro [c for c in cols if c in df.columns] evita que reviente si alguna
# columna no llegó por el merge.
sp_out  = sp_aus[[c for c in cols_sp if c in sp_aus.columns]].sort_values("taxon")
gen_out = gen_aus[[c for c in cols_gen if c in gen_aus.columns]].sort_values("taxon")

sp_out.to_csv(os.path.join(ruta_d, "indicadores_especies_ausentes.csv"), index=False, sep=";")
gen_out.to_csv(os.path.join(ruta_d, "indicadores_generos_ausentes.csv"), index=False, sep=";")

# [NUEVO] listados completos de los detectados solo por eDNA (105 especies, 97 géneros)
cols_sp_edna  = ["level", "taxon", "mean_cox1", "es_indicador", "AMBI", "NSI"]
cols_gen_edna = ["level", "taxon", "mean_cox1", "es_indicador", "AMBI", "NSI"]
sp_edna_out  = sp_edna[[c for c in cols_sp_edna if c in sp_edna.columns]].sort_values("taxon")
gen_edna_out = gen_edna[[c for c in cols_gen_edna if c in gen_edna.columns]].sort_values("taxon")
sp_edna_out.to_csv(os.path.join(ruta_d, "indicadores_especies_solo_eDNA.csv"), index=False, sep=";")
gen_edna_out.to_csv(os.path.join(ruta_d, "indicadores_generos_solo_eDNA.csv"), index=False, sep=";")

tabla5.to_csv(os.path.join(ruta_d, "tabla5_indicadores_resumen.csv"), index=False, sep=";")

print(f"\nGuardado en download_from_python/:")
print(f"  indicadores_especies_ausentes.csv   ({len(sp_out)} especies, listado completo)")
print(f"  indicadores_generos_ausentes.csv    ({len(gen_out)} géneros, listado completo)")
print(f"  indicadores_especies_solo_eDNA.csv  ({len(sp_edna_out)} especies, listado completo)")
print(f"  indicadores_generos_solo_eDNA.csv   ({len(gen_edna_out)} géneros, listado completo)")
print(f"  tabla5_indicadores_resumen.csv      (Tabla 5, con las 4 filas)")


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
