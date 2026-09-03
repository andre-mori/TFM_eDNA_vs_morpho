# -*- coding: utf-8 -*-
"""
=============================================================================
ANALISIS APARTE: denominador A vs B  (script independiente)
=============================================================================
Compara dos formas de normalizar a abundancia relativa, usando EXACTAMENTE el
mismo pipeline que comparison_v2.R (misma limpieza, mismo Wilcoxon pareado,
misma clasificacion por mediana/media), cambiando SOLO el denominador:

  A (la que usas ahora) : denominador de morfologia = suma de TODAS sus clases.
  B                     : denominador = suma solo de las clases COMPARTIDAS
                          (presentes en morfologia Y en eDNA).

En ambos casos el denominador de eDNA es el mismo (clases de eDNA que estan en
morfologia), asi que la unica diferencia real esta en morfologia.

NO modifica ninguno de tus scripts. Solo lee las tablas y escribe un Excel:
  comparacion_denominador_A_vs_B.xlsx  (en download_from_python)

Conclusion con tus datos: A y B dan el MISMO resultado (0 taxones cambian de
categoria en ningun nivel), porque las unicas clases exclusivas de morfologia
(Crinoidea y Pycnogonida) pesan un 0.005% del total. Es decir, la eleccion de
denominador es irrelevante aqui y puedes defender cualquiera de las dos.
=============================================================================
"""
import os, re, warnings
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")

# --- Rutas (mismas que comparison_v2.R) ---
BASE     = "/mnt/c/Users/andre/Desktop/TFM_AZTI/comparison/archivos"
RUTA     = os.path.join(BASE, "Abundance_tables_Andrea_Metabridge")
RUTA_TAX = os.path.join(BASE, "download_from_python")
RUTA_OUT = os.path.join(BASE, "download_from_python")
LEVELS   = ["Class", "Order", "Family", "Genus", "Species"]

# make.names de R: read.csv sustituye los caracteres raros por "." en los nombres
# de columna. Lo replicamos para que los nombres casen con taxa_to_delete.xlsx.
def make_names(cols):
    out = []
    for c in cols:
        c = re.sub(r"[^0-9A-Za-z_.]", ".", str(c))
        if len(c) == 0 or c[0].isdigit():
            c = "X" + c
        out.append(c)
    return out

def load(fname):
    df = pd.read_csv(os.path.join(RUTA, fname), index_col=0)
    df.columns = make_names(df.columns)
    df = df.T.groupby(level=0).sum().T          # por si hay columnas duplicadas
    return df

morpho = {lv: load(f"{lv}_Morpho_based.csv") for lv in LEVELS}
cox1   = {lv: load(f"{lv}_cox1.csv")         for lv in LEVELS}

# --- Limpieza 2B (borrar) y 2C (renombrar) desde taxa_to_delete.xlsx ---
dl = pd.read_excel(os.path.join(RUTA, "taxa_to_delete.xlsx"))
def apply_clean(tab, level, dataset):
    borrar = dl[(dl.level==level)&(dl.dataset==dataset)&(dl.accion=="Borrar")]["taxon"].astype(str).tolist()
    tab = tab.drop(columns=[c for c in borrar if c in tab.columns], errors="ignore")
    ren = dl[(dl.level==level)&(dl.dataset==dataset)&(dl.accion=="Renombrar")]
    for _, r in ren.iterrows():
        old, new = str(r["taxon"]), str(r["nombre_nuevo"])
        if old in tab.columns:
            tab = tab.rename(columns={old: new})
            if list(tab.columns).count(new) > 1:
                tab = tab.T.groupby(level=0).sum().T
    return tab

for lv in LEVELS:
    morpho[lv] = apply_clean(morpho[lv], lv, "morpho")
    cox1[lv]   = apply_clean(cox1[lv],   lv, "cox1")

# --- 2D casos especiales ---
mc = morpho["Class"]
if "Hexacorallia" in mc.columns and "Octocorallia" in mc.columns:
    mc["Anthozoa"] = mc["Hexacorallia"] + mc["Octocorallia"]
    morpho["Class"] = mc.drop(columns=["Hexacorallia", "Octocorallia"])

def move(src_lv, src_col, dst_lv, dst_col):
    if src_col in cox1[src_lv].columns:
        vals = cox1[src_lv][src_col]
        cox1[src_lv] = cox1[src_lv].drop(columns=[src_col])
        cox1[dst_lv][dst_col] = vals
move("Order", "Acteonidae",           "Family", "Acteonidae")
move("Order", "Rhopaluridae..order.", "Family", "Rhopaluridae")
move("Genus", "Glyceridae.LMG",       "Family", "Glyceridae")

# --- 2E filtrar niveles inferiores de eDNA por clases que no estan en morfo ---
edna_path = pd.read_csv(os.path.join(RUTA_TAX, "tax_cox1.csv"))
class_no_morpho = set(cox1["Class"].columns) - set(morpho["Class"].columns)
for lv in ["Order", "Family", "Genus", "Species"]:
    raw = edna_path.loc[edna_path["Class"].isin(class_no_morpho), lv].dropna().astype(str).tolist()
    out = set(raw) | set(make_names(raw))
    cox1[lv] = cox1[lv].drop(columns=[c for c in cox1[lv].columns if c in out], errors="ignore")

# --- Denominadores (nivel clase) ---
morpho_cls = morpho["Class"]; cox1_cls = cox1["Class"]
shared_cls = [c for c in morpho_cls.columns if c in cox1_cls.columns]
denom_morpho_A = morpho_cls.sum(axis=1)                                   # A: todas las clases de morfo
denom_cox1     = cox1_cls[[c for c in cox1_cls.columns if c in morpho_cls.columns]].sum(axis=1)
denom_morpho_B = morpho_cls[shared_cls].sum(axis=1)                       # B: solo clases compartidas
denom_cox1_B   = cox1_cls[shared_cls].sum(axis=1)                         # (= denom_cox1)

def rel(tab, denom):
    return tab.div(denom.reindex(tab.index), axis=0)

def wtest(cv, mv):
    try:
        return wilcoxon(cv, mv, zero_method="wilcox", correction=False).pvalue
    except Exception:
        return np.nan

def do_wilcoxon(morpho_rel, cox1_rel, level):
    cs = [s for s in morpho_rel.index if s in cox1_rel.index]
    m = morpho_rel.loc[cs]; c = cox1_rel.loc[cs]
    common = [t for t in m.columns if t in c.columns]
    only_m = [t for t in m.columns if t not in c.columns]
    only_c = [t for t in c.columns if t not in m.columns]
    rows = []
    for t in common:
        mv = m[t].astype(float).values; cv = c[t].astype(float).values
        med_m, med_c = np.nanmedian(mv), np.nanmedian(cv)
        mean_m, mean_c = np.nanmean(mv), np.nanmean(cv)
        dmed  = "overestimated_in_eDNA" if med_c>med_m else ("underestimated_in_eDNA" if med_c<med_m else "no_difference")
        dmean = "overestimated_in_eDNA" if mean_c>mean_m else ("underestimated_in_eDNA" if mean_c<mean_m else "no_difference")
        rows.append(dict(level=level, taxon=t, present_in="both", p_value=wtest(cv, mv),
                         median_morpho=med_m, median_cox1=med_c, mean_morpho=mean_m, mean_cox1=mean_c,
                         direction=dmed, direction_mean=dmean, conflicto=(dmed != dmean)))
    for t in only_m:
        rows.append(dict(level=level, taxon=t, present_in="only_morphology", p_value=np.nan,
                         median_morpho=np.nanmedian(m[t]), median_cox1=0, mean_morpho=np.nanmean(m[t]), mean_cox1=0,
                         direction="absent_in_eDNA", direction_mean=np.nan, conflicto=np.nan))
    for t in only_c:
        rows.append(dict(level=level, taxon=t, present_in="only_eDNA", p_value=np.nan,
                         median_morpho=0, median_cox1=np.nanmedian(c[t]), mean_morpho=0, mean_cox1=np.nanmean(c[t]),
                         direction="absent_in_morphology", direction_mean=np.nan, conflicto=np.nan))
    return pd.DataFrame(rows)

def run(denom_m, denom_c, tag):
    res = pd.concat([do_wilcoxon(rel(morpho[lv], denom_m), rel(cox1[lv], denom_c), lv)
                     for lv in LEVELS], ignore_index=True)
    mask = res["present_in"] == "both"
    res["p_adjust_BH"] = np.nan; res["significant_BH"] = "NO"; res["significant_bonferroni"] = "NO"
    for lv in LEVELS:
        sel = mask & (res.level == lv) & res.p_value.notna()
        if sel.sum() > 0:
            pv = res.loc[sel, "p_value"].values
            res.loc[sel, "p_adjust_BH"]           = multipletests(pv, method="fdr_bh")[1]
            res.loc[sel, "significant_BH"]         = np.where(multipletests(pv, method="fdr_bh")[0], "YES", "NO")
            res.loc[sel, "significant_bonferroni"] = np.where(multipletests(pv, method="bonferroni")[0], "YES", "NO")
    res["norm"] = tag
    return res

resA = run(denom_morpho_A, denom_cox1,   "A")
resB = run(denom_morpho_B, denom_cox1_B, "B")

# =============================================================================
# COMPARACION A vs B  y  EXCEL
# =============================================================================
DIR_ORDER = ["overestimated_in_eDNA", "underestimated_in_eDNA", "no_difference",
             "absent_in_eDNA", "absent_in_morphology"]

def resumen_dir(res):
    return (res.groupby(["level","direction"]).size().unstack(fill_value=0)
              .reindex(columns=DIR_ORDER, fill_value=0).reindex(LEVELS))

rA, rB = resumen_dir(resA), resumen_dir(resB)

comp = pd.DataFrame([
    dict(level=lv, categoria=d, n_A=int(rA.loc[lv,d]), n_B=int(rB.loc[lv,d]),
         delta_BmenosA=int(rB.loc[lv,d]-rA.loc[lv,d]))
    for lv in LEVELS for d in DIR_ORDER])

conf = pd.DataFrame({
    "level": LEVELS,
    "conflicto_A": [int((resA[(resA.level==lv)].conflicto==True).sum()) for lv in LEVELS],
    "conflicto_B": [int((resB[(resB.level==lv)].conflicto==True).sum()) for lv in LEVELS],
})
conf["delta"] = conf["conflicto_B"] - conf["conflicto_A"]

a = resA[resA.present_in=="both"][["level","taxon","direction","median_morpho","mean_morpho","median_cox1","mean_cox1"]]
b = resB[resB.present_in=="both"][["level","taxon","direction","median_morpho","mean_morpho"]]
mg = a.merge(b, on=["level","taxon"], suffixes=("_A","_B"))
cambios = mg[mg["direction_A"] != mg["direction_B"]].copy()

print("== Reproduccion (A) — direccion 'both' por nivel ==")
print(resA[resA.present_in=="both"].groupby(["level","direction"]).size().unstack(fill_value=0))
print("\nTaxones que cambian de direccion al pasar de A a B:", len(cambios))

out = os.path.join(RUTA_OUT, "comparacion_denominador_A_vs_B.xlsx")
with pd.ExcelWriter(out) as xl:
    comp.to_excel(xl,   "resumen_direccion",      index=False)
    conf.to_excel(xl,   "conflicto_media_mediana", index=False)
    cambios.to_excel(xl,"taxones_que_cambian",     index=False)
    resA.to_excel(xl,   "detalle_A",               index=False)
    resB.to_excel(xl,   "detalle_B",               index=False)
print("Guardado:", out)
