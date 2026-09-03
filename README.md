# Comparación entre morfotaxonomía y metabarcoding de eDNA (COI) en macrofauna bentónica

Código y resultados del Trabajo Fin de Máster de **Andrea Moriñigo Velaz**
(Máster en Bioinformática / Toxicología Ambiental, UNIR — prácticas en AZTI).

El trabajo compara la identificación morfológica clásica y el metabarcoding de ADN
ambiental (marcador COI) sobre las **mismas estaciones de muestreo** de la plataforma
continental noruega, en el contexto del monitoreo regulatorio del impacto de la
extracción petrolífera *offshore*. Los datos proceden del proyecto **MetaBridge**
(NORCE / AZTI, 2018–2022).

## Objetivos

1. Identificar taxones y grupos con diferencias significativas de abundancia entre métodos.
2. Determinar los taxones exclusivos de cada método y las lagunas de cobertura de la base
   de datos de referencia COI (MIDORI253darn).
3. Analizar los factores que explican ausencias y subestimaciones en eDNA.
4. Comparar los índices bióticos (AMBI, NSI, NQI1, ISI) calculados con ambos métodos.

## Estructura del repositorio

```
R/
  comparison_v2.R                  Curación de nombres, normalización a abundancia
                                   relativa (denominador a nivel de clase) y test de
                                   Wilcoxon pareado por nivel taxonómico, con corrección
                                   múltiple (Benjamini-Hochberg y Bonferroni).
  generate_taxonomic_summary_EN.R  Tablas resumen por nivel taxonómico.
  graficos/
    graficos_resultados_TFM_v3_reestructurado.R   Figuras de resultados.
    grafico_dispersion_niveles_clasificacion.R    Dispersión morfología vs. eDNA.

python/
  worms_tax_builder.py             Reconstrucción de las rutas taxonómicas completas
                                   vía la API REST de WoRMS.
  midori_comparison.py             Evalúa si los taxones ausentes en eDNA tienen
  midori_comparison_completo.py    secuencia de referencia en MIDORI253darn.
  indicator_comparison.py          Cruce con los catálogos de indicadores AMBI y NSI.
  denominador_A_vs_B.py            Comparación entre dos estrategias de normalización.

resultados/
  tablas/     Salidas de los scripts (csv y xlsx): resultados del test de Wilcoxon,
              resúmenes por nivel, cobertura de MIDORI, direcciones del sesgo.
  figuras/    Figuras generadas para la memoria.

anexos/       Tablas de apoyo y anexos de la memoria (especies ausentes, familias y
              órdenes sin detección en eDNA, rasgos por taxón, catálogo de indicadores).
```

## Flujo de trabajo

El procesamiento bioinformático previo (demultiplexado con cutadapt, *denoising* con
DADA2, agrupamiento en OTUs con swarm, curación con MUMU y asignación taxonómica con
CREST contra MIDORI253darn) fue realizado por el consorcio MetaBridge. **Este repositorio
parte de las tablas de abundancia ya procesadas**: una tabla de eDNA COI
(1.183 OTUs × 144 muestras) y las tablas morfológicas en cinco niveles taxonómicos
(clase, orden, familia, género y especie).

A partir de ahí:

1. `R/comparison_v2.R` — curación de nombres, eliminación de taxones no válidos,
   normalización y contraste estadístico.
2. `python/worms_tax_builder.py` — reconstrucción de la taxonomía de los OTUs de eDNA.
3. `python/midori_comparison.py` — cobertura de la base de datos de referencia.
4. `python/indicator_comparison.py` — condición de indicador biótico de cada taxón.
5. `R/graficos/` — figuras finales.

## Requisitos

- **R** ≥ 4.4 con `dplyr`, `tidyr`, `readxl`, `writexl`, `ggplot2`.
- **Python** ≥ 3.10 con `pandas` y `requests`.

## Disponibilidad de los datos

Las tablas de abundancia originales del proyecto MetaBridge **no se incluyen** en este
repositorio por no ser datos propios. Están disponibles a través del consorcio
MetaBridge (NORCE / AZTI). Los scripts esperan encontrarlas en una carpeta `datos_brutos/`
en la raíz del proyecto.

## Cita

Si utilizas este código, cita el TFM asociado.
