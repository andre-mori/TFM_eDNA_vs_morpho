# =============================================================================
# RESULTADOS TFM — Andrea Moriñigo  ·  SCRIPT REESTRUCTURADO (v3)
# Morfología clásica vs. metabarcoding de eDNA (COI) — MetaBridge
# -----------------------------------------------------------------------------
# Genera TODAS las figuras y tablas de los Resultados, ordenadas siguiendo la
# lógica de la memoria,
# en tres bloques:
#
#   BLOQUE 1 — Panorama general: cuántos taxones se comparten y cuántos no.
#              (Venn de clase · tabla de exclusivos de eDNA · barplot sin clase)
#   BLOQUE 2 — Taxones exclusivos de un método: causa de las ausencias.
#              (cobertura MIDORI por nivel · cuántos de los perdidos son AMBI)
#   BLOQUE 3 — Taxones compartidos: diferencias significativas de abundancia.
#              (tabla de significativos · scatter de clase · tabla de Miriam)
#
# Nomenclatura de salida: figX_Y_*.png / tabX_Y_*.png  (X = bloque, Y = orden)
# =============================================================================

library(ggplot2); library(dplyr); library(tidyr)
library(readr);   library(scales); library(ggrepel)
library(gt);      library(eulerr);  library(magick)

# ---- PALETA Y TEMA (idénticos a la galería previa) --------------------------
teal  <- "#035F72"; teal_p  <- "#7FB3C0"
rose  <- "#D77186"; rose_p  <- "#EBB6C1"
amber <- "#E69F00"; amber_p <- "#F2D08A"
sage  <- "#9CB380"; sage_p  <- "#C5D6B0"
grey_n<- "#B0B0B0"; ink     <- "#37474F"

tema_tfm <- theme_minimal(base_size = 13) +
  theme(
    text                = element_text(color = ink),
    plot.title          = element_text(face = "bold", size = 15, hjust = 0, margin = margin(b = 12)),
    plot.title.position = "plot",
    plot.subtitle       = element_text(size = 10.5, color = "grey45", margin = margin(b = 8)),
    axis.title          = element_text(size = 11.5),
    axis.text           = element_text(color = ink, size = 11),
    panel.grid.minor    = element_blank(),
    panel.grid.major    = element_blank(),
    axis.ticks          = element_blank(),
    legend.position     = "top",
    legend.title        = element_blank(),
    legend.text         = element_text(size = 10.5),
    plot.margin         = margin(16, 18, 14, 14)
  )

# Estilo común para las tablas gt (evita repetir el bloque en cada tabla)
estilo_gt <- function(g, color_cab = teal) {
  g %>%
    tab_style(style = cell_fill(color = color_cab),
              locations = cells_column_labels()) %>%
    tab_style(style = cell_text(color = "white", weight = "bold"),
              locations = cells_column_labels()) %>%
    tab_options(table.font.name = "Arial", table.font.size = px(13),
                heading.title.font.size = px(16), column_labels.padding = px(8),
                data_row.padding = px(6), table.border.top.color = "white",
                table_body.hlines.color = "#e5e5e5")
}

niveles     <- c("Class","Order","Family","Genus","Species")
niveles_sin_clase <- c("Order","Family","Genus","Species")
etiq_es     <- c(Class="Clase", Order="Orden", Family="Familia",
                 Genus="Género", Species="Especie")

# ---- RUTAS ------------------------------------------------------------------
ruta      <- "C:/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_r/"
ruta_py   <- "C:/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_python/"
ruta_tabs <- "C:/Users/andre/Desktop/TFM_AZTI/comparison/archivos/Abundance_tables_Andrea_Metabridge/"
ruta_bio  <- "C:/Users/andre/Desktop/TFM_AZTI/comparison/archivos/indicaciones_anders/"
ruta_d    <- "C:/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_r/graphics/"

# ---- LECTURA DE DATOS -------------------------------------------------------
sumario     <- read_delim(paste0(ruta, "summary_percentage_CLASS.csv"), delim = ";",
                          locale = locale(decimal_mark = ","), show_col_types = FALSE)
wilcox      <- read_delim(paste0(ruta, "wilcoxon_denominador_CLASS.csv"), delim = ";",
                          locale = locale(decimal_mark = ","), show_col_types = FALSE)
solo_edna   <- read_delim(paste0(ruta, "class_only_eDNA.csv"), delim = ";",
                          locale = locale(decimal_mark = ","), show_col_types = FALSE)
midori_sp   <- read.csv2(paste0(ruta_py, "compare_MIDORI_con_order.csv"),
                         sep = ";", stringsAsFactors = FALSE,
                         fileEncoding = "UTF-8")

# #############################################################################
# BLOQUE 1 — PANORAMA GENERAL: TAXONES COMPARTIDOS VS EXCLUSIVOS
# #############################################################################

# =============================================================================
# FIG 1.1 — Diagrama de Venn (Euler) a NIVEL DE CLASE
# -----------------------------------------------------------------------------
# Es la primera imagen del bloque: fija de un vistazo el panorama a nivel de clase
# (22 compartidas / 23 solo eDNA / 2 solo morfología).
# =============================================================================
euler_taxa <- function(nivel, titulo, archivo_salida, shape_euler = "ellipse") {
  # Conjuntos a partir de los datos YA PROCESADOS (wilcox), NO de las tablas
  # crudas Class_cox1.csv / Class_Morpho_based.csv. Así el diagrama coincide con
  # el resto del análisis (limpieza, renombrados, fusión Hexacorallia+Octocorallia
  # y filtrado de clases de wilcoxon_comparison.R). Leyendo las tablas crudas salían los
  # números viejos (a nivel de clase 20/39/6 en vez de 22/23/2).
  sub        <- wilcox[wilcox$level == nivel, ]
  taxa_both  <- sub$taxon[sub$present_in == "both"]
  set_cox1   <- c(taxa_both, sub$taxon[sub$present_in == "only_eDNA"])
  set_morpho <- c(taxa_both, sub$taxon[sub$present_in == "only_morphology"])
  fit <- euler(list("eDNA COI" = set_cox1, "Morfología" = set_morpho), shape = shape_euler)
  png(archivo_salida, width = 6, height = 6, units = "in", res = 300, bg = "white")
  print(plot(fit,
             fills      = list(fill = c(amber_p, teal_p), alpha = 0.75),
             edges      = list(col = c(amber, teal), lwd = 1.5),
             labels     = list(col = ink, font = 2, cex = 1),
             quantities = list(col = ink, cex = 1.1),
             main       = list(label = titulo, fontsize = 15, fontface = "bold", col = ink)))
  dev.off()
}
euler_taxa("Class", "Clase",
           paste0(ruta_d, "fig1_1_venn_clase.png"))

# ACTUALIZADO: la Figura 3 del TFM (composición clase/familia/especie) pasa a
# ser ORDEN/familia/especie. El nivel de clase ya se describe íntegramente en
# el texto justo antes de la figura (47 clases: 22 en ambos, 23 solo eDNA,
# 2 solo morfología), así que no hace falta repetirlo en la imagen; se
# sustituye por orden para dar un escalón intermedio entre clase y familia.
euler_taxa("Order",   "Orden",   paste0(ruta_d, "sup_venn_orden.png"),   shape_euler = "circle")
euler_taxa("Family",  "Familia", paste0(ruta_d, "sup_venn_familia.png"), shape_euler = "circle")
euler_taxa("Species", "Especie", paste0(ruta_d, "sup_venn_especie.png"), shape_euler = "circle")
imgs <- image_read(c(paste0(ruta_d, "sup_venn_orden.png"),
                     paste0(ruta_d, "sup_venn_familia.png"),
                     paste0(ruta_d, "sup_venn_especie.png")))
image_write(image_append(imgs), paste0(ruta_d, "sup_venn_por_nivel.png"))

print(imgs)

# =============================================================================
# TAB 1.2 — Clases EXCLUSIVAS del eDNA ("unique in DNA")  [antes faltaba como fig]
# -----------------------------------------------------------------------------
# Tabla completa (las 23 clases), no solo el top-10: es la que faltaba por poner
# como figura. Se ordena por abundancia relativa media y se añade el % de reads.
# =============================================================================
tab1_2 <- solo_edna %>%
  arrange(desc(mean_cox1)) %>%
  mutate(pct = mean_cox1 * 100) %>%
  select(taxon, pct)

gt1_2 <- tab1_2 %>%
  gt() %>%
  tab_header(title = md("**Clases detectadas únicamente por el metabarcoding (eDNA)**"),
             subtitle = paste0("n = ", nrow(tab1_2), " clases exclusivas de eDNA · grupos no considerados por la morfología")) %>%
  cols_label(taxon = "Clase", pct = "Abund. rel. media (%)") %>%
  fmt_number(columns = pct, decimals = 2) %>%
  data_color(columns = pct,
             fn = scales::col_numeric(palette = c("#fdf0f3", rose_p), domain = NULL)) %>%
  tab_style(style = cell_text(weight = "bold", style = "italic"),
            locations = cells_body(columns = taxon)) %>%
  tab_source_note(md("*Se retiran antes de calcular las abundancias relativas (denominador de clase); por eso no entran en la comparación cuantitativa del bloque 3.*")) %>%
  estilo_gt(color_cab = rose)

print(gt1_2)
gtsave(gt1_2, paste0(ruta_d, "tab1_2_clases_solo_edna.png"), expand = 10)

# =============================================================================
# FIG 1.3 — Barras apiladas: compartidos vs exclusivos  ·  SIN NIVEL DE CLASE
# -----------------------------------------------------------------------------
# El nivel de clase ya se cuenta en el Venn (Fig 1.1), así que aquí se elimina
# para no repetir y para que el gráfico arranque en ORDEN (evita la confusión
# que salió en la reunión al mezclar el dataset filtrado y sin filtrar).
# =============================================================================
df1_3 <- sumario %>%
  filter(level %in% niveles_sin_clase) %>%
  mutate(level = factor(level, levels = niveles_sin_clase),
         present_in = recode(present_in, both = "Compartidos",
                             only_morphology = "Solo morfología",
                             only_eDNA = "Solo eDNA"),
         present_in = factor(present_in,
                             levels = c("Compartidos","Solo morfología","Solo eDNA")))

fig1_3 <- ggplot(df1_3, aes(level, n, fill = present_in)) +
  geom_col(aes(color = present_in), width = 0.66, linewidth = 0.6) +
  geom_text(aes(label = ifelse(n > 1, n, "")),
            position = position_stack(vjust = 0.5),
            fontface = "bold", size = 3.4, show.legend = FALSE) +
  scale_x_discrete(labels = etiq_es) +
  scale_y_continuous(breaks = seq(0, 900, by = 100), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Compartidos"=rose_p, "Solo morfología"=teal_p, "Solo eDNA"=amber_p)) +
  scale_color_manual(values = c("Compartidos"=rose, "Solo morfología"=teal, "Solo eDNA"=amber),
                     guide = "none") +
  labs(title = "Taxones compartidos y exclusivos por nivel taxonómico",
       subtitle = "De orden a especie (el nivel de clase se muestra en el diagrama de Venn)",
       x = NULL, y = "Número de taxones") +
  tema_tfm

print(fig1_3)
ggsave(paste0(ruta_d, "fig1_3_compartidos_exclusivos.png"), fig1_3,
       width = 8.4, height = 5.2, dpi = 300, bg = "white")

# (Complementario) % de taxones compartidos según resolución (misma info, en línea)
df_pct <- sumario %>% filter(present_in == "both", level %in% niveles_sin_clase) %>%
  mutate(level = factor(level, levels = niveles_sin_clase), idx = as.integer(level))
fig1_3b <- ggplot(df_pct, aes(idx, porcentaje)) +
  geom_area(fill = sage_p, alpha = 0.35) +
  geom_line(color = sage, linewidth = 1.3) +
  geom_point(shape = 21, fill = "white", color = sage, stroke = 1.6, size = 4.5) +
  geom_text(aes(label = paste0(porcentaje, "%")), vjust = -1.4,
            fontface = "bold", size = 3.8, color = teal) +
  scale_x_continuous(breaks = 1:4, labels = etiq_es[niveles_sin_clase]) +
  ylim(0, 100) +
  labs(title = "Concordancia entre métodos según la resolución taxonómica",
       x = NULL, y = "% de taxones compartidos") +
  tema_tfm

print(fig1_3b)
ggsave(paste0(ruta_d, "sup_pct_compartidos.png"), fig1_3b,
       width = 8, height = 5, dpi = 300, bg = "white")

# #############################################################################
# BLOQUE 2 — TAXONES EXCLUSIVOS DE UN MÉTODO: MIDORI Y PÉRDIDA DE INDICADORES
# #############################################################################

# --- Preparación MIDORI (nivel más fino con referencia por especie) ----------
midori_sp <- midori_sp %>%
  mutate(across(c(en_MIDORI_especie, en_MIDORI_genero, family_en_MIDORI, order_en_MIDORI), as.logical),
         primer_nivel_con_ref = case_when(
           en_MIDORI_especie ~ "Especie",
           en_MIDORI_genero  ~ "Género",
           family_en_MIDORI  ~ "Familia",
           order_en_MIDORI   ~ "Orden",
           TRUE              ~ "Clase"),
         primer_nivel_con_ref = factor(primer_nivel_con_ref,
                                       levels = c("Especie","Género","Familia","Orden","Clase")))

resumen_midori <- midori_sp %>%
  count(primer_nivel_con_ref, name = "n", .drop = FALSE) %>%
  mutate(pct = round(n / sum(n) * 100, 1))
print(resumen_midori)   # debe sumar nrow(midori_sp) y traer las 5 categorías

# =============================================================================
# TAB 2.1 — Cobertura de MIDORI a distintos niveles taxonómicos
# -----------------------------------------------------------------------------
# Para cada especie exclusiva de morfología, nivel más fino con referencia en
# MIDORI253darn. Distingue "no hay referencia" de "hay pero solo a rango alto".
# =============================================================================
tab2_1 <- resumen_midori %>%
  rename(nivel = primer_nivel_con_ref, n_especies = n) %>%
  bind_rows(tibble(nivel = "Total",
                   n_especies = sum(resumen_midori$n), pct = 100))

gt2_1 <- tab2_1 %>%
  gt() %>%
  tab_header(title = md("**Cobertura de MIDORI253darn para las especies exclusivas de morfología**"),
             subtitle = paste0("n = ", nrow(midori_sp), " especies no recuperadas por eDNA")) %>%
  cols_label(nivel = "Nivel más fino con referencia", n_especies = "Nº especies", pct = "% del total") %>%
  fmt_number(columns = pct, decimals = 1) %>%
  data_color(columns = pct,
             fn = scales::col_numeric(palette = c("#f7fbff", teal_p), domain = c(0, 100))) %>%
  tab_style(style = cell_fill(color = "#f2f2f2"),
            locations = cells_body(rows = nivel == "Total")) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = nivel, rows = nivel == "Total")) %>%
  tab_source_note(md("*Fuente: elaboración propia (midori_comparison.py) sobre MIDORI253darn.*")) %>%
  estilo_gt(color_cab = teal)

print(gt2_1)
gtsave(gt2_1, paste0(ruta_d, "tab2_1_cobertura_midori.png"), expand = 10)

# (Complementario) misma información en barras
colores_midori <- c("Especie"=teal, "Género"=sage, "Familia"=amber, "Orden"=rose, "Clase"=grey_n)
fig2_1b <- ggplot(resumen_midori, aes(primer_nivel_con_ref, n, fill = primer_nivel_con_ref)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.6) +
  geom_text(aes(label = paste0(n, "\n(", pct, "%)")), vjust = -0.3,
            fontface = "bold", size = 3.6, color = ink) +
  scale_fill_manual(values = colores_midori, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Nivel taxonómico más fino con referencia en MIDORI",
       subtitle = paste0("Especies exclusivas de morfología (n = ", nrow(midori_sp), ")"),
       x = NULL, y = "Nº de especies") +
  tema_tfm

print(fig2_1b)
ggsave(paste0(ruta_d, "sup_cobertura_midori_barras.png"), fig2_1b,
       width = 8, height = 5.5, dpi = 300, bg = "white")

# =============================================================================
# FIG/TAB 2.2 — ¿Cuántos de los taxones PERDIDOS son indicadores AMBI? [NUEVO]
# -----------------------------------------------------------------------------
# Cruza los taxones exclusivos de morfología (table4_missing_in_eDNA_bottomup.csv)
# con la lista de taxones usados en índices bióticos (Taxa_in_biotic_indices.csv).
# Cuantifica la pérdida de valor para el biomonitoreo regulatorio: cuántas de las
# especies/géneros/familias que el eDNA NO detecta son taxones indicadores.
# =============================================================================
# ACTUALIZADO (reunión 24-07). Tres cambios respecto a la versión anterior:
#   1) La lista de indicadores ya no es Taxa_in_biotic_indices.csv (contenía
#      'Animalia' y entradas de rango alto que inflaban el recuento), sino las
#      listas LIMPIAS de Anders: actual_indicator_*_NSI_AMBI_only.csv.
#   2) El denominador ya no es table4_missing_in_eDNA_bottomup.csv (161 géneros,
#      124 familias) sino TODOS los taxones ausentes en eDNA: 632 especies y
#      399 géneros, que es el número que aparece en el bar plot.
#   3) Se elimina la fila de FAMILIA: no hay indicadores a rango de familia en
#      AMBI/NSI ("just remove it, doesn't make sense" - Anders).
# El cálculo se hace en Python (indicator_comparison.py); aquí solo se lee.
norm_txt <- function(x) tolower(trimws(gsub("\\.", " ", x)))   # se reutiliza más abajo

# ACTUALIZADO (añadido el recíproco): tabla5_indicadores_resumen.csv ahora trae
# 4 filas -> Species/Genus x only_morphology (perdidos en eDNA) y
# Species/Genus x only_eDNA (aportados en exclusiva por el eDNA, 105 especies y
# 97 géneros). Columnas nuevas: n_taxones (antes n_ausentes, nombre genérico
# porque ahora vale para las dos direcciones) y las ya existentes n_con_AMBI /
# n_con_NSI, que antes se calculaban pero no se mostraban en la tabla.
ambi_resumen <- read_delim(paste0(ruta_py, "tabla5_indicadores_resumen.csv"),
                           delim = ";", show_col_types = FALSE) %>%
  transmute(level, direccion,
            taxones     = n_taxones,
            indicadores = n_indicador,
            no_indic    = n_no_indic,
            con_AMBI    = n_con_AMBI,
            con_NSI     = n_con_NSI,
            pct_ind     = pct_indicador) %>%
  mutate(level     = factor(level, levels = c("Species","Genus")),
         direccion = factor(direccion, levels = c("only_morphology","only_eDNA"))) %>%
  arrange(direccion, level)

print(ambi_resumen)   # Species perdidas 632/527 (83,4%) · Genus perdidos 399/242 (60,7%)
                      # Species solo eDNA 105/39 (37,1%) · Genus solo eDNA 97/29 (29,9%)

# --- Tabla resumen AMBI/NSI (perdidos en eDNA + aportados solo por eDNA) ---
gt2_2 <- ambi_resumen %>%
  mutate(
    nivel_es = recode(as.character(level), Species="Especie", Genus="Género"),
    nivel    = case_when(
      direccion == "only_morphology" ~ paste0(nivel_es, " (solo morfología)"),
      direccion == "only_eDNA"       ~ paste0(nivel_es, " (solo eDNA)")
    )
  ) %>%
  select(nivel, taxones, indicadores, con_AMBI, con_NSI, pct_ind) %>%
  gt() %>%
  tab_header(title = md("**Taxones indicadores de índices bióticos: perdidos y aportados en exclusiva por cada método**"),
             subtitle = "Cuánta información para el biomonitoreo (AMBI y NSI) no recupera el eDNA, y cuánta aporta en exclusiva") %>%
  cols_label(nivel = "Nivel", taxones = "Nº de taxones",
             indicadores = "De ellos, indicadores", con_AMBI = "Con AMBI",
             con_NSI = "Con NSI", pct_ind = "% indicadores") %>%
  fmt_number(columns = pct_ind, decimals = 1) %>%
  data_color(columns = pct_ind,
             fn = scales::col_numeric(palette = c("#f2f7ec", sage_p), domain = c(0, 100))) %>%
  tab_source_note(md("*Cruce con las listas de indicadores AMBI y NSI (Anders, 24-07-2026), restringidas a entradas de especie y de género. Match exacto sobre el nombre completo. \"Con AMBI\"/\"Con NSI\": de los indicadores, cuántos tienen valor en cada índice (una especie puede tener ambos).*")) %>%
  estilo_gt(color_cab = sage)

print(gt2_2)
gtsave(gt2_2, paste0(ruta_d, "tab2_2_perdidos_ambi.png"), expand = 10)

# --- Figura: barras apiladas indicador vs no indicador, por nivel ---
# Se mantiene igual que antes: solo la dirección "perdidos en eDNA" (los
# nuevos datos de "solo eDNA" quedan en la tabla gt2_2, no en esta figura).
df2_2 <- ambi_resumen %>%
  filter(direccion == "only_morphology") %>%
  select(level, indicadores, no_indic) %>%
  pivot_longer(c(indicadores, no_indic), names_to = "clase", values_to = "n") %>%
  mutate(clase = recode(clase, indicadores = "Indicador biótico",
                        no_indic = "No indicador"),
         clase = factor(clase, levels = c("Indicador biótico","No indicador")))

fig2_2 <- ggplot(df2_2, aes(level, n, fill = clase)) +
  geom_col(aes(color = clase), width = 0.62, linewidth = 0.6) +
  geom_text(aes(label = ifelse(n > 0, n, "")),
            position = position_stack(vjust = 0.5),
            fontface = "bold", size = 3.5, show.legend = FALSE) +
  scale_x_discrete(labels = c(Species="Especie", Genus="Género", Family="Familia")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Indicador biótico"=sage_p, "No indicador"=grey_n)) +
  scale_color_manual(values = c("Indicador biótico"=sage, "No indicador"="#8a8a8a"), guide = "none") +
  labs(title = "Taxones perdidos por el eDNA que son indicadores bióticos",
       subtitle = "Por nivel taxonómico · exclusivos de morfología",
       x = NULL, y = "Nº de taxones") +
  tema_tfm

print(fig2_2)
ggsave(paste0(ruta_d, "fig2_2_perdidos_ambi.png"), fig2_2,
       width = 8, height = 5.2, dpi = 300, bg = "white")

# #############################################################################
# BLOQUE 3 — TAXONES COMPARTIDOS: DIFERENCIAS SIGNIFICATIVAS DE ABUNDANCIA
# #############################################################################

# =============================================================================
# TAB 3.1 — Diferencias significativas sobre el total de compartidos, por nivel
# -----------------------------------------------------------------------------
# Cuántos de los taxones compartidos difieren significativamente (BH y Bonferroni).
# =============================================================================
tab3_1 <- sumario %>%
  filter(present_in == "both") %>%
  mutate(level = factor(level, levels = niveles)) %>%
  arrange(level) %>%
  transmute(nivel = recode(as.character(level), !!!etiq_es),
            compartidos = n,
            sig_BH  = sig_BH,
            sig_bon = sig_bon,
            pct_sig_bon = round(sig_bon / n * 100, 1))

gt3_1 <- tab3_1 %>%
  gt() %>%
  tab_header(title = md("**Taxones compartidos con diferencias significativas de abundancia**"),
             subtitle = "Test de Wilcoxon pareado · corrección BH y Bonferroni") %>%
  cols_label(nivel = "Nivel", compartidos = "Compartidos",
             sig_BH = "Sig. BH", sig_bon = "Sig. Bonferroni", pct_sig_bon = "% sig. (Bonf.)") %>%
  fmt_number(columns = pct_sig_bon, decimals = 1) %>%
  data_color(columns = pct_sig_bon,
             fn = scales::col_numeric(palette = c("#f7fbff", teal_p), domain = c(0, 100))) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_body(columns = nivel)) %>%
  estilo_gt(color_cab = teal)

print(gt3_1)
gtsave(gt3_1, paste0(ruta_d, "tab3_1_significativos_por_nivel.png"), expand = 10)

# =============================================================================
# FIG 3.2 — Dispersión morfología vs eDNA a NIVEL DE CLASE (log-log)
# -----------------------------------------------------------------------------
# Escala logarítmica a propósito: las abundancias relativas cubren varios órdenes
# de magnitud (~1e-4 a ~0.6); en escala lineal todo se aplasta contra el origen.
# =============================================================================
df3_2 <- wilcox %>%
  filter(level == "Class", present_in == "both") %>%
  mutate(
    categoria = case_when(
      significant_bonferroni != "YES"       ~ "No significativa",
      direction == "conflict"               ~ "Conflicto media-mediana",
      direction == "underestimated_in_eDNA" ~ "Subestimado en eDNA",
      direction == "overestimated_in_eDNA"  ~ "Sobreestimado en eDNA"),
    categoria = factor(categoria, levels = c("Subestimado en eDNA","Sobreestimado en eDNA",
                                             "Conflicto media-mediana","No significativa")),
    mean_morpho = pmax(mean_morpho, 1e-5),
    mean_cox1   = pmax(mean_cox1,   1e-5))

# NOTA: la POSICIÓN de cada punto respecto a la diagonal refleja la relación entre
# MEDIAS (arriba = media eDNA mayor); el COLOR refleja la clasificación oficial, que
# combina media y mediana. Los "conflicto media-mediana" son puntos donde media y
# mediana apuntan en direcciones opuestas (típicamente mediana eDNA = 0 pero media
# eDNA mayor): por eso pueden quedar por encima de la diagonal sin ser "sobreestimados".
fig3_2 <- ggplot(df3_2, aes(mean_morpho, mean_cox1)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#90a4ae", linewidth = 0.6) +
  geom_point(aes(fill = categoria, color = categoria), shape = 21, size = 4, stroke = 1.4) +
  geom_text_repel(aes(label = taxon), size = 3, fontface = "italic",
                  color = ink, max.overlaps = 30) +
  annotate("text", x = 0.4,  y = 0.04, label = "media eDNA\n< morfología",
           color = "grey55", size = 2.8, fontface = "italic", alpha = 0.9) +
  annotate("text", x = 0.04, y = 0.4,  label = "media eDNA\n> morfología",
           color = "grey55", size = 2.8, fontface = "italic", alpha = 0.9) +
  scale_x_log10(labels = label_percent(accuracy = 0.01)) +
  scale_y_log10(labels = label_percent(accuracy = 0.01)) +
  scale_fill_manual(values = c("Subestimado en eDNA"=teal, "Sobreestimado en eDNA"=amber,
                               "Conflicto media-mediana"=grey_n, "No significativa"="white"),
                    name = NULL) +
  scale_color_manual(values = c("Subestimado en eDNA"=teal, "Sobreestimado en eDNA"=amber,
                                "Conflicto media-mediana"="#8a8a8a", "No significativa"=grey_n),
                     name = NULL) +
  labs(title = "Abundancia relativa media por clase: morfología vs eDNA",
       subtitle = "La posición muestra la MEDIA; el color, la clasificación (media + mediana)",
       x = "Abundancia relativa media — Morfología (%)",
       y = "Abundancia relativa media — eDNA COI (%)") +
  coord_equal() + tema_tfm + theme(legend.position = "bottom")

print(fig3_2)
ggsave(paste0(ruta_d, "fig3_2_dispersion_clase.png"), fig3_2,
       width = 8.4, height = 5.2, dpi = 300, bg = "white")

print(fig3_2)
ggsave(paste0(ruta_d, "fig3_2_dispersion_clase_over_under.png"), fig3_2,
       width = 8.4, height = 5.2, dpi = 300, bg = "white")
# =============================================================================
# TAB 3.3 — TABLA DE MIRIAM: direccionalidad del sesgo por nivel  [NUEVO]
# -----------------------------------------------------------------------------
# Entre los taxones SIGNIFICATIVOS (Bonferroni), se separan tres situaciones
# según la columna 'direction' (que ya combina media y mediana):
#   · eDNA sobreestima         -> "overestimated_in_eDNA"
#   · eDNA subestima           -> "underestimated_in_eDNA"
#   · conflicto media-mediana  -> "conflict" (media y mediana en direcciones
#     opuestas; los empates de mediana en 0 quedan fuera, van como over/under)
# =============================================================================
tab3_3 <- wilcox %>%
  filter(present_in == "both") %>%
  group_by(level) %>%
  summarise(
    compartidos = n(),
    sig         = sum(significant_bonferroni == "YES"),
    edna_sobre  = sum(significant_bonferroni == "YES" & direction == "overestimated_in_eDNA"),
    edna_infra  = sum(significant_bonferroni == "YES" & direction == "underestimated_in_eDNA"),
    conflicto   = sum(significant_bonferroni == "YES" & direction == "conflict"),
    .groups = "drop") %>%
  mutate(level = factor(level, levels = niveles)) %>%
  arrange(level) %>%
  mutate(nivel = recode(as.character(level), !!!etiq_es)) %>%
  select(nivel, compartidos, sig, edna_sobre, edna_infra, conflicto)

gt3_3 <- tab3_3 %>%
  gt() %>%
  tab_header(title = md("**Direccionalidad de las diferencias significativas (test de Wilcoxon)**"),
             subtitle = "Entre los taxones significativos (Bonferroni): sobreestima, subestima o conflicto media-mediana") %>%
  cols_label(
    nivel = "Nivel",
    compartidos = "Compartidos (total)",
    sig = "Dif. significativas",
    edna_sobre = "eDNA sobreestima",
    edna_infra = "eDNA subestima",
    conflicto = "Conflicto media-mediana") %>%
  data_color(columns = edna_infra,
             fn = scales::col_numeric(palette = c("#f7fbff", teal_p), domain = NULL)) %>%
  data_color(columns = edna_sobre,
             fn = scales::col_numeric(palette = c("#fff8f0", amber_p), domain = NULL)) %>%
  data_color(columns = conflicto,
             fn = scales::col_numeric(palette = c("#f5f5f5", grey_n), domain = NULL)) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_body(columns = nivel)) %>%
  tab_source_note(md("*'eDNA sobreestima/subestima' combinan media y mediana; varios de estos taxones tienen mediana = 0 en ambos métodos y su dirección la marca la media. 'Conflicto media-mediana' = media y mediana apuntan a direcciones opuestas (se comentan en el texto).*")) %>%
  estilo_gt(color_cab = teal)

print(gt3_3)
gtsave(gt3_3, paste0(ruta_d, "tab3_3_tabla_miriam_direccion.png"), expand = 10)

# (Complementario) versión gráfica de la tabla de Miriam: barras divergentes
df3_3 <- tab3_3 %>%
  select(nivel, edna_infra, edna_sobre, conflicto) %>%
  pivot_longer(-nivel, names_to = "cat", values_to = "n") %>%
  mutate(cat = recode(cat,
                      edna_infra = "eDNA subestima",
                      edna_sobre = "eDNA sobreestima",
                      conflicto  = "Conflicto media-mediana"),
         cat = factor(cat, levels = c("eDNA subestima","Conflicto media-mediana","eDNA sobreestima")),
         nivel = factor(nivel, levels = etiq_es))
# Etiquetas del eje X con la fraccion de taxones significativos sobre el total de
# compartidos en cada nivel (p. ej. Clase -> "Clase\n(13/22)").
etiq_frac <- setNames(
  paste0(tab3_3$nivel, "\n(", tab3_3$sig, "/", tab3_3$compartidos, ")"),
  tab3_3$nivel)

fig3_3 <- ggplot(df3_3, aes(nivel, n, fill = cat)) +
  geom_col(aes(color = cat), width = 0.66, linewidth = 0.6) +
  geom_text(aes(label = ifelse(n > 0, n, "")), position = position_stack(vjust = 0.5),
            fontface = "bold", size = 3.4, show.legend = FALSE) +
  scale_x_discrete(labels = etiq_frac) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("eDNA subestima"=teal_p, "Conflicto media-mediana"=grey_n,
                               "eDNA sobreestima"=amber_p)) +
  scale_color_manual(values = c("eDNA subestima"=teal, "Conflicto media-mediana"="#8a8a8a",
                                "eDNA sobreestima"=amber), guide = "none") +
  labs(title = "Direccionalidad del sesgo en los taxones significativos",
       subtitle = "Versión gráfica de la tabla de direccionalidad",
       x = NULL, y = "Nº de taxones significativos") +
  tema_tfm

print(fig3_3)
ggsave(paste0(ruta_d, "sup_direccion_sesgo.png"), fig3_3,
       width = 8.4, height = 5.2, dpi = 300, bg = "white")

# (Complementario / apéndice) Detalle de las clases significativas a nivel de clase
tab3_ap <- wilcox %>%
  filter(level == "Class", present_in == "both", significant_bonferroni == "YES") %>%
  mutate(direccion = case_when(
           direction == "conflict"               ~ "Conflicto media-mediana",
           direction == "underestimated_in_eDNA" ~ "Subestimada en eDNA",
           direction == "overestimated_in_eDNA"  ~ "Sobreestimada en eDNA")) %>%
  select(taxon, mean_morpho, mean_cox1, median_morpho, median_cox1, direccion, p_adjust_bonferroni) %>%
  arrange(desc(mean_morpho))
gt3_ap <- tab3_ap %>%
  gt() %>%
  tab_header(title = md("**Clases con diferencias significativas — detalle**"),
             subtitle = "Nivel de clase · Wilcoxon pareado (Bonferroni)") %>%
  cols_label(taxon="Clase", mean_morpho="Media morfo", mean_cox1="Media eDNA",
             median_morpho="Mediana morfo", median_cox1="Mediana eDNA",
             direccion="Dirección", p_adjust_bonferroni="p (Bonf.)") %>%
  fmt_number(columns = c(mean_morpho, mean_cox1, median_morpho, median_cox1), decimals = 4) %>%
  fmt_scientific(columns = p_adjust_bonferroni, decimals = 1) %>%
  data_color(columns = direccion,
             fn = function(x) dplyr::case_when(
               grepl("Subestim",   x) ~ teal_p,
               grepl("Sobreestim", x) ~ amber_p,
               TRUE                   ~ grey_n)) %>%
  tab_style(style = cell_text(weight = "bold", style = "italic"),
            locations = cells_body(columns = taxon)) %>%
  estilo_gt(color_cab = teal)

print(gt3_ap)
gtsave(gt3_ap, paste0(ruta_d, "sup_clases_significativas_detalle.png"), expand = 10)

# =============================================================================
cat("\n================  RESUMEN DE SALIDAS  ================\n")
cat("BLOQUE 1: fig1_1_venn_clase · tab1_2_clases_solo_edna · fig1_3_compartidos_exclusivos\n")
cat("BLOQUE 2: tab2_1_cobertura_midori · tab2_2_perdidos_ambi · fig2_2_perdidos_ambi\n")
cat("BLOQUE 3: tab3_1_significativos_por_nivel · fig3_2_dispersion_clase · tab3_3_tabla_miriam_direccion\n")
cat("Complementarios (sup_*): venn por nivel, % compartidos, barras midori, direccion sesgo, detalle clases\n")
print(ambi_resumen)
print(tab3_3)
cat("=====================================================\n")

# =============================================================================
# GRÁFICOS — Dispersión abundancia relativa media: morfología vs eDNA
# Andrea Moriñigo — TFM
# -----------------------------------------------------------------------------
# Scatter log-log de los taxones compartidos ("both") a distintos niveles.
#   · COLOR    = significancia (Bonferroni): rosa relleno = Dif. significativa ·
#                círculo abierto = No significativa
#   · POSICIÓN respecto a la diagonal = dirección del sesgo (arriba = eDNA mayor
#                = sobreestimado; abajo = eDNA menor = subestimado)
#
# Una función, tres niveles. En CLASE se etiquetan todos los taxones (son pocos);
# en ORDEN y ESPECIE se etiquetan SOLO los significativos (hay muchos más puntos).
# El suelo de la escala log y las etiquetas de dirección se adaptan al rango de
# cada nivel. print() para verlos en RStudio y ggsave() para descargar los PNG.
# =============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(ggrepel)

# ---- PALETA Y TEMA ----------------------------------------------------------
teal   <- "#035F72"; teal_p  <- "#7FB3C0"
rose   <- "#D77186"; rose_p  <- "#EBB6C1"
amber  <- "#E69F00"; amber_p <- "#F2D08A"
grey_n <- "#B0B0B0"; ink     <- "#37474F"

tema_tfm <- theme_minimal(base_size = 13) +
  theme(
    text                = element_text(color = ink),
    plot.title          = element_text(face = "bold", size = 15, hjust = 0, margin = margin(b = 12)),
    plot.title.position = "plot",
    axis.title          = element_text(size = 11.5),
    axis.text           = element_text(color = ink, size = 11),
    panel.grid.minor    = element_blank(),
    panel.grid.major    = element_blank(),
    axis.ticks          = element_blank(),
    legend.position     = "top",
    legend.title        = element_blank(),
    legend.text         = element_text(size = 10.5),
    plot.margin         = margin(16, 18, 14, 14)
  )

# ---- RUTAS ------------------------------------------------------------------
ruta   <- "C:/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_r/"
ruta_d <- "C:/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_r/graphics/"

wilcox <- read_delim(paste0(ruta, "wilcoxon_denominador_CLASS.csv"), delim = ";",
                     locale = locale(decimal_mark = ","), show_col_types = FALSE)

# ---- FUNCIÓN: scatter de un nivel -------------------------------------------
scatter_nivel <- function(datos, nivel, titulo, archivo_png, solo_sig_labels = FALSE) {
  
  df <- datos %>%
    filter(level == nivel, present_in == "both") %>%
    mutate(
      sig       = ifelse(significant_bonferroni == "YES", "Dif. significativa", "No significativa"),
      taxon_lab = gsub("\\.", " ", taxon)                 # "Laonice.sarsi" -> "Laonice sarsi"
    )
  
  # suelo just por debajo del menor valor real, para que los 0 no rompan el log10
  piso <- min(c(df$mean_morpho[df$mean_morpho > 0],
                df$mean_cox1[df$mean_cox1  > 0]), na.rm = TRUE) / 3
  df <- df %>% mutate(mean_morpho = pmax(mean_morpho, piso),
                      mean_cox1   = pmax(mean_cox1,   piso))
  
  # qué puntos llevan nombre: todos, o solo los significativos
  df_labels <- if (solo_sig_labels) filter(df, sig == "Dif. significativa") else df
  
  # posiciones de las etiquetas de dirección, adaptadas al rango de cada nivel
  lx <- range(log10(df$mean_morpho)); ly <- range(log10(df$mean_cox1))
  pos <- function(r, f) 10^(r[1] + diff(r) * f)
  
  fig <- ggplot(df, aes(mean_morpho, mean_cox1)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "#90a4ae", linewidth = 0.6) +
    geom_point(aes(fill = sig, color = sig), shape = 21, size = 4, stroke = 1.4) +
    geom_text_repel(data = df_labels, aes(label = taxon_lab), size = 3,
                    fontface = "italic", color = ink, max.overlaps = Inf, seed = 42) +
    annotate("text", x = pos(lx, 0.12), y = pos(ly, 0.90), label = "sobreestimado\nen eDNA",
             color = amber, size = 3, fontface = "italic", alpha = 0.9) +
    annotate("text", x = pos(lx, 0.90), y = pos(ly, 0.12), label = "subestimado\nen eDNA",
             color = teal,  size = 3, fontface = "italic", alpha = 0.9) +
    scale_x_log10(labels = label_percent(accuracy = 0.01)) +
    scale_y_log10(labels = label_percent(accuracy = 0.01)) +
    scale_fill_manual(values  = c("Dif. significativa" = rose, "No significativa" = "white")) +
    scale_color_manual(values = c("Dif. significativa" = rose, "No significativa" = grey_n)) +
    labs(title = titulo,
         x = "Abundancia relativa media — Morfología (%)",
         y = "Abundancia relativa media — eDNA COI (%)") +
    coord_equal() + tema_tfm +
    theme(legend.position = c(0.85, 0.12))
  
  print(fig)                                              # verlo en el panel Plots
  ggsave(archivo_png, fig, width = 8.4, height = 5.2, dpi = 300, bg = "white")
  cat("Guardado en:", archivo_png, "\n")
  invisible(fig)
}

# ---- GENERAR LOS TRES SCATTERS ----------------------------------------------
# CLASE — se etiquetan todos (22 taxones)
scatter_nivel(wilcox, "Class",
              "Abundancia relativa media por clase: morfología vs eDNA",
              paste0(ruta_d, "grafico_dispersion_clase.png"),
              solo_sig_labels = FALSE)

# ORDEN — solo se etiquetan los significativos (40 puntos, 22 significativos)
scatter_nivel(wilcox, "Order",
              "Abundancia relativa media por orden: morfología vs eDNA",
              paste0(ruta_d, "grafico_dispersion_orden.png"),
              solo_sig_labels = TRUE)

# ESPECIE — solo se etiquetan los significativos (88 puntos, 41 significativos)
scatter_nivel(wilcox, "Species",
              "Abundancia relativa media por especie: morfología vs eDNA",
              paste0(ruta_d, "grafico_dispersion_especie.png"),
              solo_sig_labels = TRUE)


# #############################################################################
# BLOQUE 4 — ACTUALIZACIÓN REUNIÓN 17-07  [NUEVO]
# -----------------------------------------------------------------------------
# Añade las figuras y tablas acordadas en la reunión del 17-07 con la MISMA
# estética del resto del script (paleta, tema_tfm, estilo_gt, naming figX/tabX).
# Reutiliza los objetos ya cargados arriba (wilcox, midori_sp, niveles, etiq_es,
# rutas). Por eso debe ejecutarse DESPUÉS del resto del script.
#
#   4.1  Datos base: tabla maestra de las 633 especies (abundancia, MIDORI, AMBI/NSI).
#   4.2  Dispersión abundancia media morfología vs eDNA a NIVEL DE ORDEN.
#   4.3  Cobertura MIDORI a 4 columnas: total y subconjunto mas abundante (>= 0,1 %).
#   4.4  Composicion taxonomica (nivel de orden) de las ausentes: total y >= 0,1 %.
#   4.5  Apéndice: tabla por especie (633) + candidatas + rangos altos ausentes.
# #############################################################################

# Tema completo (con subtítulo) por si se ejecuta desde el final del script.
tema_tfm <- theme_minimal(base_size = 13) +
  theme(
    text                = element_text(color = ink),
    plot.title          = element_text(face = "bold", size = 15, hjust = 0, margin = margin(b = 12)),
    plot.title.position = "plot",
    plot.subtitle       = element_text(size = 10.5, color = "grey45", margin = margin(b = 8)),
    axis.title          = element_text(size = 11.5),
    axis.text           = element_text(color = ink, size = 11),
    panel.grid.minor    = element_blank(),
    panel.grid.major    = element_blank(),
    axis.ticks          = element_blank(),
    legend.position     = "top",
    legend.title        = element_blank(),
    legend.text         = element_text(size = 10.5),
    plot.margin         = margin(16, 18, 14, 14)
  )

# ---- 4.1 Datos base: abundancia relativa (morfología) y lista AMBI/NSI -------
# NOTA: la direccion (subestimado / sobreestimado / conflicto media-mediana) se
# toma tal cual de la columna 'direction' del Wilcoxon: la mediana 0 en ambos ya
# viene resuelta por la media desde el pipeline y el conflicto se mantiene como
# categoria aparte (no se toca).
abund_morpho <- read_delim(paste0(ruta, "sp_morpho_compare.csv"), delim = ";",
                           locale = locale(decimal_mark = ","), show_col_types = FALSE) %>%
  filter(level == "Species") %>%
  transmute(Specie = gsub("\\.", " ", taxon), rel_ab_morpho = mean_morpho) %>%
  distinct(Specie, .keep_all = TRUE)

# ACTUALIZADO (24-07): antes se leía species_in_AMBI_or_NSI.csv, que era la lista
# sucia (mezclaba géneros, familias y hasta 'Animalia'). Ahora se usa la lista
# limpia de Anders: solo AMBI + NSI y solo entradas de especie.
# Se conservan además los VALORES de AMBI y NSI para poder llevarlos al anexo.
if (!exists("norm_txt")) norm_txt <- function(x) tolower(trimws(gsub("\\.", " ", x)))

ind_sp_new <- read_delim(paste0(ruta_bio, "actual_indicator_species_NSI_AMBI_only.csv"),
                         delim = "\t", show_col_types = FALSE) %>%
  mutate(key = norm_txt(Species)) %>%
  distinct(key, .keep_all = TRUE)          # evita duplicar filas en el join

ambi_list <- ind_sp_new$key
cat("[4.1] Lista de indicadores (especies, AMBI+NSI):", length(ambi_list), "entradas únicas\n")

# Tabla maestra: 633 especies exclusivas de morfología con todo lo necesario.
master <- midori_sp %>%
  mutate(Specie = trimws(Specie),
         # BLOQUE 4 autosuficiente: si no se ha corrido el BLOQUE 2, calculamos
         # aqui las logicas de MIDORI y el nivel mas fino con referencia.
         across(c(en_MIDORI_especie, en_MIDORI_genero, family_en_MIDORI, order_en_MIDORI), as.logical),
         primer_nivel_con_ref = factor(case_when(
           en_MIDORI_especie ~ "Especie",
           en_MIDORI_genero  ~ "Género",
           family_en_MIDORI  ~ "Familia",
           order_en_MIDORI   ~ "Orden",
           TRUE              ~ "Clase"),
           levels = c("Especie","Género","Familia","Orden","Clase"))) %>%
  left_join(abund_morpho, by = "Specie") %>%
  mutate(rel_ab_morpho = ifelse(is.na(rel_ab_morpho), 0, rel_ab_morpho),
         en_01pct = rel_ab_morpho >= 0.001,                 # subconjunto mas abundante: >= 0,1 %
         # (rel_ab_morpho es una PROPORCION 0-1, asi que 0.001 = 0,1 %, no 1 %)
         # norm_txt en los DOS lados: la lista viene con espacios y los taxones
         # de R con puntos. Sin normalizar, el match daría 0.
         AMBI_NSI = norm_txt(Specie) %in% ambi_list,
         cat_ord  = ifelse(is.na(Order) | Order == "",
                           paste0(Class, " (s/ orden)"), Order)) %>%
  # NUEVO (24-07): trae los VALORES de AMBI (grupo ecológico 1-5) y NSI de cada
  # especie indicadora, para poder mostrarlos como columnas del anexo. Las no
  # indicadoras quedan con NA, que es la información correcta.
  mutate(key = norm_txt(Specie)) %>%
  left_join(ind_sp_new %>% select(key, AMBI_val = AMBI, NSI_val = NSI), by = "key") %>%
  select(-key)

cat("\n[4.1] Especies:", nrow(master),
    "| subconjunto >= 0,1 %:", sum(master$en_01pct),
    "| indicadores AMBI/NSI:", sum(master$AMBI_NSI), "\n")
stopifnot(nrow(master) == 632)   # si el join duplicara filas, aquí saltaría

# CORRECCION: tamanos de muestra dinamicos. Antes los rotulos de titulos, cabeceras
# y notas al pie llevaban "633" y "78" escritos a mano; si wilcoxon_comparison.R cambiaba,
# las cifras de las tablas se actualizaban pero los textos seguian diciendo lo viejo.
N_TOT  <- nrow(master)                 # especies exclusivas de morfologia
N_1PCT <- sum(master$en_01pct)          # subconjunto mas abundante (>= 0,1 %)
# NUEVO (24-07): especies SIN referencia de especie en MIDORI. Es el subconjunto
# que pasa a ser la tarta (b): "donde hay que hacer el esfuerzo de barcoding".
N_SINREF <- sum(!master$en_MIDORI_especie)
cat("[4.1] Sin referencia de especie en MIDORI:", N_SINREF, "\n")

# ---- 4.2 Dispersión abundancia media morfología vs eDNA — NIVEL DE ORDEN -----
# Misma idea y estética que la Fig 3.2 (dispersión de clase), pero a nivel de
# ORDEN. COLOR = clasificación (media + mediana): subestimado / sobreestimado /
# conflicto media-mediana / no significativa. POSICIÓN respecto a la diagonal =
# relación entre medias (arriba = media eDNA mayor). Se etiquetan solo los
# significativos porque hay 40 puntos.
df4_orden <- wilcox %>%
  filter(level == "Order", present_in == "both") %>%
  mutate(
    categoria = case_when(
      significant_bonferroni != "YES"       ~ "No significativa",
      direction == "conflict"               ~ "Conflicto media-mediana",
      direction == "underestimated_in_eDNA" ~ "Subestimado en eDNA",
      direction == "overestimated_in_eDNA"  ~ "Sobreestimado en eDNA"),
    categoria = factor(categoria, levels = c("Subestimado en eDNA","Sobreestimado en eDNA",
                                             "Conflicto media-mediana","No significativa")),
    taxon_lab = gsub("\\.", " ", taxon))

# suelo bajo el menor valor real para que los 0 no rompan el log10
piso <- min(c(df4_orden$mean_morpho[df4_orden$mean_morpho > 0],
              df4_orden$mean_cox1  [df4_orden$mean_cox1   > 0]), na.rm = TRUE) / 3
df4_orden <- df4_orden %>% mutate(mean_morpho = pmax(mean_morpho, piso),
                                  mean_cox1   = pmax(mean_cox1,   piso))
df4_lab <- df4_orden %>% filter(significant_bonferroni == "YES")

# posiciones de las etiquetas de dirección, adaptadas al rango del nivel
lx <- range(log10(df4_orden$mean_morpho)); ly <- range(log10(df4_orden$mean_cox1))
pos <- function(r, f) 10^(r[1] + diff(r) * f)

fig4_1 <- ggplot(df4_orden, aes(mean_morpho, mean_cox1)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#90a4ae", linewidth = 0.6) +
  geom_point(aes(fill = categoria, color = categoria), shape = 21, size = 4, stroke = 1.4) +
  geom_text_repel(data = df4_lab, aes(label = taxon_lab), size = 3, fontface = "italic",
                  color = ink, max.overlaps = Inf, seed = 42) +
  annotate("text", x = pos(lx, 0.12), y = pos(ly, 0.90), label = "sobreestimado\nen eDNA",
           color = amber, size = 3, fontface = "italic", alpha = 0.9) +
  annotate("text", x = pos(lx, 0.90), y = pos(ly, 0.12), label = "subestimado\nen eDNA",
           color = teal,  size = 3, fontface = "italic", alpha = 0.9) +
  scale_x_log10(labels = label_percent(accuracy = 0.01)) +
  scale_y_log10(labels = label_percent(accuracy = 0.01)) +
  scale_fill_manual(values = c("Subestimado en eDNA"=teal, "Sobreestimado en eDNA"=amber,
                               "Conflicto media-mediana"=grey_n, "No significativa"="white"),
                    name = NULL, drop = FALSE) +
  scale_color_manual(values = c("Subestimado en eDNA"=teal, "Sobreestimado en eDNA"=amber,
                                "Conflicto media-mediana"="#8a8a8a", "No significativa"=grey_n),
                     name = NULL, drop = FALSE) +
  labs(title = "Abundancia relativa media por orden: morfología vs eDNA",
       subtitle = "La posición muestra la MEDIA; el color, la clasificación (media + mediana)",
       x = "Abundancia relativa media — Morfología (%)",
       y = "Abundancia relativa media — eDNA COI (%)") +
  coord_equal() + tema_tfm + theme(legend.position = "bottom")

print(fig4_1)
ggsave(paste0(ruta_d, "fig4_1_dispersion_orden.png"), fig4_1,
       width = 8.4, height = 6.2, dpi = 300, bg = "white")


# ---- 4.3 Cobertura MIDORI a 4 columnas (total y subconjunto >= 0,1 %) --------
cob633 <- master %>% count(primer_nivel_con_ref, name = "n633", .drop = FALSE) %>%
  mutate(pct633 = round(n633 / sum(n633) * 100, 1))
cob78  <- master %>% filter(en_01pct) %>%
  count(primer_nivel_con_ref, name = "n78", .drop = FALSE) %>%
  mutate(pct78 = round(n78 / sum(n78) * 100, 1))

tab4_2 <- cob633 %>% left_join(cob78, by = "primer_nivel_con_ref") %>%
  transmute(nivel = as.character(primer_nivel_con_ref), n633, pct633, n78, pct78) %>%
  bind_rows(tibble(nivel = "Total", n633 = sum(cob633$n633), pct633 = 100,
                   n78 = sum(cob78$n78), pct78 = 100))

# REESTRUCTURADA (reunión 24-07, petición de Miriam):
#   - Las cuatro columnas se agrupan en DOS bloques con cabecera común
#     ("Total" y "Más abundantes"), en vez de cuatro títulos largos sueltos.
#     Así se ve de un vistazo que las dos primeras van juntas y las dos últimas
#     también ("you see that the two first columns belong together").
#   - Dentro de cada bloque, las columnas se llaman solo "n" y "%".
#   - El símbolo % va DENTRO de la celda, para poder quitarlo del texto.
#   - Línea vertical separando los dos bloques.
gt4_2 <- tab4_2 %>%
  gt() %>%
  tab_header(title = md("**Cobertura en MIDORI253darn: total y subconjunto más abundante**"),
             subtitle = "Nivel taxonómico más fino con secuencia de referencia") %>%
  cols_label(nivel = "Nivel más fino con referencia",
             n633 = "n", pct633 = "%", n78 = "n", pct78 = "%") %>%
  tab_spanner(label = "Total",           columns = c(n633, pct633), id = "sp_total") %>%
  tab_spanner(label = "Más abundantes",  columns = c(n78,  pct78),  id = "sp_abund") %>%
  # Los dos bloques miden lo mismo de ancho: 70 + 80 px cada uno.
  cols_width(nivel ~ px(240),
             n633  ~ px(70), pct633 ~ px(80),
             n78   ~ px(70), pct78  ~ px(80)) %>%
  # El % dentro de la celda: pattern añade el símbolo sin convertir a texto,
  # así data_color sigue funcionando sobre el valor numérico.
  fmt_number(columns = c(pct633, pct78), decimals = 1, pattern = "{x} %") %>%
  data_color(columns = pct633,
             fn = scales::col_numeric(palette = c("#f7fbff", teal_p), domain = c(0, 100))) %>%
  data_color(columns = pct78,
             fn = scales::col_numeric(palette = c("#f7fbff", teal_p), domain = c(0, 100))) %>%
  cols_align(align = "center", columns = c(n633, pct633, n78, pct78)) %>%
  # --- Las cabeceras de bloque, del mismo azul que el resto de la cabecera ----
  # estilo_gt() solo colorea cells_column_labels(), no los spanners, así que hay
  # que pintarlos aquí. Va antes de estilo_gt() pero no se pisan entre sí.
  tab_style(style = list(cell_fill(color = teal),
                         cell_text(color = "white", weight = "bold")),
            locations = cells_column_spanners(spanners = everything())) %>%
  # --- Líneas verticales: blancas sobre el azul, grises en el cuerpo ----------
  tab_style(style = cell_borders(sides = "left", color = "white", weight = px(1.5)),
            locations = list(cells_column_labels(columns = c(n633, n78)),
                             cells_column_spanners(spanners = everything()))) %>%
  tab_style(style = cell_borders(sides = "left", color = "#b0bec5", weight = px(1.5)),
            locations = cells_body(columns = c(n633, n78))) %>%
  # --- Fila de totales resaltada ---
  tab_style(style = cell_fill(color = "#f2f2f2"),
            locations = cells_body(rows = nivel == "Total")) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(rows = nivel == "Total")) %>%
  tab_source_note(md(paste0("*Especies exclusivas de morfología (n = ", N_TOT,
                            "). El subconjunto más abundante son las de abundancia relativa media >= 0,1 % (n = ",
                            N_1PCT, ").*"))) %>%
  estilo_gt(color_cab = teal)
print(gt4_2)
gtsave(gt4_2, paste0(ruta_d, "tab4_2_cobertura_4col.png"), expand = 10)

# ---- 4.4 Composicion taxonomica (nivel de orden): total y >= 0,1 % ----------
# Paleta cualitativa apagada, en el mismo espacio cromático que la del TFM.
paleta_tax <- c(teal, amber, rose, sage, teal_p, "#B07AA1", "#5B8AA6",
                "#C57B57", "#8A9A5B", rose_p, amber_p, sage_p, "#9B8AA6", "#D9B382")
top_cats <- master %>% count(cat_ord, sort = TRUE) %>% head(14) %>% pull(cat_ord)
col_map  <- setNames(paleta_tax[seq_along(top_cats)], top_cats)
col_map["Otros órdenes"] <- grey_n

comp_orden <- function(dat, titulo, subtit, archivo) {
  d <- dat %>%
    mutate(cat = ifelse(cat_ord %in% top_cats, cat_ord, "Otros órdenes")) %>%
    count(cat)
  lev <- c(top_cats, "Otros órdenes"); lev <- lev[lev %in% d$cat]
  d <- d %>% mutate(cat = factor(cat, levels = lev)) %>% arrange(cat)
  tot <- sum(d$n)
  d <- d %>% mutate(pct = round(n / tot * 100, 1), etq = paste0(n, "  (", pct, "%)"))
  alt <- max(4, 0.46 * nrow(d) + 1.5)
  g <- ggplot(d, aes(n, cat, fill = cat)) +
    geom_col(width = 0.76, color = "white", linewidth = 0.5) +
    geom_text(aes(label = etq), hjust = -0.12, size = 3.4, color = ink) +
    scale_fill_manual(values = col_map, guide = "none") +
    scale_y_discrete(limits = rev(lev)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = titulo, subtitle = subtit,
         x = "Nº de especies ausentes en eDNA", y = NULL) +
    tema_tfm +
    theme(panel.grid.major.x = element_line(color = "#ECEFF1", linewidth = 0.6))
  print(g)
  ggsave(archivo, g, width = 8.8, height = alt, dpi = 300, bg = "white")
  invisible(g)
}

fig4_3 <- comp_orden(master,
                     "Composición taxonómica de las especies ausentes en eDNA",
                     "Nivel de orden · las 15 categorías mayores + 'Otros'",
                     paste0(ruta_d, "fig4_3_composicion_634.png"))
# Versión de barras del subconjunto abundante. Se mantiene como alternativa,
# pero la figura que va al TFM es la tarta (b), que ahora son las SIN referencia.
fig4_4 <- comp_orden(dplyr::filter(master, en_01pct),
                     "Composición taxonómica del subconjunto más abundante",
                     "Nivel de orden · abundancia relativa media >= 0,1 %",
                     paste0(ruta_d, "fig4_4_composicion_78.png"))

# --- Versión PIE de la composición (fig4_3b / fig4_4b) ------------------------
# Mismo contenido y colores que las barras (col_map, top_cats); útil como
# alternativa visual. Etiqueta el % dentro de las porciones grandes (>= 4 %) y
# lista categoría + n + % en la leyenda.
comp_orden_pie <- function(dat, titulo, subtit, archivo) {
  d <- dat %>%
    mutate(cat = ifelse(cat_ord %in% top_cats, cat_ord, "Otros órdenes")) %>%
    count(cat)
  lev <- c(top_cats, "Otros órdenes"); lev <- lev[lev %in% d$cat]
  d <- d %>% mutate(cat = factor(cat, levels = lev)) %>% arrange(cat) %>%
    mutate(prop = n / sum(n) * 100,
           # etiqueta en DOS lineas: el nombre del orden arriba y "n (%)" debajo.
           # Asi la etiqueta mas larga (Polychaeta (s/ orden)) ya no marca el ancho
           # de la leyenda y la figura entra mejor a lo ancho de la pagina.
           leg  = paste0(cat, "\n", n, " (", round(prop, 1), "%)"))
  g <- ggplot(d, aes(x = 1, y = prop, fill = cat)) +
    geom_col(width = 1, color = "white", linewidth = 0.7) +
    coord_polar(theta = "y", start = 0) +
    geom_text(aes(label = ifelse(prop >= 4, paste0(round(prop), "%"), "")),
              position = position_stack(vjust = 0.5),
              color = ink, size = 3.1, fontface = "bold") +
    scale_fill_manual(values = col_map, breaks = as.character(d$cat),
                      labels = d$leg, name = NULL) +
    labs(title = titulo, subtitle = subtit) +
    theme_void(base_size = 13) +
    theme(plot.title      = element_text(face = "bold", size = 15, color = ink,
                                         hjust = 0, margin = margin(b = 6)),
          plot.subtitle   = element_text(size = 10.5, color = "grey45", margin = margin(b = 10)),
          legend.text      = element_text(size = 13, color = ink),
          legend.key.size  = unit(0.8, "cm"),
          legend.spacing.y = unit(0.15, "cm"),
          legend.position  = "right",
          plot.margin     = margin(16, 18, 14, 14))
  print(g)
  ggsave(archivo, g, width = 11, height = 6.8, dpi = 300, bg = "white")
  invisible(g)
}

# CAMBIO (reunión 24-07): la tarta (b) ya NO es el subconjunto más abundante,
# sino las especies SIN referencia de especie en MIDORI. Cada tarta responde a
# una pregunta distinta:
#   (a) ¿hay una tendencia de organismos cuyo ADN no se recupera?
#   (b) ¿en qué grupos hay que concentrar el esfuerzo de barcoding?
fig4_3b <- comp_orden_pie(master,
                          "Composición taxonómica de las especies ausentes en eDNA",
                          "Nivel de orden · las 15 categorías mayores + 'Otros'",
                          paste0(ruta_d, "fig4_3b_composicion_634_pie.png"))
fig4_4b <- comp_orden_pie(dplyr::filter(master, !en_MIDORI_especie),
                          "Composición taxonómica de las especies sin referencia en MIDORI",
                          "Nivel de orden · dónde concentrar el esfuerzo de barcoding",
                          paste0(ruta_d, "fig4_4b_composicion_sinref_pie.png"))

# ============================================================================
# DOS FORMAS DE COMBINAR LAS TARTAS EN UNA SOLA FIGURA (se generan las dos)
# ============================================================================

# --- OPCION 1: VERTICAL (una encima de otra), SIN titulos, con las DOS leyendas
#     completas (cada orden con su n y %). No se pierde informacion.
#     Se generan dos tartas sin titulo (titulo = NULL, subtitulo = NULL) y se
#     apilan con magick.
comp_orden_pie(master, NULL, NULL,
               paste0(ruta_d, "_pie_a_sintitulo.png"))
comp_orden_pie(dplyr::filter(master, !en_MIDORI_especie), NULL, NULL,
               paste0(ruta_d, "_pie_b_sintitulo.png"))
pies_vert <- image_read(c(paste0(ruta_d, "_pie_a_sintitulo.png"),
                          paste0(ruta_d, "_pie_b_sintitulo.png")))
image_write(image_append(pies_vert, stack = TRUE),
            paste0(ruta_d, "fig4_34b_composicion_pies_VERTICAL.png"))
cat("  Guardada OPCION 1 (vertical, sin titulos, dos leyendas): fig4_34b_composicion_pies_VERTICAL.png\n")

# --- OPCION 2: HORIZONTAL (lado a lado) con LEYENDA COMPARTIDA ----------------
# Las categorias del subconjunto son un subconjunto de las del total y ambas usan
# la misma paleta (col_map), asi que con facet_wrap se dibujan las dos tartas una
# al lado de la otra y ggplot genera una unica leyenda para las dos. Se usa
# scale_fill_manual(drop = FALSE) para que la leyenda liste TODAS las categorias
# aunque alguna no aparezca en el subconjunto. El n y % de cada orden ya no van en
# la leyenda (serian distintos por tarta); el % se muestra dentro de las porciones
# grandes (>= 4 %) y los recuentos por orden estan en las Figuras de barras.
comp_pies_combo <- function(dat_a, dat_b, tit_a, tit_b, archivo) {
  lev <- c(top_cats, "Otros órdenes")
  prep <- function(dd, panel) {
    dd %>%
      mutate(cat = ifelse(cat_ord %in% top_cats, cat_ord, "Otros órdenes")) %>%
      count(cat) %>%
      mutate(cat = factor(cat, levels = lev)) %>%
      arrange(cat) %>%
      mutate(prop = n / sum(n) * 100, panel = panel)
  }
  d <- bind_rows(prep(dat_a, tit_a), prep(dat_b, tit_b)) %>%
    mutate(panel = factor(panel, levels = c(tit_a, tit_b)))
  g <- ggplot(d, aes(x = 1, y = prop, fill = cat)) +
    geom_col(width = 1, color = "white", linewidth = 0.6) +
    coord_polar(theta = "y", start = 0) +
    geom_text(aes(label = ifelse(prop >= 4, paste0(round(prop), "%"), "")),
              position = position_stack(vjust = 0.5),
              color = ink, size = 3.0, fontface = "bold") +
    facet_wrap(~ panel) +
    scale_fill_manual(values = col_map, drop = FALSE, name = "Orden") +
    theme_void(base_size = 13) +
    theme(strip.text      = element_text(face = "bold", size = 13, color = ink,
                                         margin = margin(b = 8)),
          legend.text     = element_text(size = 12, color = ink),
          legend.title    = element_text(size = 12.5, face = "bold", color = ink),
          legend.key.size = unit(0.6, "cm"),
          legend.position = "right",
          plot.margin     = margin(14, 16, 14, 14))
  print(g)
  ggsave(archivo, g, width = 12.5, height = 6.0, dpi = 300, bg = "white")
  invisible(g)
}
fig4_34b <- comp_pies_combo(
  master, dplyr::filter(master, !en_MIDORI_especie),
  "(a) Especies ausentes en eDNA",
  "(b) Especies sin referencia en MIDORI",
  paste0(ruta_d, "fig4_34b_composicion_pies_HORIZONTAL.png"))
cat("  Guardada OPCION 2 (horizontal, leyenda compartida): fig4_34b_composicion_pies_HORIZONTAL.png\n")

# --- OPCION 3: lado a lado, UNA LEYENDA A CADA LADO (con n y %), tartas pequenas
#     La (a) lleva su leyenda a la IZQUIERDA y la (b) a la DERECHA, cada una con
#     el n y % de sus propios ordenes (como las tartas individuales). Las tartas
#     se hacen mas pequenas para dar aire a las leyendas. Se pegan con magick, de
#     modo que las leyendas quedan en los extremos y las dos tartas, juntas en el
#     centro.
pie_lado <- function(dat, titulo, archivo, leg_pos, t_hjust) {
  d <- dat %>%
    mutate(cat = ifelse(cat_ord %in% top_cats, cat_ord, "Otros órdenes")) %>%
    count(cat)
  lev <- c(top_cats, "Otros órdenes"); lev <- lev[lev %in% d$cat]
  d <- d %>% mutate(cat = factor(cat, levels = lev)) %>% arrange(cat) %>%
    mutate(prop = n / sum(n) * 100,
           leg  = paste0(cat, " — ", n, " (", round(prop, 1), "%)"))
  g <- ggplot(d, aes(x = 1, y = prop, fill = cat)) +
    geom_col(width = 1, color = "white", linewidth = 0.7) +
    coord_polar(theta = "y", start = 0) +
    geom_text(aes(label = ifelse(prop >= 4, paste0(round(prop), "%"), "")),
              position = position_stack(vjust = 0.5),
              color = ink, size = 3.0, fontface = "bold") +
    scale_fill_manual(values = col_map, breaks = as.character(d$cat),
                      labels = d$leg, name = NULL) +
    labs(title = titulo) +
    theme_void(base_size = 13) +
    # Todo lo mas compacto posible para que, al ajustar la figura al ancho de la
    # pagina de Word, las tartas no queden pequenas: leyenda apretada, margenes
    # minimos y casi cero margen en el lado por el que se pegan las dos tartas.
    theme(plot.title       = element_text(face = "bold", size = 11, color = ink,
                                          hjust = t_hjust, margin = margin(b = 3)),
          legend.text      = element_text(size = 9, color = ink, lineheight = 0.95),
          legend.key.size  = unit(0.42, "cm"),
          legend.margin    = margin(0, 0, 0, 0),
          legend.position  = leg_pos,
          plot.margin      = if (leg_pos == "left") margin(5, 1, 5, 3)
                             else                    margin(5, 3, 5, 1))
  print(g)
  ggsave(archivo, g, width = 5.9, height = 5.4, dpi = 300, bg = "white")
  invisible(g)
}
# Titulos descriptivos, centrados sobre cada panel (hjust = 0.5). Al ir cada uno
# dentro de su propia imagen, no se pisan aunque las tartas esten muy juntas.
pie_lado(master,
         "(a) Especies ausentes en eDNA",
         paste0(ruta_d, "_pie_a_legizq.png"), "left", 0.5)
pie_lado(dplyr::filter(master, !en_MIDORI_especie),
         "(b) Especies sin referencia en MIDORI",
         paste0(ruta_d, "_pie_b_legder.png"), "right", 0.5)
pies_lados <- image_read(c(paste0(ruta_d, "_pie_a_legizq.png"),
                           paste0(ruta_d, "_pie_b_legder.png")))
image_write(image_append(pies_lados),
            paste0(ruta_d, "fig4_34b_composicion_pies_LEYENDAS_LADOS.png"))
cat("  Guardada OPCION 3 (leyendas a los lados, con n y %): fig4_34b_composicion_pies_LEYENDAS_LADOS.png\n")

# =============================================================================
# FIG 4.5 — BRECHA DE BARCODING: % de cada orden SIN referencia en MIDORI
# -----------------------------------------------------------------------------
# ALTERNATIVA a la tarta (b). Las dos tartas se parecen mucho porque (b) es un
# subconjunto del 56 % de (a): al ser tan grande, arrastra la misma composición.
# Aquí no se compara composición sino TASA: de cada orden ausente en eDNA, qué
# porcentaje carece de referencia de especie en MIDORI. Eso sí distingue los
# grupos y da una recomendación concreta de dónde barcodear.
# Se excluyen los órdenes con menos de MIN_N especies: con 2 o 3 el porcentaje
# es ruido (un 100 % sobre 2 especies no dice nada).
# =============================================================================
MIN_N <- 8

gap_orden <- master %>%
  mutate(cat = cat_ord) %>%
  group_by(cat) %>%
  summarise(total   = n(),
            sin_ref = sum(!en_MIDORI_especie),
            .groups = "drop") %>%
  filter(total >= MIN_N) %>%
  mutate(pct = round(sin_ref / total * 100, 1),
         etq = paste0(sin_ref, "/", total, "  (", pct, " %)")) %>%
  arrange(pct)

media_global <- round(sum(!master$en_MIDORI_especie) / nrow(master) * 100, 1)

# Color según esté por encima o por debajo de la media global: lo que queda a la
# derecha de la línea es donde el esfuerzo de barcoding rinde más.
gap_orden <- gap_orden %>%
  mutate(cat = factor(cat, levels = cat),
         grupo = ifelse(pct >= media_global, "Por encima de la media", "Por debajo de la media"))

fig4_5 <- ggplot(gap_orden, aes(pct, cat, fill = grupo)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.5) +
  geom_vline(xintercept = media_global, linetype = "dashed",
             color = "#546e7a", linewidth = 0.7) +
  annotate("text", x = media_global, y = 0.3,
           label = paste0("media global: ", media_global, " %"),
           hjust = -0.05, vjust = 0, size = 3.1, color = "#546e7a", fontface = "italic") +
  geom_text(aes(label = etq), hjust = -0.10, size = 3.3, color = ink) +
  # Azul (teal), que es el color asignado a "Solo morfología" en el resto de
  # figuras: estas 632 especies son exclusivas de la morfología. Tono fuerte por
  # encima de la media y tono claro por debajo, para no meter un color nuevo.
  scale_fill_manual(values = c("Por encima de la media" = teal,
                               "Por debajo de la media" = teal_p)) +
  scale_x_continuous(limits = c(0, 118), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, " %"),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Brecha de barcoding por grupo taxonómico",
       subtitle = paste0("Porcentaje de especies ausentes en eDNA que carecen de referencia de especie en MIDORI253darn\n",
                         "Solo órdenes con al menos ", MIN_N, " especies ausentes"),
       x = "Especies exclusivas de morfología sin referencia en MIDORI",
       y = "Grupo taxonómico (orden)") +
  tema_tfm +
  theme(panel.grid.major.x = element_line(color = "#ECEFF1", linewidth = 0.6),
        axis.line.y        = element_line(color = "#CFD8DC", linewidth = 0.4),
        axis.title.y       = element_text(margin = margin(r = 10)))

print(fig4_5)
ggsave(paste0(ruta_d, "fig4_5_brecha_barcoding.png"), fig4_5,
       width = 9.2, height = 0.42 * nrow(gap_orden) + 2.2, dpi = 300, bg = "white")
cat("Guardada: fig4_5_brecha_barcoding.png |", nrow(gap_orden), "órdenes | media",
    media_global, "%\n")

# Versión tabla, por si prefieres números a barras (o para el anexo).
write.csv2(gap_orden %>% arrange(desc(pct)) %>% select(Orden = cat, Ausentes = total,
                                                       Sin_referencia = sin_ref, Pct = pct),
           paste0(ruta_d, "tab4_7_brecha_barcoding.csv"),
           row.names = FALSE, fileEncoding = "UTF-8")

print(fig4_3)
print(fig4_4)
print(fig4_3b)
print(fig4_4b)


# ---- 4.5 Apéndice: tabla por especie (633) + candidatas + rangos altos ------
tabla_apendice <- master %>%
  transmute(Especie = Specie, Clase = Class,
            Orden   = ifelse(is.na(Order)  | Order  == "", "-", Order),
            Familia = ifelse(is.na(Family) | Family == "", "-", Family),
            Genero  = Genus,
            Ab_rel_morfo_pct   = round(rel_ab_morpho * 100, 4),
            Subconjunto_1pct   = ifelse(en_01pct, "Si", ""),   # nombre de columna conservado
            Nivel_MIDORI       = as.character(primer_nivel_con_ref),
            Congeneres_MIDORI  = n_especies_genero_MIDORI,
            # ACTUALIZADO (24-07): calculado con la lista limpia de Anders.
            Indicador_AMBI_NSI = ifelse(AMBI_NSI, "Si", "No"),
            # NUEVO: los valores del índice, no solo el Si/No.
            AMBI = AMBI_val,
            NSI  = round(NSI_val, 2)) %>%
  arrange(Clase, Orden, Familia, Genero, Especie)
write.csv2(tabla_apendice, paste0(ruta_d, "tab4_3_apendice_633_especies.csv"),
           row.names = FALSE, fileEncoding = "UTF-8")

# ---- NUEVO: anexo equivalente a NIVEL DE GÉNERO (399 géneros ausentes) -------
# Mismo formato que el de especies. Sale de indicator_comparison.py, que ya trae
# el Si/No y los valores de AMBI y NSI para cada uno de los 399 géneros.
anexo_generos <- read_delim(paste0(ruta_py, "indicadores_generos_ausentes.csv"),
                            delim = ";", show_col_types = FALSE) %>%
  transmute(Genero  = taxon,
            Clase   = Class,
            Orden   = ifelse(is.na(Order)  | Order  == "", "-", Order),
            Familia = ifelse(is.na(Family) | Family == "", "-", Family),
            Ab_rel_morfo_pct   = round(mean_morpho * 100, 4),
            # as.logical: Python escribe "True"/"False" como texto.
            Indicador_AMBI_NSI = ifelse(as.logical(es_indicador), "Si", "No"),
            AMBI = AMBI,
            NSI  = round(NSI, 2)) %>%
  arrange(Clase, Orden, Familia, Genero)
write.csv2(anexo_generos, paste0(ruta_d, "tab4_6_apendice_399_generos.csv"),
           row.names = FALSE, fileEncoding = "UTF-8")
cat("Anexo de géneros:", nrow(anexo_generos), "filas |",
    sum(anexo_generos$Indicador_AMBI_NSI == "Si"), "indicadores\n")

# Candidatas: poliquetos del subconjunto >= 0,1 % con referencia de especie en MIDORI.
candidatas <- master %>%
  filter(en_01pct, en_MIDORI_especie, Class == "Polychaeta") %>%
  transmute(Especie = Specie,
            Orden   = ifelse(is.na(Order) | Order == "", "-", Order),
            Familia = Family,
            Ab_rel_morfo_pct  = round(rel_ab_morpho * 100, 3),
            Congeneres_MIDORI = n_especies_genero_MIDORI,
            AMBI_NSI = ifelse(AMBI_NSI, "Si", "No")) %>%
  arrange(Orden, desc(Ab_rel_morfo_pct))

gt4_3 <- candidatas %>%
  gt() %>%
  tab_header(title = md("**Candidatas para el ejercicio de asignación (poliquetos del 0,1 %)**"),
             subtitle = paste0("Con referencia de especie en MIDORI253darn · n = ", nrow(candidatas))) %>%
  cols_label(Especie = "Especie (poliqueto)", Orden = "Orden", Familia = "Familia",
             Ab_rel_morfo_pct = "Ab. rel. morfo (%)",
             Congeneres_MIDORI = "Congéneres MIDORI", AMBI_NSI = "AMBI/NSI") %>%
  tab_style(style = cell_text(style = "italic"), locations = cells_body(columns = Especie)) %>%
  tab_source_note(md("*De esta lista se elige la decena para el BLAST/OTU falso (Miriam).*")) %>%
  estilo_gt(color_cab = teal)
print(gt4_3)
gtsave(gt4_3, paste0(ruta_d, "tab4_4_candidatas.png"), expand = 10)

# Rangos altos AUSENTES en eDNA (para descartar de golpe especies de esos grupos).
rangos_ausentes <- function(nivel) {
  wilcox %>% filter(level == nivel, present_in == "only_morphology") %>%
    pull(taxon) %>% sort()
}
clases_aus   <- rangos_ausentes("Class")
ordenes_aus  <- rangos_ausentes("Order")
familias_aus <- rangos_ausentes("Family")
write.csv2(data.frame(orden_ausente_eDNA   = ordenes_aus),
           paste0(ruta_d, "tab4_5_ordenes_ausentes.csv"),  row.names = FALSE, fileEncoding = "UTF-8")
write.csv2(data.frame(familia_ausente_eDNA = familias_aus),
           paste0(ruta_d, "tab4_5_familias_ausentes.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("\n================  BLOQUE 4 (reunión 17-07)  ================\n")
cat("Figuras: fig4_1_dispersion_orden · fig4_3_composicion_634 (+ fig4_3b pie) · fig4_4_composicion_78 (+ fig4_4b pie)\n")
cat("Tablas : tab4_2_cobertura_4col · tab4_4_candidatas\n")
cat("CSV    : tab4_3_apendice_633_especies · tab4_5_ordenes/familias_ausentes\n")
cat("Rangos altos ausentes -> clases:", length(clases_aus),
    "| órdenes:", length(ordenes_aus), "| familias:", length(familias_aus), "\n")
cat("===========================================================\n")
