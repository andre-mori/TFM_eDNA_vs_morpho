# =============================================================================
# GRÁFICOS — Dispersión morfología vs eDNA por nivel, COLOR = CLASIFICACIÓN
# Andrea Moriñigo — TFM
# -----------------------------------------------------------------------------
# Misma forma que la FIG 3.2 (nivel de clase), pero para los 5 niveles:
#   · POSICIÓN = media de abundancia relativa (log-log, en %)
#   · COLOR    = clasificación (combina media y mediana):
#        Subestimado en eDNA (teal) · Sobreestimado en eDNA (ámbar) ·
#        Conflicto media-mediana (gris relleno) · No significativa (círculo abierto)
#
# En CLASE se etiquetan todos los taxones; en ORDEN/FAMILIA/GÉNERO/ESPECIE solo
# los significativos (los coloreados), porque hay muchos puntos. El suelo de la
# escala log y las anotaciones de dirección se adaptan al rango de cada nivel.
# print() para verlos en RStudio y ggsave() para descargar los PNG.
# NOTA: en género/especie hay muchos significativos (~50-64); si las etiquetas se
# amontonan, sube width/height en ggsave o pon solo_sig_labels con menos taxones.
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
    plot.subtitle       = element_text(size = 10.5, color = "grey45", margin = margin(b = 8)),
    axis.title          = element_text(size = 11.5),
    axis.text           = element_text(color = ink, size = 11),
    panel.grid.minor    = element_blank(),
    panel.grid.major    = element_blank(),
    axis.ticks          = element_blank(),
    legend.position     = "bottom",
    legend.title        = element_blank(),
    legend.text         = element_text(size = 10.5),
    plot.margin         = margin(16, 18, 14, 14)
  )

# ---- RUTAS ------------------------------------------------------------------
ruta   <- "C:/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_r/"
ruta_d <- "C:/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_r/graphics/"

wilcox <- read_delim(paste0(ruta, "wilcoxon_denominador_CLASS.csv"), delim = ";",
                     locale = locale(decimal_mark = ","), show_col_types = FALSE)

# ---- FUNCIÓN: scatter de un nivel (color = clasificación) --------------------
scatter_clasif <- function(datos, nivel, titulo, archivo_png, solo_sig_labels = FALSE) {

  df <- datos %>%
    filter(level == nivel, present_in == "both") %>%
    mutate(
      categoria = case_when(
        significant_bonferroni != "YES"       ~ "No significativa",
        direction == "conflict"               ~ "Conflicto media-mediana",
        direction == "underestimated_in_eDNA" ~ "Subestimado en eDNA",
        direction == "overestimated_in_eDNA"  ~ "Sobreestimado en eDNA"),
      categoria = factor(categoria, levels = c("Subestimado en eDNA", "Sobreestimado en eDNA",
                                               "Conflicto media-mediana", "No significativa")),
      taxon_lab = gsub("\\.", " ", taxon))

  # suelo justo por debajo del menor valor real, para que los 0 no rompan el log10
  piso <- min(c(df$mean_morpho[df$mean_morpho > 0],
                df$mean_cox1[df$mean_cox1 > 0]), na.rm = TRUE) / 3
  df <- df %>% mutate(mean_morpho = pmax(mean_morpho, piso),
                      mean_cox1   = pmax(mean_cox1,   piso))

  # etiquetas: todas, o solo los significativos (todo lo que no sea "No significativa")
  df_labels <- if (solo_sig_labels) filter(df, categoria != "No significativa") else df

  # anotaciones de dirección adaptadas al rango de cada nivel
  lx <- range(log10(df$mean_morpho)); ly <- range(log10(df$mean_cox1))
  pos <- function(r, f) 10^(r[1] + diff(r) * f)

  fig <- ggplot(df, aes(mean_morpho, mean_cox1)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#90a4ae", linewidth = 0.6) +
    geom_point(aes(fill = categoria, color = categoria), shape = 21, size = 4, stroke = 1.4) +
    geom_text_repel(data = df_labels, aes(label = taxon_lab), size = 3, fontface = "italic",
                    color = ink, max.overlaps = Inf, seed = 42) +
    annotate("text", x = pos(lx, 0.90), y = pos(ly, 0.10), label = "media eDNA\n< morfología",
             color = "grey55", size = 2.8, fontface = "italic", alpha = 0.9) +
    annotate("text", x = pos(lx, 0.10), y = pos(ly, 0.90), label = "media eDNA\n> morfología",
             color = "grey55", size = 2.8, fontface = "italic", alpha = 0.9) +
    scale_x_log10(labels = label_percent(accuracy = 0.01)) +
    scale_y_log10(labels = label_percent(accuracy = 0.01)) +
    scale_fill_manual(values = c("Subestimado en eDNA" = teal, "Sobreestimado en eDNA" = amber,
                                 "Conflicto media-mediana" = grey_n, "No significativa" = "white"),
                      name = NULL, drop = FALSE) +
    scale_color_manual(values = c("Subestimado en eDNA" = teal, "Sobreestimado en eDNA" = amber,
                                  "Conflicto media-mediana" = "#8a8a8a", "No significativa" = grey_n),
                       name = NULL, drop = FALSE) +
    labs(title = titulo,
         subtitle = "La posición muestra la MEDIA; el color, la clasificación (media + mediana)",
         x = "Abundancia relativa media — Morfología (%)",
         y = "Abundancia relativa media — eDNA COI (%)") +
    coord_equal() + tema_tfm

  print(fig)
  ggsave(archivo_png, fig, width = 8.4, height = 5.2, dpi = 300, bg = "white")
  cat("Guardado en:", archivo_png, "\n")
  invisible(fig)
}

# ---- GENERAR LOS 5 SCATTERS --------------------------------------------------
# CLASE — se etiquetan todos (22 taxones)
scatter_clasif(wilcox, "Class",
               "Abundancia relativa media por clase: morfología vs eDNA",
               paste0(ruta_d, "disp_clasif_clase.png"),   solo_sig_labels = FALSE)

# ORDEN / FAMILIA / GÉNERO / ESPECIE — solo se etiquetan los significativos
scatter_clasif(wilcox, "Order",
               "Abundancia relativa media por orden: morfología vs eDNA",
               paste0(ruta_d, "disp_clasif_orden.png"),   solo_sig_labels = TRUE)

scatter_clasif(wilcox, "Family",
               "Abundancia relativa media por familia: morfología vs eDNA",
               paste0(ruta_d, "disp_clasif_familia.png"), solo_sig_labels = TRUE)

scatter_clasif(wilcox, "Genus",
               "Abundancia relativa media por género: morfología vs eDNA",
               paste0(ruta_d, "disp_clasif_genero.png"),  solo_sig_labels = TRUE)

scatter_clasif(wilcox, "Species",
               "Abundancia relativa media por especie: morfología vs eDNA",
               paste0(ruta_d, "disp_clasif_especie.png"), solo_sig_labels = TRUE)
