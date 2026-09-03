# =============================================================================
# ANÁLISIS WILCOXON CON NORMALIZACIÓN
# Andrea Moriñigo
# =============================================================================



# =============================================================================
# BLOQUE 0 — LIBRERÍAS
# =============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)



# =============================================================================
# BLOQUE 1 — CARGAR DATOS
# =============================================================================
ruta     <- "C:/Users/andre/Desktop/TFM_AZTI/comparison/archivos/Abundance_tables_Andrea_Metabridge/"
ruta_out <- "C:/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_r/"  
ruta_tax <- "C:/Users/andre/Desktop/TFM_AZTI/comparison/archivos/download_from_python/"

# Morfología — ejemplo resuelto para class, completa el resto tú:
morpho_class   <- read.csv(paste0(ruta, "Class_Morpho_based.csv"),  row.names = 1)
morpho_order   <- read.csv(paste0(ruta, "Order_Morpho_based.csv"),  row.names = 1)
morpho_family  <- read.csv(paste0(ruta, "Family_Morpho_based.csv"),  row.names = 1)
morpho_genus   <- read.csv(paste0(ruta, "Genus_Morpho_based.csv"),  row.names = 1)
morpho_species <- read.csv(paste0(ruta, "Species_Morpho_based.csv"),  row.names = 1)

# eDNA — mismo patrón pero los archivos se llaman _cox1.csv
cox1_class   <- read.csv(paste0(ruta, "Class_cox1.csv"),  row.names = 1)
cox1_order   <- read.csv(paste0(ruta, "Order_cox1.csv"),  row.names = 1)
cox1_family  <- read.csv(paste0(ruta, "Family_cox1.csv"),  row.names = 1)
cox1_genus   <- read.csv(paste0(ruta, "Genus_cox1.csv"),  row.names = 1)
cox1_species <- read.csv(paste0(ruta, "Species_cox1.csv"),  row.names = 1)

# Rutas taxonómicas de todos los tax
# Tabla TAX de eDNA (para filtrar por clase)
edna_path <- read.csv(paste0(ruta_tax, "tax_cox1.csv"), sep=",")
morpho_path <- read.csv(paste0(ruta, "Metabridge_shadow_report_morpho_TAX.csv"), sep = ";")


# Comprobación: ejecuta esto y mira que da 144 en los dos casos
cat("Morphological samples:", nrow(morpho_class), "\n")
cat("eDNA samples:      ", nrow(cox1_class),   "\n")



# =============================================================================
# BLOQUE 2A — LIMPIAR PARÉNTESIS EN NOMBRES DE eDNA
# =============================================================================
clean_names <- function(name) {        # name: el nombre que le doy al argumento, es decir, a lo que va a coger el valor
  name <- gsub("\\s*\\(.*?\\)", "", name)  
  name <- trimws(name)
  return(name)
}

# Aplícala a cada tabla de eDNA con colnames():
# colnames() devuelve los nombres de columna de una tabla
# Si le asignas algo, los cambia: colnames(tabla) <- nuevos_nombres
colnames(cox1_class)   <- clean_names(colnames(cox1_class))
colnames(cox1_order)   <- clean_names(colnames(cox1_order))
colnames(cox1_family)  <- clean_names(colnames(cox1_family))
colnames(cox1_genus)   <- clean_names(colnames(cox1_genus))
colnames(cox1_species) <- clean_names(colnames(cox1_species))


# =============================================================================
# BLOQUE 2B — BORRAR TAXONES INVÁLIDOS
# =============================================================================
# Cargar el Excel:
deleted_list <- read_excel(paste0(ruta, "taxa_to_delete.xlsx"))   # <<< nombre de tu Excel

# Define la función — rellena los ??? :
delete_column <- function(table, col_name) {
  if (col_name %in% colnames(table)) {              # si la columna existe...
    table <- table[, !colnames(table) %in% col_name, drop = FALSE]           # ...bórrala
    cat("  Borrado:", col_name, "\n")
  }
  return(table)   # devuelve la tabla (con o sin la columna borrada)
}


# Bucles para borrar los taxones de cada una de las tablas
#     Extraer listas 
delete_morpho_species <- deleted_list$taxon[   # De la lista "deeltred_list" quiero coger solo lo que esté en la columna taxones que coindida con: "Species" y "morpho"
  deleted_list$level   == "Species" &
    deleted_list$dataset == "morpho" &
    deleted_list$accion  == "Borrar"
]
delete_morpho_genus <- deleted_list$taxon[   
  deleted_list$level   == "Genus" &
    deleted_list$dataset == "morpho" &
    deleted_list$accion  == "Borrar"
]
delete_morpho_family <- deleted_list$taxon[   
  deleted_list$level   == "Family" &
    deleted_list$dataset == "morpho" &
    deleted_list$accion  == "Borrar"
]
delete_morpho_order <- deleted_list$taxon[   
  deleted_list$level   == "Order" &
    deleted_list$dataset == "morpho" &
    deleted_list$accion  == "Borrar"
]
delete_morpho_class <- deleted_list$taxon[   
  deleted_list$level   == "Class" &
    deleted_list$dataset == "morpho" &
    deleted_list$accion  == "Borrar"
]

delete_edna_species <- deleted_list$taxon[
  deleted_list$level   == "Species" &
    deleted_list$dataset == "cox1" &
    deleted_list$accion  == "Borrar"
]
delete_edna_genus <- deleted_list$taxon[
  deleted_list$level   == "Genus" &
    deleted_list$dataset == "cox1" &
    deleted_list$accion  == "Borrar"
]
delete_edna_family <- deleted_list$taxon[
  deleted_list$level   == "Family" &
    deleted_list$dataset == "cox1" &
    deleted_list$accion  == "Borrar"
]
delete_edna_order <- deleted_list$taxon[
  deleted_list$level   == "Order" &
    deleted_list$dataset == "cox1" &
    deleted_list$accion  == "Borrar"
]
delete_edna_class <- deleted_list$taxon[
  deleted_list$level   == "Class" &
    deleted_list$dataset == "cox1" &
    deleted_list$accion  == "Borrar"
]


#     for para eliminar columnas con la función
morpho_species_clean <- morpho_species
for (taxon in delete_morpho_species) {
  morpho_species_clean <- delete_column(morpho_species_clean, taxon)
}

morpho_genus_clean <- morpho_genus
for (taxon in delete_morpho_genus) {
  morpho_genus_clean <- delete_column(morpho_genus_clean, taxon)
}

morpho_family_clean <- morpho_family
for (taxon in delete_morpho_family) {
  morpho_family_clean <- delete_column(morpho_family_clean, taxon)
}

morpho_order_clean <- morpho_order
for (taxon in delete_morpho_order) {
  morpho_order_clean <- delete_column(morpho_order_clean, taxon)
}

morpho_class_clean <- morpho_class
for (taxon in delete_morpho_class) {
  morpho_class_clean <- delete_column(morpho_class_clean, taxon)
}

cox1_species_clean <- cox1_species
for (taxon in delete_edna_species) {
  cox1_species_clean <- delete_column(cox1_species_clean, taxon)
}

cox1_genus_clean <- cox1_genus
for (taxon in delete_edna_genus) {
  cox1_genus_clean <- delete_column(cox1_genus_clean, taxon)
}

cox1_family_clean <- cox1_family
for (taxon in delete_edna_family) {
  cox1_family_clean <- delete_column(cox1_family_clean, taxon)
}

cox1_order_clean <- cox1_order
for (taxon in delete_edna_order) {
  cox1_order_clean <- delete_column(cox1_order_clean, taxon)
}

cox1_class_clean <- cox1_class
for (taxon in delete_edna_class) {
  cox1_class_clean <- delete_column(cox1_class_clean, taxon)
}


# =============================================================================
# BLOQUE 2C — RENOMBRAR TAXONES
# =============================================================================
rename_column <- function(table, old_name, new_name) {
  if (old_name %in% colnames(table)) {
    colnames(table)[colnames(table) == old_name] <- new_name   # <<< rellena
    cat("  Renombrado:", old_name, "->", new_name, "\n")
  }
  return(table)
}

# Extraer lista de renombrados del Excel
rename_morpho_species <- deleted_list[
  deleted_list$level   == "Species"  &
    deleted_list$dataset == "morpho"   &
    deleted_list$accion  == "Renombrar",
]
rename_morpho_genus <- deleted_list[
  deleted_list$level   == "Genus"  &
    deleted_list$dataset == "morpho"   &
    deleted_list$accion  == "Renombrar",
]
rename_morpho_family <- deleted_list[
  deleted_list$level   == "Family"  &
    deleted_list$dataset == "morpho"   &
    deleted_list$accion  == "Renombrar",
]
rename_morpho_order <- deleted_list[
  deleted_list$level   == "Order"  &
    deleted_list$dataset == "morpho"   &
    deleted_list$accion  == "Renombrar",
]
rename_morpho_class <- deleted_list[
  deleted_list$level   == "Class"  &
    deleted_list$dataset == "morpho"   &
    deleted_list$accion  == "Renombrar",
]
rename_cox1_species <- deleted_list[
  deleted_list$level   == "Species"  &
    deleted_list$dataset == "cox1"   &
    deleted_list$accion  == "Renombrar",
]
rename_cox1_genus <- deleted_list[
  deleted_list$level   == "Genus"  &
    deleted_list$dataset == "cox1"   &
    deleted_list$accion  == "Renombrar",
]
rename_cox1_family <- deleted_list[
  deleted_list$level   == "Family"  &
    deleted_list$dataset == "cox1"   &
    deleted_list$accion  == "Renombrar",
]
rename_cox1_order <- deleted_list[
  deleted_list$level   == "Order"  &
    deleted_list$dataset == "cox1"   &
    deleted_list$accion  == "Renombrar",
]
rename_cox1_class <- deleted_list[
  deleted_list$level   == "Class"  &
    deleted_list$dataset == "cox1"   &
    deleted_list$accion  == "Renombrar", # este si que es una tabla porque tiene la coma
]

# Aplicar renombrados con for loop — ejemplo resuelto:
# seq_len(n) genera la secuencia 1, 2, 3, ..., n
# nrow(tabla) devuelve el número de filas de una tabla
for (i in seq_len(nrow(rename_cox1_species))) {    # seq_len coge el valor de un numero, en este caso el nº de la fila
  cox1_species_clean <- rename_column(
    cox1_species_clean,
    rename_cox1_species$taxon[i],        
    rename_cox1_species$nombre_nuevo[i]  
  )
}
for (i in seq_len(nrow(rename_cox1_genus))) {   
  cox1_genus_clean <- rename_column(
    cox1_genus_clean,
    rename_cox1_genus$taxon[i],        
    rename_cox1_genus$nombre_nuevo[i]  
  )
}
for (i in seq_len(nrow(rename_cox1_family))) {   
  cox1_family_clean <- rename_column(
    cox1_family_clean,
    rename_cox1_family$taxon[i],        
    rename_cox1_family$nombre_nuevo[i]  
  )
}
for (i in seq_len(nrow(rename_cox1_order))) {   
  cox1_order_clean <- rename_column(
    cox1_order_clean,
    rename_cox1_order$taxon[i],        
    rename_cox1_order$nombre_nuevo[i]  
  )
}
for (i in seq_len(nrow(rename_cox1_class))) {    
  cox1_class_clean <- rename_column(
    cox1_class_clean,
    rename_cox1_class$taxon[i],        
    rename_cox1_class$nombre_nuevo[i]  
  )
}
for (i in seq_len(nrow(rename_morpho_species))) {   
  morpho_species_clean <- rename_column(
    morpho_species_clean,
    rename_morpho_species$taxon[i],        
    rename_morpho_species$nombre_nuevo[i]  
  )
}
for (i in seq_len(nrow(rename_morpho_genus))) {   
  morpho_genus_clean <- rename_column(
    morpho_genus_clean,
    rename_morpho_genus$taxon[i],        
    rename_morpho_genus$nombre_nuevo[i]  
  )
}
for (i in seq_len(nrow(rename_morpho_family))) {   
  morpho_family_clean <- rename_column(
    morpho_family_clean,
    rename_morpho_family$taxon[i],        
    rename_morpho_family$nombre_nuevo[i]  
  )
}
for (i in seq_len(nrow(rename_morpho_order))) {   
  morpho_order_clean <- rename_column(
    morpho_order_clean,
    rename_morpho_order$taxon[i],        
    rename_morpho_order$nombre_nuevo[i]  
  )
}
for (i in seq_len(nrow(rename_morpho_class))) {   
  morpho_class_clean <- rename_column(
    morpho_class_clean,
    rename_morpho_class$taxon[i],        
    rename_morpho_class$nombre_nuevo[i]  
  )
}

# =============================================================================
# BLOQUE 2D — CASOS ESPECIALES
# =============================================================================
# A) JUNTAR Hexacorallia y Octocorallia en Anthozoa (en morpho_class)
if ("Hexacorallia" %in% colnames(morpho_class_clean) &
    "Octocorallia" %in% colnames(morpho_class_clean)) {
  morpho_class_clean$Anthozoa <- morpho_class_clean$Hexacorallia +
    morpho_class_clean$Octocorallia
  morpho_class_clean <- morpho_class_clean[,
                                           !colnames(morpho_class_clean) %in% c("Hexacorallia", "Octocorallia")]
  cat("  Creada columna Anthozoa (Hexacorallia + Octocorallia)\n")
}

# B) Order to Family
acteonidae_valores <- cox1_order[, "Acteonidae"]        # Guardar los valores de Acteonidae antes de borrarla de Order
cox1_order_clean <- delete_column(cox1_order_clean, "Acteonidae")  # Borrarla de Order
cox1_family_clean$Acteonidae <- acteonidae_valores            # Añadirla a Family con los mismos valores

# C) Order to Family
rhopaluridae_valores <- cox1_order_clean[, "Rhopaluridae..order."]        # Guardar los valores antes de borrarla de Order
cox1_order_clean <- delete_column(cox1_order_clean, "Rhopaluridae..order.")  # Borrarla de Order
cox1_family_clean$Rhopaluridae <- rhopaluridae_valores            # Añadirla a Family con los mismos valores

# D) Genus to Family 
glyceridae_valores <- cox1_genus_clean[, "Glyceridae.LMG"]        # Guardar los valores antes de borrarla de Order
cox1_genus_clean <- delete_column(cox1_genus_clean, "Glyceridae.LMG")  # Borrarla de Order
cox1_family_clean$Glyceridae <- glyceridae_valores            # Añadirla a Family con los mismos valores

# E) FUSIONAR "Spiophanes cf. kroyeri VVP" en "Spiophanes kroyeri" (en cox1_species)
#    Son la misma especie. En la tabla cruda del eDNA coexisten dos columnas
#    (Spiophanes.kroyeri y Spiophanes.cf..kroyeri.VVP). Antes se renombraba la
#    segunda en taxa_to_delete.xlsx, pero como el destino YA existia se generaba
#    una columna duplicada que R rebautizaba como "Spiophanes.kroyeri.1" y acababa
#    contada como una especie exclusiva del eDNA. Se ha quitado esa regla de
#    taxa_to_delete.xlsx y aqui se SUMAN las dos columnas en una sola, igual que
#    se hace con Anthozoa (Hexacorallia + Octocorallia). No se renombra ni se borra
#    en taxa_to_delete: es un caso especial de fusion.
col_cf <- "Spiophanes.cf..kroyeri.VVP"
if (col_cf %in% colnames(cox1_species_clean)) {
  if ("Spiophanes.kroyeri" %in% colnames(cox1_species_clean)) {
    cox1_species_clean$Spiophanes.kroyeri <- cox1_species_clean$Spiophanes.kroyeri +
      cox1_species_clean[, col_cf]
    cat("  Fusionada", col_cf, "en Spiophanes.kroyeri (eDNA)\n")
  } else {
    colnames(cox1_species_clean)[colnames(cox1_species_clean) == col_cf] <- "Spiophanes.kroyeri"
    cat("  ", col_cf, "-> Spiophanes.kroyeri (no existia previamente)\n")
  }
  cox1_species_clean <- delete_column(cox1_species_clean, col_cf)
}



# =============================================================================
# BLOQUE 2E — Filtrado por clases presentes en morfo
# =============================================================================
# Clases que NO están en morpho → hay que eliminar todos sus taxones inferiores
class_NO_morpho <- setdiff(colnames(cox1_class_clean), colnames(morpho_class_clean))         # setdiff: evalua las diferencias (que esté en A pero no en B)

cat("Clases solo en eDNA (se eliminan):", length(class_NO_morpho), "\n")
head(class_NO_morpho, 10)
print(class_NO_morpho)

orders_out   <- unique(na.omit(edna_path[edna_path$Class %in% class_NO_morpho, "Order"]))    # unique: elimina los duplicados
families_out <- unique(na.omit(edna_path[edna_path$Class %in% class_NO_morpho, "Family"]))
genus_out    <- unique(na.omit(edna_path[edna_path$Class %in% class_NO_morpho, "Genus"]))
species_out  <- unique(na.omit(edna_path[edna_path$Class %in% class_NO_morpho, "Species"]))

cox1_order_clean   <- cox1_order_clean[,   !colnames(cox1_order_clean)   %in% orders_out]
cox1_family_clean  <- cox1_family_clean[,  !colnames(cox1_family_clean)  %in% families_out]
cox1_genus_clean   <- cox1_genus_clean[,   !colnames(cox1_genus_clean)   %in% genus_out]
cox1_species_clean <- cox1_species_clean[, !colnames(cox1_species_clean) %in% species_out]

# Comprobación
cat("Órdenes eDNA tras filtro:",   ncol(cox1_order_clean),   "\n")
cat("Familias eDNA tras filtro:",  ncol(cox1_family_clean),  "\n")
cat("Géneros eDNA tras filtro:",   ncol(cox1_genus_clean),   "\n")
cat("Especies eDNA tras filtro:",  ncol(cox1_species_clean), "\n")



# =============================================================================
# BLOQUE 3 — CALCULAR EL DENOMINADOR
# =============================================================================
# Denominador de morfo
denom_morpho <- rowSums(morpho_class_clean,  na.rm = TRUE)                       # suma los valores de todos los taxones de esa fila


# Denominador eDNA
class_in_morpho  <- colnames(morpho_class_clean)
denom_cox1 <- rowSums(cox1_class_clean[,  colnames(cox1_class_clean)  %in% class_in_morpho],  na.rm = TRUE)

# Comprobación (ninguno debería tener ceros)
cat("¿Ceros en el denominador?", any(denom_morpho == 0), "\n")


# =============================================================================
# BLOQUE 4 — CALCULAR ABUNDANCIAS RELATIVAS
# =============================================================================
# --- ANÁLISIS 1: denominador = denom_A1 ---
morpho_class_rel   <- sweep(morpho_class_clean,   1, denom_morpho, FUN = "/")
morpho_order_rel   <- sweep(morpho_order_clean, 1, denom_morpho, FUN = "/")   
morpho_family_rel  <- sweep(morpho_family_clean,   1, denom_morpho, FUN = "/")  
morpho_genus_rel   <- sweep(morpho_genus_clean,   1, denom_morpho, FUN = "/")
morpho_species_rel <- sweep(morpho_species_clean,   1, denom_morpho, FUN = "/")

cox1_class_rel     <- sweep(cox1_class_clean, 1, denom_cox1, FUN = "/")
cox1_order_rel     <- sweep(cox1_order_clean, 1, denom_cox1, FUN = "/")
cox1_family_rel    <- sweep(cox1_family_clean, 1, denom_cox1, FUN = "/")
cox1_genus_rel     <- sweep(cox1_genus_clean, 1, denom_cox1, FUN = "/")
cox1_species_rel   <- sweep(cox1_species_clean, 1, denom_cox1, FUN = "/")

# Comprobación A1: la tabla que tiene el mismo nivel que el denominador
# (morpho_class_rel_A1) debería sumar exactamente 1 en cada fila
cat("Comprobación A1 fila 1:", round(sum(morpho_class_rel[1,], na.rm=TRUE), 4),
    "← debe ser 1\n")




# =============================================================================
# BLOQUE 5 — FUNCIÓN DO_WILCOXON
# =============================================================================
do_wilcoxon <- function(morpho_rel_table, cox1_rel_table, taxonommic_level) {
  
  cat("\n--- Nivel:", taxonommic_level, "---\n")
  
  # nos quedamos para el paired solo con las muestras comunes (row) (para asegurarme que comparo la misma muestra en ambos metodos/por si una muestra se ha perdido)
  common_samples      <- intersect(rownames(morpho_rel_table), rownames(cox1_rel_table))   # me da una lista de las MUESTRAS comunes
  morpho_rel_table    <- morpho_rel_table[common_samples, ]                                # cojo los taxones comunes de morpho y edna en filas y las columnas se quedan todas igual
  cox1_rel_table      <- cox1_rel_table[common_samples, ]
  
  # separo los taxones (col) en 3 listas difernetes: common, only morpho, only cox
  common_taxa         <- intersect(colnames(morpho_rel_table), colnames(cox1_rel_table))   # me da una lista de los TAXONES comunes
  only_morpho_tax     <- setdiff(colnames(morpho_rel_table), colnames(cox1_rel_table))
  only_cox1_tax       <- setdiff(colnames(cox1_rel_table),   colnames(morpho_rel_table))
  
  cat("  En ambos:", length(common_taxa),
      "| Solo morfo:", length(only_morpho_tax),
      "| Solo eDNA:", length(only_cox1_tax), "\n")
  
  common_results <- data.frame()                                                        # se contruye un dataframe para los resultados en común
  for (taxon in common_taxa) {
    morpho_values <- morpho_rel_table[, taxon]                                         # coges todas las filas y de todas esas filas unicamente la columna del taxon que se está analizando
    cox1_values   <- cox1_rel_table[, taxon]
    
    test <- tryCatch(
      wilcox.test(cox1_values, morpho_values, paired = TRUE, exact = FALSE,
      alternative = "two.sided"),
      error = function(e) NULL
    )
    # tryCatch: intenta ejecutar el wilcox.test, pero si por algún motivo falla, en lugar de romper todo el script devuelve NULL y sigue adelante
    
    if (!is.null(test)) {
      med_m <- median(morpho_values, na.rm = TRUE)
      med_c <- median(cox1_values,   na.rm = TRUE)
      mean_m <- mean(morpho_values,   na.rm = TRUE)
      mean_c <- mean(cox1_values,     na.rm = TRUE)

      # Dirección según la MEDIANA y según la MEDIA por separado:
      #   "over" = eDNA mayor que morfología · "under" = eDNA menor · "tie" = empate
      # === CORRECCION (direccion segun la mediana) ===============================
      # El test de Wilcoxon pareado evalua la MEDIANA DE LAS DIFERENCIAS PAREADAS,
      # no la diferencia entre las medianas de cada grupo. Usar
      # median(cox1) - median(morpho) es incorrecto: puede dar direccion cuando
      # la mediana de las diferencias es 0 (y al reves). Se calcula sobre las
      # diferencias muestra a muestra.
      med_dif    <- median(cox1_values - morpho_values, na.rm = TRUE)
      dir_median <- ifelse(med_dif > 0, "over",
                           ifelse(med_dif < 0, "under", "tie"))
      dir_mean   <- ifelse(mean_c > mean_m, "over",
                           ifelse(mean_c < mean_m, "under", "tie"))

      # Clasificación oficial combinando MEDIA Y MEDIANA:
      #   - si mediana y media coinciden en dirección -> esa dirección (over/under)
      #   - si la mediana empata (típicamente mediana = 0 en AMBOS métodos) la
      #     dirección la decide la media -> over/under. Estos casos se conservan
      #     como sobre/subestimados aunque tengan mediana = 0 (varios de ellos la
      #     tienen; se comenta en el texto).
      #   - si la media empata, la decide la mediana
      #   - si mediana y media apuntan a direcciones OPUESTAS -> "conflict"
      #     (los empates de mediana en 0 quedan excluidos: van arriba como over/under)
      #   - si ninguna marca dirección (empate total, muy raro) -> "no_difference"
      clasif <- ifelse(dir_median == dir_mean, dir_median,
                ifelse(dir_median == "tie",    dir_mean,
                ifelse(dir_mean   == "tie",    dir_median,
                       "conflict")))
      direccion_final <- switch(clasif,
                                over     = "overestimated_in_eDNA",
                                under    = "underestimated_in_eDNA",
                                tie      = "no_difference",
                                conflict = "conflict")

      common_results <- rbind(common_results, data.frame(
        level         = taxonommic_level,
        taxon         = taxon,
        present_in    = "both",
        p_value       = test$p.value,
        median_morpho = med_m,
        median_cox1   = med_c,
        median_diff   = med_dif,                                  # mediana de las diferencias pareadas (cox1 - morpho)
        mean_morpho   = mean_m,
        mean_cox1     = mean_c,
        direction               = direccion_final,                # clasificación oficial (media + mediana)
        direction_mean          = dir_mean,                       # dirección solo según la media (over/under/tie)
        conflicto_media_mediana = (direccion_final == "conflict") # TRUE solo cuando es conflicto de dirección
      ))
    }
  }
  
  only_morpho_results <- data.frame(
    level         = taxonommic_level,
    taxon         = only_morpho_tax,
    present_in    = "only_morphology",
    p_value       = 0,
    median_morpho = sapply(only_morpho_tax, function(t) median(morpho_rel_table[, t], na.rm=TRUE)),
    median_cox1   = 0,
    median_diff   = NA,
    mean_morpho   = sapply(only_morpho_tax, function(t) mean(morpho_rel_table[, t],   na.rm=TRUE)),
    mean_cox1     = 0,
    direction     = "absent_in_eDNA",
    direction_mean = NA,
    conflicto_media_mediana = NA
  )
  
  only_cox1_results <- data.frame(
    level         = taxonommic_level,
    taxon         = only_cox1_tax,
    present_in    = "only_eDNA",
    p_value       = 0,
    median_morpho = 0,
    median_cox1   = sapply(only_cox1_tax, function(t) median(cox1_rel_table[, t], na.rm=TRUE)),
    median_diff   = NA,
    mean_morpho   = 0,
    mean_cox1     = sapply(only_cox1_tax, function(t) mean(cox1_rel_table[, t],   na.rm=TRUE)),
    direction     = "absent_in_morphology",
    direction_mean = NA,
    conflicto_media_mediana = NA
  )
  
  return(rbind(common_results, only_morpho_results, only_cox1_results))
}


# =============================================================================
# BLOQUE 6 — EJECUTAR EL WILCOXON PARA CADA ANÁLISIS
# =============================================================================
# --- ANÁLISIS 1 ---
cat("\n===== ANÁLISIS 1: denominador CLASS =====\n")

class_results     <- do_wilcoxon(morpho_class_rel, cox1_class_rel, "Class")
order_results     <- do_wilcoxon(morpho_order_rel, cox1_order_rel, "Order")
family_results    <- do_wilcoxon(morpho_family_rel, cox1_family_rel, "Family")
genus_results     <- do_wilcoxon(morpho_genus_rel, cox1_genus_rel, "Genus")
species_results   <- do_wilcoxon(morpho_species_rel, cox1_species_rel, "Species")

# Unir todos los niveles en una tabla:
results <- rbind(
  class_results,
  order_results,
  family_results,
  genus_results,
  species_results
)



# =============================================================================
# BLOQUE 7 — CORRECCIÓN DE P-VALORES
# =============================================================================
# Solo corregimos los que tienen p-valor real (present_in == "both"):
idx_1 <- !is.na(results$p_value) & results$p_value != 0

results$p_adjust_BH         <- NA
results$p_adjust_bonferroni <- NA

results$p_adjust_BH[idx_1] <- p.adjust(
  results$p_value[idx_1], method = "BH")

results$p_adjust_bonferroni[idx_1] <- p.adjust(
  results$p_value[idx_1], method = "bonferroni")

results$significant_BH <- ifelse(
  !is.na(results$p_adjust_BH) & results$p_adjust_BH < 0.05,
  "YES", "NO")

results$significant_bonferroni <- ifelse(
  !is.na(results$p_adjust_bonferroni) & results$p_adjust_bonferroni < 0.05,
  "YES", "NO")

cat("A1 — significativos Bonferroni:", 
    sum(results$significant_bonferroni == "YES", na.rm=TRUE), "\n")




# =============================================================================
# BLOQUE 8 — TABLAS DE RESUMENES
# =============================================================================
order_levels <- c("Class", "Order", "Family", "Genus", "Species")

summary_per <- results %>%
  group_by(level, present_in) %>%
  summarise(n         = n(), 
            sig_BH    = sum(significant_BH == "YES", na.rm = TRUE),
            sig_bon   = sum(significant_bonferroni == "YES", na.rm = TRUE),
            .groups   = "drop") %>%                                              # summarise: todas las filas del grupo en una / n(): cuenta las filas del grupo
  
  group_by(level) %>%
  mutate(                                                                        # añade columnas nuevas sin alterar las filas
    total_level = sum(n),
    porcentaje  = round(n / total_level * 100, 1)
  ) %>%
  ungroup() %>%
  arrange(factor(level, levels = order_levels))                                  # arrange ordena las filas de la tabla. Por defecto R las ordena alfabéticamente



# Under/Over ratio
# Ahora "direction" tiene 3 categorías para los "both": overestimated_in_eDNA,
# underestimated_in_eDNA y conflict (media y mediana en direcciones opuestas).
# Los sobre/subestimados con mediana = 0 (empate de mediana desempatado por la
# media) siguen contando como over/under; solo son "conflict" los que discrepan
# en dirección con mediana distinta de 0.
direcciones <- results %>%
  filter(present_in == "both", significant_bonferroni == "YES") %>%
  group_by(level, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(level) %>%
  mutate(
    total_sig = sum(n),
    pct       = round(n / total_sig * 100, 1)
  ) %>%
  ungroup() %>%
  arrange(factor(level, levels = order_levels))

print(direcciones)
write.csv(direcciones, paste0(ruta_out, "direcciones_significativos.csv"), row.names = FALSE)


# =============================================================================
# TABLA RESUMEN — direcciones (over / under / conflict) + mediana 0
# =============================================================================
# Objetivo: una fila por nivel taxonómico con el total de taxones significativos
# (Bonferroni), cuántos están sobreestimados / subestimados en eDNA y cuántos
# son "conflict" (media y mediana discrepan en dirección). Además, de los
# sobre/subestimados, cuántos tienen mediana = 0 en AMBOS métodos: son casos
# reales de sesgo cuya dirección la marca la media (el test es significativo por
# las pocas muestras con valores, pero la mediana es 0 en la mayoría). Ese
# recuento es el que se comenta en el texto.
# library(writexl) se añade solo para exportar a Excel al final de este bloque.
library(writexl)

tabla_direcciones <- results %>%
  filter(present_in == "both", significant_bonferroni == "YES") %>%
  mutate(
    over_under   = direction %in% c("overestimated_in_eDNA", "underestimated_in_eDNA"),
    mediana_cero = (median_morpho == 0 & median_cox1 == 0)
  ) %>%
  group_by(level) %>%
  summarise(
    total                   = n(),
    sobreestimados          = sum(direction == "overestimated_in_eDNA"),
    subestimados            = sum(direction == "underestimated_in_eDNA"),
    conflicto               = sum(direction == "conflict"),
    over_under_mediana_cero = sum(over_under & mediana_cero),
    .groups = "drop"
  ) %>%
  arrange(factor(level, levels = order_levels))

# Fila de totales (todos los niveles juntos)
tabla_direcciones <- bind_rows(
  tabla_direcciones,
  tabla_direcciones %>%
    summarise(
      level                   = "TOTAL",
      total                   = sum(total),
      sobreestimados          = sum(sobreestimados),
      subestimados            = sum(subestimados),
      conflicto               = sum(conflicto),
      over_under_mediana_cero = sum(over_under_mediana_cero)
    )
)

print(tabla_direcciones)

write_xlsx(tabla_direcciones, paste0(ruta_out, "tabla_direcciones_comparacion.xlsx"))
cat("✓ Guardado: tabla_direcciones_comparacion.xlsx\n")


# Tabla de clases que están solo en eDNA
only_edna_class <- results %>%
  filter(level == "Class", present_in == "only_eDNA") %>%
  select(taxon, mean_cox1, median_cox1) %>%
  arrange(desc(mean_cox1))
print(only_edna_class)


# Tabla de 1er nivel donde empiezan a faltar taxones en eDNA
missing_by_level <- list(
  Class   = results$taxon[results$level=="Class"   & results$present_in=="only_morphology"],
  Order   = results$taxon[results$level=="Order"   & results$present_in=="only_morphology"],
  Family  = results$taxon[results$level=="Family"  & results$present_in=="only_morphology"],
  Genus   = results$taxon[results$level=="Genus"   & results$present_in=="only_morphology"],
  Species = results$taxon[results$level=="Species" & results$present_in=="only_morphology"]
)


table3 <- data.frame()

for (lv in order_levels) {
  
  # taxa missing at the current level
  taxa_at_level <- missing_by_level[[lv]]
  
  # all levels above the current one (e.g. if lv = "Family", higher_levels = Class, Order)
  position_lv    <- match(lv, order_levels)
  higher_levels  <- order_levels[seq_len(position_lv - 1)]
  
  for (tx in taxa_at_level) {
    
    # species names use dots instead of spaces,
    # we need to convert them before searching in morpho_path
    if (lv == "Species") {
      tx_to_search <- gsub("\\.", " ", tx)
      row_found <- morpho_path[morpho_path[["Species.full"]] == tx_to_search, ]
    } else {
      row_found <- morpho_path[morpho_path[[lv]] == tx, ]
    }
    
    # check, one by one, if any level above this taxon is also missing
    already_explained <- FALSE
    
    for (lv_up in higher_levels) {
      parent_taxon <- row_found[[lv_up]][1]
      
      if (!is.na(parent_taxon) && parent_taxon %in% missing_by_level[[lv_up]]) {
        already_explained <- TRUE
      }
    }
    
    # only keep it if no higher level already explains this absence
    if (already_explained == FALSE) {
      new_row <- data.frame(first_missing_level = lv, taxon = tx)
      table3  <- rbind(table3, new_row)
    }
  }
}

cat("Total rows in table3:", nrow(table3), "\n")
print(table(table3$first_missing_level))



# Tabla de 1er nivel donde empiezan a faltar taxones en MORFO
# start with the full list of missing species
table4 <- data.frame(level = "Species", taxon = missing_by_level$Species)

levels_to_check <- c("Genus", "Family", "Order", "Class")

for (current_level in levels_to_check) {
  
  # the level right below the current one
  position_current <- match(current_level, order_levels)
  lower_level       <- order_levels[position_current + 1]
  
  taxa_at_level <- missing_by_level[[current_level]]
  
  for (tx in taxa_at_level) {
    
    # find all the "children" of this missing taxon at the level below
    children <- morpho_path[morpho_path[[current_level]] == tx, lower_level]
    children <- unique(children)
    children <- children[children %in% missing_by_level[[lower_level]]]
    
    # remove those children from table4, they will be replaced by the parent
    table4 <- table4[!(table4$level == lower_level & table4$taxon %in% children), ]
    
    # add the parent instead
    new_row <- data.frame(level = current_level, taxon = tx)
    table4  <- rbind(table4, new_row)
  }
}

cat("Total rows in table4:", nrow(table4), "\n")
print(table(table4$level))



# =============================================================================
# BLOQUE 9 — GUARDAR RESULTADOS
# =============================================================================
write.csv2(results,
           paste0(ruta_out, "wilcoxon_denominador_CLASS.csv"),
           row.names = FALSE)

write.csv2(summary_per,
           paste0(ruta_out, "summary_percentage_CLASS.csv"),
           row.names = FALSE)


write.csv2(only_edna_class,
           paste0(ruta_out, "class_only_eDNA.csv"),
           row.names = FALSE)


write.csv2(table3, paste0(ruta_out, "table3_first_absence_in_eDNA.csv"), row.names = FALSE)


write.csv2(table4, paste0(ruta_out, "table4_missing_in_eDNA_bottomup.csv"), row.names = FALSE)



# Exportar especies exclusivas de morfología para midori_comparison.py
morpho_especies_para_MIDORI <- results %>%
  filter(
    level == "Species",
    present_in == "only_morphology",
  )

cat("Especies para comparar con MIDORI:", nrow(morpho_especies_para_MIDORI), "\n")

write.csv2(morpho_especies_para_MIDORI,
           paste0(ruta_out, "sp_morpho_compare.csv"),
           row.names = FALSE)


gen_morpho_compare <- results %>%
  filter(level == "Genus", present_in == "only_morphology")

cat("Géneros ausentes en eDNA:", nrow(gen_morpho_compare), "\n")   # debe dar 399

write.csv2(gen_morpho_compare,
           paste0(ruta_out, "gen_morpho_compare.csv"),
           row.names = FALSE)

# --- RECÍPROCO: especies y géneros detectados SOLO por eDNA (ausentes en morfo) ---
# Para la Tabla 4 (tab2_2_perdidos_ambi): además de cuántos taxones PIERDE el
# eDNA que son indicadores AMBI/NSI, se añade cuántos de los taxones que el
# eDNA APORTA (no vistos en morfología) son también indicadores. Mismo filtro
# de "results" que sp_morpho_compare/gen_morpho_compare pero con
# present_in == "only_eDNA" en vez de "only_morphology".
sp_edna_compare <- results %>%
  filter(level == "Species", present_in == "only_eDNA")

cat("Especies detectadas solo por eDNA:", nrow(sp_edna_compare), "\n")   # debe dar 105

write.csv2(sp_edna_compare,
           paste0(ruta_out, "sp_edna_compare.csv"),
           row.names = FALSE)

gen_edna_compare <- results %>%
  filter(level == "Genus", present_in == "only_eDNA")

cat("Géneros detectados solo por eDNA:", nrow(gen_edna_compare), "\n")   # debe dar 97

write.csv2(gen_edna_compare,
           paste0(ruta_out, "gen_edna_compare.csv"),
           row.names = FALSE)

# UNA LISTA DE CADA NIVEL TAXONÓMICO PARA PYTHON (worms_tax_builder.py):
# DESCOMENTADO: estas escrituras estaban apagadas y por eso los nombres_cox1_*.csv
# quedaban congelados (arrastrando "Spiophanes.kroyeri.1"). Ahora se regeneran a
# partir de las tablas ya curadas y fusionadas.
write.csv(data.frame(taxon = colnames(cox1_class_clean)),
          paste0(ruta_out, "nombres_cox1_class.csv"),   row.names = FALSE)
write.csv(data.frame(taxon = colnames(cox1_order_clean)),
          paste0(ruta_out, "nombres_cox1_order.csv"),   row.names = FALSE)
write.csv(data.frame(taxon = colnames(cox1_family_clean)),
          paste0(ruta_out, "nombres_cox1_family.csv"),  row.names = FALSE)
write.csv(data.frame(taxon = colnames(cox1_genus_clean)),
          paste0(ruta_out, "nombres_cox1_genus.csv"),   row.names = FALSE)
write.csv(data.frame(taxon = colnames(cox1_species_clean)),
          paste0(ruta_out, "nombres_cox1_species.csv"), row.names = FALSE)