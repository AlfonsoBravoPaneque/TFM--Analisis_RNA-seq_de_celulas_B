#### ---- Paquetes necesarios ---- ####

# Paquetes
library(DESeq2)
library(dplyr)
library(apeglm)
library(ashr)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(stringr)
library(tidyr)
library(openxlsx)
library(extrafont)


loadfonts(device = "win") # carga las fuentes


#### ---- Busqueda de genes en cada subpoblación usando los IDs de Ensembl (mejor realizarlo asi) ---- ####

# Conectarse a Ensembl versión 109 (basada en GRCm39)
ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "mmusculus_gene_ensembl",
  version = 109
)

# Obtener coordenadas de todos los genes
attributes <- c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position")
all_genes <- getBM(attributes = attributes, mart = ensembl)

# Eliminar genes sin nombre
all_genes <- all_genes[all_genes$external_gene_name != "", ]

# Renombrar columnas para su integración posterior
colnames(all_genes)[1] <- "Gene.ID"
colnames(all_genes)[2] <- "Gene.Name"

# Visualizar los primeros genes
head(all_genes)


#### ---- Se carga la matriz de cuentas y se elimina la región contaminada ---- ####

# Se carga la matriz de cuentas con las muestras KO (knock-out) y PO (wild-type)
counts_data <- read.table("//wsl.localhost/Ubuntu/home/alfonso/results/counts/gene_counts_ALL_SAMPLES.txt", header = TRUE, row.names = 1, sep = "\t")

# Se definen las regiones problemáticas del cromosoma 6 (desde 108524660 (inicio de la banda E2) hasta 139322197 (final de la banda G1))
contaminated_start <- 108524660
contaminated_end <- 139322197

# Se define la región del gen Aicda presente en el cromosoma 6
aicda_start <- 122530760
aicda_end <- 122541139

# Se identifican las columnas de muestras Knock Out (KO) y Wild Type (PO)
ko_samples <- grep("KO", colnames(counts_data), value = TRUE)
po_samples <- grep("PO", colnames(counts_data), value = TRUE)


# Se crea una copia de la matriz de cuentas original
filtered_counts_data <- counts_data

# Se la primera coordenada en caso de que haya múltiples valores separados por ";"
counts_data$Start_clean <- as.numeric(gsub(";.*", "", as.character(counts_data$Start)))
counts_data$End_clean <- as.numeric(gsub(";.*", "", as.character(counts_data$End)))
counts_data$Chr_clean <- gsub(";.*", "", as.character(counts_data$Chr))


# Aplicar el filtrado usando las columnas limpias
genes_filt <- with(counts_data, 
                   Chr_clean == "chr6" & 
                     Start_clean >= contaminated_start & End_clean <= contaminated_end & 
                     !(Start_clean >= aicda_start & End_clean <= aicda_end))


# Se visualizan los genes filtrados de las muestras KO
filtered_genes <- counts_data[genes_filt, c("Chr", "Start", "End")]
View(filtered_genes)


# Elimina completamente las filas contaminadas de la matriz
filtered_counts_data <- filtered_counts_data[!genes_filt, ]



#### ---- Se filtra el gen Aicda en las muestras Knock-out ---- ####
# Al realizar los análisis de expresión diferencial, me percaté que había expresión del gen Aicda en los ratones Knock-out. Esto podría explicarse por el hecho de que
# al hacer el alineamiento se reconocio partes del gen Aicda, ya que el Knock-out solo se inhibía la region funcional del gen y no se eliminaba la secuencaia como tal. Por tanto,
# se le intodujo un valor de 0 a las muestras Knock-out de cada subpoblación celular.

# Se extrae el ID del gen Aicda
aicda_gene_id <- "ENSMUSG00000040627.15"

# Además, se elimina las cuentas del gen Aicda en las muestras Knock-out
if (aicda_gene_id %in% rownames(filtered_counts_data)) {
  # Asignar valores de 0 a Aicda en las muestras KO
  filtered_counts_data[aicda_gene_id, ko_samples] <- 0
} else {
  print("El gen Aicda no se encontró en la matriz de cuentas en las muestras Knock-out.")
}

# Además, finalmente se verifica que Aicda solo está presente en las muestras PO
filtered_counts_data[aicda_gene_id, ]



#### ---- Filtrado de genes ribosomales, RNA ribosómicos y de transferencia---- ####

# En la mayoría de los estudios de expresión diferencial, es recomendable eliminar los genes ribosomales antes del análisis. Esto mejora la normalización, evita sesgos y 
# permite una mejor identificación de genes relevantes. Además, se eliminan los RNAs ribosómicos y de tranferencia, ya que solo nos interesan los RNA mensajeros.

# Se guardan los rownames originales antes del merge
filtered_counts_data$Original.ID <- rownames(filtered_counts_data)

# Se crea una columna Gene.ID sin el sufijo para hacer el merge correctamente
filtered_counts_data$Gene.ID <- sub("\\..*", "", filtered_counts_data$Original.ID)

# Se realiza el merge con la anotación de genes
counts_data_filt_merge <- merge(
  x = filtered_counts_data,
  y = all_genes,
  by = "Gene.ID",
  all.x = TRUE
)

# Se restaura el orden original según filtered_counts_data
counts_data_filt_merge <- counts_data_filt_merge[match(filtered_counts_data$Gene.ID, counts_data_filt_merge$Gene.ID), ]

# Se restaura los rownames originales
rownames(counts_data_filt_merge) <- filtered_counts_data$Original.ID

# Se filtran las proteínas ribosomales nucleares y mitocondriales, ya que suelen estar constitutivamente expresadas y pueden sesgar los resultados en el análisis de expresión diferencial. 
# Además, también se eliminan algunas RNAs ribosómicos (rRNA) y de transferencia (tRNA) presentes en la matriz de cuentas, que son RNAs altamente abundantes y de función estructural,
# por lo que no suelen aportar información útil para este tipo de análisis.

# Los genes ribosomales nucleares y mitocondriales suelen eliminarse previo al análisis de expresión diferencial, ya que presentan una expresión constitutiva elevada que puede introducir sesgos
# en los resultados. Asimismo, los RNAs ribosómicos (rRNA) y los RNAs de transferencia (tRNA) se excluyen por ser RNA constitutivos y no informativos para este estudio transcriptómico.


# Se emplea la función grepl para hacer el filtrado. Se utiliza \\. (evita que el . se entienda como cualquier caracter). Además, se usa el ignore.case = TRUE para
# que no se distinga entre mayúscula y minúscula para buscar coincidencias
genes_delete <- grepl("^rpl|^rps|^mrpl|^mrps|^trna|rrna|^mt-r|^rn[0-9]|^n-r[0-9]|rna5s|18s|28s|5\\.8s|5s", counts_data_filt_merge$Gene.Name, ignore.case = TRUE)
sum(genes_delete) # número de genes eliminados (597)

counts_data_filtrado <- counts_data_filt_merge[!genes_delete, ]


#### ---- Preparación de datos para DESeq2 ---- ####

# Se selecciona la columna de la cuentas, y se modifican los nombres de las muetras para favorecer el manejo
counts <- dplyr::select(counts_data_filtrado, starts_with("X.home.alfonso.results.alignment"))

# Se limpian los nombres de las muestras
colnames(counts) <- gsub("X.home.alfonso.results.alignment.", "", colnames(counts))
colnames(counts) <- gsub("_sorted.bam", "", colnames(counts))

# Se convierte el data.frame a matriz:
counts_matrix <- as.matrix(counts)

# Verificar las dimensiones de la matriz de cuentas (donde las filas son los genes y las columnas las muestras)
dim(counts_matrix)

# En DESeq2, es fundamental que la matriz de cuentas y la tabla de metadatos (colData) tengan un orden consistente de muestras, porque DESeq2 no 
# infiere automáticamente qué fila de colData corresponde a qué columna en la matriz de cuentas. Por tanto:

#   -Las columnas de la matriz de cuentas representan las muestras.
#   -Las filas de la tabla de metadatos también representan esas mismas muestras.
#   -El orden de las muestras en las columnas de la matriz de cuentas debe de coincidir exactamente
#    con el orden de la filas en la tabla de metadatos.

metadata <- data.frame(condition = rep(c("Aicda_KO", "Aicda_PO"), c(1,1)), cell_type = rep(c("Activated", "GC", "Naive"), c(8,4,8)), row.names = colnames(counts_matrix))

# Se convierte a factor las variables categóricas
metadata$condition = as.factor(metadata$condition)
metadata$cell_type = as.factor(metadata$cell_type)

# Comprobar si se ha realizado bien la conversión a factor
str(metadata)

# Verificar si el orden de las muestras en las columnas de la matriz de cuentas coincide exactamente con el orden de la filas en la tabla de metadatos.
all(rownames(metadata) %in% colnames(counts_matrix)) # comprueba si las muestras presentes en metadata están en la matriz de cuentas (debe de dar TRUE)
all(rownames(metadata) == colnames(counts_matrix)) # comprueba el orden (debe de dar TRUE)


#### ---- Análisis de expresión diferencial con DESeq2 y filtrado de datos ---- ####

# Se crea un nuevo factor que combine condición y tipo celular, ya que se quiere comparar el tipo de celular según la condición experimental
metadata$group <- factor(paste(metadata$condition, metadata$cell_type, sep = "_"))


# Se genera el objeto dds (DESeqDataSet) con el conteo y los metadatos, especificando como diseño experimental el nuevo factor generado. Es más sencillo hacer esto
# que realizando interacción entre tipo celular y condición experimental (cell_type:condition), ya que facilita los análisis posteriores
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = metadata,
                              design = ~ group)


# Debido a que en el experimento hay tres subpoblaciones celulares y en cada una de ellas hay dos condiciones (Aicda+ y Aicda-), con un número
# variable de réplicas por grupo. Para el filtrado de genes se seleccionó la subpoblación celular con el tamaño del grupo más pequeño (2 réplicas para el caso
# de las células del centro germinal)

# Definir el tamaño del grupo más pequeño
smallestGroupSize <- 2

# Filtrar genes que tengan al menos 10 cuentas en al menos "smallestGroupSize" muestras
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep, ]

# Se comprueva el número de genes eliminados (numero de filas)
dim(dds)

# Debido a que mis replicas son replicas biológicas (provenientes de distintos ratones) no se usa la función collapseReplicates para agrupar las
# diferentes réplicas

# Se ejecuta la función DESeq(), la cual lleva a cabo el análisis de expresión diferencial aplicando un modelo basado en la distribución Binomial Negativa.
# Esta función realiza los siguientes pasos:
#
#  - Estima los factores de tamaño: ajusta las diferencias en la profundidad de secuenciación entre muestras utilizando la mediana de razones geométricas.
#                                    Esto normaliza los datos sin sesgar la expresión relativa de los genes.
#
#  - Estima la dispersión: modela la variabilidad biológica de cada gen en función de su nivel de expresión medio. Para ello,
#                          ajusta una curva suavizada a la relación media-dispersión y aplica un enfoque bayesiano empírico
#                          para estabilizar las estimaciones en genes con pocos conteos.
#
#  - Ajusta un modelo GLM y realiza la prueba de Wald: emplea un modelo lineal generalizado (GLM) con distribución Binomial Negativa para estimar los cambios de expresión.
#                                                      Luego, aplica la prueba de Wald para evaluar la significancia estadística de los cambios, generando valores de log2 
#                                                      fold change y p-values.
#
#  - Ajusta los valores p por pruebas múltiples: utiliza el método de Benjamini-Hochberg para controlar la tasa de falsos descubrimientos (FDR).
#
#  - Opcionalmente, permite aplicar shrinkage en los log2 fold changes (lfcShrink()): reduce la sobreestimación de los cambios en genes de baja expresión,
#                                                                                     mejorando la interpretabilidad de los resultados.

dds <- DESeq(dds) # ya normaliza los datos directamente

# La funciónresultsNames() permite obtener los nombres de los coeficientes del modelo
resultsNames(dds)

# Seguidamente, se usa la función results(), que pemite extraer los resultados del análisis de expresión diferencial. Además, se define un alpha de 0,05


#### ---- Células Naive -> expresión diferencial ---- ####

# No hace falta definir la referencia, ya que con contrast se indican los grupos a comparar. Se comparan la Naive_Aicda_KO vs Naive_Aicda_PO
res_naive <- results(dds, contrast = c("group", "Aicda_KO_Naive", "Aicda_PO_Naive"), alpha = 0.05)
mcols(res_naive)$description # proporciona información sobre las columnas de los metadatos del objeto res_naive


#### ---- Células GC -> expresión diferencial ---- ####

# No hace falta definir la referencia, ya que con contrast se indican los grupos a comparar. Se comparan la GC_Aicda_KO vs GC_Aicda_PO
res_GC <- results(dds, contrast = c("group", "Aicda_KO_GC", "Aicda_PO_GC"), alpha = 0.05)
mcols(res_GC)$description # proporciona información sobre las columnas de los metadatod del objeto res_GC


#### ---- Células Activated -> expresión diferencial ---- ####

# No hace falta definir la referencia, ya que con contrast se indican los grupos a comparar. Se comparan la Activated_Aicda_KO vs Activated_Aicda_PO
res_activated <- results(dds, contrast = c("group", "Aicda_KO_Activated", "Aicda_PO_Activated"), alpha = 0.05)
mcols(res_activated)$description # proporciona información sobre las columnas de los metadatod del objeto res_activated


#### ---- Aplicación de lfcShrink() a cada cubpoblación celular ---- ####

# Adicionlamente, se aplicó el Shrinkage en los log2 fold changes, permitiendo mejorar la interpretabilidad de los resultados y reducir el ruido (favoreciendo
# principalmente a la interpretabilidad de los genes con baja expresión). Esto se realizó mediante la función lfcShrink(), que permite reducir la
# variabilidad en las estimaciones de log2 fold change (LFC). En este caso, se empleó el método de reducción del log2FoldChange llamado "ashr", el cual funciona
# cuando se usa contrast

res_naive_shrink <- lfcShrink(dds, contrast = c("group", "Aicda_KO_Naive", "Aicda_PO_Naive"), type = "ashr")
res_GC_shrink <- lfcShrink(dds, contrast = c("group", "Aicda_KO_GC", "Aicda_PO_GC"), type = "ashr")
res_activated_shrink <- lfcShrink(dds, contrast = c("group", "Aicda_KO_Activated", "Aicda_PO_Activated") ,type = "ashr")



#### ---- Procesado de Células Naive e introducción de genes (nombres comunes) ---- ####

# Se seleccionan los genes que poseen un p valor ajustado no igual a NA y los que poseen un p valor ajustado menor a 0.05. Además,
# se añade el ID de Ensembl presente en las filas como una nueva columna del data.frame, y se eliminan los sufijos de los IDs
res_naive_shrink_df <- as.data.frame(res_naive_shrink)
res_naive_shrink_df_filt <- res_naive_shrink_df[!is.na(res_naive_shrink_df$padj) & res_naive_shrink_df$padj < 0.05, ]
res_naive_shrink_df_filt$Gene.ID <- rownames(res_naive_shrink_df_filt)
res_naive_shrink_df_filt$Gene.ID <- sub("\\..*", "", res_naive_shrink_df_filt$Gene.ID)

# Se fusiona el data.frame filtrado (res_naive_shrink_df_filt) con las coordenadas genéticas de Ensembl. Además, se
# eliminan los genes que no poseen Gene.Name
res_naive_shrink_df_filt_merge <- merge(
  x = res_naive_shrink_df_filt,
  y = all_genes,
  by = "Gene.ID",
  all.x = TRUE)                     

res_naive_shrink_df_filt_merge <- na.omit(res_naive_shrink_df_filt_merge)



#### ---- Procesado de Células GC e introducción de genes (nombres comunes) ---- ####

# Se seleccionan los genes que poseen un p valor ajustado no igual a NA y los que poseen un p valor ajustado menor a 0.05. Además,
# se añade el ID de Ensembl presente en las filas como una nueva columna del data.frame, y se eliminan los sufijos de los IDs
res_GC_shrink_df <- as.data.frame(res_GC_shrink)
res_GC_shrink_df_filt <- res_GC_shrink_df[!is.na(res_GC_shrink_df$padj) & res_GC_shrink_df$padj < 0.05, ]
res_GC_shrink_df_filt$Gene.ID <- rownames(res_GC_shrink_df_filt)
res_GC_shrink_df_filt$Gene.ID <- sub("\\..*", "", res_GC_shrink_df_filt$Gene.ID)

# Se fusiona el data.frame filtrado (res_GC_shrink_df_filt) con las coordenadas genéticas de Ensembl. Además, se
# eliminan los genes que no poseen Gene.Name
res_GC_shrink_df_filt_merge <- merge(
  x = res_GC_shrink_df_filt,
  y = all_genes,
  by.x = "Gene.ID",
  all.x = TRUE) 

res_GC_shrink_df_filt_merge <- na.omit(res_GC_shrink_df_filt_merge)



#### ---- Procesado de Células Activated e introducción de genes (nombres comunes) ---- ####

# Se seleccionan los genes que poseen un p valor ajustado no igual a NA y los que poseen un p valor ajustado menor a 0.05. Además,
# se añade el ID de Ensembl presente en las filas como una nueva columna del data.frame, y se eliminan los sufijos de los IDs
res_activated_shrink_df <- as.data.frame(res_activated_shrink)
res_activated_shrink_df_filt <- res_activated_shrink_df[!is.na(res_activated_shrink_df$padj) & res_activated_shrink_df$padj < 0.05, ]
res_activated_shrink_df_filt$Gene.ID <- rownames(res_activated_shrink_df_filt)
res_activated_shrink_df_filt$Gene.ID <- sub("\\..*", "", res_activated_shrink_df_filt$Gene.ID)

# Se fusiona el data.frame filtrado (res_activated_shrink_df_filt) con las coordenadas genéticas de Ensembl. Además, se
# eliminan los genes que no poseen Gene.Name
res_activated_shrink_df_filt_merge <- merge(
  x = res_activated_shrink_df_filt,
  y = all_genes,
  by = "Gene.ID",
  all.x = TRUE) 

res_activated_shrink_df_filt_merge <- na.omit(res_activated_shrink_df_filt_merge)


#### ---- Volcano plots ---- ####

# Se define el threshold del pvalue y del log2FoldChange
pvalue_threshold = 0.05
FC_threshold = 0.5


#### ---- Células GC -> Volcano plot ---- ####

# Se introducen los nombres de los genes al data.frame de expresión diferencial SIN FILTRAR. Después, se eliminan los sufijos
# de los IDs, y se eliminan los valores NA (se eliminan las filas que poseen valores NA.
res_GC_shrink_df$Gene.ID <- rownames(res_GC_shrink_df)
res_GC_shrink_df$Gene.ID <- sub("\\..*", "", res_GC_shrink_df$Gene.ID)

markers_GC_names <- merge(
  x = res_GC_shrink_df,
  y = all_genes,
  by = "Gene.ID",
  all.x = TRUE) 

markers_GC_names <- na.omit(markers_GC_names) # se eliminan tanto genes que no han sido asignados al hacer
# la búsqueda de su nombre en Ensembl como los genes con un p valor de NA

# Además, se conservan los genes que poseen un p valor no igual a 0, ya que los genes con p valor igual a 0 ocasionan problemas de graficado
# (un p valor de 0 da error al transformarlos con el log10(p_value))
markers_GC_filt <- markers_GC_names %>% filter(pvalue != 0)

# Se agrega una columna llamada State, y se pone Up o down según el log2FoldChange
markers_GC_filt <- markers_GC_filt %>%
  mutate(State = case_when(
    log2FoldChange > 0 ~ "up",
    log2FoldChange < 0 ~ "down"
  ))

# Se crea un data.frame para los upregulated
upregulated_GC <- markers_GC_filt %>% filter(State == "up", pvalue < pvalue_threshold, log2FoldChange > FC_threshold)

# Me quedo con los mejores 5 genes más expresados o genes que yo quiera seleccionar
top_upregulated_GC <- upregulated_GC %>% arrange(-log2FoldChange) %>% head( n = 5)
top_upregulated_GC <- markers_GC_filt[markers_GC_filt$Gene.Name %in% c("H2ac1", "Fcrl5", "Wnt10a", "H1f4", "Ighm", "Ighd", "Adora2a", "Hck"), ]

# Se crea un data.frame para los downregulated
downregulated_GC <- markers_GC_filt %>% filter(State == "down", pvalue < pvalue_threshold, log2FoldChange < -FC_threshold)

# Me quedo con los mejores 5 genes menos expresados o genes que yo quiera seleccionar
top_down_GC <- downregulated_GC %>% arrange(-log2FoldChange) %>% tail( n = 5)
top_down_GC <- markers_GC_filt[markers_GC_filt$Gene.Name %in% c("Aicda", "Ighg2c", "Ighe", "Ighg2b", "Ptpn3", "Ccnb1", "Idh2"), ]

# Se crea el plot base sin color
volcano1 <- ggplot(data = markers_GC_filt, mapping = aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(alpha = 0.1, shape = 1, color = "#2E294E") + 
  theme(text = element_text(family = "Arial", size = 12),  # Cambiar la fuente a Arial y tamaño 12
        axis.title.y = element_text(family = "Arial", size = 12, margin = margin(r = 10)),  # Margen derecho para el título del eje Y
        axis.title.x = element_text(family = "Arial", size = 12, margin = margin(t = 10)),  # Margen superior para el título del eje X
        axis.text = element_text(family = "Arial", size = 12))  # Fuente Arial y tamaño 12 para texto de ejes

# Se agrega líneas para marcar los thresholds del pvalue
volcano2 <- volcano1 + geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", linewidth = 1)

# Se generar um vector para marcar el threshold del log2FoldChange
vertical_lines <- c(-FC_threshold, FC_threshold)
volcano3 <- volcano2 + geom_vline(xintercept = vertical_lines, linetype = "dashed", linewidth = 1)

# Se agregan capas de puntos con color para los genes sobre expresados y menos expresados
volcano4 <- volcano3 + 
  geom_point(data = upregulated_GC, mapping = aes(color = "Aicda_KO"), shape = 16, size = 1.5, alpha = 0.7) +
  geom_point(data = downregulated_GC, mapping = aes(color = "Aicda_PO"), shape = 16, size = 1.5, alpha = 0.7)

# Se agregan capas de etiquetas para genes sobre expresados
volcano5 <- volcano4 + 
  geom_label_repel(
    data = top_upregulated_GC, 
    mapping = aes(label = Gene.Name),
    max.overlaps = 50,
    size = 4,  # Ajusta el tamaño del texto
    box.padding = 1.5, # Espacio alrededor de la etiqueta
    point.padding = 0.1, # Espacio entre la etiqueta y el punto
    segment.color = 'grey50', # Color de la línea que conecta la etiqueta con el punto
    segment.size = 0.5, # Grosor de la línea
    fill = alpha("white", 0.75), # Fondo blanco con transparencia
    color = "black",  # Color del texto
    force = 4,
    family = "Calibri"  # Fuente Calibri
  )

# Se agregan capas de etiquetas para genes menos expresados
volcano6 <- volcano5 + 
  geom_label_repel(
    data = top_down_GC, 
    mapping = aes(label = Gene.Name),
    max.overlaps = 50,
    size = 4,  # Ajusta el tamaño del texto
    box.padding = 1.7, # Espacio alrededor de la etiqueta
    point.padding = 0.1, # Espacio entre la etiqueta y el punto
    segment.color = 'grey50', # Color de la línea que conecta la etiqueta con el punto
    segment.size = 0.5, # Grosor de la línea
    fill = alpha("white", 0.75), # Fondo blanco con transparencia
    color = "black",  # Color del texto
    force = 4,
    family = "Calibri"  # Fuente Calibri
  )

# Se mejora el eje x, definiendo el valor absoluto máximo de log2FoldChange presente en nuestro dataframe
absolute_max_x <- markers_GC_filt$log2FoldChange %>% abs() %>% max() %>% ceiling()

# Se agrega otra capa
volcano7 <- volcano6 + 
  scale_x_continuous(limits = c(-absolute_max_x, absolute_max_x), 
                     breaks = seq(floor(-absolute_max_x/2)*2, ceiling(absolute_max_x/2)*2, by = 2)) 

# Se limpia el plot y se agrega el título y la leyenda
volcano8 <- volcano7 + 
  labs(title = "Expresión diferencial de GC según la condición", color = "Condición") + 
  scale_color_manual(values = c("Aicda_KO" = "#D55E00", "Aicda_PO" = "#56B4E9")) +
  theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5, family = "Calibri", size = 14),
    legend.title = element_text(family = "Calibri", size = 12,  face = "bold"),  # Cambiar la fuente y tamaño de la leyenda
    legend.text = element_text(family = "Calibri", size = 12),  # Cambiar la fuente y tamaño del texto de la leyenda
    text = element_text(family = "Calibri", size = 12),  # Cambiar la fuente a Arial y tamaño 12
    axis.title.y = element_text(family = "Calibri", size = 12, margin = margin(r = 10)),  # Margen derecho para el título del eje Y
    axis.title.x = element_text(family = "Calibri", size = 12, margin = margin(t = 10)),  # Margen superior para el título del eje X
    axis.text = element_text(family = "Calibri", size = 12)
  )

# Se visualiza el plot final (698x523)
volcano8


# Se exportan los genes sobre expresados y menos expresados (sin genes ribosomales) con p valor menor a 0.05, y con log2FoldChange mayor o menor a 0.05 respectivamente. Se ordenan los genes por Log2FoldChange
upregulated_GC_sorted <- upregulated_GC %>% arrange(-log2FoldChange)
write.xlsx(upregulated_GC_sorted, file = "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Volcano plots/GC/Genes_sobreexpresados_GC_Aicda_KO.xlsx", rowNames = FALSE)

downregulated_GC_sorted <- downregulated_GC %>% arrange(log2FoldChange)
write.xlsx(downregulated_GC_sorted , file = "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Volcano plots/GC/Genes_sobreexpresados_GC_Aicda_PO.xlsx", rowNames = FALSE)



#### ---- Células Activated -> Volcano plot---- ####

# Se introducen los nombres de los genes al data.frame de expresión diferencial SIN FILTRAR. Después, se eliminan los sufijos
# de los IDs, y se eliminan los valores NA (se eliminan las filas que poseen valores NA.
res_activated_shrink_df$Gene.ID <- rownames(res_activated_shrink_df)
res_activated_shrink_df$Gene.ID <- sub("\\..*", "", res_activated_shrink_df$Gene.ID)

markers_activated_names <- merge(
  x = res_activated_shrink_df,
  y = all_genes,
  by = "Gene.ID",
  all.x = TRUE) 

markers_activated_names <- na.omit(markers_activated_names) # se eliminan tanto genes que no han sido asignados al hacer
# la búsqueda de su nombre en Ensembl como los genes con un p valor de NA


# Además, se conservan los genes que poseen un p valor no igual a 0, ya que los genes con p valor igual a 0 ocasionan problemas de graficado
# (un p valor de 0 da error al transformarlos con el log10(p_value))
markers_activated_filt <- markers_activated_names %>% filter(pvalue != 0)


# Se agrega una columna llamada State, y se pone Up o down según el log2FoldChange
markers_activated_filt <- markers_activated_filt %>%
  mutate(State = case_when(
    log2FoldChange > 0 ~ "up",
    log2FoldChange < 0 ~ "down"
  ))


# Se crea un data.frame para los upregulated
upregulated_activated <- markers_activated_filt %>% filter(State == "up", pvalue < pvalue_threshold, log2FoldChange > FC_threshold)

# Me quedo con los mejores 5 genes más expresados o genes que yo quiera seleccionar
top_upregulated_activated <- upregulated_activated %>% arrange(-log2FoldChange) %>% head( n = 5)
top_upregulated_activated <- markers_activated_filt[markers_activated_filt$Gene.Name %in% c("H2ac1", "Wdfy1", "Ighd", "Eps8l1", "Prlr"), ]

# Se crea un data.frame para los downregulated
downregulated_activated <- markers_activated_filt %>% filter(State == "down", pvalue < pvalue_threshold, log2FoldChange < -FC_threshold)

# Me quedo con los mejores 5 genes menos expresados o genes que yo quiera seleccionar
top_down_activated <- downregulated_activated %>% arrange(-log2FoldChange) %>% tail( n = 5)
top_down_activated <- markers_activated_filt[markers_activated_filt$Gene.Name %in% c("Aicda", "Ifi208", "Ighg2b", "Ighg3", "Pirb", "Pik3cg", "Irs2"), ]


# Se crea el plot base sin color
volcano1 <- ggplot(data = markers_activated_filt, mapping = aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(alpha = 0.1, shape = 1, color = "#2E294E") + 
  theme(text = element_text(family = "Arial", size = 12),  # Cambiar la fuente a Arial y tamaño 12
        axis.title.y = element_text(family = "Arial", size = 12, margin = margin(r = 10)),  # Margen derecho para el título del eje Y
        axis.title.x = element_text(family = "Arial", size = 12, margin = margin(t = 10)),  # Margen superior para el título del eje X
        axis.text = element_text(family = "Arial", size = 12))  # Fuente Arial y tamaño 12 para texto de ejes

# Se agrega líneas para marcar los thresholds del pvalue
volcano2 <- volcano1 + geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", linewidth = 1)

# Se generar um vector para marcar el threshold del log2FoldChange
vertical_lines <- c(-FC_threshold, FC_threshold)
volcano3 <- volcano2 + geom_vline(xintercept = vertical_lines, linetype = "dashed", linewidth = 1)

# Se agregan capas de puntos con color para los genes sobre expresados y menos expresados
volcano4 <- volcano3 + 
  geom_point(data = upregulated_activated, mapping = aes(color = "Aicda_KO"), shape = 16, size = 1.5, alpha = 0.7) +
  geom_point(data = downregulated_activated, mapping = aes(color = "Aicda_PO"), shape = 16, size = 1.5, alpha = 0.7)

# Se agregan capas de etiquetas para genes sobre expresados
volcano5 <- volcano4 + 
  geom_label_repel(
    data = top_upregulated_activated, 
    mapping = aes(label = Gene.Name),
    max.overlaps = 50,
    size = 4,  # Ajusta el tamaño del texto
    box.padding = 2.2, # Espacio alrededor de la etiqueta
    point.padding = 0.1, # Espacio entre la etiqueta y el punto
    segment.color = 'grey50', # Color de la línea que conecta la etiqueta con el punto
    segment.size = 0.5, # Grosor de la línea
    fill = alpha("white", 0.75), # Fondo blanco con transparencia
    color = "black",  # Color del texto
    force = 4,
    family = "Calibri"  # Fuente Calibri
  )

# Se agregan capas de etiquetas para genes menos expresados
volcano6 <- volcano5 + 
  geom_label_repel(
    data = top_down_activated, 
    mapping = aes(label = Gene.Name),
    max.overlaps = 50,
    size = 4,  # Ajusta el tamaño del texto
    box.padding = 2.2, # Espacio alrededor de la etiqueta
    point.padding = 0.1, # Espacio entre la etiqueta y el punto
    segment.color = 'grey50', # Color de la línea que conecta la etiqueta con el punto
    segment.size = 0.5, # Grosor de la línea
    fill = alpha("white", 0.75), # Fondo blanco con transparencia
    color = "black",  # Color del texto
    force = 4.3,
    family = "Calibri"  # Fuente Calibri
  )

# Se mejora el eje x, definiendo el valor absoluto máximo de log2FoldChange presente en nuestro dataframe
absolute_max_x <- markers_activated_filt$log2FoldChange %>% abs() %>% max() %>% ceiling()

# Se agrega otra capa
volcano7 <- volcano6 + 
  scale_x_continuous(limits = c(-absolute_max_x, absolute_max_x), 
                     breaks = seq(-absolute_max_x, absolute_max_x, by = 2))

# Se limpia el plot y se agrega el título y la leyenda
volcano8 <- volcano7 + 
  labs(title = "Expresión diferencial de Activated según la condición", color = "Condición") + 
  scale_color_manual(values = c("Aicda_KO" = "#D55E00", "Aicda_PO" = "#56B4E9")) +
  theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5, family = "Calibri", size = 14),
    legend.title = element_text(family = "Calibri", size = 12,  face = "bold"),  # Cambiar la fuente y tamaño de la leyenda
    legend.text = element_text(family = "Calibri", size = 12),  # Cambiar la fuente y tamaño del texto de la leyenda
    text = element_text(family = "Calibri", size = 12),  # Cambiar la fuente a Arial y tamaño 12
    axis.title.y = element_text(family = "Calibri", size = 12, margin = margin(r = 10)),  # Margen derecho para el título del eje Y
    axis.title.x = element_text(family = "Calibri", size = 12, margin = margin(t = 10)),  # Margen superior para el título del eje X
    axis.text = element_text(family = "Calibri", size = 12)
  )


# Se visualiza el plot final (700x460)
volcano8


# Se exportan los genes sobre expresados y menos expresados (sin genes ribosomales) con p valor menor a 0.05, y con log2FoldChange mayor o menor a 0.05 respectivamente. Se ordenan los genes por Log2FoldChange
upregulated_activated_sorted <- upregulated_activated %>% arrange(-log2FoldChange)
write.xlsx(upregulated_activated_sorted, file = "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Volcano plots/Activated/Genes_sobreexpresados_Activated_Aicda_KO.xlsx", rowNames = FALSE)

downregulated_activated_sorted <- downregulated_activated %>% arrange(log2FoldChange)
write.xlsx(downregulated_activated_sorted , file = "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Volcano plots/Activated/Genes_sobreexpresados_Activated_Aicda_PO.xlsx", rowNames = FALSE)


