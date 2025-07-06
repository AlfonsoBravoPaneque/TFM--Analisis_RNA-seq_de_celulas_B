#### ---- Paquetes necesarios ---- ####

# Paquetes
library(DESeq2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(extrafont)


loadfonts(device = "win") # carga las fuentes


#### ---- Heatmaps ---- ####

# En el análisis de expresión diferencial se empleo le función DESeq que normaliza internamente los valores calculando los factores de normalización por muestra
# y asignando estos factores al objeto dds y normalizando mediante estos factores las cuentas de cada gen. Sin embargo, en el heatmap es recomendable emplear
# otras normalizaciones/tranformaciones que permiten comparar mejor los valores entre si. En este caso se empleó la normalización vst.


#### ---- Heatmap con trasnformación vst: genes más expresados con rowMeans por condición y tipo celular ---- ####

# Se aplica la transformación VST para estabilizar la varianza (permite que los datos tenga una varianza homogenea)
vsd <- vst(dds, blind=FALSE)
View(assay(vsd)) # se visualiza la transformación
mat <- assay(vsd) # se obtiene la matriz transformada

# Seguidamente, se usa la función scale (Z.score), normalizando por filas (gen)
mat_scaled <- t(scale(t(mat)))

# Se crea un data.frame con las anotaciones
df <- as.data.frame(colData(dds)[, c("condition", "cell_type")])

# Se reodena la matriz
df <- df[order(df$cell_type, df$condition), ]
mat_scaled <- mat_scaled[, rownames(df)]

# Se filtran los genes sin nombre antes de seleccionar los más expresados
clean_gene_ids <- gsub("\\..*", "", rownames(mat_scaled))  # Limpiar IDs Ensembl
gene_map <- setNames(all_genes$Gene.Name, all_genes$Gene.ID)  # Mapeo ID a Nombre
selected_gene_names <- gene_map[clean_gene_ids]  # Obtener nombres de genes

# Se filtran solo los genes que tienen un nombre real
genes_con_nombre <- !is.na(selected_gene_names)  # Vector lógico TRUE/FALSE
mat_scaled <- mat_scaled[genes_con_nombre, ]  # Mantener solo genes con nombre
selected_gene_names <- selected_gene_names[genes_con_nombre]  # Mantener nombres correcto

# Se asignan los nombres a la matriz
rownames(mat_scaled) <- selected_gene_names

# Después, se crea un vector vacio para alamcenar los resultados del bucle for. Basicamente, se selecciona los genes más expresados en cada grupo
# (teniendo en cuenta la condición y tipo celular: hay 6 grupos) calculando la media por filas de la matriz mat_scaled
top_genes_per_group <- c()

for (group in unique(paste(df$cell_type, df$condition))) {
  samples_in_group <- rownames(df)[paste(df$cell_type, df$condition) == group]
  group_means <- rowMeans(mat_scaled[, samples_in_group, drop = FALSE])
  top_genes <- names(sort(group_means, decreasing = TRUE))[1:100]  # Se seleccionan los 100 genes más expresados por grupo
  top_genes_per_group <- c(top_genes_per_group, top_genes)
}

# Observar si hay más de un gen que se expresa por grupo
gene_counts <- table(top_genes_per_group)
genes_recurrentes <- names(gene_counts[gene_counts > 1])
print(genes_recurrentes)  # Lista de genes presentes en más de un grupo

length(top_genes_per_group) # número de genes presentes en el heatmap

# Se crea la nueva matriz con los genes seleccionado
mat_scaled_selected <- mat_scaled[top_genes_per_group, ]

# Se personaliza las leyendas del heatmap
df <- setNames(as.data.frame(colData(dds)[, c("condition", "cell_type")]), 
               c("Condición Experimental", "Tipo Celular"))

# Es normal que los genes más expresasos por condición y tipo celular en el heatmap no coincidan con los genes más en el volcano plot. Esto es debido a que
# en el heatmap se grafican los genes mas expresados teniendo en cuanta la expresión absoluta (empleando la función rowMeans). En cambio, en el volcano plot
# se representan los genes diferencialmente expresados por condición y tipo celular, teniendo en cuenta tanto el log2FC y el p valor.


df <- df[match(colnames(mat_scaled_selected), rownames(df)), ]  

# Definir colores para las anotaciones. Además, se oculta la anotación
col_annotation <- HeatmapAnnotation(
  "Condición Experimental" = df$`Condición Experimental`, 
  "Tipo Celular" = df$`Tipo Celular`,
  col = list(
    "Condición Experimental" = c(Aicda_PO = "#56B4E9", Aicda_KO = "#D55E00"),
    "Tipo Celular" = c(Naive = "cyan", Activated = "#FFD700", GC = "#32CD32")
  ),
  gp = gpar(col = "black"),  # Agregar bordes a cada caja de color
  show_legend = FALSE,  # Ocultar leyendas de las anotaciones
  annotation_name_gp = gpar(col = NA)
)


# Crear heatmap sin mostrar la leyenda del heatmap
ht <- Heatmap(mat_scaled_selected, 
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              show_row_names = FALSE,
              column_names_rot = 50,  # Rotar nombres de columnas
              col = colorRampPalette(c("blue", "white", "red"))(50),
              top_annotation = col_annotation,  # Mantener anotaciones arriba
              show_heatmap_legend = FALSE,  # Ocultar la leyenda del heatmap
              column_names_gp = gpar(fontsize = 12, fontfamily = "Calibri") 
)


# Crear leyendas manuales
ht_legend <- Legend(title = "Nivel de Expresión", 
                    col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")), 
                    title_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = "Calibri"),
                    labels_gp = gpar(fontsize = 9, fontfamily = "Calibri"))

legend_condicion <- Legend(title = "Condición Experimental", 
                           labels = c("Aicda_PO", "Aicda_KO"), 
                           legend_gp = gpar(fill = c("#56B4E9", "#D55E00")),
                           title_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = "Calibri"),
                           labels_gp = gpar(fontsize = 9, fontfamily = "Calibri"))

legend_tipo <- Legend(title = "Tipo Celular", 
                      labels = c("Naive", "Activated", "GC"), 
                      legend_gp = gpar(fill = c("cyan", "#FFD700", "#32CD32")),
                      title_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = "Calibri"),
                      labels_gp = gpar(fontsize = 9, fontfamily = "Calibri"))


# Dibujar heatmap sin leyendas automáticas
draw(ht, padding = unit(c(40, 30, 10, 10), "mm")) 


# 1250x1250 (mejor este tamaño). Cuadra perfecto con fuente calibri
draw(ht_legend, x = unit(290.5, "mm"), y = unit(290, "mm"))   # Leyenda del heatmap en
draw(legend_condicion, x = unit(290.7, "mm"), y = unit(261.1, "mm"))  # Condición Experimental
draw(legend_tipo, x = unit(291, "mm"), y = unit(236, "mm"))  # Tipo Celular

