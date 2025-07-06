#### ---- Redes empleando GENIE3 (más adecuado cuando hay pocas muestras) ---- ####

# Paquetes
library(DESeq2)
library(dplyr)
library(tidyr)
library(openxlsx)
library(GENIE3)
library(RCy3)

loadfonts(device = "win") # carga las fuentes


# Es importante fijar una semilla de aleatorización que permite la reproducibilidad
set.seed(1234)

# Además, es mejor correr el código poco a poco ya que se se corrre todo el script directamente pueden aparecer
# fallos en la aplicación cytoscape asociadas al estilo de la red

#### ---- Matriz de expresión transformada y limpiada---- ####

# Transformación VST
vsd <- vst(dds, blind = FALSE)
expr_matrix <- assay(vsd)

# Se limpian los IDs de Ensmbl de la matriz de expresión, y también de all_genes (data.frame obtenido de Biomart que asocia cada ID de
# Ensmbl con un nombre del gen)
rownames(expr_matrix) <- sub("\\..*", "", rownames(expr_matrix))
all_genes$Gene.ID <- sub("\\..*", "", all_genes$Gene.ID)


#### ---- # Cargar los TFs desde AnimalTFDB ---- ####

# Me bajé una lista más completa que la proporcionadas por JASPAR2022 para descrubrir los genes que codifican a TF
# en mis datos de expresión diferencial
tfs_df <- read.delim("C:/Users/Alfonso/Desktop/Asignaturas máster/Segundo cuatrimestre/TFM/TFs/Mus_musculus_TF.txt", header = TRUE, sep = "\t")
tfs_symbols <- unique(na.omit(tfs_df$Symbol))


#### ---- Datos DEGS GC con filtro de padj ---- ####

# Se filtra solo con el p valor adj para aumentar el número de genes presentes en la red
res_GC_shrink_df_padj <- res_GC_shrink_df[!is.na(res_GC_shrink_df$padj) & 
                                                  res_GC_shrink_df$padj < 0.1, ]

res_GC_shrink_df_padj$Gene.ID <- rownames(res_GC_shrink_df_padj)
res_GC_shrink_df_padj$Gene.ID <- sub("\\..*", "", res_GC_shrink_df_padj$Gene.ID)


# Se fusiona el data.frame (res_GC_shrink_df_ppi) con las coordenadas genéticas de Ensembl. Además, se
# eliminan los genes que no poseen Gene.Name
res_GC_shrink_df_padj_merge <- merge(
  x = res_GC_shrink_df_padj,
  y = all_genes,
  by = "Gene.ID",
  all.x = TRUE)  

res_GC_shrink_df_padj_merge <- na.omit(res_GC_shrink_df_padj_merge)





#### ---- Red GENIE3 GC con los genes Aicda_PO expresados ---- ####

# Se seleccionan las muestras GC Aicda_PO de la matrix trasnformada
muestras_GC_PO <- colnames(expr_matrix)[grepl("GC.*Aicda_PO", colnames(expr_matrix))]

# Se seleccionan los DEGs deseados: solo con log2FC <  0 (genes downregulados AICDA_PO) y padj < 0.1 (filtrado anteriormente)
res_down_GC_PO <- res_GC_shrink_df_padj_merge %>%
  filter(log2FoldChange < 0)

genes_down_GC_PO <- res_down_GC_PO$Gene.Name # se listan los genes (nombre)

gene_ids_down_GC_PO <- all_genes$Gene.ID[all_genes$Gene.Name %in% genes_down_GC_PO] # se listan los genes (ID de Ensmbl)

# Se filtra seleccionan de la matriz de expresión los genes DEGs y las muestras Aicda_PO
expr_deg_down_GC_PO <- expr_matrix[rownames(expr_matrix) %in% gene_ids_down_GC_PO, muestras_GC_PO]

# MAPEAR SYMBOLS PARA LA MATRIZ
expr_annotated_GC_PO <- merge(
  data.frame(Gene.ID = rownames(expr_deg_down_GC_PO), expr_deg_down_GC_PO),
  all_genes[, c("Gene.ID", "Gene.Name")],
  by = "Gene.ID"
)

rownames(expr_annotated_GC_PO) <- expr_annotated_GC_PO$Gene.Name
expr_named_GC_PO <- expr_annotated_GC_PO[, !(colnames(expr_annotated_GC_PO) %in% c("Gene.ID", "Gene.Name"))]
expr_named_GC_PO <- as.matrix(expr_named_GC_PO)

# Sebuscan los genes reguladores (TFs) presentes en la matriz
regulators_GC_PO <- intersect(rownames(expr_named_GC_PO), tfs_symbols)

if (length(regulators_GC_PO) < 2) {
  stop("No hay suficientes TFs para correr GENIE3.")
}

# Se ejecuta GENIE3
weightMatrix_GC_PO <- GENIE3(expr_named_GC_PO, regulators = regulators_GC_PO, nCores = 2)
linkList_GC_PO <- getLinkList(weightMatrix_GC_PO)
colnames(linkList_GC_PO) <- c("regulatoryGene", "targetGene", "weight")

# Se filtra la linklist: peso > 0.01 y top 5 por TF (5 conexiones por TF)
filtered_links_GC_PO <- linkList_GC_PO %>%
  filter(weight > 0.01) %>%
  group_by(regulatoryGene) %>%
  slice_max(order_by = weight, n = 5) %>%
  ungroup()

# Se exporta el resultado a Cytoscape para su visualización
top_links_cyto_GC_PO <- filtered_links_GC_PO %>%
  rename(source = regulatoryGene, target = targetGene) %>%
  mutate(interaction = "regulates", name = paste(source, "regulates", target))

# Se guarda el resultado
write.xlsx(top_links_cyto_GC_PO, "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Redes nuevas/GC/GENIE3_GC_PO_DEGs_Tabla.xlsx")

createNetworkFromDataFrames(
  edges = top_links_cyto_GC_PO,
  title = "GENIE3_GC_PO_DEGs",
  collection = "GENIE3_PO"
)

# Se anotan los nodos para aplicar después el color a cada uno de ellos
all_nodes_GC_PO <- unique(c(top_links_cyto_GC_PO$source, top_links_cyto_GC_PO$target))

node_type_df_GC_PO <- data.frame(
  id = all_nodes_GC_PO,
  type = ifelse(all_nodes_GC_PO %in% top_links_cyto_GC_PO$source, "TF", "Target")
)

loadTableData(node_type_df_GC_PO, data.key.column = "id", table.key.column = "name")


# Se modifica el estilo visual de la red y se aplica color a los nodos
style_name <- "GENIE3_Purple_Yellow"

createVisualStyle(style_name,
                  defaults = list(
                    NODE_SHAPE = "ellipse",
                    NODE_BORDER_WIDTH = 1,
                    EDGE_STROKE_UNSELECTED_PAINT = "#AAAAAA",
                    NETWORK_BACKGROUND_PAINT = "#FFFFFF"
                  )
)

setNodeColorMapping("type",
                    table.column.values = c("TF", "Target"),
                    colors = c("gold", "#55B1FD"),
                    mapping.type = "d",
                    style.name = style_name
)

setNodeLabelMapping("name", style.name = style_name)


setCurrentNetwork("GENIE3_GC_PO_DEGs") # se setea la red antes de aplicar el estilo

setVisualStyle(style_name)




#### ---- Red GENIE3 GC con los genes Aicda_KO expresados: en este red vamos a tener que seleccionar con detenimiento los TF relevantes, ya que hay bastantes ---- ####

# Se seleccionan las muestras GC Aicda_PO de la matrix trasnformada
muestras_GC_KO <- colnames(expr_matrix)[grepl("GC.*Aicda_KO", colnames(expr_matrix))]

# Se seleccionan los DEGs deseados: solo con log2FC > 0 (genes upregulados AICDA_KO) y padj < 0.1 (filtrado anteriormente)
res_up_GC_KO <- res_GC_shrink_df_padj_merge %>%
  filter(log2FoldChange > 0)

genes_up_GC_KO <- res_up_GC_KO$Gene.Name # se listan los genes (nombre)

gene_ids_up_GC_KO <- all_genes$Gene.ID[all_genes$Gene.Name %in% genes_up_GC_KO] # se listan los genes (ID de Ensmbl)

# Se filtra seleccionan de la matriz de expresión los genes DEGs y las muestras Aicda_KO
expr_deg_up_GC_KO <- expr_matrix[rownames(expr_matrix) %in% gene_ids_up_GC_KO, muestras_GC_KO]

# MAPEAR SYMBOLS PARA LA MATRIZ
expr_annotated_GC_KO <- merge(
  data.frame(Gene.ID = rownames(expr_deg_up_GC_KO), expr_deg_up_GC_KO),
  all_genes[, c("Gene.ID", "Gene.Name")],
  by = "Gene.ID"
)

rownames(expr_annotated_GC_KO) <- expr_annotated_GC_KO$Gene.Name
expr_named_GC_KO <- expr_annotated_GC_KO[, !(colnames(expr_annotated_GC_KO) %in% c("Gene.ID", "Gene.Name"))]
expr_named_GC_KO <- as.matrix(expr_named_GC_KO)

# Sebuscan los genes reguladores (TFs) presentes en la matriz
regulators_GC_KO <- intersect(rownames(expr_named_GC_KO), tfs_symbols)

if (length(regulators_GC_KO) < 2) {
  stop("No hay suficientes TFs para correr GENIE3.")
}

# Se ejecuta GENIE3
weightMatrix_GC_KO <- GENIE3(expr_named_GC_KO, regulators = regulators_GC_KO, nCores = 2)
linkList_GC_KO <- getLinkList(weightMatrix_GC_KO)
colnames(linkList_GC_KO) <- c("regulatoryGene", "targetGene", "weight")

# Se filtra la linklist: peso > 0.01 y top 5 por TF (5 conexiones por TF)
filtered_links_GC_KO <- linkList_GC_KO %>%
  filter(weight > 0.01) %>%
  group_by(regulatoryGene) %>%
  slice_max(order_by = weight, n = 5) %>%
  ungroup()

# Se exporta el resultado a Cytoscape para su visualización
top_links_cyto_GC_KO <- filtered_links_GC_KO %>%
  rename(source = regulatoryGene, target = targetGene) %>%
  mutate(interaction = "regulates", name = paste(source, "regulates", target))


# Se guarda el resultado
write.xlsx(top_links_cyto_GC_KO, "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Redes nuevas/GC/GENIE3_GC_KO_DEGs_Tabla.xlsx")

##

#Para seleccionar los TFs concretos:
selected_TFs <- c("Egr1", "Irf4", "Batf", "Myc", "Rel", "Relb", "Nfkb1", "Stat3", "Foxp1", "Klf9", "Bcl11a", "Zbtb32", "Egr3", "Stat6")


filtered_links_selected <- filtered_links_GC_KO %>%
  filter(regulatoryGene %in% selected_TFs) %>%
  rename(source = regulatoryGene, target = targetGene) %>%
  mutate(interaction = "regulates",
         name = paste(source, "regulates", target))

# Se exporta el resultado a Cytoscape para su visualización
top_links_cyto_GC_KO <- filtered_links_GC_KO %>%
  filter(regulatoryGene %in% selected_TFs) %>%
  rename(source = regulatoryGene, target = targetGene) %>%
  mutate(interaction = "regulates",
         name = paste(source, "regulates", target))

# Se guarda el resultado
write.xlsx(top_links_cyto_GC_KO, "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Redes nuevas/GC/GENIE3_GC_KO_DEGs_Tabla_2.xlsx")

##


createNetworkFromDataFrames(
  edges = top_links_cyto_GC_KO,
  title = "GENIE3_GC_KO_DEGs",
  collection = "GENIE3_KO"
)

# Se anotan los nodos para aplicar después el color a cada uno de ellos
all_nodes_GC_KO <- unique(c(top_links_cyto_GC_KO$source, top_links_cyto_GC_KO$target))

node_type_df_GC_KO <- data.frame(
  id = all_nodes_GC_KO,
  type = ifelse(all_nodes_GC_KO %in% top_links_cyto_GC_KO$source, "TF", "Target")
)

loadTableData(node_type_df_GC_KO, data.key.column = "id", table.key.column = "name")


# Se modifica el estilo visual de la red y se aplica color a los nodos
style_name <- "GENIE3_Purple_Yellow"

createVisualStyle(style_name,
                  defaults = list(
                    NODE_SHAPE = "ellipse",
                    NODE_BORDER_WIDTH = 1,
                    EDGE_STROKE_UNSELECTED_PAINT = "#AAAAAA",
                    NETWORK_BACKGROUND_PAINT = "#FFFFFF"
                  )
)

setNodeColorMapping("type",
                    table.column.values = c("TF", "Target"),
                    colors = c("gold", "#55B1FD"),
                    mapping.type = "d",
                    style.name = style_name
)

setNodeLabelMapping("name", style.name = style_name)

setCurrentNetwork("GENIE3_GC_KO_DEGs")

setVisualStyle(style_name)





#### ---- Datos DEGS Activated con filtro de padj ---- ####

# Se filtra solo con el p valor adj para aumentar el número de genes presentes en la red
res_activated_shrink_df_padj <- res_activated_shrink_df[!is.na(res_activated_shrink_df$padj) & 
                                                  res_activated_shrink_df$padj < 0.1, ]

res_activated_shrink_df_padj$Gene.ID <- rownames(res_activated_shrink_df_padj)
res_activated_shrink_df_padj$Gene.ID <- sub("\\..*", "", res_activated_shrink_df_padj$Gene.ID)


# Se fusiona el data.frame (res_activated_shrink_df_ppi) con las coordenadas genéticas de Ensembl. Además, se
# eliminan los genes que no poseen Gene.Name
res_activated_shrink_df_padj_merge <- merge(
  x = res_activated_shrink_df_padj,
  y = all_genes,
  by = "Gene.ID",
  all.x = TRUE)  

res_activated_shrink_df_padj_merge <- na.omit(res_activated_shrink_df_padj_merge)




#### ---- Red GENIE3 Activated con los genes Aicda_PO expresados ---- ####

# Se seleccionan las muestras Activated Aicda_PO de la matrix trasnformada
muestras_activated_PO <- colnames(expr_matrix)[grepl("Activated.*Aicda_PO", colnames(expr_matrix))]

# Se seleccionan los DEGs deseados: solo con log2FC <  0 (genes downregulados AICDA_PO) y padj < 0.1 (filtrado anteriormente)
res_down_activated_PO <- res_activated_shrink_df_padj_merge %>%
  filter(log2FoldChange < 0)

genes_down_activated_PO <- res_down_activated_PO$Gene.Name # se listan los genes (nombre)

gene_ids_down_activated_PO <- all_genes$Gene.ID[all_genes$Gene.Name %in% genes_down_activated_PO] # se listan los genes (ID de Ensmbl)

# Se filtra seleccionan de la matriz de expresión los genes DEGs y las muestras Aicda_PO
expr_deg_down_activated_PO <- expr_matrix[rownames(expr_matrix) %in% gene_ids_down_activated_PO, muestras_activated_PO]

# MAPEAR SYMBOLS PARA LA MATRIZ
expr_annotated_activated_PO <- merge(
  data.frame(Gene.ID = rownames(expr_deg_down_activated_PO), expr_deg_down_activated_PO),
  all_genes[, c("Gene.ID", "Gene.Name")],
  by = "Gene.ID"
)

rownames(expr_annotated_activated_PO) <- expr_annotated_activated_PO$Gene.Name
expr_named_activated_PO <- expr_annotated_activated_PO[, !(colnames(expr_annotated_activated_PO) %in% c("Gene.ID", "Gene.Name"))]
expr_named_activated_PO <- as.matrix(expr_named_activated_PO)

# Sebuscan los genes reguladores (TFs) presentes en la matriz
regulators_activated_PO <- intersect(rownames(expr_named_activated_PO), tfs_symbols)

if (length(regulators_activated_PO) < 2) {
  stop("No hay suficientes TFs para correr GENIE3.")
}

# Se ejecuta GENIE3
weightMatrix_activated_PO <- GENIE3(expr_named_activated_PO, regulators = regulators_activated_PO, nCores = 2)
linkList_activated_PO <- getLinkList(weightMatrix_activated_PO)
colnames(linkList_activated_PO) <- c("regulatoryGene", "targetGene", "weight")

# Se filtra la linklist: peso > 0.01 y top 5 por TF (5 conexiones por TF)
filtered_links_activated_PO <- linkList_activated_PO %>%
  filter(weight > 0.01) %>%
  group_by(regulatoryGene) %>%
  slice_max(order_by = weight, n = 5) %>%
  ungroup()

# Se exporta el resultado a Cytoscape para su visualización
top_links_cyto_activated_PO <- filtered_links_activated_PO %>%
  rename(source = regulatoryGene, target = targetGene) %>%
  mutate(interaction = "regulates", name = paste(source, "regulates", target))

# Se guarda el resultado
write.xlsx(top_links_cyto_activated_PO, "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Redes nuevas/Activated/GENIE3_Activated_PO_DEGs_Tabla.xlsx")


createNetworkFromDataFrames(
  edges = top_links_cyto_activated_PO,
  title = "GENIE3_Activated_PO_DEGs",
  collection = "GENIE3_PO"
)

# Se anotan los nodos para aplicar después el color a cada uno de ellos
all_nodes_activated_PO <- unique(c(top_links_cyto_activated_PO$source, top_links_cyto_activated_PO$target))

node_type_df_activated_PO <- data.frame(
  id = all_nodes_activated_PO,
  type = ifelse(all_nodes_activated_PO %in% top_links_cyto_activated_PO$source, "TF", "Target")
)

loadTableData(node_type_df_activated_PO, data.key.column = "id", table.key.column = "name")


# Se modifica el estilo visual de la red y se aplica color a los nodos
style_name <- "GENIE3_Purple_Yellow"

createVisualStyle(style_name,
                  defaults = list(
                    NODE_SHAPE = "ellipse",
                    NODE_BORDER_WIDTH = 1,
                    EDGE_STROKE_UNSELECTED_PAINT = "#AAAAAA",
                    NETWORK_BACKGROUND_PAINT = "#FFFFFF"
                  )
)

setNodeColorMapping("type",
                    table.column.values = c("TF", "Target"),
                    colors = c("gold", "#55B1FD"),
                    mapping.type = "d",
                    style.name = style_name
)

setNodeLabelMapping("name", style.name = style_name)


setCurrentNetwork("GENIE3_Activated_PO_DEGs") # se setea la red antes de aplicar el estilo

setVisualStyle(style_name)





#### ---- Red GENIE3 Activated con los genes Aicda_KO expresados ---- ####

# Se seleccionan las muestras Activated Aicda_PO de la matrix trasnformada
muestras_activated_KO <- colnames(expr_matrix)[grepl("Activated.*Aicda_KO", colnames(expr_matrix))]

# Se seleccionan los DEGs deseados: solo con log2FC > 0 (genes upregulados AICDA_KO) y padj < 0.1 (filtrado anteriormente)
res_up_activated_KO <- res_activated_shrink_df_padj_merge %>%
  filter(log2FoldChange > 0)

genes_up_activated_KO <- res_up_activated_KO$Gene.Name # se listan los genes (nombre)

gene_ids_up_activated_KO <- all_genes$Gene.ID[all_genes$Gene.Name %in% genes_up_activated_KO] # se listan los genes (ID de Ensmbl)

# Se filtra seleccionan de la matriz de expresión los genes DEGs y las muestras Aicda_KO
expr_deg_up_activated_KO <- expr_matrix[rownames(expr_matrix) %in% gene_ids_up_activated_KO, muestras_activated_KO]

# MAPEAR SYMBOLS PARA LA MATRIZ
expr_annotated_activated_KO <- merge(
  data.frame(Gene.ID = rownames(expr_deg_up_activated_KO), expr_deg_up_activated_KO),
  all_genes[, c("Gene.ID", "Gene.Name")],
  by = "Gene.ID"
)

rownames(expr_annotated_activated_KO) <- expr_annotated_activated_KO$Gene.Name
expr_named_activated_KO <- expr_annotated_activated_KO[, !(colnames(expr_annotated_activated_KO) %in% c("Gene.ID", "Gene.Name"))]
expr_named_activated_KO <- as.matrix(expr_named_activated_KO)

# Se buscan los genes reguladores (TFs) presentes en la matriz
regulators_activated_KO <- intersect(rownames(expr_named_activated_KO), tfs_symbols)

if (length(regulators_activated_KO) < 2) {
  stop("No hay suficientes TFs para correr GENIE3.")
}

# Se ejecuta GENIE3
weightMatrix_activated_KO <- GENIE3(expr_named_activated_KO, regulators = regulators_activated_KO, nCores = 2)
linkList_activated_KO <- getLinkList(weightMatrix_activated_KO)
colnames(linkList_activated_KO) <- c("regulatoryGene", "targetGene", "weight")

# Se filtra la linklist: peso > 0.01 y top 5 por TF (5 conexiones por TF)
filtered_links_activated_KO <- linkList_activated_KO %>%
  filter(weight > 0.01) %>%
  group_by(regulatoryGene) %>%
  slice_max(order_by = weight, n = 5) %>%
  ungroup()

# Se exporta el resultado a Cytoscape para su visualización
top_links_cyto_activated_KO <- filtered_links_activated_KO %>%
  rename(source = regulatoryGene, target = targetGene) %>%
  mutate(interaction = "regulates", name = paste(source, "regulates", target))


# Se guarda el resultado
write.xlsx(top_links_cyto_activated_KO, "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Redes nuevas/Activated/GENIE3_Activated_KO_DEGs_Tabla.xlsx")



createNetworkFromDataFrames(
  edges = top_links_cyto_activated_KO,
  title = "GENIE3_Activated_KO_DEGs",
  collection = "GENIE3_KO"
)

# Se anotan los nodos para aplicar después el color a cada uno de ellos
all_nodes_activated_KO <- unique(c(top_links_cyto_activated_KO$source, top_links_cyto_activated_KO$target))

node_type_df_activated_KO <- data.frame(
  id = all_nodes_activated_KO,
  type = ifelse(all_nodes_activated_KO %in% top_links_cyto_activated_KO$source, "TF", "Target")
)

loadTableData(node_type_df_activated_KO, data.key.column = "id", table.key.column = "name")


# Se modifica el estilo visual de la red y se aplica color a los nodos
style_name <- "GENIE3_Purple_Yellow"

createVisualStyle(style_name,
                  defaults = list(
                    NODE_SHAPE = "ellipse",
                    NODE_BORDER_WIDTH = 1,
                    EDGE_STROKE_UNSELECTED_PAINT = "#AAAAAA",
                    NETWORK_BACKGROUND_PAINT = "#FFFFFF"
                  )
)

setNodeColorMapping("type",
                    table.column.values = c("TF", "Target"),
                    colors = c("gold", "#55B1FD"),
                    mapping.type = "d",
                    style.name = style_name
)

setNodeLabelMapping("name", style.name = style_name)

setCurrentNetwork("GENIE3_Activated_KO_DEGs")

setVisualStyle(style_name)


