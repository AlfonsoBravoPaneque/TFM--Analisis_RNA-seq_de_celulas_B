#### ---- Paquetes necesarios ---- ####

# Paquetes
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(openxlsx)
library(extrafont)
library(enrichplot)

loadfonts(device = "win") # carga las fuentes


#### ---- Enriquecimiento funcional - GSEA (Gene Set Enrichment Analysis) ---- ####

# Para llevar a cabo el GSEA primero se ordenan los genes de mayor a menor Log2FoldChange (se hace el ranking en función del foldchange).
# Para el enriquecimiento se emplea la ontología Biological Processes de Gen ontology (GO). En este caso, no hace falta eliminar los valores NA
# presentes en las columnas de p value y padj, ya que no afecta al GSEA (afectarías si el NA estuviera en la columan de log2FoldChange).

# Una gran diferencia con respecto al ORA es que el GSEA usa la función gseaGO(), la cual emplea todos los genes (incluso los no mapeados o identificados). En estos
# analisis se selecciona todo el pull de genes y se introducen en la función gseaGO(), ya que no afecta a los terminos GO obtenidos sino más bien afectan ligeramente al p value
# y al NES score puediéndose mantener los genes no mapeados que actúan como un background general.

set.seed(1234) # se filla semilla para favorecer la reproducibilidad


#### ---- Células GC -> Enriquecimiento - GSEA ---- ####

# Se calcula el número de NA por columna
na_counts <- colSums(is.na(res_GC_shrink_df))

# Se muestrar sólo las columnas que tienen al menos un NA. Hay 5 NA en la columna pvalue y 3967 en padj (p valor ajustado):
na_counts[na_counts > 0]

# Se seleccionan los genes para realizar el enriquecimiento (todos los genes) y se ordena el data.frame de mayor a menor log2FoldChange. 
all_genes_GC <- arrange(res_GC_shrink_df, -log2FoldChange)
genes_list_all_genes_GC <- all_genes_GC$Gene.ID


# Se comprueba que los genes presentes en genes_list_all_genes_GC sean reconozidos en su mayoría por org.Mm.eg.db. Para ello, primero
# se obtienen todos los IDs de ensembl disponibles en org.Mm.eg.db, y después se comprueba que la mayoría están presentes en Ensmbl
all_ensembl_ids <- keys(org.Mm.eg.db, keytype = "ENSEMBL")

# 13534 genes encontrados de 18577 (sin embargo, aunque algunos de estos genes no sean mapeados, para el calculo estadístico
# la función gseGO tiene en cuenta a todos los genes del análisis, por lo que no hace falta seleccionar solo los 18577 ya que un 70% de los genes
# mapeados se considera un buen valor. En cambio, al hacer el enrichGO la propia función descarta los genes no mapeados. Es verdad que modifica un poco el NES score
# y el p valor pero de forma mínima, y se pueden mantener los genes no mapeados auctuando como un background general)
valid_ids_all_genes_GC <- genes_list_all_genes_GC[genes_list_all_genes_GC %in% all_ensembl_ids]
cat("Genes encontrados en org.Mm.eg.db:", length(valid_ids_all_genes_GC), "de", length(genes_list_all_genes_GC))

# Se prepara el cojunto total de genes, generando un vector numérico ordenado (de mayor a menor log2FoldChange) que posea tanto los valores
# de log2FoldChange como los ID de Ensmbl
geneList_all_genes_GC <- all_genes_GC$log2FoldChange
names(geneList_all_genes_GC) <- all_genes_GC$Gene.ID

# Se lleva a cabo el enriquecimiento GSEA
ego_gsea_GC <- gseGO(geneList = geneList_all_genes_GC,
                        OrgDb = org.Mm.eg.db,
                        ont = "BP",
                        keyType = "ENSEMBL",
                        minGSSize = 50,
                        maxGSSize = 1000,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        nPermSimple = 10000, # permite afinar el p valor y que sea mas robusto (permutaciones aleatorias)
                        verbose = TRUE)


# Se exporta el resultado en un data.frame
go_results_gsea_GC <- as.data.frame(ego_gsea_GC)
View(go_results_gsea_GC)

# Se separan los ID de los genes del data.frame con los resultados (1 ID por fila)
tabla_expanded_GC_GSEA <- tidyr::separate_rows(go_results_gsea_GC, core_enrichment, sep = "/")

# Se convierte la lista a data.frame para poder hacer el merge. Además, se hace un merge, permitiendo asociar el nombre de los genes con su ID
gene_list_GC_df <- data.frame(Gene.ID = names(geneList_all_genes_GC))


gene_list_GC_df_all <- merge(gene_list_GC_df, all_genes[, c("Gene.ID", "Gene.Name")], 
                                by = "Gene.ID", all.x = TRUE)

colnames(gene_list_GC_df_all)[1] <- "core_enrichment"

# Se hace un merge de "tabla_expanded_GC_GSEA" con "gene_list_GC_df_all" (presenta también los ID de ensembl)
tabla_merged_GC <- merge(tabla_expanded_GC_GSEA,
                            gene_list_GC_df_all,
                            by = "core_enrichment",
                            all.x = TRUE)


# Se agrupan los datos para que cada término GO tenga una lista de genes asociada
tabla_final_GC <- tabla_merged_GC %>%
  group_by(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, rank, leading_edge) %>%
  summarise(core_enrichment = paste(core_enrichment, collapse = "/"), 
            Gene.Name = paste(Gene.Name, collapse = "/"), .groups = "drop")


# Se ordena la tabla de mayor a menor p valor ajustado
tabla_final_GC_sorted <- tabla_final_GC %>%
  arrange(p.adjust)


# Se exporta el data.frame tabla_final_GC_sorted a formato excel
write.xlsx(tabla_final_GC_sorted, "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Enriquecimiento funcional GSEA/GC")


# Se plotea el resultado realizando un GSEA plot (elegir los que tengan sentido biológico, y que tenga relación con los datos)

# Extraer el término
term_row <- ego_gsea_GC@result[ego_gsea_GC@result$ID == "GO:0045621", ]
term_row_2 <- ego_gsea_GC@result[ego_gsea_GC@result$ID == "GO:0002683", ]
term_row_3 <- ego_gsea_GC@result[ego_gsea_GC@result$ID == "GO:0044282", ]


# Se convierte el p.adjust a notación científica base 10 como texto
p_value_text <- formatC(term_row$p.adjust, format = "e", digits = 2)
p_value_text_2 <- formatC(term_row_2$p.adjust, format = "e", digits = 2)
p_value_text_3 <- formatC(term_row_3$p.adjust, format = "e", digits = 2)


# Crear el título como texto plano con esa notación
gsea_title <- paste0(
  "Enriquecimiento GSEA: Regulación positiva de la diferenciación linfocitaria",
  "\n",
  "NES = ", round(term_row$NES, 2),
  ", p.adjust = ", p_value_text
)


gsea_title_2 <- paste0(
  "Enriquecimiento GSEA: Regulación negativa de procesos del sistema inmunitario",
  "\n",
  "NES = ", round(term_row_2$NES, 2),
  ", p.adjust = ", p_value_text_2
)

gsea_title_3 <- paste0(
  "Enriquecimiento GSEA: Proceso catabólico de pequeñas moléculas",
  "\n",
  "NES = ", round(term_row_3$NES, 2),
  ", p.adjust = ", p_value_text_3
)



# Se grafica

p1 <- gseaplot2(ego_gsea_GC, geneSetID = "GO:0045621", title = gsea_title) # REGULACIÓN POSITIA DE LA DIERENCIACIÓN LINFOCITARIA
p1

p2 <- gseaplot2(ego_gsea_GC, geneSetID = "GO:0002683", title = gsea_title_2) # REGULACIÓN NEGATIVA DE LOS PROCESOS DEL SISTEMA INMUNE
p2

p3 <- gseaplot2(ego_gsea_GC, geneSetID = "GO:0044282", title = gsea_title_3) # PROCESO CATABÓLICO DE PEQUEÑAS MOLÉCULAS
p3


#### ---- Células Activated -> Enriquecimiento - GSEA ---- ####

# Se calcula el número de NA por columna
na_counts <- colSums(is.na(res_activated_shrink_df))

# Se muestrar sólo las columnas que tienen al menos un NA. Hay 5 NA en la columna pvalue y otros 5 en padj (p valor ajustado)
na_counts[na_counts > 0]

# Se seleccionan los genes para realizar el enriquecimiento (todos los genes) y se ordena el data.frame de mayor a menor log2FoldChange. 
all_genes_activated <- arrange(res_activated_shrink_df, -log2FoldChange)
genes_list_all_genes_activated <- all_genes_activated$Gene.ID

# Se comprueba que los genes presentes en genes_list_all_genes_activated sean reconozidos en su mayoría por org.Mm.eg.db. Para ello, primero
# se obtienen todos los IDs de ensembl disponibles en org.Mm.eg.db, y después se comprueba que la mayoría están presentes en Ensmbl
all_ensembl_ids <- keys(org.Mm.eg.db, keytype = "ENSEMBL")

# 13534 genes encontrados de 18577 (sin embargo, aunque algunos de estos genes no sean mapeados, para el calculo estadístico
# la función gseGO tiene en cuenta a todos los genes del análisis, por lo que no hace falta seleccionar solo los 13534 ya que un 70% de los genes
# mapeados se considera un buen valor. En cambio, al hacer el enrichGO la propia función descarta los genes no mapeados. Es verdad que modifica un poco el NES score
# y el p valor pero de forma mínima, y se pueden mantener los genes no mapeados auctuando como un background general)
valid_ids_all_genes_activated <- genes_list_all_genes_activated[genes_list_all_genes_activated %in% all_ensembl_ids]
cat("Genes encontrados en org.Mm.eg.db:", length(valid_ids_all_genes_activated), "de", length(genes_list_all_genes_activated))

# Se prepara el cojunto total de genes, generando un vector numérico ordenado (de mayor a menor log2FoldChange) que posea tanto los valores
# de log2FoldChange como los ID de Ensmbl
geneList_all_genes_activated <- all_genes_activated$log2FoldChange
names(geneList_all_genes_activated) <- all_genes_activated$Gene.ID

# Se lleva a cabo el enriquecimiento GSEA
ego_gsea_activated <- gseGO(geneList = geneList_all_genes_activated,
                        OrgDb = org.Mm.eg.db,
                        ont = "BP",
                        keyType = "ENSEMBL",
                        minGSSize = 10, # número mínimos de genes asociados al término GO (se disminutye el valor ya que salen pocos resultados con 50)
                        maxGSSize = 8000, # número máximo de genes asociados al término GO (se aumenta el valor ya que salen pocos resultados con 1000)
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        nPermSimple = 10000, # permite afinar el p valor y que sea mas robusto (permutaciones aleatorias)
                        verbose = TRUE)

# Se exporta el resultado en un data.frame
go_results_gsea_activated <- as.data.frame(ego_gsea_activated)
View(go_results_gsea_activated)


# Se separan los ID de los genes del data.frame con los resultados (1 ID por fila)
tabla_expanded_activated_GSEA <- tidyr::separate_rows(go_results_gsea_activated, core_enrichment, sep = "/")

# Se convierte la lista a data.frame para poder hacer el merge. Además, se hace un merge, permitiendo asociar el nombre de los genes con su ID
gene_list_activated_df <- data.frame(Gene.ID = names(geneList_all_genes_activated))

gene_list_activated_df_all <- merge(gene_list_activated_df, all_genes[, c("Gene.ID", "Gene.Name")], 
                                by = "Gene.ID", all.x = TRUE)

colnames(gene_list_activated_df_all)[1] <- "core_enrichment"


# Se hace un merge de "tabla_expanded_activated_GSEA" con "gene_list_activated_df_all" (presenta también los ID de ensembl)
tabla_merged_activated <- merge(tabla_expanded_activated_GSEA,
                                gene_list_activated_df_all,
                                by = "core_enrichment",
                                all.x = TRUE)

# Se agrupan los datos para que cada término GO tenga una lista de genes asociada
tabla_final_activated <- tabla_merged_activated %>%
  group_by(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, rank, leading_edge) %>%
  summarise(core_enrichment = paste(core_enrichment, collapse = "/"), 
            Gene.Name = paste(Gene.Name, collapse = "/"), .groups = "drop")


# Se ordena la tabla de mayor a menor p valor ajustado
tabla_final_activated_sorted <- tabla_final_activated %>%
  arrange(p.adjust)


# Se exporta el data.frame tabla_final_activated_sorted a formato excel
write.xlsx(tabla_final_activated_sorted, "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Enriquecimiento funcional GSEA/Activated")


# Se plotea el resultado realizando un GSEA plot (elegir los que tengan sentido biológico, y que tenga relación con los datos)

# Extraer el término
term_row <- ego_gsea_activated@result[ego_gsea_activated@result$ID == "GO:1902275", ]
term_row_2 <- ego_gsea_activated@result[ego_gsea_activated@result$ID == "GO:0080111", ]
term_row_3 <- ego_gsea_activated@result[ego_gsea_activated@result$ID == "GO:0040029", ]


# Se convierte el p.adjust a notación científica base 10 como texto
p_value_text <- formatC(term_row$p.adjust, format = "e", digits = 2)
p_value_text_2 <- formatC(term_row_2$p.adjust, format = "e", digits = 2)
p_value_text_3 <- formatC(term_row_3$p.adjust, format = "e", digits = 2)


# Crear el título como texto plano con esa notación
gsea_title <- paste0(
  "Enriquecimiento GSEA: Regulación de la organización de la cromatina",
  "\n",
  "NES = ", round(term_row$NES, 2),
  ", p.adjust = ", p_value_text
)


gsea_title_2 <- paste0(
  "Enriquecimiento GSEA: Desmetilación del ADN",
  "\n",
  "NES = ", round(term_row_2$NES, 2),
  ", p.adjust = ", p_value_text_2
)

gsea_title_3 <- paste0(
  "Enriquecimiento GSEA: Regulación epigenética de la expresión génica",
  "\n",
  "NES = ", round(term_row_3$NES, 2),
  ", p.adjust = ", p_value_text_3
)



# Se grafica

p1 <- gseaplot2(ego_gsea_activated, geneSetID = "GO:1902275", title = gsea_title) # REGULACIÓN DE LA ORGANIZACIÓN DE LA CROMATINA
p1

p2 <- gseaplot2(ego_gsea_activated, geneSetID = "GO:0080111", title = gsea_title_2) # DESMETILACIÓN DEL ADN
p2

p3 <- gseaplot2(ego_gsea_activated, geneSetID = "GO:0040029", title = gsea_title_3) # REGULACIÓN EPIGENÉTICA DE LA EXPRESIÓN GÉNICA
p3

