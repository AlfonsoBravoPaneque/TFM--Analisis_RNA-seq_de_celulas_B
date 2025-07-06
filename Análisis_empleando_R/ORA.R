#### ---- Paquetes necesarios ---- ####

# Paquetes
library(dplyr)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyr)
library(openxlsx)
library(extrafont)


loadfonts(device = "win") # carga las fuentes


### ---- Enriquecimiento funcional - ORA (Over Representation Analysis) ---- ####

# Para llevar a cabo el ORA se seleccionan los genes diferencialmente expresados en ambas condiciones y tipos celulares. No se emplea un threshold log2foldchange
# y si se emplea un threshold de p valor ajustado (no de pa valor ajustado). Para realizar el enriquecimiento funcional se emplea la ontología de Biological
# Processes de Gene Ontology (GO).

# Una gran diferencia entre el ORA y el GSEA es la función que usan. En este caso, para el ORA se usa la función enrichGO(), que analiza únicamente los genes que se mapean a términos GO en org.Mm.eg.db
# (por eso es importante comprobar que la mayoría de los genes en cada ORA realizado están presentes en org.Mm.eg.db)


#### ---- Células GC -> Enriquecimiento - ORA ---- ####

# Se eligen los genes para realizar el enriquecimiento. Se eligen los genes sobre expresados y menos expreados (eliminando con anterioridad los genes ribosomales)
# con un threshold del p value ajustado menor a 0.1
upregulated_GC_pvalTH <- subset(markers_GC_filt, padj < 0.1 & State == "up")
gene_list_upregulated_GC <- upregulated_GC_pvalTH$Gene.ID

downregulated_GC_pvalTH <- subset(markers_GC_filt, padj < 0.1 & State == "down")
gene_list_downregulated_GC <- downregulated_GC_pvalTH$Gene.ID

# Se comprueba que los genes presentes en gene_list_upregulated_GC y gene_list_downregulated_GC sean reconozidos en su mayoría por org.Mm.eg.db. Para ello, primero
# se obtienen todos los IDs de ensembl disponibles en org.Mm.eg.db, y después se determina el número de genes disponibles
all_ensembl_ids <- keys(org.Mm.eg.db, keytype = "ENSEMBL")

# 616 genes encontrados de 651
valid_ids_upregulated_GC <- gene_list_upregulated_GC[gene_list_upregulated_GC %in% all_ensembl_ids]
cat("Genes encontrados en org.Mm.eg.db:", length(valid_ids_upregulated_GC), "de", length(gene_list_upregulated_GC))

# 265 genes encontrados de 290
valid_ids_downregulated_GC <- gene_list_downregulated_GC[gene_list_downregulated_GC %in% all_ensembl_ids]
cat("Genes encontrados en org.Mm.eg.db:", length(valid_ids_downregulated_GC), "de", length(gene_list_downregulated_GC))


#### ---- Enriquecimiento (ORA) de genes sobreexpresados (Aicda_KO) en células GC ---- ####

# Se realiza un enriquecimiento GO para los genes más expresados con un threshold de p ajustado menor a 0.1


# 999 resultados (se puede usar la función dropGO() para filtrar los terminos más generales del enriquecimiento)
ego_GC_Aicda_KO <- enrichGO(gene = gene_list_upregulated_GC,
                            OrgDb = org.Mm.eg.db,
                            keyType = "ENSEMBL",
                            ont = "BP",  # Se puede usar "BP" (procesos biológicos), "MF" (funciones moleculares), o "CC" (componentes celulares)
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)


# Se extrae el data.frame con los resultados del enriquecimiento
go_results_GC_Aicda_KO <- as.data.frame(ego_GC_Aicda_KO)

# Se crea una nueva columna log_p.adjust con los valores log transformados de p.adjust
go_results_GC_Aicda_KO$log_p.adjust <- -log10(go_results_GC_Aicda_KO$p.adjust)


# Seleccionar los términos que se desean mostrar
lista = c("GO:0030098", "GO:0002440", "GO:0042113", "GO:0046631", "GO:0002764", "GO:0051250", "GO:0002698", "GO:0002923", "GO:0002639", "GO:0070228")

go_results_GC_Aicda_KO_filt_term <- dplyr::filter(go_results_GC_Aicda_KO, ID %in% lista)

ego_GC_Aicda_KO@result <-  go_results_GC_Aicda_KO_filt_term


# Se grafica el dot plot. Se eligen los terminos GO con mejor p.adjust, ordenando estos 10 terminos en el plot en función de
# su Gene ratio
p1 <- dotplot(ego_GC_Aicda_KO, showCategory = 10) +
  ggtitle("Enriquecimiento ORA - Aicda_KO") +
  aes(fill = log_p.adjust) +  # Mapear el relleno usando la nueva columna log_p.adjust
  scale_fill_gradient(low = "#FF69B4", high = "#00CED1", name = "-log10(p.adjust)") +  # Colores similares a los originales
  guides(color = "none") +  # Eliminar la leyenda de color original
  theme(
    text = element_text(family = "Calibri", size = 12),
    axis.title.y = element_text(family = "Calibri", size = 12),
    axis.text.x = element_text(family = "Calibri", size = 12),
    legend.title = element_text(family = "Calibri", size = 12),
    legend.text = element_text(family = "Calibri", size = 10),
    plot.title = element_text(family = "Calibri", size = 14, hjust = 0.5)
  )

# Resulado del dot plot (700 x 700)
p1


# Se separan los geneID del data.frame con los resultados (1 geneID por fila)
tabla_expanded_GC_KO <- tidyr::separate_rows(go_results_GC_Aicda_KO, geneID, sep = "/")

# Se convierte la lista a data.frame para poder hacer el merge. Además, se hace un merge, permitiendo asociar el nombre de los genes con su ID
gene_list_upregulated_GC_df <- data.frame(Gene.ID = gene_list_upregulated_GC)

gene_list_upregulated_GC_df_all <- merge(gene_list_upregulated_GC_df, all_genes[, c("Gene.ID", "Gene.Name")], 
                                         by = "Gene.ID", all.x = TRUE)


# Se hace un merge de "tabla_expanded_GC_KO" con "gene_list_upregulated_GC_df_all" (presenta también los ID de ensembl)
tabla_merged_GC_KO <- merge(tabla_expanded_GC_KO,
                            gene_list_upregulated_GC_df_all,
                            by.x = "geneID",
                            by.y = "Gene.ID",
                            all.x = TRUE)

# Se agrupan los datos para que cada término GO tenga una lista de genes asociada
tabla_final_GC_KO <- tabla_merged_GC_KO %>%
  group_by(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count) %>%
  summarise(geneID = paste(geneID, collapse = "/"), 
            geneName = paste(Gene.Name, collapse = "/"), .groups = "drop")


# Se ordena la tabla de mayor a menor p valor ajustado
tabla_final_GC_KO_sorted <- tabla_final_GC_KO %>%
  arrange(p.adjust)

# Se exporta el resultado en un excel
write.xlsx(tabla_final_GC_KO_sorted, file = "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Enriquecimiento funcional con 0.1 en ORA/ORA/GC/Knock out", rowNames = FALSE)



#### ---- Enriquecimiento (ORA) de genes menos expresados (Aicda_PO) en células GC ---- ####

# Se realiza un enriquecimiento GO para los genes menos expresados con un threshold de p value ajustado menor a 0.1

# 159 resultados (se puede usar la función dropGO() para filtrar los terminos más generales del enriquecimiento)
ego_GC_Aicda_PO <- enrichGO(gene = gene_list_downregulated_GC,
                            OrgDb = org.Mm.eg.db,
                            keyType = "ENSEMBL",
                            ont = "BP",  # Se puede usar (procesos biológicos), "MF" (funciones moleculares), o "CC" (componentes celulares)
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)


# Se extrae el data.frame con los resultados del enriquecimiento
go_results_GC_Aicda_PO <- as.data.frame(ego_GC_Aicda_PO)

# Se crea una nueva columna log_p.adjust con los valores log transformados de p.adjust
go_results_GC_Aicda_PO$log_p.adjust <- -log10(go_results_GC_Aicda_PO$p.adjust)


# Seleccionar los términos que se desean mostrar
lista = c("GO:0045787", "GO:0051304","GO:0007059", "GO:0000280", "GO:1905820", "GO:0045786", "GO:0044839", "GO:0009060", "GO:0043467", "GO:0045333")

go_results_GC_Aicda_PO_filt_term <- dplyr::filter(go_results_GC_Aicda_PO, ID %in% lista)

ego_GC_Aicda_PO@result <-  go_results_GC_Aicda_PO_filt_term


# Se grafica el dot plot. Se eligen los terminos GO con mejor p.adjust, ordenando estos 10 terminos en el plot en función de
# su Gene ratio
p1 <- dotplot(ego_GC_Aicda_PO, showCategory = 10) +
  ggtitle("Enriquecimiento ORA - Aicda_PO") +
  aes(fill = log_p.adjust) +  # Mapear el relleno usando la nueva columna log_p.adjust
  scale_fill_gradient(low = "#FF69B4", high = "#00CED1", name = "-log10(p.adjust)") +  # Colores similares a los originales
  guides(color = "none") +  # Eliminar la leyenda de color original
  theme(
    text = element_text(family = "Calibri", size = 12),
    axis.title.y = element_text(family = "Calibri", size = 12),
    axis.text.x = element_text(family = "Calibri", size = 12),
    legend.title = element_text(family = "Calibri", size = 12),
    legend.text = element_text(family = "Calibri", size = 10),
    plot.title = element_text(family = "Calibri", size = 14, hjust = 0.5)
  )

# Resulado del dot plot (700 x 750)
p1


# Se separan los geneID del data.frame con los resultados (1 geneID por fila)
tabla_expanded_GC_PO <- tidyr::separate_rows(go_results_GC_Aicda_PO, geneID, sep = "/")

# Se convierte la lista a data.frame para poder hacer el merge. Además, se hace un merge, permitiendo asociar el nombre de los genes con su ID
gene_list_downregulated_GC_df <- data.frame(Gene.ID = gene_list_downregulated_GC)

gene_list_downregulated_GC_df_all <- merge(gene_list_downregulated_GC_df, all_genes[, c("Gene.ID", "Gene.Name")], 
                                           by = "Gene.ID", all.x = TRUE)

# Se hace un merge de "tabla_expanded_GC_PO" con "gene_list_downregulated_GC_df_all" (presenta también los ID de ensembl)
tabla_merged_GC_PO <- merge(tabla_expanded_GC_PO,
                            gene_list_downregulated_GC_df_all,
                            by.x = "geneID",
                            by.y = "Gene.ID",
                            all.x = TRUE)

# Se agrupan los datos para que cada término GO tenga una lista de genes asociada
tabla_final_GC_PO <- tabla_merged_GC_PO %>%
  group_by(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count) %>%
  summarise(geneID = paste(geneID, collapse = "/"), 
            geneName = paste(Gene.Name, collapse = "/"), .groups = "drop")

# Se ordena la tabla de mayor a menor p valor ajustado
tabla_final_GC_PO_sorted <- tabla_final_GC_PO %>%
  arrange(p.adjust)

# Se exporta el resultado en un excel
write.xlsx(tabla_final_GC_PO_sorted, file = "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Enriquecimiento funcional con 0.1 en ORA/ORA/GC/Positive", rowNames = FALSE)



#### ---- Células Activated -> Enriquecimiento - ORA ---- ####

# Se eligen los genes para realizar el enriquecimiento. Se eligen los genes sobre expresados y menos expreados (eliminando con anterioridad los genes ribosomales)
# con un threshold del p value ajustado menor a 0.1
upregulated_activated_pvalTH <- subset(markers_activated_filt, padj < 0.1 & State == "up")
gene_list_upregulated_activated <- upregulated_activated_pvalTH$Gene.ID

downregulated_activated_pvalTH <- subset(markers_activated_filt, padj < 0.1 & State == "down")
gene_list_downregulated_activated <- downregulated_activated_pvalTH$Gene.ID

# Se comprueba que los genes presentes en gene_list_upregulated_activated y gene_list_downregulated_activated sean reconozidos en su mayoría por org.Mm.eg.db. Para ello, primero
# se obtienen todos los IDs de ensembl disponibles en org.Mm.eg.db, y después se determina el número de genes disponibles
all_ensembl_ids <- keys(org.Mm.eg.db, keytype = "ENSEMBL")

# 33 genes encontrados de 34
valid_ids_upregulated_activated <- gene_list_upregulated_activated[gene_list_upregulated_activated %in% all_ensembl_ids]
cat("Genes encontrados en org.Mm.eg.db:", length(valid_ids_upregulated_activated), "de", length(gene_list_upregulated_activated))

# 50 genes encontrados de 51
valid_ids_downregulated_activated <- gene_list_downregulated_activated[gene_list_downregulated_activated %in% all_ensembl_ids]
cat("Genes encontrados en org.Mm.eg.db:", length(valid_ids_downregulated_activated), "de", length(gene_list_downregulated_activated))


#### ---- Enriquecimiento (ORA) de genes sobreexpresados (Aicda_KO) en células Activated ---- ####

# Se realiza un enriquecimiento GO para los genes más expresados con un threshold de p value ajustado menor a 0.1

# 2 resultados (no es necesario usar la función dropGO() ya que hay pocos resultados)
ego_activated_Aicda_KO <- enrichGO(gene = gene_list_upregulated_activated,
                                   OrgDb = org.Mm.eg.db,
                                   keyType = "ENSEMBL",
                                   ont = "BP",  # Se puede usar "BP" (procesos biológicos), "MF" (funciones moleculares), o "CC" (componentes celulares)
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.05)


# Se extrae el data.frame con los resultados del enriquecimiento
go_results_activated_Aicda_KO <- as.data.frame(ego_activated_Aicda_KO)

# Se crea una nueva columna log_p.adjust con los valores log transformados de p.adjust
go_results_activated_Aicda_KO$log_p.adjust <- -log10(go_results_activated_Aicda_KO$p.adjust)


# Seleccionar los términos que se desean mostrar
lista = c("GO:0016601", "GO:0050847")

go_results_activated_Aicda_KO_filt_term <- dplyr::filter(go_results_activated_Aicda_KO, ID %in% lista)

ego_activated_Aicda_KO@result <-  go_results_activated_Aicda_KO_filt_term


# Se grafica el dot plot. Se eligen los terminos GO con mejor p.adjust, ordenando estos 10 terminos en el plot en función de
# su Gene ratio
p1 <- dotplot(ego_activated_Aicda_KO, showCategory = 2) +
  ggtitle("Enriquecimiento ORA - Aicda_KO") +
  aes(fill = log_p.adjust) +  # Mapear el relleno usando la nueva columna log_p.adjust
  scale_fill_gradient(low = "#FF69B4", high = "#00CED1", name = "-log10(p.adjust)") +  # Colores similares a los originales
  guides(color = "none") +  # Eliminar la leyenda de color original
  theme(
    text = element_text(family = "Calibri", size = 12),
    axis.title.y = element_text(family = "Calibri", size = 12),
    axis.text.x = element_text(family = "Calibri", size = 12),
    legend.title = element_text(family = "Calibri", size = 12),
    legend.text = element_text(family = "Calibri", size = 10),
    plot.title = element_text(family = "Calibri", size = 14, hjust = 0.5)
  )

# Resulado del dot plot (700 x 600)
p1


# Se separan los geneID del data.frame con los resultados (1 geneID por fila)
tabla_expanded_activated_KO <- tidyr::separate_rows(go_results_activated_Aicda_KO, geneID, sep = "/")

# Se convierte la lista a data.frame para poder hacer el merge. Además, se hace un merge, permitiendo asociar el nombre de los genes con su ID
gene_list_upregulated_activated_df <- data.frame(Gene.ID = gene_list_upregulated_activated)

gene_list_upregulated_activated_df_all <- merge(gene_list_upregulated_activated_df, all_genes[, c("Gene.ID", "Gene.Name")], 
                                                by = "Gene.ID", all.x = TRUE)

# Se hace un merge de "tabla_expanded_activated_KO" con "gene_list_upregulated_activated_df_all" (presenta también los ID de ensembl)
tabla_merged_activated_KO <- merge(tabla_expanded_activated_KO,
                                   gene_list_upregulated_activated_df_all,
                                   by.x = "geneID",
                                   by.y = "Gene.ID",
                                   all.x = TRUE)

# Se agrupan los datos para que cada término GO tenga una lista de genes asociada
tabla_final_activated_KO <- tabla_merged_activated_KO %>%
  group_by(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count) %>%
  summarise(geneID = paste(geneID, collapse = "/"), 
            geneName = paste(Gene.Name, collapse = "/"), .groups = "drop")

# Se ordena la tabla de mayor a menor p valor ajustado
tabla_final_activated_KO_sorted <- tabla_final_activated_KO %>%
  arrange(p.adjust)

# Se exporta el resultado en un excel
write.xlsx(tabla_final_activated_KO_sorted, file = "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Enriquecimiento funcional con 0.1 en ORA/ORA/Activated/Knock out", rowNames = FALSE)




#### ---- Enriquecimiento (ORA) de genes menos expresados (Aicda_PO) en células Activated ---- ####

# Se realiza un enriquecimiento GO para los genes más expresados con un threshold de p value ajustado de 0.1

# 86 resultados (se puede usar la función dropGO() para filtrar los terminos más generales del enriquecimiento)
ego_activated_Aicda_PO <- enrichGO(gene = gene_list_downregulated_activated,
                                   OrgDb = org.Mm.eg.db,
                                   keyType = "ENSEMBL",
                                   ont = "BP",  # Se puede usar (procesos biológicos), "MF" (funciones moleculares), o "CC" (componentes celulares)
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.05)


# Se extrae el data.frame con los resultados del enriquecimiento
go_results_activated_Aicda_PO <- as.data.frame(ego_activated_Aicda_PO)

# Se crea una nueva columna log_p.adjust con los valores log transformados de p.adjust
go_results_activated_Aicda_PO$log_p.adjust <- -log10(go_results_activated_Aicda_PO$p.adjust)


# Seleccionar los términos que se desean mostrar
lista = c("GO:0002455", "GO:0050853", "GO:0019724", "GO:0042113", "GO:0016064", "GO:0002429", "GO:0002697", "GO:0030888", "GO:1902105", "GO:0051251")

go_results_activated_Aicda_PO_filt_term <- dplyr::filter(go_results_activated_Aicda_PO, ID %in% lista)

ego_activated_Aicda_PO@result <-  go_results_activated_Aicda_PO_filt_term


# Se grafica el dot plot. Se eligen los terminos GO con mejor p.adjust, ordenando estos 10 terminos en el plot en función de
# su Gene ratio
p1 <- dotplot(ego_activated_Aicda_PO, showCategory = 10) +
  ggtitle("Enriquecimiento ORA - Aicda_PO") +
  aes(fill = log_p.adjust) +  # Mapear el relleno usando la nueva columna log_p.adjust
  scale_fill_gradient(low = "#FF69B4", high = "#00CED1", name = "-log10(p.adjust)") +  # Colores similares a los originales
  guides(color = "none") +  # Eliminar la leyenda de color original
  theme(
    text = element_text(family = "Calibri", size = 12),
    axis.title.y = element_text(family = "Calibri", size = 12),
    axis.text.x = element_text(family = "Calibri", size = 12),
    legend.title = element_text(family = "Calibri", size = 12),
    legend.text = element_text(family = "Calibri", size = 10),
    plot.title = element_text(family = "Calibri", size = 14, hjust = 0.5)
  )

# Resulado del dot plot (700 x 750)
p1


# Se separan los geneID del data.frame con los resultados (1 geneID por fila)
tabla_expanded_activated_PO <- tidyr::separate_rows(go_results_activated_Aicda_PO, geneID, sep = "/")

# Se convierte la lista a data.frame para poder hacer el merge. Además, se hace un merge, permitiendo asociar el nombre de los genes con su ID
gene_list_downregulated_activated_df <- data.frame(Gene.ID = gene_list_downregulated_activated)

gene_list_downregulated_activated_df_all <- merge(gene_list_downregulated_activated_df, all_genes[, c("Gene.ID", "Gene.Name")], 
                                                  by = "Gene.ID", all.x = TRUE)

# Se hace un merge de "tabla_expanded_activated_PO" con "gene_list_downregulated_activated_df_all" (presenta también los ID de ensembl)
tabla_merged_activated_PO <- merge(tabla_expanded_activated_PO, 
                                   gene_list_downregulated_activated_df_all, 
                                   by.x = "geneID", 
                                   by.y = "Gene.ID", 
                                   all.x = TRUE)

# Se agrupan los datos para que cada término GO tenga una lista de genes asociada
tabla_final_activated_PO <- tabla_merged_activated_PO %>%
  group_by(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count) %>%
  summarise(geneID = paste(geneID, collapse = "/"), 
            geneName = paste(Gene.Name, collapse = "/"), .groups = "drop")

# Se ordena la tabla de mayor a menor p valor ajustado
tabla_final_activated_PO_sorted <- tabla_final_activated_PO %>%
  arrange(p.adjust)

# Se exporta el resultado en un excel
write.xlsx(tabla_final_activated_PO_sorted, file = "C:/Users/Alfonso/Desktop/Eliminando regiones contaminantes/Enriquecimiento funcional con 0.1 en ORA/ORA/Activated/Positive", rowNames = FALSE)

