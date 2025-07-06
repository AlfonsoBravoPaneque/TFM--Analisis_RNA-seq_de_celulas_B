#### ---- Paquetes necesarios ---- ####

# Paquetes
library(DESeq2)
library(ggplot2)
library(extrafont)

loadfonts(device = "win") # carga las fuentes


#### ---- PCA ---- ####

# Para observar como las diferentes muestras se agrupan, se realizó un PCA.

pcaData <- plotPCA(vsd, intgroup=c("cell_type", "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

colores_cell_type <- c(
  "Naive" = "cyan",
  "Activated" = "#FFD700",
  "GC" = "#32CD32"
)


ggplot(pcaData, aes(PC1, PC2, color = cell_type, shape = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% varianza")) +
  ylab(paste0("PC2: ",percentVar[2],"% varianza")) +
  labs(color = "Tipo Celular", shape = "Condición Experimental") +
  scale_color_manual(values = colores_cell_type) +
  theme_classic() +
  theme(axis.title.x = element_text(family = "Calibri", size = 12, margin = margin(5,0,0,0)),
        axis.title.y = element_text(family = "Calibri", size = 12, margin = margin(0,5,0,0)),
        axis.text = element_text(family = "Calibri", size = 10),
        legend.title = element_text(family = "Calibri", size = 11, face = "bold"),
        legend.text = element_text(family = "Calibri", size = 9))
coord_fixed()