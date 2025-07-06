| Autor: | Alfonso Bravo Paneque                                                                          |
:--------------------------------------------------------------------------------------------------------|:-

# **Repositorio TFM: Análisis de datos de expresión de células B**

### **1. Propósito**

<div align="justify">
Este repositorio fue creado para mostrar el código que desarrollé durante la elaboración de mi Trabajo de Fin de Máster (TFM), con el objetivo principal de garantizar la reproducibilidad del procesamiento y análisis de los datos. Dichos datos provienen del experimento realizado por Hogenbirk y colaboradores (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47703" target="_blank">GSE47703</a>). En este trabajo, llevé a cabo análisis de expresión diferencial, análisis de enriquecimiento funcional (ORA y GSEA) e inferencia de redes de regulación génica en subpoblaciones de células B, bajo condiciones de tipo silvestre (WT) y con el gen <i>Aicda</i> inactivado (<i>knockout</i>).
</div>

<br> 

### **2. Estructura de las carpetas**

📁 **Procesado_datos_crudos_empleado_WSL**

<div align="justify">

La carpeta <a href="https://github.com/AlfonsoBravoPaneque/TFM--Analisis_RNA-seq_de_celulas_B/tree/main/An%C3%A1lisis_empleando_R" target="_blank">Procesado_datos_crudos_empleado_WSL</a> contiene el <i>script</i> empleado para procesar los datos de secuenciación (archivos FASTQ), hasta obtener la matrtiz de cuenta, que se empleará para loa análisis posteriores reflejados en la carpeta <a href="https://github.com/AlfonsoBravoPaneque/TFM--Analisis_RNA-seq_de_celulas_B/tree/main/An%C3%A1lisis_empleando_R" target="_blank">Análisis_empleando_R</a>.

</div>

<br> 

📁 **Análisis_empleando_R**

<div align="justify">

La carpeta <a href="https://github.com/AlfonsoBravoPaneque/TFM--Analisis_RNA-seq_de_celulas_B/tree/main/An%C3%A1lisis_empleando_R" target="_blank">Análisis_empleando_R</a> contiene los diferentes <i>scripts</i> realizados en R:
- Procesamiento_datos_y_expresión_diferencial.R: procesado de la matriz de cuentas y análisis de expresión diferencial, en función de la condición experimental, de las subpoblaciones de células B naive, activadas y del centro germinal. Además, tambien contiene el script empleado para realizar los *volcano plots*.
- ORA.R: análisis de enriquecimiento funcional de los diferentes genes sobreexpresados, en función de la condición experimental, en las subpoblaciones de células B del centro germinal y activadas.
- GSEA.R: análisis de enriquecimiento funcional de todos los genes diferencialmente expresados, en función de la condición experimental, en las subpoblaciones de células B del centro germinal y activadas. Previamente, se ordenó estos genes en función del Log2FC.
- Heatmap.R: mapa de calor de los 100 genes con mayor expresión promedio por tipo celular y condición experimental.
- PCA.R: PCA de las diferentes muestras empeladas (20 en total), mostrando los dos primeros componentes principales.
- Redes.R: inferencia de redes génicas, seleccionando los genes diferencialmente expresados en células B del centro geminal y activadas en función de la condición experimental.

</div> 

