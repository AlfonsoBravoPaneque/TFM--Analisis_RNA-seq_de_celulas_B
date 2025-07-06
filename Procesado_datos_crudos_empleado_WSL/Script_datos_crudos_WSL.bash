## Procesamiento final de datos crudos ##

### PASO 1: CREAR EL ENTORNO Y DESCARGAR LOS SOFTWARES ###

# Crear un nuevo entorno en Conda para el análisis
conda create --name rna_seq_pipeline python=3.9 -y

# Activar el entorno
conda activate rna_seq_pipeline

# Instalar herramientas bioinformáticas esenciales
conda install -c bioconda fastqc trimmomatic hisat2 samtools subread -y # featureCounts es un paquete dentro de subread (paquete de anaconda que contiene a featureCounts y otros paquetes realacionados).
																	    # En bioconda no está el paquete featureCounts por separado, y por eso se instala con subread.

## Versiones de softwares ##
# fastqc v0.12.1
# hisat2 v2.2.1
# samtools v1.13
# featureCounts v2.0.1
# fastp v0.23.2


### PASO 2: CREAR DIRECTORIOS ###
# Se crean los directorios de forma manual con rutas absolutas
mkdir -p /home/alfonso/genome
mkdir -p /home/alfonso/data
mkdir -p /home/alfonso/results
mkdir -p /home/alfonso/results/quality_control
mkdir -p /home/alfonso/results/alignment
mkdir -p /home/alfonso/results/counts


### PASO 3: DESCARGAR EL GENOMA DE REFERENCIA Y LA ANOTACIÓN ###
# Se descarga el genoma de referencia y su anotación en /home/alfonso/genome
cd /home/alfonso/genome
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.primary_assembly.annotation.gtf.gz

# Descomprimir los archivos
gunzip /home/alfonso/genome/GRCm39.primary_assembly.genome.fa.gz
gunzip /home/alfonso/genome/gencode.vM36.primary_assembly.annotation.gtf.gz


### PASO 4: CONTROL DE CALIDAD CON FASTQC ### (UNA VEZ POR MUESTRA)
# Evaluar calidad de la muestra (se ejecuta para cada muestra individualmente)
fastqc /home/alfonso/data/SAMPLE.fastq.gz -o /home/alfonso/results/quality_control


### PASO 5: FILTRADO Y RECORTE (TRIMMING) CON FASTP ### (UNA VEZ POR MUESTRA)

# Se puede recortar las primeras bases usando fastp
# Argumentos de la función fastp: #
# -i: indica el archivo de entrada
# -o: indica el archivo de salida
# --trim_front1: recortar los nucleotidos indicados de cada lectura
# -a: indica la secuencia adaptadora a eliminar
# -l: indica el número mínimo de bases que tiene que tener una lectura. Si tiene menos la descarta
fastp -i /home/alfonso/data/SAMPLE.fastq.gz -o /home/alfonso/results/quality_control/SAMPLE_trimmed.fastq.gz -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -l 25 --trim_front1 8


# Para eliminar el adaptador se usa fastp. La detección del adaptador no siempre es tan efectiva para lecturas de un solo extremo, por lo que se recomienda proporcionar 
# la secuencia del adaptador. Además, se decicidó eliminar de todas las muestra el adptador, ya que había ocasiones donde el HTML porporcionado por FASTQC daba información
# de que no habia adpatadores en las lecturas, y después al hacer el filtrado con fastp se eliminaban algunas lectura que si poseían adaptadores. Además, también se eliminaron
# las lecturas que poseían menos de 25 pb, y las 8 priemras bases de cada lectura (esto se realizó ya que al prepara las librerías se introducian random hexamers, y esto perjudicaba a las primeras
# bases de cada lecturas de todas las muestras).

# Después se vuelve a comprobar la calidad(se ejecuta para cada muestra individualmente)
fastqc /home/alfonso/results/quality_control/SAMPLE_trimmed.fastq.gz -o /home/alfonso/results/quality_control

### PASO 6: INDEXAR Y ALINEAR LAS LECTURAS CON HISAT2 ###
# Indexación del genoma (SOLO SE HACE UNA VEZ. ESTE UN PASO MUY COSTOSO COMPUTACIONALMENTE)
hisat2-build /home/alfonso/genome/GRCm39.primary_assembly.genome.fa /home/alfonso/genome/hisat2_index/GRCm39 # hisat2-build generó 8 archivos de índice (.ht2)

# Alineación de la muestra (UNA VEZ POR MUESTRA)

# Argumentos de la función hisat2: # 
# -x: se especifica la indexacion anteriormente realizada
# -U: indica la muestra (fastq) que se quiere ALINEAR
# -S: indica el output se salida en formato SAM
# -p: indica el número de hilos empleandos en la alineación. Basicamente, permite paralelizar (indicando el numero de núcleos del CPU), permitiendo que el tiempo de computo sea menor

hisat2 -p 8 -x /home/alfonso/genome/hisat2_index/GRCm39 -U /home/alfonso/results/quality_control/SAMPLE_trimmed.fastq.gz -S /home/alfonso/results/alignment/SAMPLE.sam

### PASO 7: CONVERTIR SAM A BAM Y ORDENAR ### (UNA VEZ POR MUESTRA)

# Argumento de la función samtools view: #
# -b: indica que la salida debe ser un archivo BAM
#-@: indica el número de hilos empleados (disminuye el tiempo de cómputo)
samtools view -@ 8 -b /home/alfonso/results/alignment/SAMPLE.sam > /home/alfonso/results/alignment/SAMPLE.bam

# Argumento de la función samtools sort (ordena un archivo BAM, para después facilitar la cuantificaión): #
# -o: permite redirigir la salida ordenada a otro archivo
#-@: indica el número de hilos empleados (disminuye el tiempo de cómputo)
samtools sort -@ 8 /home/alfonso/results/alignment/SAMPLE.bam -o /home/alfonso/results/alignment/SAMPLE_sorted.bam

# Eliminar el archivo SAM para ahorrar espacio (UNA VEZ POR MUESTRA). Este paso es opcional, pero se puede eliminar ya que el archivo pesa en torno a 10 Gb
rm /home/alfonso/results/alignment/SAMPLE.sam

# Eliminar el archivo BAM para ahorrar espacio (UNA VEZ POR MUESTRA), ya que nos interesa solo el _sorted.bam
rm /home/alfonso/results/alignment/SAMPLE.bam

### PASO 8: CONTAR LAS LECTURAS CON FEATURECOUNTS ### (UNA VEZ POR MUESTRA). Para ello, se emplea el archivo GTF original.
# Argumento de la función featureCounts: el archivo de salida suele ser un archivo tabulado que se puede abrir en Excel#
# -T: indica el número de hilos empleados (disminuye el tiempo de cómputo)
# -t: indica la caracteristica que se va a cuantificar. En mi caso uso los exones, ya que los fragmentos de RNA se alinean sobre los exones y se agrupan por genes. CUENTA LAS LECUTAS QUE CAEN EN EXONES
# -a: nombre del archivo de anotación (archivo GTF). En Gencode vM36, el identificador estándar de los genes es gene_id, por lo que este es el atributo más preciso para obtener cuentas por gen. AGRUPA LAS LECTURAS POR EL ID DEL GEN EN ELA RCHUVO GTF
# -o: salida de la matriz de cuentas (en formato .txt)
# -O: si una lectura cae en varios genes, se cuenta en todos ellos
# --largestOverlap (modifica al argumento -O): si una lectura solapa con múltiples genes, se asigna solo al gen con mayor solapamiento.
# -M: cuenta la lecturas multipareadas (NO ES ÓPTIMO PARA MI RNA-SEQ, YA QUE DESPUÉS VOY A REALIAR ANÁLISIS DE EXPRESIÓN DIFERENCIAL).
# --fraction: divide las cuentas proporcionalmente

# Debido a que mi objetivo es hacer un análisis de expresión diferencial, creo que el enfoque más óptimo es:
# -contar lecturas a nivel de exón, asignándolas al gen con mayor solapamiento
# -evitar sesgos exluyendo el multi-mapeos o distribuir proporcionalmente las lecturas (no usar ni -M ni --fraction)
# -si una lectura solapa con más de un gen se cuenta la lectura al gen con más solapamiento

featureCounts -T 8 --largestOverlap -t exon -g gene_id -a /home/alfonso/genome/gencode.vM36.primary_assembly.annotation.gtf -o /home/alfonso/results/counts/gene_counts_SAMPLE.txt /home/alfonso/results/alignment/SAMPLE_sorted.bam

# Se pueden agrupar todas las muestras en la misma matriz de cuentas:
featureCounts -T 8 --largestOverlap -t exon -g gene_id -a /home/alfonso/genome/gencode.vM36.primary_assembly.annotation.gtf -o /home/alfonso/results/counts/gene_counts_ALL_SAMPLES.txt /home/alfonso/results/alignment/*_sorted.bam
