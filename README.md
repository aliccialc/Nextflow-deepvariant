# DeepVariant en Nextflow

El flujo de trabajo está dividido en 4 módulos:
### FORMAT
*Input* es el directorio donde tengas los archivos bam localizados.

El objetivo de este módulo es crear una lista donde se encuentren los nombres   de los achivos bam, sin la terminación del formato del archivo.

### DEEPVARIANT
*Input* :
val(sample) = Buscará de la lista creada cada documento llamado de ese modo. path(bam) = el directorio donde se encuentran los bam.
path(ref2) = El documento fasta con el genoma de referencia.
path(ref_fai2) = El docuemento fai del genoma de referencia.
path(ref_dict) = Necesario ya que se encuentra el HEADER del documento fasta
path(ref_gzi) = Formato comprimido e indexado

El objetivo del módulo es la detección de variantes utilizando DeepVariant. Los archivos generados son en formato .vcf

### FILTER_VCF
*Input*:
Tupla con los val(sample) y el procesado de cada fichero del directorio vcf.

El objetivo del módulo es filtrar aquellas variantes que han pasado el filtro de calidad  (FILTER = PASS) y que se encuentren en los cromosomas seleccionados.

### CLINVAR
*Input* :
En este caso, encontramos más de un archivo Clinvar, ya que se estaba realizando una comparativa entre dos fechas.


