



Pipeline de trabajo
================================================================================

- Desde github clonar la rama "devel"
- Desde esta version clonada de github la copiais a vuestra maquina (`git clone`)
- Trabajais en local
- Subis un commit con los cambios a vuestra rama de github
- Desde github haceis un "pull request" que llegara a la rama devel de Genometra
- Actualizo la rama devel de Genometra y mas adelante la rama master


Trabajo en local
--------------------------------------------------------------------------------

El script `make.r` deberia compilar la libreria, vignetas incluidas y demas.

Una vez se ejecuta el script `make.r` aparecera un fichero `check/agilent.Rcheck/00check.log`.
En este fichero no debe haber ningun error y cuantos menos warnings o notas mejor.

La documentacion de cada una de las funciones se edita en el script R de la propia funcion.
NO en el directorio `man`.
Los ficheros del directorio man se generan automaticamente en el `make.r`
con la libreria `roxygen2`
<https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html>


Cosas para hacer por orden de prioridad
================================================================================

[x] Revisar ejemplos. Algunos fallan cuando se ejecutan. Esto se ve en los logs de la compilacion del paquete.

--------------------------------------------------------------------------------

[x] Hacer que detecte automaticamente la columna de foreground.

Esta puede ser:

- "gMeanSignal"
- "gMedianSignal"

Si existen las dos usar la de "gMedianSignal"

Creo que habra que cambiar solamente la funcion `readAgilent.r`

Esto incluye la documnetacion. Mas o menos con lo que hay hecho deberiais poder reescribirla...
pero por si acaso este es el paquete que se utiliza para procesar la documentacion <https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html>

--------------------------------------------------------------------------------

[x] Quitar los warnings.

Cuando escribi el paquete puse muchos mensajes y warnings (funciones de R `message` y `warning`)
la mayoria de ellos deberian cambiarse por funciones `print` o mejor `cat`.

--------------------------------------------------------------------------------

Vigneta

--------------------------------------------------------------------------------

Revisar si es necesaria la dependencia del paquete `genefilter`.
No se por que la puse en su momento... creo que por alguna prueba o dato adicional para la documentacion.

--------------------------------------------------------------------------------

Revisar el codigo y ayuda para que 

--------------------------------------------------------------------------------

<!--
Use the function colMedians form pakage [matrixStats](http://cran.fhcrc.org/web/packages/matrixStats/index.html)
within the function averageDuplicatedRows
-->
