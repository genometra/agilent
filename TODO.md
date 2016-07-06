

Cosas para hacer por orden de prioridad
================================================================================


Hacer que detecte automaticamente la coumna de foreground.
Esta puede ser:

- "gMeanSignal"
- "gMedianSignal"

Si existen las dos usar la de "gMedianSignal"

Creo que habra que cambiar solamente la funcion `readAgilent.r`

--------------------------------------------------------------------------------







Quitar los warnings cuando no encuentra algunas columnas.... pasarlo a mensajes.



Incluir dependencia de affy

falla en estas funciones

normexp.fit
normexp.signal




Use the function colMedians form pakage [matrixStats](http://cran.fhcrc.org/web/packages/matrixStats/index.html)
within the function averageDuplicatedRows
