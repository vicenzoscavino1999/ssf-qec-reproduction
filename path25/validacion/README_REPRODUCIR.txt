REPRODUCCIÓN DEL PILOTO PATH25

Artefactos principales:
  INFORME_PATH25.txt           resultados, interpretación y decisión;
  LEMA_DENSIDAD_ORDEN.txt       prueba de la cota uniforme para el orden fijado;
  REVISION_RESULTADOS.json     recepción posterior sin repetir el decoder;
  attempt_01/                  tablas, espectros, catálogo, trazas y costes;
  inputs/                     diseño, preflight y código P1 congelados;
  AUDIT.json                  hashes y alcance, producido al cerrar la unidad.

No ejecutar supervisor.py encima del intento recibido. Para una reproducción
nueva, crear una carpeta vacía aparte y copiar allí worker.py, supervisor.py e
inputs/ completos. Usar Python con biblioteca estándar en Windows y ejecutar:
  python -B -X utf8 supervisor.py
Se creará attempt_01/; si ya existe, se rechazará la ejecución. No hay reintentos
automáticos. El código importa funciones de la copia congelada de P1 sin ejecutar
su main ni escribir en sus carpetas originales.

El supervisor necesita Windows para consultar memoria y lanzar un proceso oculto.
Su límite de 60 s puede detener una reproducción en una máquina más lenta; una
parada es un resultado de recursos, no evidencia matemática negativa. Conservar
los archivos del intento y revisar el coste antes de cualquier nueva ejecución.

Resultados exactos esperados, independientes de las mediciones de tiempo:
  4096 síndromes comparados, 8192 clases, 12438 consultas locales;
  376 decoraciones en 1664 soportes interrogados;
  565 fuentes revisadas para recomposición, 123 multicomponente, 111 adaptativas;
  candidatos 103, huellas distintas 43, eta=1, densidad refinada 1/3;
  primeros coeficientes de fallo total: [0,4,115]; silencioso: [0,0,14].

Para la recepción matemática, copiar también revisar_resultados.py a la carpeta
de reproducción y ejecutarlo después del supervisor. No vuelve a ejecutar el
decoder. Los hashes de tablas binarias uint64 dependen del orden nativo de bytes;
los espectros enteros y tablas JSON permiten la comparación semántica.

preparar.py y el registro de controles pertenecen al workspace original y NO son
necesarios para reproducir la ciencia en una carpeta independiente. El manuscrito,
las entregas F17 y los censos N11/N12 no se necesitan ni se modifican.
