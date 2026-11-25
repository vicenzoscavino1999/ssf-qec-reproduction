# Reproducción de *Silent logical failure in small-set-flip decoding*

Esta edición acompaña al manuscrito de 34 páginas incluido en `manuscrito/`. Un solo archivo contiene los datos, las semillas, los programas, los 1.412 fragmentos del catálogo de tiras y las referencias numéricas. La numeración de esta guía corresponde a la versión limpia del artículo.

## Primera comprobación

Use Python 3.12 de 64 bits. La versión utilizada en Windows fue 3.12.10. La ruta rápida utiliza la biblioteca estándar. Reserve 15 GB si va a conservar el ZIP junto a su contenido extraído. Desde PowerShell, en la carpeta en la que descargó el archivo, ejecute lo siguiente. Puede sustituir `D:\SSF` por otra carpeta corta y vacía.

```powershell
py -3.12 -m zipfile -e "REPRODUCCION_SSF_QEC_20261124.zip" "D:\SSF"
Set-Location -LiteralPath "D:\SSF"
py -3.12 reproducir.py rapido
```

En un entorno con Python 3.12 como intérprete predeterminado, `python` sustituye a `py -3.12`. También puede extraer el archivo con una utilidad que admita ZIP64. Abra la carpeta que contiene `reproducir.py`. Las pequeñas bibliotecas ZIP incluidas dentro de algunos módulos se leen directamente desde Python.

El último renglón muestra `OK` y la ruta a `resultado.json`. En la carpeta de salida encontrará los registros de cada cálculo. La acción `rapido` comprueba las huellas SHA-256, vuelve a calcular los intervalos Bernoulli a partir de los censos guardados y genera las filas de las tiras en `tablas_tiras/`. También verifica los cortes de F17 y de (3,4), el cero negativo, el ejemplo por fibras, los diagnósticos de las cadenas y 62 valores de las figuras 4–6 contra sus datos. El paso posterior para volver a llamar al decodificador se elige según la tabla siguiente.

Dos resultados que permiten reconocer la ejecución son los intervalos deterministas de la tabla 2. A `p = 10^-6` y `q = 10^-4`, el intervalo de 225 cúbits queda en `[3.24009, 3.27319] × 10^-14` y el de 900 en `[8.63039, 9.69644] × 10^-14`, con redondeo hacia fuera. El fichero `01-calculos_del_articulo.log` conserva más cifras.

**Quick start in English.** Extract the entire ZIP64 archive to a short path. From the directory containing this guide, run `python reproducir.py rapido` with 64-bit Python 3.12. The last line gives the result and output path. Run `python reproducir.py sectores` to rebuild Table 2 from the two seeds. A fresh output directory is used for each invocation.

## Del artículo al cálculo

Use `py -3.12 reproducir.py ACCION`, reemplazando `ACCION` por una de las entradas siguientes. La ruta predeterminada es una carpeta nueva dentro de `resultados_ssf`, junto al directorio extraído. Para fijarla, añada, por ejemplo, `--salida "D:\salidas_ssf\sectores-01"`. El programa indica qué registro consultar si falla un paso.

| Resultado en la versión limpia | Acción | Qué se vuelve a obtener |
| --- | --- | --- |
| §5, tabla 1 y apéndices B y C | `tasas` | Filas racionales de las tres esquinas (3,4), los dos órdenes del circuito, el traslado directo por fibras y los certificados Bernoulli y circuital de (17,18). Incluye los controles locales del primer contacto. |
| §5, cero negativo y cortes de precisión | `rapido` | Dos conteos del polinomio de grado cuatro, signos racionales que encierran el cero en `(-0.097080,-0.097079)` y los cortes suficientes `6482/6575` y `3028/3057` en 225/900 cúbits. |
| §6 y apéndice D, tabla 4 | `f17` | Construye las máscaras H1, reproduce la comparación implícita H3 y genera y comprueba las trazas de 1.214 configuraciones del sector afín. Verifica los cortes suficientes 253.000 y 126.000 y la contribución de circuito mayor que `1.02 × 10^-499` a las tasas del manuscrito. |
| §7.1, tabla 2, figura 4(b) y apéndice E.1 | `sectores` | Compila el núcleo C++17, construye las dos semillas desde cero, ejecuta los censos por órbitas y coteja los histogramas de clase y los intervalos con los censos directos guardados. Los conteos de (2,1) son 324 y 864. |
| §7.1, tabla 5, figura 4(a) y apéndice E.2 | `piloto` | Recalcula las tres filas Bernoulli y las ocho filas de circuito en 225 cúbits. La segunda implementación comprueba los intervalos, calendarios y conteos de las 16 muestras. La cota de riesgo conjunta de esas muestras es `16 exp(-7)`. |
| §7.2, tablas 3, 6 y 7, figuras 5(a) y 6 | `rapido` | Regenera las filas de error a partir de los extremos de las tiras y recupera los recursos medidos. En 523 cúbits, con tasas iguales a `10^-5` y corte cuatro, las cotas redondeadas son 1,221 % por componentes y 99,513 % con corte global. |
| §7, tablas 8 y 9 y apéndice E.4 | `comparaciones` | Reanaliza los registros guardados. En 225 reconstruye el transporte de Bennett y sus diagnósticos. En 61 el verificador original comprueba el expediente completo de la corrida. Se cotejan también la frecuencia exacta y el intervalo del estrato DSS. |
| Apéndice E.3, PATH25 | `path25` | Repite la comparación de dos decodificadores sobre 4.096 síndromes y el cálculo de 524.288 asignaciones. El estimador de circuito alcanza resto cero en el corte ocho. El verificador compara los ficheros científicos deterministas con los de referencia. |
| Figuras 1–6 y artículo completo | `pdf` | Ejecuta pdfLaTeX tres veces en la carpeta de salida. Las figuras geométricas usan TikZ y los gráficos usan los valores incorporados al TeX. El PDF distribuido conserva los bytes de la versión limpia. |

Los extremos uniformes para (17,18) valen para toda la familia. El extremo circuital `r_c = 5 × 10^-33, r_m = 2 × 10^-11` y la contribución del sector de la tabla 4 utilizan el F17 afín concreto. Las dos instancias de la tabla 2 tienen distancia seis. El producto de distancia diez citado en la discusión es una comparación bibliográfica.

## Repetir el muestreo o reconstruir los catálogos

| Acción | Cálculo y recursos |
| --- | --- |
| `splitting-225` | Compila el motor, construye el modelo VV primero, repite el ancla y los 28 niveles hasta `r_c = 2 × 10^-7, r_m = 10^-5`, reanaliza las cadenas y vuelve a medir el estimador de referencia. La corrida publicada de splitting tardó 161,91 s. |
| `dss-225` | Vuelve a construir las tablas del modelo y compila el núcleo. Ejecuta un millón de muestras con la semilla registrada y censa las 2.721.600 configuraciones del estrato (2,1). Se esperan 119 aciertos en la muestra y 324 en el censo. |
| `splitting-61` | Ejecuta las semillas, el ancla y los niveles del protocolo original para la tira de 61 cúbits. La corrida publicada duró 2,49 h. La acción `comparaciones` permite examinar antes los registros incluidos. |
| `tiras-1` y `tiras-2` | Reconstruyen los resultados de producción a partir de los catálogos incluidos. La primera corrida tardó históricamente 3,32 h, con presupuesto de 6,5 h. La segunda incluye 523 cúbits y dispone de 9 h. |
| `catalogos` | Regenera los 1.412 fragmentos y compara sus huellas con el catálogo distribuido. El presupuesto original fue de 8 h, 24 GiB de memoria agregada y hasta 60 GB de salida. |

Para repetir las dos corridas de las tiras, desde la raíz del paquete:

```powershell
py -3.12 reproducir.py tiras-1 --salida "D:\salidas_ssf\tiras-1"
py -3.12 reproducir.py tiras-2 --previa "D:\salidas_ssf\tiras-1" --salida "D:\salidas_ssf\tiras-2"
```

El segundo comando lee el resultado y el sello de la primera corrida. Los programas conservan sus presupuestos originales y dejan las salidas parciales para inspección si alcanzan un límite. Para estas corridas se recomiendan al menos 32 GiB de RAM y espacio adicional para las nuevas salidas. Las duraciones publicadas corresponden al equipo del artículo. El directorio `tiras/` conserva las especificaciones de producción y sus límites.

## Entorno e interpretación

- La reproducción integral se preparó para Windows x64 con Python 3.12.10. La comprobación rápida y los cálculos racionales usan la biblioteca estándar. Mantenga activadas las aserciones de Python.
- `sectores`, `splitting-225` y `dss-225` compilan código C++17. En Windows use Visual Studio Build Tools con el componente C++ x64 y una CPU compatible con AVX2/BMI2. `sectores` y `dss-225` también incluyen la compilación con GNU C++17 en x86-64. La plataforma ensayada para esta edición fue Windows.
- `pdf` requiere pdfLaTeX con los paquetes indicados en la cabecera del TeX, entre ellos Latin Modern, TikZ, PGFPlots, microtype y flushend. Las figuras y la bibliografía están incorporadas. Las versiones de TeX y los metadatos de fecha pueden cambiar la huella del PDF recompilado.
- Los conteos enteros y las desigualdades racionales se comprueban con la precisión de sus verificadores. Los gráficos se cotejan a la precisión de sus rótulos numéricos. Para los intervalos híbridos y las cadenas se conservan semillas, órdenes, recuentos y diagnósticos. Un `OK` en splitting acredita la reproducción del resultado publicado, incluido su diagnóstico de mezcla. El intervalo de referencia en 61 es determinista. El de circuito en 225 tiene la cobertura conjunta del piloto. Los tiempos y la memoria se miden de nuevo cuando corresponde.
- Un informe numérico se guarda por acción. Para repetir una acción basta ejecutarla otra vez con la salida automática o dar otro nombre a `--salida`.

## Comprobación de esta edición

En una extracción nueva en Windows se reconstruyeron los dos censos por órbitas, los 11 intervalos del piloto, el muestreo DSS, el splitting de 225 cúbits y los controles PATH25. El PDF recompilado tiene 34 páginas y su texto extraído coincide página por página con el PDF distribuido. Para las corridas de tiras de varias horas y el splitting de 61 cúbits se verificaron los resultados archivados. Los comandos de producción anteriores permiten repetirlos desde sus entradas.

Para ejecutar por separado el verificador original de la producción archivada de las tiras, use `py -3.12 tiras/reproduce.py verify-recorded`. Comprueba las huellas del catálogo y las 128 filas de resultados.

## Correspondencia con el apéndice F

| Paquete citado en el artículo | Contenido disponible aquí |
| --- | --- |
| `SOURCE_CUMULANT_CERTIFICATE.zip` | `certificados/f17/` |
| `FIBRE_CERTIFICATES.zip` | `certificados/fibras/` |
| `REPRODUCCION_PAPER_20261117_v3.zip` | `certificados/supplement/h1/` y `h3/`, los dos cálculos citados |
| `F17_CIRCUIT_SECTOR_20261107.zip` | `certificados/sector_f17/`, geometría y trazas. `calculos_del_articulo.py` evalúa las tasas de la versión limpia. |
| `SSF_SECTORES_225_900_20261124.zip` | `sectores/` |
| `HGP225_PILOT.zip` | `comparaciones/`, generador y verificador independiente de los 11 intervalos |
| `QEC_RARE_EVENT_COMPARISONS_20261124.zip` | `comparaciones/`, más los programas originales y los registros por cadena de 61 cúbits en `tiras/workspace/T2_2_*/` |
| `PATH25_VALIDATION_20261106.zip` | `path25/validacion/` |
| `PATH25_CERTIFIED_ESTIMATORS_20261112.zip` | `path25/estimadores/` |
| `STRIP_COMPONENT_CERTIFICATES_20261116.zip` | `tiras/tables/` |
| `T2_1_REPRODUCCION_COMPLETA_20261116.zip` | `tiras/`, incluidos los insumos de la producción completa |

Las carpetas conservan los nombres y las referencias internas de los códigos originales para mantener comprobables sus huellas. La tabla de resultados anterior identifica cada cálculo con la numeración vigente. El listado `SHA256.json` cubre todos los ficheros distribuidos salvo el propio listado.

