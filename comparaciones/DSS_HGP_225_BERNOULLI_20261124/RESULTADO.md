# Muestreo por subconjuntos en el código de 225 cúbits

24 de noviembre de 2026. Resultados de la tabla 13 del manuscrito y comparación con la tabla 2 de las tiras.

## Un millón de llamadas y censo completo

Una ronda de ruido Bernoulli independiente en el HGP afín \([[225,9,6]]\), con \(p=10^{-6}\) en los 225 datos y \(q=10^{-4}\) en las 108 lecturas. Se usa el mismo SSF/FGL global del piloto, \(\beta=1/100\), con su desempate fijo. Un fallo silencioso requiere síndrome físico residual cero y clase lógica residual no trivial. El estrato \((a,b)=(2,1)\) condiciona a que haya exactamente dos datos erróneos y una lectura invertida.

| Evaluación del estrato | Llamadas al evaluador SSF | Fallos | Frecuencia | Tiempo del proceso |
|---|---:|---:|---:|---:|
| Muestreo uniforme con reemplazo | 1.000.000 | 119 | \(1,19000\times10^{-4}\) | 0,345 s |
| Censo de todos los pares y lecturas | 2.721.600 | 324 | \(1/8400=1,190476\ldots\times10^{-4}\) | 0,795 s |

El intervalo binomial de Clopper–Pearson bilateral al 95 % de la muestra es \([0,98582,\,1,42400]\times10^{-4}\), redondeado hacia fuera. El punto obtenido quedó un 0,04 % por debajo del censo. Los márgenes del intervalo respecto del punto son 17,16 % y 19,67 % hacia abajo y hacia arriba. El censo reproduce las 324 configuraciones silenciosas del análisis anterior de pares.

Se fijaron un millón de llamadas y la semilla entera 202611242151 antes de ejecutar el muestreo. En cada llamada se sortean dos datos distintos y un índice de lectura. Cada par no ordenado tiene dos órdenes equiprobables. El tamaño de muestra es fijo y la extracción se repite con reemplazo entre llamadas. El censo recorre cada par no ordenado y cada lectura una sola vez en otro proceso. Ambos combinan las etiquetas físicas e invocan directamente la misma función C++ que decide el fallo en el experimento de *splitting* de 225 cúbits. Ninguno de los evaluadores usa la lista de las 324 configuraciones.

La generación de las tablas compartidas tomó 1,535 s y reprodujo el archivo del piloto, con 1.584 casos de comprobación del decodificador. Los dos tiempos de proceso incluyen carga y autocomprobación. Los bucles de evaluación duraron 0,308 y 0,764 s. La compilación y el análisis estadístico posterior quedan fuera de esos tiempos. Windows x64, Intel Core i5-14400F y MSVC con optimización y BMI2. [Registros numéricos](resultado.json) · [Preparación del modelo](preparacion.json)

## Qué aporta el estrato a la probabilidad total

Escribamos \(w_{ab}\) para la probabilidad del estrato y \(f_{ab}\) para su probabilidad condicional de fallo silencioso. Suponiendo sorteos independientes con reemplazo y tamaños fijos \(N_{ab}\), la varianza de la suma ponderada de los estratos muestreados es exactamente \(\sum_{a,b} w_{ab}^{2} f_{ab}(1-f_{ab})/N_{ab}\). Para este estrato,

\[
S=\binom{225}{2}108,\qquad
w_{21}=S p^2(1-p)^{223}q(1-q)^{107},\qquad
w_{21} f_{21}\simeq 3,20480035\times10^{-14}.
\]

El intervalo determinista del piloto para la probabilidad total es \([3,24009158,\,3,26238669]\times10^{-14}\), redondeado hacia fuera. Por tanto, este estrato aporta entre el **98,23 % y el 98,92 %** del fallo total. La muestra estima su contribución en \(3,20352\times10^{-14}\), con intervalo al 95 % de \([2,65387,\,3,83345]\times10^{-14}\). Las demás combinaciones de fallos se incluyen en la referencia total del piloto.

La varianza relativa de la frecuencia muestreada es \(8399/N\). Para una semiamplitud relativa nominal \(\epsilon\), la aproximación normal bilateral al 95 % da \(N\simeq 1,959964^2\,8399/\epsilon^2\). Al redondear al entero superior se obtienen 3.226.442 llamadas al 10 % y 322.644.127 al 1 %. Son presupuestos aproximados para muestreo con reemplazo del estrato. Para muestreo uniforme sin reemplazo se multiplica la varianza por \((S-N)/(S-1)\). El censo ya obtiene la frecuencia exacta con \(S=2.721.600\) llamadas.

## Cómo se comparan las tiras y el splitting

En las tiras de 61, 127, 259 y 523 cúbits, \(W\) cuenta sitios de respuesta efectivos no nulos. Para el mismo circuito, ruido, decodificador y filtro, el muestreo bajo la ley condicional correcta en cada peso \(w\le4\) converge a \(L=\sum_{w\le4}\Pr(W=w)\Pr(F\mid W=w)\). Si toda la masa del resto se deja sin resolver, el intervalo límite es \([L,L+\Pr(W>4)]\). En la tira de 523 cúbits con \(r_c=r_m=10^{-5}\), el margen relativo del punto óptimo es como máximo 99,513 % con este corte global y 1,221 % al incluir familias de componentes compatibles. Un cambio de unidad de fallo o un tratamiento adicional de la cola requiere su propio cálculo. [Heußen, Winter, Rispler y Müller, § I.2, ecuaciones 5–9](https://arxiv.org/html/2309.12774v1)

| Instancia de circuito | Comparación del punto de splitting con la referencia | Tiempo registrado |
|---|---|---|
| Tira de 61 cúbits | A \(r_c=r_m=10^{-5}\), un 34,46–36,07 % por encima del valor compatible con el intervalo determinista. A \(10^{-6}\), un 36,10–36,13 % por encima. | 2,49 h por las cuatro tasas en una escalera compartida con su verificación. |
| Afín de 225 cúbits, VV primero, \(r_c=2\times10^{-7}\) y \(r_m=10^{-5}\) | \(4,91234\times10^{-16}\) con splitting y \(3,34223\times10^{-16}\) por sectores. Si cubre el intervalo del piloto, el exceso de splitting está entre 46,85 % y 47,11 %. La cobertura simultánea de las ocho filas de circuito supera el 98,54 %. | Splitting en 161,91 s. Sectores en 8,31 s, mediana de tres corridas con preparación. |

Estas dos filas miden el fallo total de circuito. La fila Bernoulli de la tabla 13 mide \(f_{21}\). El control a \(p=3,75\times10^{-7},q=10^{-3}\) con medianas de 1,075 y 1,065 s queda como nota de pequeño volumen en el apéndice N.

## Reproducir

En el paquete **QEC_RARE_EVENT_COMPARISONS_20261124.zip**, extraer todos los archivos y entrar en la carpeta **DSS_HGP_225_BERNOULLI_20261124**. Python 3.12 y su biblioteca estándar bastan para el análisis. El binario suministrado es para Windows x64 con AVX2/BMI2.

~~~powershell
python -I -B comprobar_tabla.py
python -I -B ejecutar.py --output repeticion.json
~~~

El primer comando contrasta la fila de las tiras, las dos filas de splitting y la del estrato a partir de los registros. Comprueba también las colas binomiales del intervalo con una suma independiente de masas. El segundo reproduce los sorteos y el censo con la semilla publicada en dos procesos nuevos. Los recuentos y el intervalo son reproducibles. Cada ejecución mide sus propios tiempos.

Para recompilar en Windows, ejecutar **compilar.cmd** desde la misma carpeta con las herramientas de Visual C++ 2022 instaladas. En sistemas x86-64 con GCC, la orden equivalente es **g++ -O3 -std=c++17 -mbmi2 -mpopcnt estrato_225.cpp -o estrato_225**. Los resultados publicados se obtuvieron con MSVC. El paquete incluye el código fuente del decodificador. La orden siguiente reconstruye las tablas desde la geometría.

~~~powershell
python -I -B -c "from runpy import run_path; run_path('../SPLITTING_HGP_225_VV_20261124/modelo_225.py')['build']('.')"
~~~

La carpeta de la tira de 61 cúbits contiene los registros de la corrida previa y su dictamen. El paquete independiente de *splitting* en 225 cúbits incluye sus cadenas, histogramas y productor completos. El nuevo experimento es el muestreo Bernoulli y su censo.
