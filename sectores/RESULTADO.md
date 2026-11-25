# Fallo silencioso Bernoulli en los códigos conectados de 225 y 900 cúbits

Una ronda del SSF global, con \(p=10^{-6}\), \(q=10^{-4}\) y \(\beta=1/100\). El suceso es \(H_XR=0\) y \([R]\ne0\). Las cotas incluyen cualquier número de errores. La estimación publicada es la media armónica de los extremos exactos. Todos los pesos y las cotas del resto se calculan con fracciones.

| Código | Intervalo para la probabilidad total | Estimación | Error relativo máximo | Tiempo |
| --- | --- | --- | --- | --- |
| \([[225,9,6]]\) | \([3.24009,3.27319]\,10^{-14}\) | \(3.25656\,10^{-14}\) | 0,51 % | 1,10 s |
| \([[900,36,6]]\) | \([8.63039,9.69644]\,10^{-14}\) | \(9.13241\,10^{-14}\) | 5,82 % | 32,79 s |

En el Core i5-14400F, el pico de memoria del proceso de preparación fue de 22,4 y 55,8 MiB, y el del ejecutable de recuento fue de 8,8 y 22,8 MiB. La suma por sectores hizo 1.579.785 y 36.739.110 consultas a configuraciones u órbitas, respectivamente. Los controles exhaustivos se midieron por separado.

El certificado publicado previamente en 225 cúbits, \([3.24009158171737,3.26238668184343]\,10^{-14}\) a la precisión archivada, queda dentro del nuevo intervalo. La tabla 6 conserva esa cota más ajustada. En el caso de 900 se cuadruplica el volumen con la distancia fijada en seis.

## Reproducción

Ejecutar `python -B reproducir.py` desde este directorio. Requiere Python 3.10 o posterior y un compilador C++17. En Windows usa las herramientas de compilación de Visual Studio. `--sheets 2` selecciona 900 cúbits. `run/` recibe matrices, archivos de órbitas, histogramas por clase y los dos JSON de resultados. La carpeta `referencia/` del paquete guarda los registros medidos. `python -B contrastar.py run referencia` compara los resultados independientes de la máquina, incluidos todos los histogramas. Para el control directo opcional, ejecutar `censo run/semilla_2.bin --direct` y `censo run/semilla_2.bin --allpairs1` (usar `.\censo.exe` en PowerShell).

La medición parte de la semilla, construye las tablas y las órbitas, ejecuta el censo y calcula el intervalo. El binario compilado es la única preparación reutilizada. Se registra el pico de conjunto de trabajo por separado para Python y para el proceso C++. Las pruebas aleatorias de comparación de síndromes del generador tienen semilla fija y no intervienen en las probabilidades.

## Semilla, clases y desempate

La semilla de una hoja tiene las doce rectas del plano afín de orden tres y nueve puntos. Se numeran como en el piloto. Con \(t\) hojas se reemplaza la incidencia del punto 0 de la recta 6 de cada hoja por el punto 0 de la hoja siguiente, cíclicamente. El grafo es conexo. La arista sustituida pertenece al ciclo punto 0, recta 0, punto 3, recta 10, punto 5, recta 6. Al retirar la columna 6, las restantes tienen rango nueve en cada hoja, lo que da rango completo para el levantamiento. Las rectas de las dos primeras clases paralelas de una hoja dan una palabra \(c\) de peso seis. La fila VV cero con columnas \(c\) da un lógico de peso seis porque \(e_0\cdot c=1\). Para \(t=1,2\) se enumeran los \(2^{3t}\) elementos del núcleo y se verifica distancia clásica seis. La proyección \(a^\top U\) de cualquier residual físico de síndrome cero pertenece al núcleo para todo \(a\in\ker H\). La base lógica del autoproducto de rango completo prueba que alguna proyección es no nula para todo residual no trivial. Esto da distancia \(X\) seis, y la simetría por transposición da la misma distancia en \(Z\).

El orden de las coordenadas físicas es VV y luego CC, por filas. En cada paso gana el máximo descenso por unidad de huella, luego el mayor número de cúbits volteados y, finalmente, la menor máscara entera física. Las marcas se calculan en una base del núcleo y sus coordenadas de información. Las representaciones de una misma huella tienen la misma marca, comprobada durante la construcción. El síndrome físico y la marca de la corrección se contrastan con un recorrido directo del catálogo.

## Reducción de dos datos

Sea \(D\) el conjunto de lecturas alteradas. En un fallo silencioso el síndrome observado termina en \(D\), que debe ser estable. Una franja VV de tres checks contiene a lo sumo un bit de \(D\), pues con dos luces el singleton sería elegible. El descenso total desde la entrada es \(B_0=|H_XE|-2|H_XE\cap D|\). Cada paso desciende un entero positivo y, en cada proyección lógica \(a^\top U\), cambia como máximo una coordenada, la columna del generador.

Con \(v\) datos VV y un total \(a_0\le1\) de datos, \(B_0\le4a_0-v<6-v\), insuficiente para una palabra final no nula. Para dos datos, sean \(o\) el tamaño de la intersección de las huellas y \(h=|H_XE\cap D|\). Respecto al mínimo de \(6-v\) pasos, el exceso disponible para la suma de los descensos menos uno es a lo sumo \(2-2(o+h)\).

Si \(o=0\), se tiene \(h\le1\) y queda encendida por completo una huella original. El primer ganador es perfecto, tiene área al menos tres y consume todo el exceso. Por tanto \(h=0\), el primer paso tiene área tres y todos los posteriores tienen descenso uno. Si \(o>0\), necesariamente \(o=1,h=0\), con descenso uno en todos los pasos. En ambos casos hay exactamente \(6-v\) pasos. La igualdad en la proyección exige que las columnas de los \(v\) datos VV pertenezcan a una misma palabra de peso seis, sean distintas y no se volteen en esos pasos.

En el primer caso la franja ganadora de área tres está en una columna nueva. Por ausencia de ciclos de cuatro, cada franja original aporta a lo sumo un punto a esa franja, y \(D\) aporta a lo sumo uno. Ambas originales deben tocarla. Un CC original aún tendría al menos tres de sus cuatro luces tras este paso y superaría la puntuación de cualquier paso de descenso uno (densidad a lo sumo \(2/3\)). En el segundo caso un CC original ya tendría tres luces al comienzo y da la misma contradicción. Quedan dos VV en la misma fila, sea por el contacto con la nueva franja, sea por su propio solapamiento. Se obtienen 792 o 1 728 pares marcados, respectivamente.

Para cada par, el censo con \(D=\varnothing\) da cero. Si el par pudiera fallar con \(D\ne\varnothing\), debe haber una lectura en su fila, o dos lecturas en algún rectángulo que toque esa fila. En efecto, para franjas originales disjuntas el primer paso perfecto exige una luz de \(D\) en la fila. En el caso con solapamiento todos los descensos deben ser uno. Mientras haya un singleton elegible en la fila, un ganador con ese descenso tiene área tres y sigue el mismo desempate en la fila que con \(D=\varnothing\). Si ya no hay tal singleton, cada línea de la fila contiene a lo sumo una luz. Con una sola lectura externa por rectángulo que toca la fila habría como máximo dos luces en cualquier huella relevante, insuficientes salvo en un singleton VV. Fuera de la fila el estado leído es el estable \(D\). Por inducción una ejecución silenciosa en estas condiciones también habría fallado con lectura perfecta. Este argumento permite múltiples lecturas en una fila del levantamiento.

Para \(b=|D|\ge3\) la unión sobre una lectura en la fila más otras dos, o un par en un rectángulo de la fila más otra lectura, da por par el factor
\[
 q^3\{B {m-1\choose2}+3A{12\choose2}(m-2)\},
\]
con \(A=12t,B=9t,m=AB\). Se multiplica por el peso de cada par concreto de datos y por su número.

## Tres y cuatro datos

Los censos exactos de los sectores \((a,b)=(2,0),(2,1),(2,2),(3,0)\) dan \(0,324,32076,36\) en 225 y \(0,864,364176,144\) en 900. Cada ejecución aporta una vez a la clase terminal no nula correspondiente. Para los tríos sin medidas se descartan las órbitas donde ninguna trayectoria sobre los empates puede fallar; en las restantes se ejecuta el SSF con el orden físico para cada miembro. Los dos censos directos sin esta reducción recuperan también los histogramas.

Con tres datos y una medida se permite cualquier empate después de maximizar puntuación y tamaño en cada paso. Las simetrías son independientes en los dos ejes. En la base se usan las 432 transformaciones afines; en el levantamiento se permiten los dos estabilizadores del punto 0 y la recta 6 y el intercambio de hojas, de orden total 288. El programa verifica cada permutación en la matriz. Se divide el producto de multiplicidades por el estabilizador de las permutaciones internas de los datos del mismo tipo. La suma recupera \({n\choose3}m\) en cada código. Así se evita suponer que una simetría respeta la prioridad entre máscaras físicas.

Si el primer ganador elimina un subconjunto de los datos originales módulo estabilizadores, la ejecución posterior equivale a la del subconjunto restante con las mismas lecturas. La igualdad se comprueba por síndrome y marca. A lo sumo \((n-2)C_{21}\) tríos con una lectura fallan por esta vía, donde \(C_{21}\) es el censo físico de parejas. La exploración de todos los empates del resto aporta, como techo, 93 312 o 1 260 288 configuraciones, respectivamente. Cada arista de esa exploración reduce el peso del síndrome. Se aceptan exclusivamente hojas con síndrome observado \(D\) y marca lógica no nula.

Para \((3,2),(3,3),(4,0),(4,1)\) se usan dos condiciones necesarias para que el primer paso no sea una eliminación genuina. O bien ninguna franja original está completamente encendida, o bien hay una acción perfecta que no elimina un subconjunto original. Para la primera condición cada dato debe tocar otro dato o una lectura. Un dato toca como máximo otros 24 y cuatro lecturas. En un trío con algún par que se toca, el tercer dato toca a ese par o a las lecturas. Sin tal par los tres pertenecen a los como máximo \(7b\) datos que tocan las lecturas. Bajo la ley uniforme del estrato, la probabilidad queda acotada por
\[
  3\,\frac{24}{n-1}\left(\frac{48}{n-2}+\frac{4b}{m}\right)
    +\frac{\binom{7b}{3}}{\binom n3}.
\]
Con cuatro datos basta imponer el contacto del primero, de probabilidad a lo sumo \(72/(n-1)+4b/m\).

Para la segunda condición se fija una huella candidata dentro de su rectángulo de doce celdas. Cada celda recibe dos datos del generador y cinco datos exteriores. Cada dato exterior afecta a lo sumo una celda del rectángulo por ausencia de ciclos de cuatro. Se enumeran los siete bits interiores y se descarta un patrón si incluye una representación de la acción o su complemento en el estabilizador del generador. Para las \(d\) celdas de la huella se impone paridad iluminada exacta con cinco posibles datos exteriores y una lectura por celda. Los restantes \(n-7-5d\) datos y \(m-d\) lecturas completan el estrato binomialmente. Una huella en una sola fila se repite en tres rectángulos y en una sola columna se repite en cuatro. Una huella que ocupa al menos dos filas y dos columnas determina el rectángulo. La suma sobre las huellas mayoriza las incidencias espurias. Es el cálculo entero de `false_perfect` en `mayorante.py`.

Se añaden los fallos posibles tras una eliminación genuina, multiplicando cada fallo de dos o tres datos por el número de formas de completarlo al peso inicial. Finalmente se carga toda la masa restante de tres datos con al menos cuatro lecturas, de cuatro datos con al menos dos lecturas y de al menos cinco datos. La suma de los cuatro censos positivos da la cota inferior. Las sumas superiores y los límites de unión incluyen cualquier número de fuentes. Los JSON desglosan cada término.

## Alcance del contraste con el compilador uniforme

El mayor número de fuentes en un estrato sometido a recuento físico u orbital es cuatro. Algunas cotas combinatorias usan hasta seis fuentes y las colas cubren los demás pesos. En los dos tamaños el testigo de peso seis da \(\ell=(1+999999^6)^{-1}\). Para el error del 10 %, la tolerancia de numerador es \(\zeta=\ell/1280\). Reservando media tolerancia logarítmica, la cota uniforme \(N e^{-0.015(M+1)}\) da los órdenes suficientes 6 436 y 6 529. El programa verifica ambos con una suma de Taylor racional inferior para la exponencial. Estos órdenes son los de esa receta de prueba. Los tiempos medidos corresponden a las sumas positivas por sectores. El coste práctico al crecer la distancia y la evaluación directa de la serie uniforme siguen abiertos.

