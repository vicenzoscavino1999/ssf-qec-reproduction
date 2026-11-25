"""H5-A en grados (3,4). Un límite a la expansión en torno a ruido cero.
18 de noviembre de 2026. Reproducir con Python 3.10 o posterior, sin -O.

Para el producto de una semilla binaria H₀ consigo misma, fijar la primera
coordenada i de un check X y tomar el subhost J de las 18 mediciones (i,r).
Las fuentes en J sólo pueden ser errores de medición D. El síndrome observado
empieza en D. Si pasa el filtro de síndrome físico residual cero también
termina en D. Como cada paso reduce estrictamente su peso, pasa el filtro
exactamente cuando el decodificador se detiene sin dar pasos. En ese caso su
soporte es D y su marca es +1 para cualquier máscara.

Los únicos checks encendidos que un generador puede ver están en una fila
de su rectángulo 4×3, en los checks vecinos de una columna de H₀. Dos o más
luces en esa fila habilitan la acción singleton VV de huella 3, con puntuación
≥1/3. Con una o cero luces en cada rectángulo ninguna acción es elegible,
pues la huella mínima es 3. Todo par de checks de J que comparta columna de
H₀ dispone de un generador que contiene su fila (elegir cualquier vecino de i).
Por tanto D pasa si y sólo si ocupa a lo sumo un check por columna de H₀.
Esta reducción usa el SSF global del manuscrito y β=1/100.

La semilla es un recubrimiento doble de la incidencia punto-recta de F₃².
Para cada incidencia r∈j y t∈{0,1} conectar (j,t) con (r,t+rj mod 2),
usando las etiquetas de seed_columns(). La matriz binaria tiene 18 filas,
24 columnas, grados (3,4) y rango fila pleno. La partición filtrada es

    Ξ(J) = I(z) = 1 + 18z + 81z² + 60z³ + 2z⁴,       z=q/(1−q).

Con q=−1/10 resulta I(−1/11)=−174/14641, mientras I(0)=1. El programa
localiza un cero en −0.097080 < q < −0.097079. El logaritmo de la partición
es singular en ese cero, luego su serie de Taylor centrada en q=0 tiene
radio <0.097080. La normalización previa que agrupa las fuentes individuales
del manuscrito da (1−q)^18 I(q/(1−q)), con el mismo cero. En este subhost
b_c=q para todo check, de modo que la expansión absoluta en bloques de
fuentes de ese logaritmo tampoco puede converger en q=+0.1. Tal convergencia
absoluta la haría analítica en todo el disco abierto |q|<0.1. El desarrollo
en las razones z alrededor de cero también encuentra el cero negativo antes
del valor z=1/9 que corresponde a q=+0.1.

H5-A exige acceso en todos los subhosts inducidos. Este ejemplo excluye
llegar uniformemente al 10 % mediante esas expansiones directas. El
polinomio I es positivo para todo q físico en [0,1). Su evaluación directa
ilustra que el cero negativo limita una representación. La posibilidad de
otras formas de estimación y el umbral de funcionamiento del código son
cuestiones distintas. El ejemplo tampoco determina el mejor radio uniforme.

Base local: apéndices A y B y definición de acceso filtrado de
BORRADOR_INTRODUCCION_Y_RESULTADOS_EN_v2.tex. El rectángulo (3,4) ya derivado
en H5A_grados_3_4_cotas_y_obstruccion.txt es independiente de este cálculo.
"""

from fractions import Fraction
from functools import lru_cache


def seed_columns():
    columns = [0] * 24
    for x in range(3):
        for y in range(3):
            point = 3*x + y
            # La recta de pendiente m tiene etiqueta 3m+b; las verticales, 9+x.
            lines = [3*m + (y - m*x) % 3 for m in range(3)] + [9+x]
            for line in lines:
                for layer in range(2):
                    check = 2*point + (layer + point*line) % 2
                    columns[2*line + layer] |= 1 << check
    return columns


def validate_seed(columns):
    rows = [sum(1 << j for j, col in enumerate(columns) if col & (1 << i))
            for i in range(18)]
    assert len(columns) == 24 and all(col.bit_count() == 3 for col in columns)
    assert all(row.bit_count() == 4 for row in rows)
    pivots = {}
    for row in rows:
        while row:
            lead = row.bit_length() - 1
            if lead not in pivots:
                pivots[lead] = row
                break
            row ^= pivots[lead]
    assert len(pivots) == 18


def check_local_rule():
    # Construir todas las acciones por XOR físico y probar los cuatro lugares
    # posibles de la fila horizontal ocupada y sus ocho configuraciones.
    stripes = [7 << (3*row) for row in range(4)]
    stripes += [sum(1 << (3*row + col) for row in range(4)) for col in range(3)]
    candidates = []
    for mask in range(1, 1 << 7):
        footprint = 0
        for site, stripe in enumerate(stripes):
            if mask & (1 << site):
                footprint ^= stripe
        size = footprint.bit_count()
        if size and 2*size >= 3*mask.bit_count():
            candidates.append((footprint, size))
    assert len(candidates) == 101 and min(size for _, size in candidates) == 3
    for row in range(4):
        for pattern in range(8):
            syndrome = pattern << (3*row)
            eligible = any(100*(2*(syndrome & f).bit_count()-size) >= size
                           for f, size in candidates)
            assert eligible == (pattern.bit_count() >= 2)


def count_fault_sets(columns):
    counts = [0] * 19
    for faults in range(1 << 18):
        if all((faults & col).bit_count() <= 1 for col in columns):
            counts[faults.bit_count()] += 1
    return counts


def graph_polynomial(columns):
    # Cálculo alternativo: conjuntos independientes, excluyendo o incluyendo
    # un check y eliminando en este último caso sus vecinos.
    closed = []
    for i in range(18):
        sites = 1 << i
        for col in columns:
            if col & (1 << i):
                sites |= col
        closed.append(sites)

    @lru_cache(maxsize=None)
    def polynomial(mask):
        if not mask:
            return (1,)
        v = mask.bit_length() - 1
        without_v = polynomial(mask ^ (1 << v))
        with_v = (0,) + polynomial(mask & ~closed[v])
        size = max(len(without_v), len(with_v))
        return tuple(x+y for x, y in zip(without_v+(0,)*(size-len(without_v)),
                                         with_v+(0,)*(size-len(with_v))))

    return polynomial((1 << 18) - 1)


def main():
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para conservar las comprobaciones.")
    columns = seed_columns()
    validate_seed(columns)
    check_local_rule()
    counts = count_fault_sets(columns)
    poly = graph_polynomial(columns)
    assert tuple(counts[:len(poly)]) == poly == (1, 18, 81, 60, 2)
    assert all(n == 0 for n in counts[len(poly):])

    def partition(q):
        z = q/(1-q)
        return sum(coef * z**degree for degree, coef in enumerate(poly))

    assert partition(Fraction(0)) == 1
    assert partition(Fraction(-1, 10)) == Fraction(-174, 14641)
    assert partition(Fraction("-0.097080")) < 0 < partition(Fraction("-0.097079"))
    print("Semilla binaria (3,4) de rango 18. Regla SSF verificada.")
    print("Enumeración de 262144 subconjuntos y recurrencia coinciden:", poly)
    print("Cero real entre -0.097080 y -0.097079. Radio del logaritmo < 0.097080.")


if __name__ == "__main__":
    main()
