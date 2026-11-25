"""Comprobaciones de H5A_comparacion_de_ejecuciones.txt. Python 3.10+."""

from fractions import Fraction as F
from math import factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import verificar_H5A_resto_global as base

prior = base.prior
P0, P, Q = base.P, F('1e-30'), base.Q
odds = prior.actividad
exp_upper = prior.exp_arriba


def exp_lower(x, order=8):
    return sum(x**j / factorial(j) for j in range(order + 1))


def geometry():
    matrix, _, _ = prior.geometria()
    vv = max(sum(row[:2]) for row in matrix[:2])
    vc = max(row[3] for row in matrix[:2])
    cv, cc = sum(matrix[3][:2]), matrix[3][3]
    return matrix, vv, vc, cv, cc


def first_difference():
    _, vv, _, cv, cc = geometry()
    assert cc == 83810
    sv, sm = (odds(r)*exp_upper(F('.015')) for r in (P, Q))
    assert 2*sv*vv < F('9.52e-23') and 2*sv*cv < F('5.62e-24')
    assert 2*sm*cc*909/8 < F('.000387') and 2*sm*cc/8 < F('4.26e-7')


def two_data():
    # Volver a comprobar el punto de partida y estrechar la fila exterior
    # de medida usando las mismas cotas para núcleo y ramas interiores.
    base.verify()
    _, vv, vc, cv, cc = geometry()
    lv, lc = F('5e-36'), F('2.05065e-11')
    s = base.SOURCE_AND_ACTION
    core = max(s[0]+s[3]+17*s[5], s[1]+s[4]+18*s[5])
    check = s[2]+35*core
    tilt = F('.025')
    yp, yq = odds(P0), odds(Q)
    assert lv >= yp*exp_upper(tilt+vv*(lv+core)+vc*(lc+check))
    assert lc >= yq*exp_upper(tilt+cv*(lv+core)+cc*(lc+check))
    assert base.XI >= yq*exp_upper(tilt+cv*lv+cc*lc)
    assert exp_upper(tilt+135*(vv*lv+vc*lc)) < F('1.03')
    rv, rm = core+18*lv, check+lc
    ke, kd = 135*(vv*rv+vc*rm), cv*rv+cc*rm
    assert ke < F('.01') and kd < F('.000002')
    # Quitar la decoración de una medida con todo su presupuesto antiguo.
    # El clúster anclado en esa decoración sólo precisa el exponente kd.
    mixed_m = rm-yq*exp_lower(tilt) + yq*exp_upper(F('.015'))*(exp_upper(kd)-1)
    assert 0 < mixed_m < F('1.25e-14')
    scale = (odds(P)/yp)**2
    return rv*scale, mixed_m*scale


def tree_rows(x, xi, y, full, planted, av, am):
    """Filas de H5A_p_general con presupuestos propios de dato y medida."""
    matrix, vv, vc, cv, cc = geometry()
    order = (0, 1, 3, 0, 1, 2)
    m = [[matrix[i][j] for j in order] for i in order]
    width, readers = (17,18,1,17,18,306), (18,17,1,18,17,306)
    h = [[w*l for l in readers] for w in width]
    need = (9,10,17)
    factors = [v/k for v,k in zip(planted,need,strict=True)]
    for i, activity in enumerate((x,x,xi,y,y,y)):
        exponent = prior.escalar(m[i],full) + prior.escalar(h[i][3:],factors)
        cap = full[i] if i < 3 else planted[i-3]
        assert activity*exp_upper(exponent) < cap
        if i >= 3:
            assert planted[i-3]*prior.escalar(h[i],full)/need[i-3] < full[i]
    rv = max(full[0]+full[3]+17*full[5],full[1]+full[4]+18*full[5])
    rm = full[2]
    assert 135*(vv*rv+vc*rm) < av and cv*rv+cc*rm < am
    return rv, rm


def activity_counts():
    x, xi, y = map(F, ('1e-11', '2e-8', '2e-5'))
    full = tuple(map(F, ('1.055476e-11','1.058684e-11','2.006808e-8',
                         '1.133360e-12','1.083308e-12','2.596133e-11')))
    planted = tuple(map(F, ('2.110952e-5','2.117368e-5','5.074245e-5')))
    av, am = F(7), F('.0031')
    rv, rm = tree_rows(x, xi, y, full, planted, av, am)
    assert rv < F('4.80e-10') and rm < F('2.007e-8')
    # H <= (10/3)e+(5/8)d. Tras la comparación quedan factores
    # menores que 0.52 por dato y 0.90 por medida.
    price_v = odds(P)*exp_upper((av+F('.015'))/16)**16 / x
    price_m = odds(Q)*exp_upper(am+F('.015')) / xi
    assert price_v**3 < F('.52')**3*y**10
    assert price_m**8 < F('.90')**8*y**5
    assert F(38,3) <= F(10,3)+15*F(5,8)

    # Otra inclinación de las mismas filas cubre H <= 4e aun con d=0.
    x0, xi0, y0 = map(F, ('2.02e-11','5.55e-11','5.6e-5'))
    full0 = tuple(map(F, ('2.13e-11','2.14e-11','5.6e-11',
                          '3.7e-13','3.6e-13','7.4e-12')))
    planted0 = tuple(map(F, ('5.9e-5','5.92e-5','1.325e-4')))
    rv0, rm0 = tree_rows(x0, xi0, y0, full0, planted0, 1, 1)
    tilt0 = exp_upper(F('1.015'))
    assert odds(P)*tilt0 < F('.015')*x0*y0**4
    assert odds(Q)*tilt0 < xi0
    # Se permite sumar ambas cotas: el solapamiento sólo las aumenta.
    return rv*F('.52')**3 + rv0*F('.015')**3, rm*F('.52')**3 + rm0*F('.015')**3


def quantitative_bounds():
    first_difference()
    data_low, measure_low = two_data()
    data_high, measure_high = activity_counts()
    ref = F('1.907e-6')*odds(P)/odds(F('1e-6'))
    assert ref < F('1.907e-30')
    bd, bm = data_low+data_high+ref, measure_low+measure_high+ref
    assert bd < F('1.94e-5') and bm < F('.00070')
    assert F('.00070')+F('.004298') == F('.004998')
    print(f'Mayorantes en la meta (valores aproximados): dato {float(bd):.7g}, '
          f'medida {float(bm):.7g}. Cotas racionales: 1.94e-5 y 0.00070.')


def physical_cycle(k):
    # En la copia activa, el tamaño de las celdas alterna 9 y 8.
    start = 0
    cells = []
    for j in range(2*k):
        size = 9 if j % 2 == 0 else 8
        cells.append(((1 << size)-1) << start)
        start += size
    windows = [cell | cells[(j+1) % len(cells)] for j,cell in enumerate(cells)]
    assert all(w.bit_count() == 17 for w in windows)
    assert all(sum((w >> i) & 1 for w in windows) == 2 for i in range(start))
    # Cada check tiene otras 16 columnas, cada una con la misma coordenada
    # en 17 copias distintas. Usar duplicados como columnas independientes.
    on_copy = (1 << start)-1
    full_cols = [w << (c*start) for c in range(17) for w in windows]
    for i in range(start):
        filler = sum(1 << (c*start+i) for c in range(17))
        full_cols.extend([filler]*16)
    assert all(w.bit_count() == 17 for w in full_cols)
    # Para una entrada en la copia activa, las columnas de relleno leen
    # a lo sumo una luz y las otras copias ninguna.
    assert all((w & on_copy).bit_count() <= 1 for w in full_cols[len(windows):])
    # Auditoría de todos los grados de la matriz, acumulando incidencias.
    row_degrees = [0]*(17*start)
    for w in full_cols:
        while w:
            bit = w & -w
            row_degrees[bit.bit_length()-1] += 1
            w ^= bit
    assert min(row_degrees) == max(row_degrees) == 18
    return windows


def winner(syndrome, windows):
    hits = [(w & syndrome).bit_count() for w in windows]
    h = max(hits)
    if 100*(2*h-17) < 17:
        return None
    j = hits.index(h)
    return (h, -j), j


def run(syndrome, windows):
    steps = []
    while (choice := winner(syndrome, windows)) is not None:
        j = choice[1]
        steps.append(j)
        syndrome ^= windows[j]
    return steps, syndrome


def compare(s, t, windows):
    sides = [s,t]
    traces = [[],[]]
    largest = (s ^ t).bit_count()
    while True:
        choices = [winner(state, windows) for state in sides]
        if choices == [None,None]:
            return traces, sides, largest
        if all(choices) and choices[0][1] == choices[1][1]:
            sides_to_update = (0,1)
        else:
            firing = max((i for i in (0,1) if choices[i] is not None),
                         key=lambda i: choices[i][0])
            sides_to_update = (firing,)
            a, b = sides[firing], sides[1-firing]
            w = windows[choices[firing][1]]
            u, v, common = ((w & mask).bit_count() for mask in (a & ~b,b & ~a,a & b))
            delta = 2*(w & a).bit_count()-w.bit_count()
            assert u > v and ((a ^ w) ^ b).bit_count()-(a ^ b).bit_count() == 2*(common-v)-delta
        for i in sides_to_update:
            j = choices[i][1]
            traces[i].append(j)
            sides[i] ^= windows[j]
        largest = max(largest,(sides[0] ^ sides[1]).bit_count())


def cycles():
    tests = 0
    for k in (3,4,5,6,7,8,12,32):
        windows = physical_cycle(k)
        all_data = (1 << k)-1
        midway = min(range(1, 2*k, 2), key=lambda j: abs(j-k))
        subsets = range(1 << k) if k <= 8 else (0,1,3,all_data >> 1, all_data)
        coefficients = [0,0]
        for subset in subsets:
            tests += 1
            selected = [2*j for j in range(k) if subset >> j & 1]
            initial = sum(windows[j] for j in selected)
            base_steps, end = run(initial, windows)
            assert end == 0 and base_steps == selected
            disturbed, end = run(initial ^ 1, windows)
            if subset == all_data:
                assert disturbed == list(range(1,2*k,2)) and end == 1
            passed = int(end == 1)
            for mask, c in enumerate((1,-1)):
                marked = passed * (-1 if mask and midway in disturbed else 1)
                predicted = 1-int(subset & 3 == 3)+c*int(subset == all_data)
                assert marked == predicted
                coefficients[mask] += (-1)**(k-subset.bit_count())*marked
            if k <= 8 or subset == all_data:
                tr, final_states, largest = compare(initial, initial ^ 1, windows)
                assert tr == [base_steps,disturbed] and final_states == [0,end]
                if subset == all_data:
                    assert largest >= 16
        if k <= 8:
            assert coefficients == [1,-1]
    print(f'Ciclos físicos: {tests} parejas de entradas. Fórmula de ambas máscaras '
          'y contraste de la comparación asíncrona correctos.')


if __name__ == '__main__':
    quantitative_bounds()
    cycles()
    print('PASS. El resto con e >= 3 y d < 15 e sigue abierto en p=1e-30.')
