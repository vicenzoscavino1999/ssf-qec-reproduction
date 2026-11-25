"""Cota uniforme del resto físico a p=4.24e-36 y q=2e-11. Python 3.10+.

Usa el certificado H5A_CUMULANTES_Y_COMPONENTES_20261117.zip del proyecto,
en la misma carpeta. Todas las desigualdades decisorias usan fracciones.
"""

from fractions import Fraction as F
from math import ceil, factorial, floor
from pathlib import Path
import sys


archive = Path(__file__).with_name('H5A_CUMULANTES_Y_COMPONENTES_20261117.zip')
if not archive.is_file():
    raise FileNotFoundError(archive)
sys.path.insert(0, str(archive))
import verificar_H5A_cumulantes as prior


TOP = 17
SCALE = 10**55
ZERO = (F(0),)*(TOP+1)
T, X, XI = F('.00665'), F('1.60e-8'), F('2.07e-11')
P, Q = F('4.24e-36'), F('2e-11')
Z_CAP = tuple(map(F, ('.00002', '.00002', '.000002')))
SOURCE_AND_ACTION = tuple(map(F, (
    '3.036918191e-20', '1.677545077e-19', '8.353195291e-24',
    '2.014322932e-18', '1.91953247e-18', '1.921077002e-17')))
PLANTED = tuple(map(F, ('5.686882931e-12', '5.395221995e-12', '5.773761201e-11')))
N3, N4, N5, NMAX = 343944, 3683016, 146395296, 306*2**34

# Supersolución para cada grado de datos del grupo. El último valor agrega
# todos los grados >=17. Se verifica cada coordenada contra su recurrencia.
U = [tuple(map(F, row.split())) for row in (
    '0 2.19275468e-18 8.20378007e-13 9.21774031e-13 4.87715245e-13 '
    '1.66298455e-13 4.20418014e-14 8.73076968e-15 1.68366295e-15 '
    '3.51931378e-16 8.63544255e-17 2.31680690e-17 6.17285844e-18 '
    '1.57779067e-18 3.91781576e-19 9.73287305e-20 2.47232993e-20 8.73511707e-21',
    '0 1.33282157e-17 3.34717209e-12 4.38694345e-12 2.52202491e-12 '
    '9.05927409e-13 2.37182117e-13 5.03430140e-14 9.77321912e-15 '
    '2.02524637e-15 4.91618729e-16 1.31707150e-16 3.52696751e-17 '
    '9.06745566e-18 2.25920184e-18 5.61471071e-19 1.42456140e-19 5.02986647e-20',
    '0 3.10418943e-16 2.58035813e-16 1.81626072e-16 8.97271097e-17 '
    '3.08162181e-17 7.95219291e-18 1.68550326e-18 3.29673527e-19 '
    '6.90629457e-20 1.68454766e-20 4.50290030e-21 1.20140648e-21 '
    '3.08297107e-22 7.68454789e-23 1.91301059e-23 4.86102969e-24 1.71715088e-24')]


def up(x):
    return F(ceil(SCALE*x), SCALE)


def constant(x):
    return (F(x),) + ZERO[1:]


def add(*rows):
    return tuple(map(sum, zip(*rows, strict=True)))


def times(x, a):
    return tuple(up(x*v) for v in a)


def mul(a, b):
    out = list(ZERO)
    for i, v in enumerate(a):
        if v:
            for j, w in enumerate(b):
                out[min(i+j, TOP)] += v*w
    return tuple(map(up, out))


def positive_series(a, *, geometric=False):
    # Escribir a=a0+b. Cada término de b tiene grado al menos uno.
    b = (F(0), *a[1:])
    if geometric:
        assert a[0] < 1
        factor = 1/(1-a[0])
        b = times(factor, b)
    else:
        factor = up(prior.exp_arriba(a[0]))
    mass = sum(b)
    term = total = constant(1)
    for j in range(1, TOP):
        term = mul(term, b)
        if not geometric:
            term = times(F(1, j), term)
        total = add(total, term)
    if geometric:
        assert mass < 1
        tail = mass**TOP/(1-mass)
    else:
        assert mass < TOP+1
        tail = mass**TOP/factorial(TOP)/(1-mass/(TOP+1))
    total = (*total[:-1], total[-1]+up(tail))
    return times(factor, total)


def linear(weights, rows):
    return add(*(times(x, a) for x, a in zip(weights, rows, strict=True)))


def group_rows(n):
    data = (F(0), F(1), *ZERO[2:])
    source = (times(X*F('.16452'), data), times(X, data), constant(XI))
    a = [mul(s, positive_series(add(constant(z), linear(row, U))))
         for s, z, row in zip(source, Z_CAP, n, strict=True)]
    v = [mul(a[i], positive_series(times(n[i][2], a[2]))) for i in (0, 1)]
    w = [mul(vertex, positive_series(linear(row[:2], v)))
         for vertex, row in zip((*v, a[2]), n, strict=True)]
    readers = linear((18, 17), w[:2])
    y = add(readers, w[2])
    g = positive_series(y, geometric=True)
    cover = add(constant(N5), times(NMAX, mul(y, g)))
    cover = add(constant(N4), mul(y, cover))
    cover = mul(mul(y, y), add(constant(N3), mul(y, cover)))
    last = add(times(5, g), mul(y, mul(g, g)))
    derivative = add(constant(4*N5), times(NMAX, mul(y, last)))
    derivative = add(constant(3*N4), mul(y, derivative))
    derivative = mul(y, add(constant(2*N3), mul(y, derivative)))
    connections = ((10700, 17600), (17600, 10100))
    k = [mul(w[i], add(linear(connections[i], w[:2]), times(40, w[2]), times(18, cover)))
         for i in (0, 1)]
    k.append(mul(w[2], add(times(800, add(*w[:2])), mul(readers, derivative))))
    d = [add(k[i], mul(v[i], linear(n[i], k))) for i in (0, 1)]
    d.append(add(k[2], mul(a[2], linear(n[2][:2], d))))
    for claimed, required in zip(U, d, strict=True):
        assert len(claimed) == TOP+1
        assert all(c >= r for c, r in zip(claimed, required, strict=True))

    # Descuento por número de datos, uniforme sobre sus dos tipos.
    discounts = tuple(map(F, (
        '0', '0', '6.43e-13', '3.23e-9', '7.18e-10', '1.10e-9', '3.93e-9',
        '1.782e-7', '1.427e-6', '4.96e-7', '7.01e-7', '9.03e-6', '3.824e-5',
        '1.194e-4', '.001618', '.01071', '.1730', '1')))
    for e in range(2, TOP):
        delta = min(F('12.66')*e-F('.36')*i-prior.perfil_grupo(i, e-i) for i in range(e+1))
        assert delta >= 0 and discounts[e]**1000 >= T**floor(1000*delta)
    assert F('.16452')**25 >= T**9 and F('1.0754e-8')**50 >= T**183
    d0 = XI*prior.exp_arriba(Z_CAP[2])
    v0 = [X*F('1.0754e-8')*prior.exp_arriba(Z_CAP[i]+n[i][2]*d0) for i in (0, 1)]
    b0 = 40*d0+18*prior.cubiertas(d0)
    solo = [u*b0 for u in v0]
    solo.append(d0*(800*max(v0)+prior.escalar(n[2][:2], v0)*b0
                    +(18*v0[0]+17*v0[1])*prior.derivada(d0)))
    costs = [s+prior.escalar(discounts, u) for s, u in zip(solo, U, strict=True)]
    assert all(c < claimed for c, claimed in zip(costs, SOURCE_AND_ACTION[:3], strict=True))
    return costs


def verify():
    # Las filas de acción, los grados y el perfil cronológico son los del
    # certificado precedente. Sólo se sustituye la suma de grupos fuente.
    lam, ke, kd, rv, rm, _ = prior.general()
    base, _, _ = prior.geometria()
    order = (0, 1, 3, 0, 1, 2)
    m = [[base[i][j] for j in order] for i in order]
    n = [row[:3] for row in m[:3]]
    n[2][2] = 0
    for i, width in enumerate((17, 18, 1)):
        z = prior.escalar(m[i], SOURCE_AND_ACTION)
        z += width*prior.escalar((18, 17, 306), PLANTED)/(1-lam)
        assert z < Z_CAP[i]
    costs = group_rows(n)

    # s_i=r_i exp(.015)/(1-2r_i), más .01 para las decoraciones físicas.
    tilt = F('.025')
    vv, vc = max(sum(row[:2]) for row in base[:2]), max(row[3] for row in base[:2])
    cv, cc = sum(base[3][:2]), base[3][3]
    lv, lc = F('5e-36'), F('2.055e-11')
    s = SOURCE_AND_ACTION
    core = max(s[0]+s[3]+17*s[5], s[1]+s[4]+18*s[5])
    check = s[2]+35*core
    yp, yq = prior.actividad(P), prior.actividad(Q)
    exp = prior.exp_arriba
    assert lv >= yp*exp(tilt+vv*(lv+core)+vc*(lc+check))
    assert lc >= yq*exp(tilt+cv*(lv+core)+cc*(lc+check))
    assert XI >= yq*exp(tilt+cv*lv+cc*lc)
    assert exp(tilt+135*(vv*lv+vc*lc)) < F('1.03')
    assert (F('1.03')*yp)**50 < X**50*T**633
    assert 0 < max(ke, kd) < F('.01') and max(rv, rm) < F('.005')
    assert rv == core+18*lv and rm == check+lc
    assert rv < F('3.47882e-16') and rm < F('2.05622e-11')

    # A p=0 todos los bloques físicos de dos o más componentes de medidas
    # se anulan. En el residuo cada soporte que subsiste contiene un dato.
    ratio = yp/prior.actividad(F('1e-6'))
    reference = ratio*F('1.907e-6')
    assert reference < F('8.1e-36')
    # El cluster físico de una medida sola pesa s_m y tiene grado uno.
    # La linealización inferior exp(x)>=1+x hace conservadora la resta.
    residual_v = rv + reference
    residual_m = rm - yq*(1+F('.015')) + reference
    assert 0 < residual_v < F('3.48e-16')
    assert 0 < residual_m < F('2.63e-13') < F('.004998')
    print('Grupos físicos por tipo de raíz', *(f'{float(c):.8g}' for c in costs))
    print('Resto completo por dato < 3.48e-16 y por medida < 2.63e-13.')
    print(f'PASS. Rectángulo general p <= {float(P):.3g}, q <= {float(Q):.3g}.')


if __name__ == '__main__':
    verify()
