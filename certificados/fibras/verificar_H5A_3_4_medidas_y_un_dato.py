"""H5-A (3,4), q <= 1e-3. Certificado reproducible sin dependencias externas."""

from fractions import Fraction as Q


def exp_upper(x, degree=24):
    assert 0 <= x < degree + 2
    total = term = Q(1)
    for k in range(1, degree + 1):
        term *= x/k
        total += term
    return total + term*x/(degree+1)/(1-x/(degree+2))


def local_profiles():
    rows = [7 << (3*i) for i in range(4)]
    cols = [sum(1 << (3*i+j) for i in range(4)) for j in range(3)]
    stripes = rows + cols
    # Construcción independiente en el producto de las dos bolas del árbol.
    # Las primeras 4 y 3 posiciones son comunes; las otras 6 por rama son privadas.
    cell_balls = []
    for i in range(4):
        for j in range(3):
            left = [*range(4), *range(4+6*i, 10+6*i)]
            right = [*range(3), *range(3+6*j, 9+6*j)]
            cell_balls.append(sum(1 << (21*r+c) for r in left for c in right))
    sizes = [0] * 4096
    for mask in range(1, 4096):
        cover = 0
        for i, ball in enumerate(cell_balls):
            if mask & (1 << i):
                cover |= ball
        sizes[mask] = cover.bit_count()
        formula = (12 + 18*sum(bool(mask & r) for r in rows)
                   + 24*sum(bool(mask & c) for c in cols) + 36*mask.bit_count())
        assert sizes[mask] == formula

    candidates = []
    for mask in range(1, 128):
        footprint = raw = 0
        for j, stripe in enumerate(stripes):
            if mask & (1 << j):
                footprint ^= stripe
                raw |= stripe
        size, sites = footprint.bit_count(), mask.bit_count()
        if size and 2*size >= 3*sites:
            candidates.append((footprint, sites, size, raw))
            if size <= 4:
                assert sites == 1 and footprint in (rows if size == 3 else cols)
    assert len(candidates) == 101
    signatures = set()
    for state in range(4096):
        best_h, best_v = -1, 1
        winners = []
        for footprint, sites, size, raw in candidates:
            lit = footprint & state
            h = lit.bit_count()
            cmp = h*best_v - best_h*size
            if cmp > 0:
                best_h, best_v, winners = h, size, []
            if cmp >= 0:
                winners.append((sites, size, h, sizes[raw]-sizes[lit]))
        if 100*(2*best_h-best_v) >= best_v:
            signatures.update(winners)
    return signatures


def budgets(signatures, cost):
    items = [(size, h, cost(f, rho)) for f, size, h, rho in signatures]
    K = max(Q(b, 2*h-size) for size, h, b in items)
    return tuple(tuple(max((K*(size-h)+b)/h for size, h, b in items if n*h >= k*size)
                       for k in range(1, n+1)) for n in (3, 4))


def first_extinctions(signatures):
    steps = budgets(signatures, lambda f, rho: 1)
    sites = budgets(signatures, lambda f, rho: f)
    increment = budgets(signatures, lambda f, rho: rho)
    assert tuple(map(sum, steps)) == (Q(7, 3), 3)
    assert tuple(1+sum(row) for row in sites) == (Q(17, 3), Q(103, 15))
    assert increment == ((150, 144, 24), (150, 150, Q(348, 5), 24))
    return tuple(initial+sum(row) for initial, row in zip((210, 252), increment))


def geometric_slots(y):
    # Cada bola con l generadores omitidos tiene tamaño A*6**l + B.
    even_v = (Q(54, 5), -Q(4, 5))
    even_c = (Q(48, 5), -Q(3, 5))
    odd_v = (Q(18, 5), -Q(3, 5))  # desde una variable
    odd_c = (Q(24, 5), -Q(4, 5))  # desde un check

    def mass(first, second):
        a, b = first
        c, d = second
        return (a*c*(1-y)/(1-36*y)
                + (a*d+b*c)*(1-y)/(1-6*y) + b*d)

    assert 0 < 36*y < 1
    caps = (55, 65, 73, 49, 155)
    actual = (mass(odd_v, even_v), mass(even_c, odd_c),
              mass(odd_c, even_v), mass(even_c, odd_v),
              mass(even_c, even_v)-1)
    assert all(0 < j < cap for j, cap in zip(actual, caps))
    return caps


def physical_rectangle(coverage):
    p, q, tau = map(Q, ("1e-12", ".001", ".01"))
    U, P, y, s = map(Q, (".00004", ".00135", ".011", ".00130"))
    x1, x2 = map(Q, ("2.31e-7", "6.85e-6"))
    u1, u2, ug, v, c = map(Q, ("2.476e-7", "7.423e-6", "7.255e-6",
                                  ".0007319", "1.1023e-6"))
    zp, zq = p/(1-p), q/(1-q)
    B, BD = 1297*U, 469*U+90*P
    D1, D2 = coverage
    assert (D1, D2) == (528, Q(3228, 5))
    tests = {
        "medidas pasivas": zq*exp_upper(tau+2*BD)/P,
        "medidas activas": zq*exp_upper(tau+180*P)/s,
        "fuente VV (cubo)": (zp*exp_upper(tau+Q(17, 3)*B+2*P*D1)/x1)**3/y**7,
        "fuente CC": zp*exp_upper(tau+Q(103, 15)*B+2*P*D2)/(y**3*x2),
    }
    j1, j2, k1, k2, jg = geometric_slots(y)
    W1 = 640*u1+441*u2+1344*ug
    W2 = 784*u1+513*u2+1596*ug
    Wg = 1792*u1+1197*u2+3648*ug
    L = jg*v
    G = (1+s)**12*exp_upper(Wg+L)
    H = k1*u1+k2*u2+jg*c
    tests.update({
        "raíz VV": x1*(1+s)**3*exp_upper(W1+j1*v)/u1,
        "raíz CC": x2*(1+s)**4*exp_upper(W2+j2*v)/u2,
        "generador alcanzado desde su causa": y*(G-1-L)/v,
        "generador alcanzado desde un hijo": y*H*(G-1)/c,
        "generador inicial o de puente": y*H*G/ug,
        "ancla de datos": max(u1+3*ug, u2+4*ug)/U,
    })
    assert max(tests.values()) < Q(".995"), [(key, float(val)) for key, val in tests.items()]
    assert max(U, P) < 1 and max(Q(17, 3), Q(103, 15), 1) < 7
    print(f"Diez desigualdades racionales, cociente máximo {float(max(tests.values())):.7f}.")
    print("Rectángulo general p <= 1e-12, q <= 1e-3. Mejora de p por un factor de 250.")


def check_one_data(signatures):
    # Un primer ganador compatible con el saldo de una entrada aceptada debe
    # agotar ese saldo. Una VV perfecta frente a un único dato CC es imposible
    # porque en la huella VV hay a lo sumo una luz de D y una del dato.
    for n in (3, 4):
        for t in range(n//2+1):
            balance = n-2*t
            admissible = set()
            for f, size, h, _ in signatures:
                delta = 2*h-size
                if delta > balance or delta*n < balance*size:
                    continue
                if n == 4 and size == 3 and h > 2:
                    continue
                admissible.add((f, size, delta))
            if balance:
                assert admissible == {(1, n, balance)}
            else:
                assert not admissible


def source_rows(p, v_limit, m_limit):
    q = Q(".001")
    tilt = exp_upper(Q(".015"))
    b, m = p/(1-2*p)*tilt, q*tilt
    M = Q(".00128")
    assert m*exp_upper(178*M) < M and 178*M < 1
    T = b*exp_upper(504*M)
    one_at_data = b*(exp_upper(504*M)-1)
    one_at_measurement = 938*M*T/(1-178*M)
    pure_measurements = m*(exp_upper(178*M)-1)
    two_data = 2*1296*b*b
    assert one_at_data+two_data < v_limit
    assert pure_measurements+one_at_measurement < m_limit
    print(f"p={float(p):g}, q=0.001: filas parciales < {float(v_limit):g} por dato "
          f"y < {float(m_limit):g} por medida.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para conservar las comprobaciones.")
    local = local_profiles()
    physical_rectangle(first_extinctions(local))
    check_one_data(local)
    source_rows(Q(".001"), Q(".003604"), Q(".003275"))
    source_rows(Q("1e-6"), Q("9.23e-7"), Q(".000263"))
    print("Verificadas las 4096 configuraciones locales y las cotas de todos los órdenes indicados.")
