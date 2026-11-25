"""Certificado uniforme (17,18), beta=1/100. Python >= 3.10, biblioteca estándar.

Ejecutar: python verificar_p_general.py
La demostración que permite usar estas desigualdades está en H5A_p_general.txt.
Todas las decisiones usan fracciones o enteros; los flotantes solo se imprimen.
"""
from fractions import Fraction as F
from math import ceil


def exigir(condicion, motivo):
    if not condicion:
        raise ArithmeticError(motivo)


# Cada fila tiene un umbral y su multiplicador dual. La primera, sin umbral,
# corresponde al descenso total. Los dos certificados se usan por separado.
PESOS_PASOS = {
    None: F(1, 17), F(1): F(1, 17), F(8, 9): F(19, 272),
    F(7, 9): F(19, 240), F(2, 3): F(19, 210), F(5, 9): F(5, 182),
    F(9, 17): F(1, 13), F(4, 9): F(1, 26), F(7, 17): F(1, 12),
    F(1, 3): F(7, 132), F(5, 17): F(1, 11), F(2, 9): F(19, 110),
    F(1, 9): F(1, 10),
}
PESOS_SITIOS = {
    None: F(2, 17), F(1): F(100, 153), F(8, 9): F(53, 72),
    F(7, 9): F(15301, 18360), F(2, 3): F(2188, 2295),
    F(5, 9): F(2186, 1989), F(4, 9): F(5105, 3978),
    F(1, 3): F(5101, 3366), F(2, 9): F(15314, 8415),
    F(1, 9): F(169, 85),
}


def verificar_perfil():
    # Para cada huella, el mayor soporte entre TODAS las formas admisibles.
    max_sitios = {}
    for x in range(19):
        for y in range(18):
            sitios = x + y
            huella = 17*x + 18*y - 2*x*y
            if sitios and 2*huella >= 17*sitios:
                max_sitios[huella] = max(max_sitios.get(huella, 0), sitios)
                if sitios >= 2:
                    exigir(huella >= 33, "huella de acción múltiple")
    exigir(min(max_sitios) == 17, "huella mínima")
    cotas = []
    for num, pesos in enumerate((PESOS_PASOS, PESOS_SITIOS)):
        exigir(all(v > 0 for v in pesos.values()), "positividad del dual")
        controles = [(t, v) for t, v in pesos.items() if t is not None]
        exigir(all(F(1, 100) < t <= 1 for t, _ in controles), "niveles del perfil")
        cota = 18*pesos[None] + sum(v*(ceil(9*(1+t))-1) for t, v in controles)
        for huella, sitios in max_sitios.items():
            for descenso in range(huella % 2, huella+1, 2):
                if 100*descenso < huella:
                    continue
                creaciones = (huella-descenso)//2
                precio = pesos[None]*descenso + sum(
                    v*(descenso if descenso*t.denominator < huella*t.numerator
                       else -creaciones) for t, v in controles
                )
                objetivo = 1 if num == 0 else sitios
                exigir(precio >= objetivo, f"dual {num}, huella {huella}, descenso {descenso}")
        cotas.append(cota)
    exigir(cotas[0] < F(38, 3) and 1 + cotas[1] < 135, "presupuestos por fuente")
    return cotas


def exp_superior(t):
    exigir(0 <= t < 22, "dominio del resto de Taylor")
    termino = suma = F(1)
    for j in range(1, 21):
        termino *= t / j
        suma += termino
    return suma + termino*t / (21 * (1 - t/22))


def producto(fila, columna):
    return sum(a*b for a, b in zip(fila, columna, strict=True))


def verificar_arboles():
    A, B = 17, 18
    ramificacion = (A-1)*(B-1)
    v2, c2 = 1+A*(B-1), 1+B*(A-1)
    v3, c3 = B*(1+ramificacion), A*(1+ramificacion)
    v4, c4 = 1+(v2-1)*(1+ramificacion), 1+(c2-1)*(1+ramificacion)
    matriz_base = [
        [v2*v4, c3*c3, c3*v4, v2*c3],
        [v3*v3, c4*c2, c4*v3, v3*c2],
        [v3*v4, c4*c3, c4*v4, v3*c3],
        [v2*v3, c3*c2, c3*v3, v2*c2],
    ]
    # Orden: fuentes VV, CC, medida; acciones singleton VV, singleton CC, múltiple.
    geometria = (0, 1, 3, 0, 1, 2)
    anchos = (A, B, 1, A, B, A*B)
    lectores = (B, A, 1, B, A, A*B)
    H = [[ancho*n for n in lectores] for ancho in anchos]
    M = [[matriz_base[i][j] for j in geometria] for i in geometria]
    k = (9, 10, 17)
    x, x_d, y = F("2.02e-11"), F("5.55e-11"), F("5.6e-5")
    completos = list(map(F, ("2.13e-11", "2.14e-11", "5.6e-11",
                             "3.7e-13", "3.6e-13", "7.4e-12")))
    plantados = list(map(F, ("5.9e-5", "5.92e-5", "1.325e-4")))
    factores_hijos = [z/divisor for z, divisor in zip(plantados, k, strict=True)]
    actividades = (x, x, x_d, y, y, y)
    for i in range(6):
        exponente = producto(H[i][3:], factores_hijos) + producto(M[i], completos)
        cubierta = completos[i] if i < 3 else plantados[i-3]
        exigir(actividades[i]*exp_superior(exponente) < cubierta, f"ramas ordinarias {i}")
        if i >= 3:
            causa = plantados[i-3] * producto(H[i], completos) / k[i-3]
            exigir(causa < completos[i], f"rama con causa distinguida {i}")

    r_v = max(completos[0]+completos[3]+A*completos[5],
              completos[1]+completos[4]+B*completos[5])
    r_c = completos[2]
    vecinos_v = max(sum(matriz_base[0][:2]), sum(matriz_base[1][:2]))
    vecinos_c = max(matriz_base[0][3], matriz_base[1][3])
    k_e = 135*(vecinos_v*r_v + vecinos_c*r_c)
    k_d = sum(matriz_base[3][:2])*r_v + matriz_base[3][3]*r_c
    exigir(0 < r_c < 1 and 0 < r_v < 1 and 0 < k_e < 1 and 0 < k_d < 1, "cargas KP")

    p, q = F(1, 10**65), F(2, 10**11)
    zp, zq = p/(1-p), q/(1-q)
    inclinacion = F(11, 4)
    exigir(exp_superior(F(101, 100)) < inclinacion, "inclinación")
    exigir(0 < y < 1 and (inclinacion*zp)**3 < x**3 * y**38, "tasa p")
    exigir(inclinacion*zq < x_d, "tasa q")
    exigir(p/F("5e-76") == 20000000000 and q == F("2e-11"), "comparación conjunta")

    # Techo para el MISMO sistema exponencial con el precio y^(38/3).
    # exp(1) > 2.7, por los primeros seis términos de su serie.
    exigir(sum(F(1, j) for j in (1, 1, 2, 6, 24, 120)) > F(27, 10), "cota inferior de e")
    exigir(F(17, 1)/(H[5][5]*F(27, 10)) < F("7e-5"), "techo de la rama múltiple")
    exigir(F("7e-5")**38 < (F("4e-63")*135*vecinos_v)**3, "techo de esta parametrización")
    return r_v, r_c, k_e, k_d


if __name__ == "__main__":
    pasos, sitios = verificar_perfil()
    rv, rc, ke, kd = verificar_arboles()
    print(f"Pasos por dato <= {pasos}; sitios volteados con repetición por dato <= {sitios}")
    print(f"R_V={rv}, R_C={rc}, K_E={ke} ({float(ke):.12g}), K_D={kd} ({float(kd):.12g})")
    print("PASS. p <= 1e-65, q <= 2e-11. Mejora en p: factor 20000000000.")
