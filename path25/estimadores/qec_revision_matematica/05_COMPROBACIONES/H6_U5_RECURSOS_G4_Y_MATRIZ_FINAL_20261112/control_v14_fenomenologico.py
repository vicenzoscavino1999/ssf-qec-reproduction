#!/usr/bin/env python3
"""H6-U5, parte C. Control complementario de V1.4 en el contrato
fenomenologico de PATH25: para todas las historias (E,D) de peso 1 a 4
sobre los 37 sitios, ejecucion global frente a restriccion y recomposicion
canonica por componentes del host K, y recuento de cancelaciones entre
componentes que han ejecutado acciones. Usa las funciones del verificador
sellado de la unidad 2. Salida determinista.

Uso: python -I -B control_v14_fenomenologico.py --output <nuevo.json>
"""
import argparse
import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as F
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
U2 = os.path.join(os.path.dirname(HERE), "H6_U2_CIRCUITO_PATH25_MULTICOMPONENTE_20261112", "verificar_h6_u2.py")
U2_SHA = "66ab5ab9e275c2a3daca1246b128a531ed9d89ee7d40fdeb274210d085da28d0"
N, M, U0 = 25, 12, 4369
CHECKS = []


class Stop(Exception):
    pass


def check(name, cond, detail=None):
    row = {"name": name, "pass": bool(cond)}
    if detail is not None:
        row["detail"] = detail
    CHECKS.append(row)
    if not cond:
        raise Stop(name)


def need(cond, name, **detail):
    if not cond:
        CHECKS.append({"name": name, "pass": False, "detail": detail})
        raise Stop(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    if os.path.exists(args.output):
        raise SystemExit("La salida ya existe; use un nombre nuevo.")
    out = {"unit": "H6-U5", "programa": "control V1.4 fenomenologico"}
    status = "PASS"
    try:
        run(out)
    except Stop as ex:
        status = "FAIL"
        out["stopped_at"] = str(ex)
    out["status"] = status
    out["checks"] = CHECKS
    with open(args.output, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(out, ensure_ascii=False, indent=1, sort_keys=True) + "\n")
    print("%s control V1.4 controles=%d fallidos=%d" % (status, len(CHECKS), sum(1 for c in CHECKS if not c["pass"])))
    if status != "PASS":
        print("DETENIDO EN:", out.get("stopped_at"))
        sys.exit(1)


def run(out):
    check("P0 verificador de la unidad 2 sin cambios", hashlib.sha256(open(U2, "rb").read()).hexdigest() == U2_SHA)
    spec = importlib.util.spec_from_file_location("u2_sellado", U2)
    u2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(u2)
    hx, hz, cols, actions, foot = u2.build_masks()
    hxs, hzs, colsets, actionsets, footsets, enc = u2.build_sets()
    T = u2.transitions_bits(foot)
    COR, END, TR, TOUCH = u2.resolve(T, actions)
    adjK = u2.host_K(actionsets, footsets, colsets)
    bz = u2.rank_basis(hz)
    psi = cols + [1 << c for c in range(M)]
    DATA = (1 << N) - 1
    wid = [bin(p).count("1") for p in foot]

    def score(s, j):
        return F(2 * bin(s & foot[j]).count("1") - wid[j], wid[j])

    by_size = {}
    cancel_examples = []
    target = None
    for size in range(1, 5):
        st = {"historias": 0, "por_componentes": {}, "varias_componentes": 0, "dos_o_mas_trazas_activas": 0,
              "cancelacion_tras_acciones": 0, "pasan_filtro_fisico": 0, "fallo_silencioso": 0}
        for S in combinations(range(N + M), size):
            sm = s = 0
            for i in S:
                sm |= 1 << i
                s ^= psi[i]
            E = sm & DATA
            tr = TR[s]
            Y = sm | TOUCH[s]
            R = E ^ COR[s]
            parts = u2.connected_parts(Y, adjK)
            hall = xr = 0
            assigned = 0
            local = []
            for P in parts:
                sp = P & sm
                need(sp != 0, "componente sin fuentes", S=S, P=P)
                sP = 0
                for i in u2.bits(sp):
                    sP ^= psi[i]
                exp = [j for j in tr if actions[j] & P]
                need(all(actions[j] & ~P == 0 for j in exp), "accion que cruza componentes", S=S, P=P)
                need((sp | TOUCH[sP]) == P, "restriccion: el soporte aislado no es la componente", S=S, P=P)
                need(list(TR[sP]) == exp, "restriccion: la traza local no es la subsucesion global", S=S, P=P)
                rP = (sp & DATA) ^ COR[sP]
                hP = u2.lin(cols, rP)
                need(hP & hall == 0, "huellas fisicas de dos componentes se cortan", S=S, P=P)
                hall |= hP
                xr ^= rP
                assigned += len(exp)
                local.append((sP, rP, hP, len(exp)))
            need(assigned == len(tr), "una accion no pertenece a exactamente una componente", S=S)
            need(xr == R, "el residuo no es el XOR de los residuos locales", S=S)
            need(hall == u2.lin(cols, R), "el sindrome fisico no es la union de los locales", S=S)
            ss = [x[0] for x in local]
            merged = []
            while True:
                best = None
                for hidx, sl in enumerate(ss):
                    j, nxt = T[sl]
                    if j >= 0:
                        kk = (score(sl, j), -j)
                        if best is None or kk > best[0]:
                            best = (kk, hidx, j, nxt)
                if best is None:
                    break
                merged.append(best[2])
                ss[best[1]] = best[3]
                need(len(merged) <= 12, "la recomposicion no termina", S=S)
            need(merged == list(tr), "recomposicion canonica distinta de la traza global", S=S)
            zero = hall == 0
            in_s = zero and u2.in_span(R, bz)
            bad_local = sum(1 for x in local if x[3] and x[2] == 0 and bin(U0 & x[1]).count("1") & 1)
            cancel = bad_local >= 2 and in_s
            st["historias"] += 1
            nc = str(len(parts))
            st["por_componentes"][nc] = st["por_componentes"].get(nc, 0) + 1
            if len(parts) > 1:
                st["varias_componentes"] += 1
            if sum(1 for x in local if x[3]) >= 2:
                st["dos_o_mas_trazas_activas"] += 1
            if zero:
                st["pasan_filtro_fisico"] += 1
                if not in_s:
                    st["fallo_silencioso"] += 1
            if cancel:
                st["cancelacion_tras_acciones"] += 1
                if len(cancel_examples) < 6:
                    cancel_examples.append({"fuentes": list(S), "traza": [actions[j] for j in tr],
                                            "residuos_locales": sorted(x[1] for x in local), "residuo": R})
            if S == (1, 2, 9, 10):
                target = {"traza": [actions[j] for j in tr], "componentes": sorted(parts),
                          "residuos_locales": sorted(x[1] for x in local), "residuo": R, "fisico": hall,
                          "en_S": in_s,
                          "locales_fuera_de_S_con_fisico_nulo": [x[2] == 0 and not u2.in_span(x[1], bz) for x in local],
                          "cancelacion": cancel}
        by_size[str(size)] = st
    low = {"1": 0, "2": 0, "3": 0}
    act = 0
    for size in ("1", "2", "3"):
        for k, v in by_size[size]["por_componentes"].items():
            low[k] = low.get(k, 0) + v
        act += by_size[size]["dos_o_mas_trazas_activas"]
    check("C0 restriccion, trazas locales, XOR de residuos, huellas fisicas disjuntas y recomposicion canonica en "
          "todas las historias de peso 1 a 4", True, {"historias": sum(v["historias"] for v in by_size.values())})
    check("C1 peso <=3: 4705 historias con una componente, 3234 con dos, 534 con tres; 3131 con dos o mas trazas "
          "activas (cifras publicadas)", low == {"1": 4705, "2": 3234, "3": 534} and act == 3131,
          {"por_componentes": low, "dos_o_mas_trazas_activas": act})
    check("C2 las 66045 historias de peso 4 pasan todas las comprobaciones", by_size["4"]["historias"] == 66045)
    check("C3 historia E={1,2,9,10}: traza 0,3,8,11; componentes {0,1,2,3} y {8,9,10,11}; residuos locales 15 y "
          "3840 fuera de S_X con fisico nulo; residuo 3855 en S_X; cancelacion",
          target is not None and target["traza"] == [1, 8, 256, 2048] and target["componentes"] == [15, 3840]
          and target["residuos_locales"] == [15, 3840] and target["residuo"] == 3855 and target["fisico"] == 0
          and target["en_S"] and all(target["locales_fuera_de_S_con_fisico_nulo"]) and target["cancelacion"], target)
    check("C4 existe al menos una cancelacion tras acciones en peso 4", by_size["4"]["cancelacion_tras_acciones"] >= 1,
          by_size["4"]["cancelacion_tras_acciones"])
    out["por_peso"] = by_size
    out["ejemplos_de_cancelacion"] = cancel_examples
    out["historia_predicha"] = target
    out["alcance"] = ("PATH25, contrato fenomenologico, historias de peso 1 a 4 (74518); contrato generico no "
                      "regular; no es un estimador")


if __name__ == "__main__":
    main()
