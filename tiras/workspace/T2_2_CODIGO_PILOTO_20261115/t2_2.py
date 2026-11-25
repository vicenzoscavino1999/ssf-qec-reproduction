#!/usr/bin/env python3
"""T2-2, piloto m = 6 (BORRADOR_PLAN_T2_2.txt v4): conductor del splitting de Bravyi y Vargo en la forma de Fail
fast para P_sil en TIRA_6x6, con r_m = r_c y las tasas 1e-3, 1e-4, 1e-5 y 1e-6.

Codigo preparado para revision (version 2: R1 a R3 del dictamen del codigo v1, DICTAMEN_CODIGO_T2_2_V1.txt).
NO SE HA EJECUTADO NINGUNA CORRIDA, ni siquiera parcial: solo las
autopruebas y las pruebas de humo de t2s_pruebas.py, con modelos de juguete y sin decoder, DP, ancla fisica ni
cadenas fisicas. Para correr hacen falta el plan sellado (SELLO_T2_2.json con los SHA-256 del plan, del codigo
y de ESCALERA_Y_CLAVES.json) y una orden del usuario; la corrida se lanza bajo t2s_supervisor.py.

Etapas, en el orden del plan (6), dentro de un solo techo global de reloj (8 h desde el inicio del supervisor,
--verificar incluido):
  preparacion (sello, dependencias, modelo, controles previos, escalera y claves recalculadas = selladas,
  referencias); ancla (<= 3e8 llamadas, <= 4 h); c_tiempo (<= 6e4 pasos, <= 15 min; proyeccion); cadenas
  (<= 4e8 llamadas, <= 3e9 pasos, fin a las 7 h 15); estimacion y diagnosticos (<= 15 min, fin a las 7 h 30);
  p_conex: NO_EJECUTADO (opcional en el plan; no se implementa en esta version del codigo).
Cada etapa queda PASS, FALLA (control fallido) o NO_EVALUADO (techo, error, cierre anormal), con el motivo y el
coste consumido. ANCLA_INCOMPLETA, ANCLA_ABORTADA (un bloque interrumpido por techo o error, o un techo del
proceso principal, aunque el prefijo ya tenga el exito s: R1), una proyeccion de C-tiempo que no cabe y
DATOS_NO_CABEN detienen la unidad (7.1). Nunca hay estimacion parcial: sin todas las etapas PASS no hay
decisiones. El plazo se consulta durante la estimacion y otra vez al terminarla, antes de aceptarla (R3).

La estimacion se recalcula SIEMPRE desde los registros escritos en disco (ancla/bloques.jsonl, cadenas/*.bin y
cadenas/registro_*.jsonl) con la misma funcion `evaluate` que usa --verificar. --verificar (R2) comprueba ademas
la configuracion (reconstruida desde el sello en el modo fisico), el inventario, las etapas, los controles
obligatorios, las decisiones de primer nivel y RECURSOS frente a RESULTADO y a los registros; su plazo cuenta
desde que empieza (R3). Devuelve 0 solo si la entrega esta completa, todas las etapas son PASS y todo coincide;
2 si es coherente pero incompleta (inspeccion de una corrida incompleta); 1 ante cualquier discrepancia o error;
4 si alcanza su techo (verificacion no concluida); 3 si no se ejecuta porque el informe ya existe o su destino
no vale (nunca sobrescribe: un informe posterior va con --salida a una ruta nueva fuera de la carpeta).

Uso (desde D:\\paper ia cuantica\\qec_revision_matematica\\05_COMPROBACIONES):
  python -I -B T2_2_CODIGO_PILOTO_20261115/t2s_supervisor.py --corrida --dir <carpeta nueva>
  python -I -B T2_2_CODIGO_PILOTO_20261115/t2_2.py --verificar <carpeta de la corrida> [--salida <nuevo json>]
  python -I -B T2_2_CODIGO_PILOTO_20261115/t2_2.py --preparar-sello --salida <nuevo json>   (tras la revision)
  (t2_2.py --corrida solo lo lanza el supervisor, con --t0.)
"""
import argparse
import json
import math
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from concurrent.futures.process import BrokenProcessPool
from decimal import Decimal
from fractions import Fraction as F

UDIR = "T2_2_CODIGO_PILOTO_20261115/"
sys.path.insert(0, os.path.abspath(UDIR))
sys.set_int_max_str_digits(0)
import t2s_mat as tm  # noqa: E402
import t2s_motor as tmo  # noqa: E402
import t2s_reg as tr  # noqa: E402

PLAN = "T2_2_BORRADOR_PLAN_SPLITTING_20261115/BORRADOR_PLAN_T2_2.txt"
SELLO = UDIR + "SELLO_T2_2.json"
ESCALERA = UDIR + "ESCALERA_Y_CLAVES.json"
CODE = ("t2_2.py", "t2s_mat.py", "t2s_reg.py", "t2s_motor.py", "t2s_fisica.py", "t2s_juguete.py",
        "t2s_supervisor.py")
CP_PATH = "T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_cp.py"
CP_SHA = "7e1a913b84ddb11941bcd8c09d26a7a37b23ced7173b51ccd0f9cd0c098f9565"
RES_T21 = "T2_1_CORRIDA_2_EJECUCION_20261115/corrida_2/RESULTADO_T2_1.json"
RES_T21_SHA = "9769c76eae7d54a5377b0397d3ad0e0ee6f68207cef8f3c70ec4971fe9f6df69"
RES_DMC = "T2_1_CORRIDA_1_EJECUCION_20261114/corrida_1/RESULTADO_DMC.json"
RES_DMC_SHA = "15c71efe0690a858310ec4953c417ab86182dc48d032a2c45a310faf96f03325"
GIB = 1 << 30
POLL = 5.0
MARGEN_ESCRITURA = 120.0     # s: el conductor termina antes del plazo de aceptacion de la fase 1 del supervisor
RESULT = "RESULTADO_T2_2.json"
RECURSOS = "RECURSOS_T2_2.json"
MANIF = "MANIFIESTO_T2_2.json"
POSTERIORES = ("CONSOLA_CONDUCTOR.txt", "SUPERVISION_T2_2.json", "VERIFICACION_T2_2.json",
               "CONSOLA_VERIFICAR.txt")
PLANNED = ("preparacion", "ancla", "c_tiempo", "cadenas", "estimacion")
ALCANCE_ANCLA = ("V6 describe el procedimiento inverso sin truncacion con historias Bernoulli(p) "
                 "independientes y un generador ideal: E[(s-1)/(T-1)] = p y Pr(l <= p <= u) >= 0.975. "
                 "No es una cobertura condicionada a terminar bajo el techo ni a pasar los filtros de "
                 "ampliacion, ni una insesgadez entre las corridas aceptadas; ANCLA_INCOMPLETA no tiene "
                 "intervalo (C3).")


class CheckFail(Exception):
    """Control fallido: la etapa queda FALLA y la unidad se detiene."""


CHECKS = []


def check(desc, ok, detail=None):
    CHECKS.append({"control": desc, "pasa": bool(ok), "detalle": detail})
    print("  [%s] %s" % ("PASS" if ok else "FALLA", desc), flush=True)
    if not ok:
        raise CheckFail(desc)


# ------------------------------------------------------------------ configuracion real (2.12, 6)
def real_config(seal_info, esc):
    return {
        "unidad": "T2-2 piloto TIRA_6x6", "tipo": "corrida",
        "modelo": {"tipo": "fisico", "r0": "1/100"},
        "semilla_maestra": tm.SEMILLA_MAESTRA,
        "escalera": {"tasas": esc["escalera"]["tasas"], "shas": esc["shas"],
                     "extremos": [tm.fstr(x) for x in tm.ENDPOINTS]},
        "ancla": {"s": tm.S_ANCHOR, "B_a": tm.B_ANCHOR, "semillas": tm.L_SEEDS, "alfa_por_cola": tm.HALF_ANCHOR,
                  "objetivo_cp": tm.TARGET_CP, "en_vuelo": 14},
        "cadenas": {"T_init": tm.T_INIT, "lambda_ini": tm.LAMBDA_INI, "ampliaciones_max": tm.MAX_EXT,
                    "eps": tm.fstr(tm.EPS), "B": tm.BLOCK, "L": tm.L_SEEDS, "M": tm.M_REPS},
        "c_tiempo": {"pasos": 10000, "repeticiones": 2},
        "procesos": 14,
        "techos": {"global_s": 8 * 3600.0, "verificar_s": 1800.0,
                   "ancla": {"llamadas": 300_000_000, "segundos": 4 * 3600.0},
                   "c_tiempo": {"llamadas": 60_000, "pasos": 60_000, "segundos": 900.0},
                   "cadenas": {"llamadas": 400_000_000, "pasos": 3_000_000_000, "fin_s": 7.25 * 3600.0},
                   "estimacion": {"segundos": 900.0, "fin_s": 7.5 * 3600.0 - MARGEN_ESCRITURA},
                   "memoria_proceso": GIB, "memoria_total": 8 * GIB,
                   "disco_bytes": 2_000_000_000, "datos_bytes": 1_800_000_000},
        "sello": seal_info,
    }


# ------------------------------------------------------------------ referencias (4.1; C2: extremos completos)
def real_refs():
    with open(RES_T21, "rb") as fh:
        raw = fh.read()
    check("T2-2 referencias: SHA-256 de RESULTADO_T2_1.json", tr.hashlib.sha256(raw).hexdigest() == RES_T21_SHA)
    R = json.loads(raw.decode("utf-8"))
    I = R["instancias"]["TIRA_6x6"]
    J, Ling, src = {}, {}, {}
    for r in tm.TARGETS:
        key = "%s|%s" % (tm.fstr(r), tm.fstr(r))
        cell = I["celdas"][key]
        row4 = cell["por_k"]["4"]
        lo, hi = (F(x) for x in row4["intervalo_P_sil_exacto"])
        check("T2-2 referencias %s: intervalo exacto con 0 < L <= U, vista coherente y e_comp(4) coherente" % key,
              0 < lo <= hi and F(Decimal(row4["intervalo_P_sil"][0])) <= lo
              and hi <= F(Decimal(row4["intervalo_P_sil"][1])) * (1 + F(1, 10 ** 20))
              and abs(tm.e_comp_of((lo, hi)) - F(Decimal(row4["e_comp"]))) <= F(1, 10 ** 20))
        J[tm.fstr(r)] = (lo, hi)
        Ling[tm.fstr(r)] = {}
        for k in tm.D5_K:
            row = cell["por_k"][str(k)]
            check("T2-2 referencias %s k = %d: H_comp COMPLETO (L_ing exacto)" % (key, k), row["H_comp"] == "COMPLETO")
            Ling[tm.fstr(r)][k] = tm.floor25_bracket(row["L_ing"])
        src[tm.fstr(r)] = ["instancias", "TIRA_6x6", "celdas", key, "por_k", "4", "intervalo_P_sil_exacto"]
    with open(RES_DMC, "rb") as fh:
        raw = fh.read()
    check("T2-2 referencias: SHA-256 de RESULTADO_DMC.json", tr.hashlib.sha256(raw).hexdigest() == RES_DMC_SHA)
    M = json.loads(raw.decode("utf-8"))
    cp = M["clopper_pearson"]
    check("T2-2 referencias: D-mc COMPATIBLE con x = 477 en N = 7e9 y 99.9 %",
          M["D_mc"]["decision"] == "COMPATIBLE" and cp["x"] == 477 and cp["N"] == 7_000_000_000
          and cp["alfa"] == "0.001")
    dmc = (F(Decimal(cp["p_inf"])), F(Decimal(cp["p_sup"])))
    return {"J": J, "L_ing": Ling, "dmc": {tm.fstr(tm.TARGETS[0]): dmc},
            "fuentes": {"RESULTADO_T2_1.json": RES_T21_SHA, "RESULTADO_DMC.json": RES_DMC_SHA,
                        "filas": src, "D_mc": ["clopper_pearson", "p_inf", "p_sup"],
                        "L_ing": "vista de 25 cifras redondeada hacia abajo -> [v, v + ulp]"}}


def toy_refs(cfg, model):
    """Referencias sinteticas del juguete: P_sil y L_ing exactos por enumeracion; J y D-mc fijados por una regla
    declarada en cfg['referencias_juguete'] (no son resultados de ningun metodo)."""
    rates = [F(x) for x in cfg["escalera"]["extremos"][1:]]
    ex, _ = model.exact(rates)
    spec = cfg["referencias_juguete"]
    J, Ling = {}, {}
    for r, (dl, dh) in zip(rates, spec["anchuras"]):
        P = ex[r]["P"]
        J[tm.fstr(r)] = (P * (1 - F(dl)), P * (1 + F(dh)))
        Ling[tm.fstr(r)] = {k: (ex[r]["L_ing"][k], ex[r]["L_ing"][k]) for k in tm.D5_K}
    P1 = ex[rates[0]]["P"]
    d1, d2 = (F(x) for x in spec["dmc"])
    return {"J": J, "L_ing": Ling, "dmc": {tm.fstr(rates[0]): (P1 * d1, P1 * d2)},
            "fuentes": {"juguete": "enumeracion exacta", "regla": spec},
            "P_exacto": {tm.fstr(r): ex[r]["P"] for r in rates}}


def refs_json(refs):
    def iv(x):
        return [tm.fstr(x[0]), tm.fstr(x[1])]
    out = {"J": {k: iv(v) for k, v in refs["J"].items()},
           "L_ing": {k: {str(kk): iv(vv) for kk, vv in v.items()} for k, v in refs["L_ing"].items()},
           "dmc": {k: iv(v) for k, v in refs["dmc"].items()}, "fuentes": refs["fuentes"]}
    if "P_exacto" in refs:
        out["P_exacto"] = {k: tm.fstr(v) for k, v in refs["P_exacto"].items()}
    return out


# ------------------------------------------------------------------ sello
def code_shas():
    return {c: tr.sha_file(os.path.join(UDIR, c)) for c in CODE}


def check_seal():
    if not os.path.exists(SELLO):
        raise CheckFail("plan sin sellar: falta %s" % SELLO)
    S = json.load(open(SELLO, encoding="ascii"))
    got = code_shas()
    ok = (S.get("plan") == PLAN and S.get("sha256") == tr.sha_file(PLAN) and S.get("codigo") == got
          and S.get("escalera", {}).get("archivo") == ESCALERA
          and S.get("escalera", {}).get("sha256") == tr.sha_file(ESCALERA) and not S.get("sintetico"))
    check("T2-2 sello: SHA-256 del plan, de los %d programas y de ESCALERA_Y_CLAVES.json = SELLO_T2_2.json"
          % len(CODE), ok, {"sello": S.get("codigo"), "actual": got})
    return {"plan": S["plan"], "sha256": S["sha256"], "codigo": S["codigo"], "escalera": S["escalera"],
            "sellado": S.get("sellado"), "sello_sha256": tr.sha_file(SELLO)}


# ------------------------------------------------------------------ etapas
class Stage:
    def __init__(self, out, name):
        self.out, self.name = out, name

    def __enter__(self):
        self.t = time.time()
        print("== etapa %s" % self.name, flush=True)
        return self

    def __exit__(self, et, ev, tb):
        rec = {"etapa": self.name, "segundos": round(time.time() - self.t, 3)}
        if et is None:
            rec["estado"] = "PASS"
        elif issubclass(et, CheckFail):
            rec.update(estado="FALLA", motivo="control fallido: %s" % ev)
        elif issubclass(et, (tmo.Agotado, tr.DatosNoCaben)):
            rec.update(estado="NO_EVALUADO", motivo="techo alcanzado: %s" % ev)
        elif issubclass(et, Parada):
            rec.update(estado="NO_EVALUADO", motivo=str(ev))
        elif issubclass(et, BrokenProcessPool):
            rec.update(estado="NO_EVALUADO", motivo="cierre anormal de un proceso de trabajo: %s" % ev)
        else:
            rec.update(estado="NO_EVALUADO", motivo="error: %s: %s" % (et.__name__, ev))
        self.out["etapas"].append(rec)
        if et is not None:
            print("  etapa %s: %s (%s)" % (self.name, rec["estado"], rec.get("motivo")), flush=True)
        return False


class Parada(Exception):
    """Parada de 7.1 con nombre propio (ANCLA_INCOMPLETA, C-tiempo no cabe, DATOS_NO_CABEN)."""


def close_stages(out):
    done = {s["etapa"] for s in out["etapas"]}
    bad = next((s for s in out["etapas"] if s["estado"] != "PASS"), None)
    for nm in PLANNED:
        if nm not in done:
            out["etapas"].append({"etapa": nm, "estado": "NO_EVALUADO", "motivo": "no alcanzada" + (
                " (%s: %s)" % (bad["etapa"], bad.get("motivo")) if bad else "")})
    out["etapas"].append({"etapa": "p_conex", "estado": "NO_EJECUTADO",
                          "motivo": "opcional (plan 3.4); no implementado en esta version del codigo; sin efecto "
                                    "en las decisiones"})
    st = [s["estado"] for s in out["etapas"] if s["etapa"] in PLANNED]
    out["estado"] = "PASS" if all(x == "PASS" for x in st) else ("FALLA" if "FALLA" in st else "NO_EVALUADO")
    if out["estado"] != "PASS":
        out["decisiones"] = {"veredicto": "PILOTO_INCOMPLETO", "motivo": "corrida no completa: " + (
            "%s: %s" % (bad["etapa"], bad.get("motivo")) if bad else "etapas no alcanzadas"),
            "VENTAJA_TEMPORAL": "NO_ACREDITADA"}


class Guard:
    """Techos del proceso principal: plazo de la etapa, memoria del proceso y disco de la corrida."""

    def __init__(self, cfg, run_dir, t0, clock=None):
        self.cfg, self.run_dir, self.t0 = cfg, run_dir, t0
        self.clock = clock or time.time            # inyectable en las pruebas (reloj sintetico, sin esperas)
        # fase 1 del supervisor: plazo duro t0 + global - verificar, aceptacion 60 s antes; se deja MARGEN_ESCRITURA
        self.deadline = t0 + cfg["techos"]["global_s"] - cfg["techos"]["verificar_s"] - MARGEN_ESCRITURA
        self.mem_peak = 0
        self.disk_peak = 0
        self.last_disk = 0.0

    def __call__(self, where):
        if self.clock() > self.deadline:
            raise tmo.Agotado("techo de tiempo (%s)" % where)
        m = tmo.mem_now()
        if m is None:
            raise tmo.Agotado("fallo de la medida de memoria del proceso principal (%s)" % where)
        self.mem_peak = max(self.mem_peak, m)
        if m > self.cfg["techos"]["memoria_proceso"]:
            raise tmo.Agotado("techo de memoria del proceso principal (%s)" % where)
        now = self.clock()
        if now - self.last_disk >= POLL:
            self.last_disk = now
            d = tr.dir_bytes(self.run_dir)
            self.disk_peak = max(self.disk_peak, d)
            if d > self.cfg["techos"]["disco_bytes"]:
                raise tr.DatosNoCaben("techo de disco: %d bytes en la carpeta de la corrida (%s)" % (d, where))


def ladder_for(cfg, model):
    lam_of = (lambda r: tm.lam_from_laws(model.laws(r)))
    ends = [F(x) for x in cfg["escalera"]["extremos"]]
    rates = tm.build_ladder(lam_of, ends)
    rec = tm.ladder_record(rates, lam_of, ends)
    return rates, rec


def stage_prep(out, cfg, model_main, refs_fn, guard):
    rates, rec = ladder_for(cfg, model_main)
    check("T2-2 escalera: recalculada desde las leyes exactas = sellada (%d niveles; t = %d)" % (
        rec["niveles"], rec["t"]), rec["tasas"] == cfg["escalera"]["tasas"])
    check("T2-2 escalera: ningun paso mayor que la regla; los extremos son niveles",
          rec["pasos_no_mayores_que_la_regla"] and all(v is not None for v in rec["indices_de_los_extremos"].values()))
    lvd = tmo.LevelData(model_main, rates)
    check("T2-2 claves y logaritmos: SHA-256 = sellados", lvd.shas() == cfg["escalera"]["shas"], lvd.shas())
    check("T2-2 formato: NS <= 255 (recuento <= 255 y d <= 255) y claves < 65536",
          model_main.NS <= 255 and lvd.kt.nkeys <= 65536)
    refs = refs_fn(cfg, model_main)
    out["escalera"] = rec
    out["claves"] = {"numero": lvd.kt.nkeys, "sha256": lvd.shas()}
    out["referencias"] = refs_json(refs)
    guard("preparacion")
    return lvd, refs


# ------------------------------------------------------------------ ancla
LECTURA_ABORTADA = ("prefijo estadistico disponible dentro de una etapa abortada: no es un ancla aceptada; sin "
                    "intervalo, semillas ni estimacion (R1)")


def anchor_cost(C, cfg, recs, done, nxt, t_start):
    """Coste del ancla: contador compartido y lo que declaran los registros de los bloques completos e
    interrumpidos (llamadas, evitadas y picos de memoria). Un bloque cuyo proceso fallo no trae recuentos."""
    known = [r for r in recs if r.get("llamadas") is not None]
    return {"llamadas": C["ancla_llamadas"].value, "techo_llamadas": cfg["techos"]["ancla"]["llamadas"],
            "llamadas_en_registros": sum(r["llamadas"] for r in known),
            "evitadas_en_registros": sum(r["evitadas"] for r in known),
            "registros_sin_recuentos": len(recs) - len(known),
            "bloques_simulados": len(done), "bloques_interrumpidos": len(recs) - len(done),
            "bloques_lanzados": nxt, "segundos": round(time.time() - t_start, 3),
            "orden_de_llegada": [r["b"] for r in recs],
            "memoria_pico_procesos": max([r.get("memoria_pico") or 0 for r in recs] or [0])}


def stage_anchor(out, cfg, ex, C, run_dir, guard, cp):
    """Bloques del ancla en vuelo (como mucho en_vuelo a la vez), registro en orden de llegada y parada exacta
    en el exito s sobre el prefijo contiguo en orden logico. El coste queda registrado aunque la etapa se corte.
    R1: cualquier bloque interrumpido (techo de tiempo, llamadas o memoria, medida fallida o error del proceso)
    y cualquier techo del proceso principal dejan la etapa NO_EVALUADA, tambien si el prefijo ya tiene el exito s:
    sin s en el prefijo es ANCLA_INCOMPLETA; con s es ANCLA_ABORTADA (prefijo disponible, no aceptado). Los
    bloques en vuelo se terminan y se registran; solo un techo del proceso principal los corta, con el motivo
    propio de la parada cooperativa."""
    A = cfg["ancla"]
    s, Ba, nfl = A["s"], A["B_a"], A["en_vuelo"]
    t_start = time.time()
    C["ancla_plazo"].value = min(t_start + cfg["techos"]["ancla"]["segundos"],
                                 guard.t0 + cfg["techos"]["cadenas"]["fin_s"])
    path = os.path.join(run_dir, "ancla", "bloques.jsonl")
    os.makedirs(os.path.dirname(path))
    done, inflight, recs, cuts = {}, {}, [], []
    nxt, gfail, stop = 0, None, False
    out["ancla"] = {"estado": "INTERRUMPIDA", "registro": "ancla/bloques.jsonl"}
    try:
        while True:
            while not stop and len(inflight) < nfl:
                inflight[ex.submit(tmo.w_anchor, nxt)] = nxt
                nxt += 1
            if not inflight:
                break
            fin, _ = wait(list(inflight), timeout=POLL, return_when=FIRST_COMPLETED)
            for f in fin:
                b = inflight.pop(f)
                try:
                    r = f.result()
                except Exception as exw:          # error del proceso de trabajo: se registra como interrupcion
                    r = {"b": b, "interrumpido": "error del proceso de trabajo: %s: %s" % (type(exw).__name__, exw),
                         "error": True, "llamadas": None, "evitadas": None, "memoria_pico": None,
                         "flujo": tmo.label_anchor(cfg["semilla_maestra"], b)}
                tr.append_jsonl(path, r)
                recs.append(r)
                if "interrumpido" in r:
                    cuts.append("bloque %d: %s" % (b, r["interrumpido"]))
                    stop = True
                else:
                    done[b] = r
            if not stop and tm.anchor_prefix(done, s, Ba)["T"] is not None:
                stop = True                       # deja de repartir; los bloques en vuelo se terminan
            try:
                guard("ancla")
            except Exception as ex2:
                if gfail is None:
                    gfail = "proceso principal: %s" % ex2
                    stop = True
                    C["ancla_plazo"].value = tmo.PARADA    # parada cooperativa, con su motivo propio
    finally:
        out["ancla"]["coste"] = anchor_cost(C, cfg, recs, done, nxt, t_start)
    pref = tm.anchor_prefix(done, s, Ba)
    causes = ([gfail] if gfail else []) + cuts
    if pref["T"] is None:
        why = "; ".join(causes) or "sin exito s"
        out["ancla"].update(estado="ANCLA_INCOMPLETA", observado=pref, motivo=why, alcance=ALCANCE_ANCLA)
        raise Parada("ANCLA_INCOMPLETA: %s; no se calcula intervalo y no se lanzan las cadenas" % why)
    if causes:
        why = "; ".join(causes)
        out["ancla"].update(estado="ANCLA_ABORTADA", prefijo_disponible=pref, motivo=why,
                            lectura=LECTURA_ABORTADA, alcance=ALCANCE_ANCLA)
        raise Parada("ANCLA_ABORTADA: %s; el prefijo tiene el exito s pero la etapa no se acepta; no se lanzan "
                     "C-tiempo ni las cadenas" % why)
    out["ancla"].update(anchor_summary(done, pref, cfg, cp))
    return out["ancla"]


def anchor_summary(done, pref, cfg, cp):
    """Desde los bloques del registro (orden logico): T, s, sobrante, semillas, estimador e intervalo."""
    A = cfg["ancla"]
    s, Ba = A["s"], A["B_a"]
    T, bstop = pref["T"], pref["bloque_de_parada"]
    prefix = {b: done[b] for b in range(bstop + 1)}
    seeds = tm.anchor_seeds(prefix, A["semillas"], Ba)
    for sd in seeds:
        r = done[sd["bloque"]]
        sd["configuracion"] = r["configuraciones"][r["posiciones"].index(sd["posicion"])]
    left = sorted(b for b in done if b > bstop)
    iv = tm.anchor_interval(s, T, cp, A["alfa_por_cola"], A["objetivo_cp"])
    ph = F(s - 1, T - 1)
    return {"estado": "COMPLETA", "T": T, "s": s, "Phat_0": tm.fstr(ph), "Phat_0_vista": tm.dstr(ph, 12),
            "intervalo_97_5": [tm.fstr(iv["l"]), tm.fstr(iv["u"])], "intervalo_vista": iv["vista"],
            "intervalo_detalle": {k: v for k, v in iv.items() if k not in ("l", "u")},
            "prefijo": pref, "recuento_por_bloque_del_prefijo": [len(done[b]["posiciones"]) for b in range(bstop + 1)],
            "sobrante": {"bloques": left, "historias_en_bloques_posteriores": len(left) * Ba,
                         "silenciosas_en_bloques_posteriores": sum(len(done[b]["posiciones"]) for b in left),
                         "historias_del_bloque_de_parada_tras_T": pref["historias_del_bloque_de_parada_tras_T"],
                         "lectura": "trabajo fisico cobrado y apartado: no entra en T, en la estimacion ni en "
                                    "las semillas"},
            "semillas": seeds, "alcance": ALCANCE_ANCLA}


def anchor_view(rec):
    """Parte del registro del ancla que se recalcula desde ancla/bloques.jsonl (sin coste)."""
    return {k: v for k, v in rec.items() if k not in ("coste", "registro")}


def read_anchor(run_dir, cfg, cp):
    """Lee ancla/bloques.jsonl. Devuelve los registros, los completos por bloque, los interrumpidos, el prefijo
    recalculado y el resumen del ancla (solo si el prefijo tiene el exito s y ningun bloque se interrumpio: R1)."""
    recs = tr.read_jsonl(os.path.join(run_dir, "ancla", "bloques.jsonl"))
    done, cuts, seen = {}, [], set()
    for r in recs:
        if r["b"] in seen:
            raise ValueError("bloque del ancla repetido: %d" % r["b"])
        seen.add(r["b"])
        if "interrumpido" in r:
            cuts.append(r)
        else:
            done[r["b"]] = r
    pref = tm.anchor_prefix(done, cfg["ancla"]["s"], cfg["ancla"]["B_a"])
    summ = anchor_summary(done, pref, cfg, cp) if pref["T"] is not None and not cuts else None
    return {"registros": recs, "completos": done, "cortes": cuts, "prefijo": pref, "resumen": summ}


# ------------------------------------------------------------------ C-tiempo
def stage_ctiempo(out, cfg, ex, C, lvd, anchor, guard):
    t = lvd.t
    ct = cfg["c_tiempo"]
    t_start = time.time()
    C["c_tiempo_plazo"].value = min(t_start + cfg["techos"]["c_tiempo"]["segundos"],
                                    guard.t0 + cfg["techos"]["cadenas"]["fin_s"])
    lv = sorted({0, t // 2, t})
    seed1 = anchor["semillas"][0]["configuracion"]
    tasks = [(rep, j, seed1, ct["pasos"]) for rep in range(1, ct["repeticiones"] + 1) for j in lv]
    res = list(ex.map(tmo.w_ctiempo, tasks))
    out["c_tiempo"] = {"medidas": res, "niveles": lv, "lectura": "inicializacion instrumental; valores descartados"}
    bad = [r for r in res if "interrumpido" in r]
    if bad:
        raise tmo.Agotado("C-tiempo interrumpido: %s" % bad[0]["interrumpido"])
    proj = project(res, lv, t, cfg, C, guard.t0 + cfg["techos"]["cadenas"]["fin_s"])
    out["c_tiempo"]["proyeccion"] = proj
    guard("c_tiempo")
    if not proj["cabe"]:
        raise Parada("C-tiempo proyecta que las cadenas no caben: %s" % proj["motivo"])


def project(res, lv, t, cfg, C, chains_end):
    per = {}
    for j in lv:
        rs = [r for r in res if r["nivel"] == j]
        n = sum(r["pasos"] for r in rs)
        per[j] = (sum(r["segundos"] for r in rs) / n, sum(r["llamadas"] for r in rs) / n,
                  sum(r["bytes"] for r in rs) / n)

    def interp(j, k):
        lo = max(x for x in lv if x <= j)
        hi = min(x for x in lv if x >= j)
        if lo == hi:
            return per[lo][k]
        w = (j - lo) / (hi - lo)
        return per[lo][k] * (1 - w) + per[hi][k] * w
    ch = cfg["cadenas"]
    nsteps = ch["L"] * ch["M"] * ch["T_init"]
    secs = sum(nsteps * interp(j, 0) for j in range(t + 1)) / cfg["procesos"]
    calls = sum(nsteps * interp(j, 1) for j in range(t + 1))
    nbytes = sum(nsteps * interp(j, 2) for j in range(t + 1))
    steps = nsteps * (t + 1)
    left = chains_end - time.time()
    caps = cfg["techos"]["cadenas"]
    why = []
    if secs > left:
        why.append("tiempo %.0f s > %.0f s restantes" % (secs, left))
    if calls > caps["llamadas"]:
        why.append("llamadas %.3g > %d" % (calls, caps["llamadas"]))
    if steps > caps["pasos"]:
        why.append("pasos %d > %d" % (steps, caps["pasos"]))
    if nbytes > cfg["techos"]["datos_bytes"] - C["bytes"].value:
        why.append("datos %.3g bytes > techo" % nbytes)
    return {"segundos": round(secs, 1), "llamadas": round(calls), "pasos": steps, "bytes": round(nbytes),
            "segundos_restantes": round(left, 1), "cabe": not why, "motivo": "; ".join(why) or None,
            "lectura": "PROYECCION con T = T_init en los t + 1 niveles, %d procesos ideales; no es una cota"
                       % cfg["procesos"]}


# ------------------------------------------------------------------ cadenas
def stage_chains(out, cfg, ex, C, anchor, run_dir, guard):
    ch = cfg["cadenas"]
    os.makedirs(os.path.join(run_dir, "cadenas"))
    C["cadenas_plazo"].value = guard.t0 + cfg["techos"]["cadenas"]["fin_s"]
    seeds = anchor["semillas"]
    tasks = [(l, i, seeds[l - 1]["configuracion"], run_dir) for i in range(1, ch["M"] + 1)
             for l in range(1, ch["L"] + 1)]
    t_start = time.time()
    futs = [ex.submit(tmo.w_pipeline, task) for task in tasks]
    res, cut, nofit = [], None, False
    pending = set(futs)
    while pending:
        fin, pending = wait(pending, timeout=POLL, return_when=FIRST_COMPLETED)
        for f in fin:
            r = f.result()
            res.append(r)
            if "interrumpido" in r:
                cut = cut or r["interrumpido"]
                nofit = nofit or r.get("datos_no_caben", False)
                C["cadenas_plazo"].value = tmo.PARADA   # parada cooperativa del resto, con su motivo propio
        try:
            guard("cadenas")
        except Exception as ex2:
            cut = cut or str(ex2)
            C["cadenas_plazo"].value = tmo.PARADA
    out["cadenas"] = {"repeticiones": sorted([[r["semilla"], r["repeticion"], r["niveles"]] for r in res]),
                      "coste": {"llamadas": C["cadenas_llamadas"].value, "pasos": C["cadenas_pasos"].value,
                                "bytes_registro": C["bytes"].value, "segundos": round(time.time() - t_start, 3),
                                "memoria_pico_procesos": max([r.get("memoria_pico", 0) for r in res] or [0])}}
    if cut:
        if nofit:
            raise Parada("DATOS_NO_CABEN: %s; se conserva todo lo escrito" % cut)
        raise tmo.Agotado(cut)


# ------------------------------------------------------------------ evaluacion desde los registros
def chain_identity_problems(regs, master):
    """Registros de cadena cuyo flujo, semilla, repeticion o nivel no son los declarados (semilla maestra)."""
    bad = []
    for (l, i), rr in sorted(regs.items()):
        for j, r in enumerate(rr):
            if (r.get("flujo") != tmo.label_chain(master, l, i, j) or r.get("semilla") != l
                    or r.get("repeticion") != i or r.get("nivel") != j):
                bad.append("%d.%d nivel %d" % (l, i, j))
    return bad


def evaluate(run_dir, cfg, lvd, refs, anchor, guard):
    """Recalcula desde disco: cadenas (contadores, continuidad, criterio de FF en cada ampliacion), Bennett,
    estimacion, diagnosticos y decisiones. La salida no contiene tiempos: es reproducible bit a bit.
    R3: el plazo se consulta al empezar cada nivel, tras cada cadena (lectura, criterios y diagnosticos), tras
    R-hat y cada Bennett, y al terminar los diagnosticos y decisiones; una operacion que vuelve tarde no pasa.
    Riesgo de memoria (no demostrado; punto 6 del dictamen del codigo v1): cada archivo de cadena se lee entero
    y se guardan los histogramas agregados de las L x M cadenas de dos niveles consecutivos."""
    ch = cfg["cadenas"]
    L, M, B = ch["L"], ch["M"], ch["B"]
    t = lvd.t
    lengths = tm.ff_lengths(ch["T_init"], ch["lambda_ini"], ch["ampliaciones_max"])
    eps = F(ch["eps"])
    reps = [(l, i) for l in range(1, L + 1) for i in range(1, M + 1)]
    regs = {}
    for l, i in reps:
        rr = tr.read_jsonl(tmo.registry_path(run_dir, l, i))
        if [r["nivel"] for r in rr] != list(range(t + 1)):
            raise CheckFail("registro de la semilla %d, repeticion %d incompleto" % (l, i))
        if any(r["estado"] not in ("CRITERIO_FF_CUMPLE", "CRITERIO_FF_NO_CUMPLE") for r in rr):
            raise CheckFail("cadena no evaluada en la semilla %d, repeticion %d" % (l, i))
        regs[(l, i)] = rr
    bad = chain_identity_problems(regs, cfg["semilla_maestra"])
    if bad:
        raise CheckFail("flujo o identidad de los registros de cadena distintos de los declarados: %s" % bad[:3])
    targets = [F(x) for x in cfg["escalera"]["extremos"][1:]]
    jt_of = {tm.fstr(r): cfg["escalera"]["tasas"].index(tm.fstr(r)) for r in targets}
    lnq = {k: [] for k in reps}
    ratios, levels = [], []
    occ = {tm.fstr(r): {k: [[None] * M for _ in range(L)] for k in tm.D5_K} for r in targets}
    prevH, prevN = None, None
    for j in range(t + 1):
        guard("estimacion, nivel %d" % j)
        H, N, halves, d3, ffst, null = {}, {}, [], {}, {"cumple": 0, "no_cumple": 0}, []
        n_half = min(regs[k][j]["T"] // 2 for k in reps)
        cnt = {"pasos": 0, "aceptaciones": 0, "llamadas": 0, "evitadas": 0, "transiciones": 0}
        sigmas = {}
        for (l, i) in reps:
            rec = regs[(l, i)][j]
            path = tmo.chain_path(run_dir, j, l, i)
            with open(path, "rb") as fh:
                raw = fh.read()
            if (tr.hashlib.sha256(raw).hexdigest() != rec["archivo"]["sha256"]
                    or len(raw) != rec["archivo"]["bytes"]):
                raise CheckFail("archivo de cadena distinto de su registro: %s" % os.path.basename(path))
            bl = tr.decode_file(raw, lvd.kt.nkeys, (j, l, i))
            hists = [h for _, h in bl]
            if len(hists) * B != rec["T"] or any(m["muestras"] != B for m, _ in bl):
                raise CheckFail("numero de bloques o muestras distinto de T en %s" % os.path.basename(path))
            tot = {k: sum(m[k] for m, _ in bl) for k in ("aceptaciones", "llamadas", "evitadas", "transiciones")}
            tot["pasos"] = rec["T"]
            if tot != rec["contadores"]:
                raise CheckFail("contadores distintos del registro en %s" % os.path.basename(path))
            for k in cnt:
                cnt[k] += tot[k]
            start = anchor["semillas"][l - 1]["configuracion"] if j == 0 else regs[(l, i)][j - 1]["final"]
            if rec["E0"] != sorted(start):
                raise CheckFail("E_0 del nivel %d no es el estado final del nivel anterior (%d.%d)" % (j, l, i))
            crits = rec["criterios"]
            for k, crit in enumerate(crits):
                again = tm.ff_criterion(hists[:lengths[k] // B], j, t, lvd.lt, eps)
                if json.loads(json.dumps(again)) != crit:
                    raise CheckFail("criterio de FF recalculado distinto en %s, comprobacion %d"
                                    % (os.path.basename(path), k))
            okseq = all(not c["cumple"] for c in crits[:-1]) and len(crits) <= len(lengths) and crits
            last = crits[-1]["cumple"]
            if not okseq or (rec["estado"] == "CRITERIO_FF_CUMPLE") != last or (
                    not last and len(crits) != len(lengths)) or rec["T"] != lengths[len(crits) - 1]:
                raise CheckFail("secuencia de ampliaciones incoherente en %s" % os.path.basename(path))
            fsig = tm.sig_from_classes(_cls_of(rec["final"], lvd), lvd.kt.key_of)
            if [list(x) for x in fsig] != rec["firma_final"] or fsig not in hists[-1]:
                raise CheckFail("firma del estado final ausente del ultimo bloque en %s" % os.path.basename(path))
            ffst["cumple" if last else "no_cumple"] += 1
            if crits[-1].get("varianza_nula"):
                null.append("%d.%d: %s" % (l, i, ",".join(crits[-1]["varianza_nula"])))
            sigmas[(l, i)] = crits[-1].get("sigma")
            H[(l, i)] = tm.merge(hists)
            N[(l, i)] = rec["T"]
            d3["%d.%d" % (l, i)] = tm.batch_means(hists, j, t, lvd.lt)
            nb = n_half // B
            lastb = hists[len(hists) - 2 * nb:]
            halves.append(tm.merge(lastb[:nb]))
            halves.append(tm.merge(lastb[nb:]))
            for r in targets:
                if jt_of[tm.fstr(r)] == j:
                    for k in tm.D5_K:
                        occ[tm.fstr(r)][k][l - 1][i - 1] = tm.occupancy_fraction(hists, k)
            guard("estimacion, nivel %d, cadena %d.%d" % (j, l, i))
        rh = {side: tm.rhat_series(halves, n_half, j, t, lvd.lt, side) for side in tm.ff_sides(j, t)}
        guard("estimacion, nivel %d, R-hat" % j)
        d3s = summarize_d3(d3)
        levels.append({"nivel": j, "tasa": cfg["escalera"]["tasas"][j], "R_hat": rh, "d3": d3s,
                       "d4": {"pasos": cnt["pasos"],
                              "tasa_de_transicion": tm.dstr(F(cnt["transiciones"], cnt["pasos"]), 8),
                              "llamadas_por_paso": tm.dstr(F(cnt["llamadas"], cnt["pasos"]), 8),
                              "evitadas_por_paso": tm.dstr(F(cnt["evitadas"], cnt["pasos"]), 8),
                              "contadores": cnt},
                       "criterio_FF": ffst, "varianza_nula": null,
                       "sigma_final": {"%d.%d" % k: v for k, v in sorted(sigmas.items())}})
        if j > 0:
            rq = {"cociente": j - 1, "ln_q": {}, "d8": {}}
            for key in reps:
                bn = tm.bennett(tm.bennett_items(prevH[key], lvd.lt, j - 1), prevN[key],
                                tm.bennett_items(H[key], lvd.lt, j - 1), N[key])
                if bn["subdesbordamientos"]:
                    raise tm.ErrorNumerico("subdesbordamiento en Bennett (cociente %d, %d.%d)" % (j - 1, *key))
                lnq[key].append(bn["ln_q"])
                rq["ln_q"]["%d.%d" % key] = repr(bn["ln_q"])
                rq["d8"]["%d.%d" % key] = overlap(prevH[key], H[key], lvd.lt, j - 1, bn["ln_q"])
                guard("estimacion, cociente %d, Bennett %d.%d" % (j - 1, *key))
            ratios.append(rq)
        prevH, prevN = H, N
    res = finish_estimation(cfg, refs, anchor, lnq, levels, ratios, occ, jt_of, regs)
    guard("estimacion: al terminar los diagnosticos y las decisiones")
    return res


def _cls_of(pairs, lvd):
    cls = [0] * lvd.kt.NS
    for s, c in pairs:
        cls[s] = c
    return cls


def summarize_d3(d3):
    out = {}
    for key, rec in sorted(d3.items()):
        for side, v in rec.items():
            o = out.setdefault(side, {"calculados": 0, "indeterminados": [], "tau": {}, "ESS": {}})
            if v["estado"] == "CALCULADO":
                o["calculados"] += 1
                o["tau"][key] = v["tau"]
                o["ESS"][key] = v["ESS"]
            else:
                o["indeterminados"].append("%s: %s" % (key, v["motivo"]))
    return out


def overlap(Hn, Hd, lt, j, lnC):
    """d8: fraccion de muestras con |ln(C rho)| <= ln 10 a cada lado (C final). Informativo."""
    l10 = math.log(10.0)
    a = sum(m for sig, m in Hn.items() if abs(lnC + lt.z(sig, j)) <= l10)
    b = sum(m for sig, m in Hd.items() if abs(lnC + lt.z(sig, j)) <= l10)
    return [tm.dstr(F(a, tm.hist_n(Hn)), 6), tm.dstr(F(b, tm.hist_n(Hd)), 6)]


def finish_estimation(cfg, refs, anchor, lnq, levels, ratios, occ, jt_of, regs):
    iv = {"l": F(anchor["intervalo_97_5"][0]), "u": F(anchor["intervalo_97_5"][1])}
    ph0 = F(anchor["Phat_0"])
    tq = tm.student_quantile(tm.STUDENT_P, cfg["cadenas"]["L"] - 1)
    per, dec = {}, {}
    for rk, jt in jt_of.items():
        est = tm.estimate_rate(lnq, jt, iv, ph0, tq, cfg["cadenas"]["L"], cfg["cadenas"]["M"])
        J = refs["J"][rk]
        a1 = [lv["nivel"] for lv in levels[:jt + 1] if any(v.get("alerta") for v in lv["R_hat"].values())]
        ind = []
        for lv in levels[:jt + 1]:
            for side, v in lv["R_hat"].items():
                if v["estado"] == tm.IND:
                    ind.append("R-hat nivel %d lado %s: %s" % (lv["nivel"], side, v["motivo"]))
            for x in lv["varianza_nula"]:
                ind.append("varianza nula en el criterio de FF, nivel %d, cadena %s" % (lv["nivel"], x))
        if est["s_Q_nula"]:
            ind.append("s_Q = 0 en la tasa %s" % rk)
        d5 = {}
        for k in tm.D5_K:
            lo_i, hi_i = refs["L_ing"][rk][k]
            Lj, Uj = J
            d5[str(k)] = tm.d5_occupancy(occ[rk][k], lo_i / Uj, min(F(1), hi_i / Lj))
        a2 = any(v["alerta"] for v in d5.values())
        ff = all(regs[key][j]["estado"] == "CRITERIO_FF_CUMPLE" for key in regs for j in range(jt + 1))
        flags = {"ff_cumple": ff, "a1": bool(a1), "a2": a2, "indeterminados": ind}
        D, num = tm.decide_rate(est, J, flags, refs["dmc"].get(rk))
        dec[rk] = D
        seed_rel = (tm.dstr(tm.NE.divide(tm.NE.sqrt(tm.d_of(est["s_Q2"])), tm.d_of(est["Qbar"])), 8)
                    if est["Qbar"] > 0 else None)
        decl = []
        for key in sorted(regs):
            s2 = sum((Decimal(lv["sigma_final"]["%d.%d" % key]) ** 2 for lv in levels[:jt + 1]
                      if lv["sigma_final"].get("%d.%d" % key) is not None), Decimal(0))
            decl.append(s2.sqrt())
        per[rk] = {
            "indice_de_nivel": jt, "x": tm.fstr(est["x"]), "x_vista": tm.dstr(est["x"], 12),
            "I": [tm.fstr(est["I"][0]), tm.fstr(est["I"][1])],
            "I_vista": [tm.dstr(est["I"][0], 12), tm.dstr(est["I"][1], 12)],
            "Qbar": tm.fstr(est["Qbar"]), "s_Q2": tm.fstr(est["s_Q2"]), "Q_lo": tm.fstr(est["Q_lo"]),
            "Q_hi": tm.fstr(est["Q_hi"]), "cuantil_student": str(tq),
            "Qbar_por_semilla": [tm.dstr(x, 17) for x in est["Qbar_por_semilla"]],
            "J": [tm.fstr(J[0]), tm.fstr(J[1])],
            "medidas": {k: (tm.fstr(v) if isinstance(v, F) else v) for k, v in num.items()},
            "medidas_vista": {k: (tm.dstr(v, 8) if isinstance(v, F) else v) for k, v in num.items()},
            "alertas": {"a1_niveles": a1, "a2": a2, "indeterminados": ind, "d5": d5},
            "d7": {"dispersion_relativa_entre_semillas": seed_rel,
                   "sigma_declarado_medio_por_repeticion": tm.dstr(sum(decl, Decimal(0)) / len(decl), 8)},
            "decisiones": D}
    ver = tm.verdict(dec)
    ext = tm.extension_rule(ver, dec)
    return {"por_tasa": per, "niveles": levels, "cocientes": ratios,
            "decisiones": {"veredicto": ver, "ampliacion": ext,
                           "lectura": "nominal por tasa; sin cobertura simultanea; CRITERIO_FF_CUMPLE y "
                                      "SIN_ALERTA no prueban mezcla"}}


# ------------------------------------------------------------------ corrida
def make_counters():
    ctx = mp.get_context("spawn")
    C = {k: ctx.Value("q", 0) for k in ("ancla_llamadas", "c_tiempo_llamadas", "c_tiempo_pasos",
                                        "cadenas_llamadas", "cadenas_pasos", "bytes")}
    C.update({k: ctx.Value("d", 0.0) for k in ("ancla_plazo", "c_tiempo_plazo", "cadenas_plazo")})
    return C


def default_executor(cfg, C):
    return ProcessPoolExecutor(cfg["procesos"], mp_context=mp.get_context("spawn"), initializer=tmo.w_init,
                               initargs=(cfg, C))


def jnorm(x):
    return json.loads(json.dumps(x))


def run(cfg, run_dir, model_fn, refs_fn, t0, prechecks=None, make_executor=None, clock=None):
    """Corrida completa (fisica o de juguete). Escribe RESULTADO, RECURSOS y MANIFIESTO en run_dir.
    model_fn() crea el modelo dentro de la etapa de preparacion (un fallo de dependencias queda registrado).
    make_executor(cfg, C) y clock solo se sustituyen en las pruebas sinteticas (resultados guionizados y reloj
    controlado); la corrida usa el ProcessPoolExecutor y el reloj de pared."""
    CHECKS.clear()
    out = {"unidad": cfg["unidad"], "tipo": cfg["tipo"], "programa_sha256": code_shas(), "configuracion": cfg,
           "etapas": [], "t0": t0}
    C = make_counters()
    guard = Guard(cfg, run_dir, t0, clock)
    ex = None
    res_cost = {}
    try:
        with Stage(out, "preparacion"):
            if cfg["modelo"]["tipo"] == "fisico":      # el sello queda entre los controles de la corrida (R2)
                check("T2-2 sello: el vigente = el de la configuracion de la corrida", check_seal() == cfg["sello"])
            cp = tm.load_cp(CP_PATH, CP_SHA)
            model_main = model_fn()
            if prechecks:
                prechecks(model_main, check)
            lvd, refs = stage_prep(out, cfg, model_main, refs_fn, guard)
        with Stage(out, "ancla"):
            ex = (make_executor or default_executor)(cfg, C)
            anchor = stage_anchor(out, cfg, ex, C, run_dir, guard, cp)
        with Stage(out, "c_tiempo"):
            stage_ctiempo(out, cfg, ex, C, lvd, anchor, guard)
        with Stage(out, "cadenas"):
            stage_chains(out, cfg, ex, C, anchor, run_dir, guard)
        ex.shutdown(wait=True)
        ex = None
        with Stage(out, "estimacion"):
            guard.deadline = min(guard.deadline, t0 + cfg["techos"]["estimacion"]["fin_s"],
                                 guard.clock() + cfg["techos"]["estimacion"]["segundos"])
            rd = read_anchor(run_dir, cfg, cp)
            check("T2-2 ancla: recalculada desde ancla/bloques.jsonl = la de la etapa, sin bloques interrumpidos",
                  rd["resumen"] is not None and not rd["cortes"]
                  and jnorm(rd["resumen"]) == jnorm(anchor_view(out["ancla"])))
            est = evaluate(run_dir, cfg, lvd, refs, rd["resumen"], guard)
            guard("estimacion: resultado completo, antes de aceptar la etapa")      # R3: nunca PASS tarde
            out["estimacion"] = est
            out["decisiones"] = est["decisiones"]
    except (CheckFail, tmo.Agotado, tr.DatosNoCaben, Parada, BrokenProcessPool, Exception):
        pass
    finally:
        if ex is not None:
            for k in ("ancla_plazo", "c_tiempo_plazo", "cadenas_plazo"):
                C[k].value = tmo.PARADA
            ex.shutdown(wait=True, cancel_futures=True)
    close_stages(out)
    out["controles"] = list(CHECKS)
    res_cost.update({k: C[k].value for k in C if not k.endswith("_plazo")})
    res_cost.update(memoria_pico_principal=guard.mem_peak, disco_pico=guard.disk_peak,
                    segundos=round(time.time() - t0, 3), techos=cfg["techos"])
    out["coste_resumen"] = res_cost
    tr.write_json(os.path.join(run_dir, RECURSOS), {"unidad": cfg["unidad"], "recursos": res_cost,
                                                    "etapas": out["etapas"]})
    tr.write_json(os.path.join(run_dir, RESULT), out)
    man = tr.manifest(run_dir, exclude=(MANIF,) + POSTERIORES)
    tr.write_json(os.path.join(run_dir, MANIF), {"archivos": man})
    print("T2-2: estado %s; veredicto %s" % (out["estado"], out.get("decisiones", {}).get("veredicto")), flush=True)
    return out


# ------------------------------------------------------------------ verificacion (R2, R3)
VERIFICAR_MAX_S = 1800.0      # techo de --verificar del plan (6); la configuracion registrada no puede ampliarlo
VERIFICAR_MEM = GIB    # plan 6: <= 1 GiB por proceso; el techo de memoria registrado limita solo a la corrida
LECTURAS = {0: "ENTREGA_VERIFICADA", 1: "DISCREPANCIA_O_ERROR",
            2: "INSPECCION_DE_UNA_CORRIDA_INCOMPLETA: coherente, pero sin todas las etapas PASS",
            4: "VERIFICACION_NO_CONCLUIDA: techo de --verificar alcanzado"}


def mandatory_controls(cfg):
    """Controles obligatorios: (subcadena del texto del control, numero exacto de apariciones) de las etapas
    preparacion y estimacion, segun el modo. Una entrega completa los tiene todos y ninguno falla."""
    prep = [("T2-2 escalera: recalculada desde las leyes exactas = sellada", 1),
            ("T2-2 escalera: ningun paso mayor que la regla", 1),
            ("T2-2 claves y logaritmos: SHA-256 = sellados", 1),
            ("T2-2 formato: NS <= 255", 1)]
    if cfg["modelo"]["tipo"] == "fisico":
        nt, nk = len(tm.TARGETS), len(tm.D5_K)
        prep += [("T2-2 sello: SHA-256 del plan", 1), ("T2-2 sello: el vigente = el de la configuracion", 1),
                 ("matrices iguales a las de B4", 1), ("ubicaciones (plan:", 1), ("clases no nulas (plan:", 1),
                 (": marginales X = formas cerradas", nt + 1),
                 (": leyes de sitio desde las marginales = law_of", nt + 1),
                 ("propagate de cada fallo unico = resp_formula", 1), ("T2-2 V4: propagate = XOR", 1),
                 ("T2-2 V5 (T2-0)", 1), ("T2-2 V5 (T2-1)", 1),
                 ("T2-2 referencias: SHA-256 de RESULTADO_T2_1.json", 1),
                 ("intervalo exacto con 0 < L <= U", nt), ("H_comp COMPLETO (L_ing exacto)", nt * nk),
                 ("T2-2 referencias: SHA-256 de RESULTADO_DMC.json", 1), ("T2-2 referencias: D-mc COMPATIBLE", 1)]
    return {"preparacion": prep, "estimacion": [("T2-2 ancla: recalculada desde ancla/bloques.jsonl", 1)]}


def contract_files(cfg, t):
    """Archivos que la carpeta de una corrida puede contener (fuera de MANIFIESTO y de los posteriores)."""
    L, M = cfg["cadenas"]["L"], cfg["cadenas"]["M"]
    files = {RESULT, RECURSOS, "ancla/bloques.jsonl"}
    for l in range(1, L + 1):
        for i in range(1, M + 1):
            files.add("cadenas/registro_s%02d_r%d.jsonl" % (l, i))
            files.update("cadenas/n%03d_s%02d_r%d.bin" % (j, l, i) for j in range(t + 1))
    return files


def record_sums(run_dir, cfg, rd, R):
    """Sumas de lo que declaran los registros en disco (ancla, cadenas y C-tiempo). No miden nada de nuevo."""
    out = {}
    if rd is not None:
        recs = rd["registros"]
        out["ancla_llamadas"] = sum(r["llamadas"] for r in recs if r.get("llamadas") is not None)
        out["ancla_sin_recuentos"] = sum(1 for r in recs if r.get("llamadas") is None)
    cl = cp_ = cb = 0
    regs = {}
    for l in range(1, cfg["cadenas"]["L"] + 1):
        for i in range(1, cfg["cadenas"]["M"] + 1):
            p = tmo.registry_path(run_dir, l, i)
            if os.path.exists(p):
                regs[(l, i)] = rr = tr.read_jsonl(p)
                for r in rr:
                    cl += r["contadores"]["llamadas"] + r.get("bloque_parcial_no_registrado", {}).get("llamadas", 0)
                    cp_ += r["T"]
                    cb += r["archivo"]["bytes"]
    cdir = os.path.join(run_dir, "cadenas")
    out.update(cadenas_llamadas=cl, cadenas_pasos=cp_, bytes=cb, registros_de_cadena=regs,
               bytes_en_disco=sum(os.path.getsize(os.path.join(cdir, f)) for f in os.listdir(cdir)
                                  if f.endswith(".bin")) if os.path.isdir(cdir) else 0)
    med = [m for m in R.get("c_tiempo", {}).get("medidas", []) if "interrumpido" not in m]
    out.update(c_tiempo_pasos=sum(m["pasos"] for m in med), c_tiempo_llamadas=sum(m["llamadas"] for m in med),
               c_tiempo_interrumpidas=len(R.get("c_tiempo", {}).get("medidas", [])) - len(med))
    return out


def verify(run_dir, model_factory=None, refs_factory=None, config_factory=None, clock=None):
    """--verificar sobre la entrega. Devuelve (codigo, informe): 0, 1, 2 o 4 (ver la cabecera).
    El plazo (min(verificar_s, 30 min)) cuenta desde la entrada, con lectura, hashes y reconstruccion (R3).
    config_factory(R) reconstruye la configuracion esperada: en el modo fisico, por defecto, desde el sello,
    ESCALERA_Y_CLAVES.json y real_config; en el sintetico solo si la prueba declara su funcion (R2).
    Las comprobaciones de recursos contrastan lo registrado entre si y con los techos; no vuelven a medir
    tiempos ni memoria historicos. El proceso que verifica mide su propia memoria contra 1 GiB (plan 6)."""
    clock = clock or time.time
    t_in = clock()
    rep = {"carpeta": run_dir, "fallos": [], "comprobaciones": []}

    def ok(desc, cond, detail=None):
        rep["comprobaciones"].append({"comprobacion": desc, "pasa": bool(cond), "detalle": detail})
        if not cond:
            rep["fallos"].append(desc)
        return bool(cond)
    complete = False
    try:
        R = json.load(open(os.path.join(run_dir, RESULT), encoding="ascii"))
        Mf = json.load(open(os.path.join(run_dir, MANIF), encoding="ascii"))["archivos"]
        cfg = R["configuracion"]
        # el proceso de --verificar tiene su techo del plan (1 GiB), no el de memoria registrado de la corrida
        g = Guard({"techos": dict(cfg["techos"], memoria_proceso=VERIFICAR_MEM)}, run_dir, t_in, clock)
        g.deadline = t_in + min(float(cfg["techos"]["verificar_s"]), VERIFICAR_MAX_S)
        rep["plazo_s"] = g.deadline - t_in
        g("verificar: lectura de RESULTADO y MANIFIESTO")
        # --- manifiesto y programas
        listed = {p for p, _, _ in Mf}
        actual = set(x for x, _, _ in tr.manifest(run_dir, exclude=(MANIF,) + POSTERIORES))
        ok("manifiesto: mismos archivos (ni faltan ni sobran)", listed == actual,
           {"faltan": sorted(listed - actual), "sobran": sorted(actual - listed)})
        bad = [p for p, h, n in Mf if p in actual and (tr.sha_file(os.path.join(run_dir, p)) != h
                                                       or os.path.getsize(os.path.join(run_dir, p)) != n)]
        ok("manifiesto: SHA-256 y bytes de cada archivo", not bad, bad[:5])
        g("verificar: manifiesto")
        ok("programa: SHA-256 actuales = los de la corrida", R["programa_sha256"] == code_shas())
        ok("unidad y tipo = los de la configuracion", R.get("unidad") == cfg.get("unidad")
           and R.get("tipo") == cfg.get("tipo"))
        # --- configuracion (R2-d)
        if cfg["modelo"]["tipo"] == "fisico":
            seal = check_seal()
            ok("sello vigente = el de la corrida", seal == cfg["sello"])
            fac = config_factory or (lambda R_: real_config(seal, json.load(open(ESCALERA, encoding="ascii"))))
            exp = jnorm(fac(R))
            ok("configuracion = reconstruida desde el sello, ESCALERA_Y_CLAVES.json y real_config", exp == cfg,
               sorted(k for k in set(exp) | set(cfg) if exp.get(k) != cfg.get(k)))
        else:
            ok("contrato sintetico de las pruebas: tipo humo, modelo de juguete y sello sintetico (no hay "
               "configuracion fisica que reconstruir)", cfg.get("tipo") == "humo" and cfg["modelo"]["tipo"] ==
               "juguete" and cfg.get("sello", {}).get("sintetico") is True)
            if config_factory is not None:
                exp = jnorm(config_factory(R))
                ok("configuracion = reconstruida por la funcion declarada de la prueba", exp == cfg,
                   sorted(k for k in set(exp) | set(cfg) if exp.get(k) != cfg.get(k)))
        g("verificar: configuracion")
        # --- modelo, escalera, claves y referencias
        model = (model_factory or tmo.make_model)(cfg["modelo"])
        rates, rec = ladder_for(cfg, model)
        ok("escalera recalculada = configurada = registrada", rec["tasas"] == cfg["escalera"]["tasas"]
           and jnorm(rec) == R.get("escalera"))
        lvd = tmo.LevelData(model, rates)
        ok("claves y logaritmos = sellados", lvd.shas() == cfg["escalera"]["shas"])
        ok("claves registradas = recalculadas", R.get("claves") == {"numero": lvd.kt.nkeys, "sha256": lvd.shas()})
        refs = (refs_factory or default_refs)(cfg, model)
        ok("referencias recalculadas desde sus fuentes = registradas", refs_json(refs) == R.get("referencias"))
        cp = tm.load_cp(CP_PATH, CP_SHA)
        g("verificar: modelo, escalera y referencias")
        # --- etapas (sin duplicados, orden, estado global)
        names = [s.get("etapa") for s in R["etapas"]]
        ok("etapas: cada una exactamente una vez y en el orden del plan, mas p_conex", names ==
           list(PLANNED) + ["p_conex"], names)
        st = {s["etapa"]: s for s in R["etapas"]}
        states = [st.get(nm, {}).get("estado") for nm in PLANNED]
        k = next((n for n, x in enumerate(states) if x != "PASS"), len(PLANNED))
        later = PLANNED[k + 1:] if k < len(PLANNED) else ()
        ok("etapas: PASS hasta la primera que no pasa (FALLA o NO_EVALUADO); las siguientes NO_EVALUADO no "
           "alcanzadas; p_conex NO_EJECUTADO",
           (k == len(PLANNED) or states[k] in ("FALLA", "NO_EVALUADO"))
           and all(st.get(nm, {}).get("estado") == "NO_EVALUADO"
                   and str(st.get(nm, {}).get("motivo", "")).startswith("no alcanzada") for nm in later)
           and st.get("p_conex", {}).get("estado") == "NO_EJECUTADO", states)
        want = "PASS" if k == len(PLANNED) else ("FALLA" if "FALLA" in states else "NO_EVALUADO")
        ok("estado global = el que dan las etapas", R.get("estado") == want, [R.get("estado"), want])
        complete = k == len(PLANNED) and R.get("estado") == "PASS"
        # --- controles (R2-b)
        ctr = R.get("controles")
        shape = isinstance(ctr, list) and all(isinstance(c, dict) and set(c) == {"control", "pasa", "detalle"}
                                              for c in ctr)
        ok("controles: lista de registros con control, pasa y detalle", shape)
        ctr = ctr if shape else []
        failed = [c for c in ctr if not c["pasa"]]
        if complete:
            ok("controles: ninguno fallido en una entrega completa", not failed, [c["control"] for c in failed])
        else:
            fall = [s for s in R["etapas"] if s.get("estado") == "FALLA"]
            if fall:
                ok("controles: el unico fallido es el ultimo y es el motivo de la etapa FALLA", len(failed) == 1
                   and ctr[-1] is failed[0] and fall[0].get("motivo") == "control fallido: " + failed[0]["control"])
            else:
                ok("controles: ninguno fallido si ninguna etapa esta en FALLA", not failed)
        mand = mandatory_controls(cfg)
        need = (mand["preparacion"] if states[0] == "PASS" else []) + (
            mand["estimacion"] if states[4] == "PASS" else [])
        miss = [[sub, n, sum(sub in c["control"] for c in ctr)] for sub, n in need
                if sum(sub in c["control"] for c in ctr) != n]
        ok("controles obligatorios de las etapas que pasaron: presentes con su numero exacto", not miss, miss)
        # --- decisiones de una corrida incompleta
        if not complete:
            badst = next((s for s in R["etapas"] if s.get("estado") != "PASS"), None)
            exp_dec = {"veredicto": "PILOTO_INCOMPLETO", "motivo": "corrida no completa: " + (
                "%s: %s" % (badst["etapa"], badst.get("motivo")) if badst else "etapas no alcanzadas"),
                "VENTAJA_TEMPORAL": "NO_ACREDITADA"}
            ok("corrida incompleta: decisiones de primer nivel = PILOTO_INCOMPLETO con el motivo de la etapa",
               R.get("decisiones") == exp_dec)
            ok("corrida incompleta: sin estimacion parcial", "estimacion" not in R)
        # --- inventario (R2-c)
        allowed = contract_files(cfg, lvd.t)
        required = allowed if complete else {RESULT, RECURSOS}
        ok("inventario: archivos obligatorios presentes", required <= actual, sorted(required - actual)[:5])
        ok("inventario: ningun archivo fuera del contrato de la corrida", actual <= allowed,
           sorted(actual - allowed)[:5])
        g("verificar: etapas, controles e inventario")
        # --- ancla (R1 incorporado)
        RA = R.get("ancla", {})
        sa = st.get("ancla", {}).get("estado")
        rd = None
        if not os.path.exists(os.path.join(run_dir, "ancla", "bloques.jsonl")):
            ok("ancla: sin registro solo si la etapa ancla no paso", sa != "PASS" and RA.get("estado") != "COMPLETA")
        else:
            rd = read_anchor(run_dir, cfg, cp)
            cuts = rd["cortes"]
            badflow = [r["b"] for r in rd["registros"]
                       if r.get("flujo") != tmo.label_anchor(cfg["semilla_maestra"], r["b"])]
            ok("ancla: flujo de cada bloque = etiqueta de la semilla maestra declarada", not badflow, badflow[:5])
            if cuts:
                ok("ancla: con bloques interrumpidos la etapa no se acepta ni el ancla queda COMPLETA (R1)",
                   sa != "PASS" and RA.get("estado") in ("ANCLA_ABORTADA", "ANCLA_INCOMPLETA"),
                   [c["b"] for c in cuts])
            ea = RA.get("estado")
            if ea == "COMPLETA":
                ok("ancla recalculada = registrada; COMPLETA solo con la etapa PASS y sin interrupciones",
                   rd["resumen"] is not None and jnorm(rd["resumen"]) == anchor_view(RA) and sa == "PASS")
            elif ea == "ANCLA_ABORTADA":
                ok("ancla ABORTADA: prefijo recalculado con el exito s = prefijo_disponible, sin intervalo ni "
                   "semillas, cada interrupcion en el motivo y la etapa no aceptada",
                   rd["prefijo"]["T"] is not None and jnorm(rd["prefijo"]) == RA.get("prefijo_disponible")
                   and not any(x in RA for x in ("intervalo_97_5", "semillas", "Phat_0")) and sa != "PASS"
                   and all("bloque %d: %s" % (c["b"], c["interrumpido"]) in RA.get("motivo", "") for c in cuts))
            elif ea == "ANCLA_INCOMPLETA":
                ok("ancla INCOMPLETA: el prefijo recalculado no tiene el exito s, = observado, y la etapa no paso",
                   rd["prefijo"]["T"] is None and jnorm(rd["prefijo"]) == RA.get("observado") and sa != "PASS"
                   and "intervalo_97_5" not in RA)
            else:
                ok("ancla interrumpida solo si la etapa ancla no paso", sa != "PASS")
            if "coste" in RA:
                cst = RA["coste"]
                recs = rd["registros"]
                ok("ancla: coste registrado = el de sus registros (bloques completos e interrumpidos, picos)",
                   cst.get("llamadas_en_registros") == sum(r["llamadas"] for r in recs if r.get("llamadas") is not None)
                   and cst.get("bloques_simulados") == len(rd["completos"])
                   and cst.get("bloques_interrumpidos") == len(cuts)
                   and cst.get("orden_de_llegada") == [r["b"] for r in recs]
                   and cst.get("memoria_pico_procesos") == max([r.get("memoria_pico") or 0 for r in recs] or [0]))
        if complete:
            ok("corrida completa: ancla COMPLETA", RA.get("estado") == "COMPLETA" and rd is not None
               and rd["resumen"] is not None)
        g("verificar: ancla")
        # --- cadenas, estimacion y decisiones (R2-a)
        sums = record_sums(run_dir, cfg, rd, R)
        if complete and rd is not None and rd["resumen"] is not None:
            regs = sums["registros_de_cadena"]
            badid = chain_identity_problems(regs, cfg["semilla_maestra"])
            ok("cadenas: flujo e identidad de cada nivel = los de la semilla maestra declarada", not badid, badid[:5])
            if not badid:
                est = evaluate(run_dir, cfg, lvd, refs, rd["resumen"], g)
                g("verificar: tras recalcular la estimacion")
                ok("estimacion, diagnosticos y decisiones recalculados = registrados",
                   jnorm(est) == R.get("estimacion"))
                ok("decisiones de primer nivel = las recalculadas", R.get("decisiones") == jnorm(est["decisiones"]))
                ok("decisiones: tasas exactas", sorted(est["por_tasa"]) == sorted(tm.fstr(F(x)) for x in
                                                                                 cfg["escalera"]["extremos"][1:]))
                ok("decisiones: claves exactas por tasa", all(
                    set(v["decisiones"]) == {"D_s1", "D_s2", "D_s3", "compatible", "D_s4", "D_s5", "D_s7",
                                             "estimacion_solo_diagnostico", "complementa_nominal"}
                    | ({"D_s6"} if kk in refs["dmc"] else set()) for kk, v in est["por_tasa"].items()))
            want_rep = sorted([[l, i, [{"nivel": r["nivel"], "estado": r["estado"], "T": r["T"]} for r in rr]]
                               for (l, i), rr in regs.items()])
            ok("cadenas: repeticiones registradas = las de sus registros", R.get("cadenas", {}).get("repeticiones")
               == want_rep)
            CT = R.get("c_tiempo", {})
            nct = cfg["c_tiempo"]["repeticiones"] * len(sorted({0, lvd.t // 2, lvd.t}))
            ok("c_tiempo: niveles, medidas completas y proyeccion que cabe",
               CT.get("niveles") == sorted({0, lvd.t // 2, lvd.t}) and len(CT.get("medidas", [])) == nct
               and sums["c_tiempo_interrumpidas"] == 0
               and all(m["pasos"] == cfg["c_tiempo"]["pasos"] for m in CT.get("medidas", []))
               and CT.get("proyeccion", {}).get("cabe") is True)
        # --- RECURSOS y coste (R2-c; no se vuelve a medir nada)
        RS = None
        if RECURSOS in actual:
            try:
                RS = json.load(open(os.path.join(run_dir, RECURSOS), encoding="ascii"))
            except ValueError:
                RS = None
        ok("RECURSOS: presente y legible", RS is not None)
        cs = R.get("coste_resumen") or {}
        if RS is not None:
            ok("RECURSOS: unidad y etapas = las de RESULTADO", set(RS) == {"unidad", "recursos", "etapas"}
               and RS.get("unidad") == cfg["unidad"] and RS.get("etapas") == R["etapas"])
            ok("RECURSOS: recursos = coste_resumen de RESULTADO", RS.get("recursos") == cs)
        ok("coste_resumen: techos = los de la configuracion", cs.get("techos") == cfg["techos"])
        T_ = cfg["techos"]
        caps = [("ancla_llamadas", T_["ancla"]["llamadas"]), ("c_tiempo_llamadas", T_["c_tiempo"]["llamadas"]),
                ("c_tiempo_pasos", T_["c_tiempo"]["pasos"]), ("cadenas_llamadas", T_["cadenas"]["llamadas"]),
                ("cadenas_pasos", T_["cadenas"]["pasos"]), ("bytes", T_["datos_bytes"])]
        over = [kk for kk, cap in caps if not (isinstance(cs.get(kk), int) and 0 <= cs[kk] <= cap)]
        ok("contadores registrados dentro de sus techos (reserva estricta; coherencia de lo registrado, no una "
           "nueva medida)", not over, over)
        if complete:
            meas = [("disco_pico", T_["disco_bytes"]), ("memoria_pico_principal", T_["memoria_proceso"]),
                    ("segundos", T_["global_s"] - T_["verificar_s"])]
            over = [kk for kk, cap in meas if not (isinstance(cs.get(kk), (int, float)) and cs[kk] <= cap)]
            ok("entrega completa: picos y segundos registrados dentro de los techos (lo registrado; no se vuelve "
               "a medir)", not over, over)
        pairs = [("ancla_llamadas", "ancla_llamadas"), ("cadenas_llamadas", "cadenas_llamadas"),
                 ("cadenas_pasos", "cadenas_pasos"), ("bytes", "bytes"), ("c_tiempo_pasos", "c_tiempo_pasos"),
                 ("c_tiempo_llamadas", "c_tiempo_llamadas")]
        if rd is None or sums.get("ancla_sin_recuentos"):
            pairs = pairs[1:]
        if complete:
            diff = [[kk, cs.get(kk), sums[sk]] for kk, sk in pairs if cs.get(kk) != sums[sk]]
            diff += [] if sums["bytes"] == sums["bytes_en_disco"] else [["bytes_en_disco", sums["bytes"],
                                                                         sums["bytes_en_disco"]]]
            cc = R.get("cadenas", {}).get("coste", {})
            diff += [] if (cc.get("llamadas"), cc.get("pasos"), cc.get("bytes_registro")) == (
                cs.get("cadenas_llamadas"), cs.get("cadenas_pasos"), cs.get("bytes")) else [["cadenas.coste", cc]]
            diff += [] if RA.get("coste", {}).get("llamadas") == cs.get("ancla_llamadas") else [["ancla.coste", 0]]
            ok("RECURSOS y contadores = sumas de los registros de ancla, C-tiempo y cadenas (entrega completa)",
               not diff, diff)
        else:
            diff = [[kk, cs.get(kk), sums[sk]] for kk, sk in pairs if not (isinstance(cs.get(kk), int)
                                                                          and cs[kk] >= sums[sk])]
            ok("RECURSOS y contadores >= sumas de los registros escritos (corrida incompleta)", not diff, diff)
        g("verificar: antes de dar el resultado")
    except (tmo.Agotado, tr.DatosNoCaben) as ex:
        rep["fallos"].append("techo de --verificar: %s" % ex)
        rep["lectura"] = LECTURAS[4]
        rep["segundos"] = round(clock() - t_in, 3)
        return 4, rep
    except Exception as ex:
        rep["fallos"].append("error: %s: %s" % (type(ex).__name__, ex))
        rep["lectura"] = LECTURAS[1]
        return 1, rep
    rep["segundos"] = round(clock() - t_in, 3)
    if rep["fallos"]:
        rep["lectura"] = LECTURAS[1]
        return 1, rep
    if not complete:
        rep["lectura"] = LECTURAS[2]
        return 2, rep
    rep["lectura"] = LECTURAS[0]
    return 0, rep


def verify_dest_problem(run_dir, dest, explicit):
    """--verificar no sobrescribe: el primer informe va a la carpeta (VERIFICACION_T2_2.json, fuera del
    manifiesto); uno posterior va con --salida a una ruta nueva fuera de la carpeta (dentro romperia el
    manifiesto). Devuelve el motivo para no ejecutar, o None."""
    if os.path.exists(dest):
        return "ya existe %s; use --salida <ruta nueva fuera de la carpeta de la corrida>" % dest
    if not os.path.isdir(os.path.dirname(os.path.abspath(dest))):
        return "no existe la carpeta de destino de %s" % dest
    if explicit:
        rd = os.path.abspath(run_dir)
        try:
            inside = os.path.commonpath([os.path.abspath(dest), rd]) == rd
        except ValueError:                       # otra unidad de disco: fuera
            inside = False
        if inside:
            return "--salida debe quedar fuera de la carpeta de la corrida"
    return None


def default_refs(cfg, model):
    if cfg["modelo"]["tipo"] == "fisico":
        return real_refs()
    return toy_refs(cfg, model)


# ------------------------------------------------------------------ preparacion del sello (no se ejecuta aqui)
def prepare_seal(path):
    if os.path.exists(path):
        raise SystemExit("la salida ya existe")
    import t2s_fisica as tf
    model = tf.make({"tipo": "fisico", "r0": "1/100"})
    lam_of = (lambda r: tm.lam_from_laws(model.laws(r)))
    rates = tm.build_ladder(lam_of)
    rec = tm.ladder_record(rates, lam_of)
    lvd = tmo.LevelData(model, rates)
    tr.write_json(path, {"unidad": "T2-2 escalera y claves", "escalera": rec, "shas": lvd.shas(),
                         "claves": lvd.kt.record(), "logaritmos": lvd.lt.record(),
                         "programa_sha256": code_shas()})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--corrida", action="store_true")
    ap.add_argument("--dir")
    ap.add_argument("--t0", type=float)
    ap.add_argument("--verificar")
    ap.add_argument("--preparar-sello", action="store_true")
    ap.add_argument("--salida")
    a = ap.parse_args()
    if a.verificar:
        dest = a.salida or os.path.join(a.verificar, "VERIFICACION_T2_2.json")
        why = verify_dest_problem(a.verificar, dest, a.salida is not None)
        if why:                                  # antes de calcular nada; 3 no es un resultado de verificacion
            print("verificar: no ejecutado: %s" % why, file=sys.stderr, flush=True)
            raise SystemExit(3)
        code, rep = verify(a.verificar)
        tr.write_json(dest, dict(rep, codigo=code))
        print("verificar: codigo %d (%s)" % (code, rep.get("lectura")), flush=True)
        for f in rep["fallos"][:5]:
            print("  fallo: %s" % f, flush=True)
        raise SystemExit(code)
    if a.preparar_sello:
        if not a.salida:
            raise SystemExit("--preparar-sello necesita --salida <nuevo json>")
        prepare_seal(a.salida)
        raise SystemExit(0)
    if not a.corrida or not a.dir or a.t0 is None:
        raise SystemExit("uso: --corrida --dir <nueva> --t0 <epoch> | --verificar <dir> [--salida <json>] "
                         "| --preparar-sello --salida <json>")
    try:
        seal = check_seal()                     # antes de crear nada ni importar el modelo fisico
    except CheckFail as ex:
        raise SystemExit("T2-2: %s" % ex)
    if os.path.exists(a.dir) and set(os.listdir(a.dir)) - set(POSTERIORES):
        raise SystemExit("--dir debe ser una carpeta nueva (solo puede contener la consola del supervisor)")
    os.makedirs(a.dir, exist_ok=True)
    esc = json.load(open(ESCALERA, encoding="ascii"))
    cfg = real_config(seal, esc)
    out = run(cfg, a.dir, lambda: tmo.make_model(cfg["modelo"]), lambda c, m: real_refs(), a.t0,
              prechecks=lambda model, chk: model.prechecks([tm.ENDPOINTS[0]] + list(tm.TARGETS),
                                                           cfg["semilla_maestra"], chk))
    raise SystemExit(0 if out["estado"] == "PASS" else 1)


if __name__ == "__main__":
    main()
