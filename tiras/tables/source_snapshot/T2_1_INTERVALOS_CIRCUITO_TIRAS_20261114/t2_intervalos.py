#!/usr/bin/env python3
"""T2-1 (BORRADOR_PLAN_T2_1.txt, N3): conductor de T2-1 en las tiras HGP(CAMINO_m, CAMINO_6), con k_log = 1.

Codigo preparado para revision (version 4). Version 2: D = 1e16 en la regla de suma de T2-1 (decision U6 del
usuario). Version 3, tras la revision de N1 v1 y N3 v2: R2 (los techos se consultan tambien en la lectura de N1,
en sus bucles largos y en la evaluacion; techo de disco), R3 (el Monte Carlo corre en un arbol de procesos propio,
un Job Object de Windows, con memoria total medida, plazo de trabajo y margen de cierre dentro de las 6.5 h y
cierre del arbol entero) y R4 (--verificar exige la entrega completa del contrato y recalcula D-t4 desde N2; una
inspeccion de parciales se etiqueta como tal y nunca devuelve exito). Version 4, tras la revision de N1 v2 y N3
v3: cierre de R3 (un fallo de una consulta o de una orden a Windows nunca se lee como cero procesos ni memoria
nula; el plazo se comprueba tambien cuando t2_mc.py sale por si mismo; un cierre forzado, no confirmado o fuera de
plazo deja D-mc NO_EVALUADO con la evidencia conservada). NO se ha ejecutado ninguna corrida. Solo se ejecutaron la
prueba de humo (--humo) y la prueba de cierre (--prueba-cierre), sin decoder ni DP. Para correr hace falta el
plan sellado (SELLO_PLAN.json con el SHA-256 del plan y del codigo) y una orden del usuario.

Contrato (N3 del borrador):
- Sin cambios en t2_motor.py, t2_calibracion.py ni b4_tiras.py (se importan tras comprobar su SHA-256). N1
  (t2_rival.py) da eps_ing y L_ing; N2 (t2_solo_u.py) se reutiliza sin repetir su DP. El resultado de la prueba
  acotada de SOLO_U se comprueba con su SHA-256, su verificador sin DP y el SHA-256 sellado de t2_solo_u.py.
- Sello: SELLO_PLAN.json debe dar el SHA-256 de PLAN_T2_1.txt y el de los cinco programas: t2_intervalos.py,
  t2_rival.py, t2_solo_u.py, t2_mc.py y t2_cp.py. Sin sello no hay corrida.
- Catalogos: los de T2-1M (trabajo/). Antes de usarlos se comprueban sus 1412 SHA-256 contra MANIFIESTO.json,
  cuyo propio SHA-256 tambien se comprueba. La DP vuelve a comprobar cada grupo al leerlo.
- DP: un trabajo de t2_motor.dp_job por tira, corte k = 1..4 y celda (modo celdas, 8 celdas con r_c en GRID_P y
  r_m en {r_c, 0}), con los odds del circuito (t2_motor.odds_circ). Regla de suma con k = 4 en m = 24 y 48
  (modo suma), evaluada en los dos extremos del intervalo de U, con D = D_SUMA = 1e16 en toda la comprobacion
  (U6): odds 1/D de la DP, X0 = D^NS U, T = prod(D + n_b), umbral D^(NS-k), cola, peso, margen y cota de
  C-anch. Con D = 1e14 la condicion previa es falsa en esas dos tiras. Las formulas de t2_motor.odds_sum y
  sum_interval y de t2_calibracion.criterion se copian aqui con D explicito (los modulos sellados leen D de
  t2_calibracion y no cambian); un control comprueba en cada corrida que la copia, con D = 1e14, da lo mismo
  que las funciones selladas. La condicion previa con D = 1e16 se comprueba antes de lanzar la DP.
  Anchura acumulada: los 384 bits fijan el redondeo de cada operacion, no la anchura final, que crece con el
  numero de terminos. Por eso se mide la anchura de U que devuelve la DP y se exige <= 1e-12 D^-k (control que
  detiene). Se informan su cociente con el umbral D^-k, la anchura de X0 sobre el umbral, la holgura
  1 - cola/umbral y una prevision antes de la DP con la anchura por termino de T2-1M.
  m = 6 y 12 conservan la procedencia de T2-0 (D = 1e14).
- Corrida 1 (m = 6, 12 y 24) y corrida 2 (m = 48), del mismo programa y con los mismos parametros. La corrida 2
  exige la lectura de la corrida 1 con la regla de paso LANZAR_CORRIDA_2: todos los controles PASS y D-mc
  COMPATIBLE. La corrida 2 vuelve a calcular, desde la evidencia guardada, las filas de la corrida 1 (deben
  coincidir) y decide D-t1, D-t2, D-t3 y D-t5 con las cuatro tiras.
- Certificado de D-mc (corrida 1):
  - se toma del resultado ya escrito en disco, de la fila TIRA_6x6, celda 1/1000|1/1000, k = 4;
  - se comprueban sus metadatos (trabajo, modo, corte, celda, completo, sin error, precision);
  - se contrasta con su origen: los extremos se recalculan desde el punto fijo y B guardados, y B desde las leyes;
  - el archivo lleva unidad, celda, k, los SHA-256 del decoder y de la propagacion, L_lo, H_hi y el SHA-256 del
    resultado de origen;
  - se valida con t2_mc.load_cert antes de lanzar el Monte Carlo.
  t2_mc.py corre despues como proceso aparte, con su techo de 1.5 h (plazo de trabajo), en un Job Object propio
  (R3): el proceso se crea suspendido, entra en el trabajo antes de reanudarse y todos sus descendientes quedan
  dentro. Cada 5 s como mucho se miden la memoria total (principal + todos los procesos del trabajo), la del
  principal y el disco; una medida que falla detiene el arbol, nunca cuenta como 0. Plazo duro del arbol: inicio +
  1.5 h + margen de cierre de 10 min, y nunca despues de las 6.5 h de la corrida; solo se lanza si caben. Plazo de
  aceptacion: el plazo duro menos una reserva de limpieza de 60 s, dentro de ese margen; ninguna espera pasa de el.
  Un techo, el plazo de aceptacion o un fallo de supervision terminan el trabajo entero (D-mc NO_EVALUADO). D-mc
  solo se acepta si la salida se observa dentro del plazo de aceptacion, el arbol queda vacio sin forzarlo (una
  consulta valida de 0 procesos activos), ninguna orden ni consulta del cierre falla y el cierre acaba antes del
  plazo duro; si no, NO_EVALUADO con el motivo, la salida y el recuento conservados y pendiente de inspeccion, sin
  reintento. Nunca se toca otro proceso.
- Paradas: toda parada (techo de tiempo, de memoria o de disco, error de un trabajo, error inesperado) deja la
  etapa en curso y las pendientes como NO_EVALUADO, con el motivo y el coste consumido (trabajos de la DP y avance
  de la lectura de N1 recibidos hasta la parada). Un control fallido deja FALLA. Nunca hay PASS parcial: sin todas
  las etapas PASS no hay decisiones ni regla de paso.
- --verificar (R4): sin DP ni catalogos. Verificacion integral: contrato de la entrega (unidad, corrida, programa
  y sello vigentes, dependencias, precision y D, estado y etapas PASS, controles, conjuntos exactos de instancias,
  celdas, cortes, trabajos, reglas de suma y decisiones; 96 filas en la corrida 1 y 128 en el resultado de la
  corrida 2), recalculo de todas las filas, N2 recalculado desde SOLO_U y decisiones con ese D-t4. Codigo 0 solo
  si todo cumple. --parcial: inspeccion de las filas presentes, etiquetada INSPECCION_PARCIAL, codigo 2.

Controles (detienen):
- SHA-256 del sello, del plan, del codigo y de las dependencias; trabajo/ contra su manifiesto.
- N1 (identidad con B4 en las 4 tiras fenomenologicas y eps_ing de S2), N2 (verificador de SOLO_U).
- Estructura de la DP: igual en las 8 celdas de cada (m, k) e igual a T2-0 o T2-1M donde existe.
- C-anc en todos los trabajos; C-anch: anchura relativa <= 1e-30 (celdas) y anchura acumulada de U medida
  <= 1e-12 D^-k (regla de suma, D = 1e16).
- Regla de suma: condicion previa cierta con D = 1e16 antes de la DP; copia con D explicito = funciones selladas
  con D = 1e14.
- Desigualdades de B4 certificadas con los extremos:
  - 0 <= eps_comp <= eps_ing;
  - L_ing <= L_comp (con L_cat como cota, solo se comprueba L_cat <= L_hi);
  - e_comp <= e_ing (solo con L_ing exacto);
  - monotonia en k de los dos rivales.

Metricas y decisiones: las del borrador, ya fijadas (A1-A4 y D-t1 a D-t5; k* con NO_DETERMINADO; presupuestos
de Monte Carlo por aproximacion normal; igual corte y proyeccion a igual tiempo, separadas).

Uso (desde D:\\paper ia cuantica\\qec_revision_matematica\\05_COMPROBACIONES):
  python -I -B T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_intervalos.py --corrida 1 --dir <carpeta nueva>
  python -I -B T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_intervalos.py --corrida 2 --dir <carpeta nueva>
      --previa <carpeta de la corrida 1>
  python -I -B T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_intervalos.py --verificar <RESULTADO json> [--parcial]
  python -I -B T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_intervalos.py --humo --dir <carpeta nueva>
  python -I -B T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_intervalos.py --prueba-cierre --dir <carpeta nueva>
"""
import argparse
import copy
import ctypes
import hashlib
import json
import math
import multiprocessing as mp
import os
import platform
import subprocess
import sys
import time
from fractions import Fraction as F

sys.set_int_max_str_digits(0)
T0 = time.perf_counter()
T0_WALL = time.time()
UDIR = "T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/"
EXE = "T2_1_PRUEBAS_ACOTADAS_EJECUCION_20261114/"
PLAN = UDIR + "PLAN_T2_1.txt"
SELLO = UDIR + "SELLO_PLAN.json"
SOLOU = EXE + "RESULTADO_SOLO_U.json"
MC_PATH = UDIR + "t2_mc.py"
CODE = ("t2_intervalos.py", "t2_rival.py", "t2_solo_u.py", "t2_mc.py", "t2_cp.py")
SHA = {SOLOU: "9c49caa981a927fcf9c68be250dad678fef8bfe83dd0b5379a3671ea3e9cd7e4"}


def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 24), b""):
            h.update(blk)
    return h.hexdigest()


for _p, _s in SHA.items():
    if sha_file(_p) != _s:
        raise SystemExit("SHA-256 distinto del registrado: " + _p)
sys.path.insert(0, os.path.abspath(UDIR))
import t2_rival as riv  # noqa: E402  (comprueba el motor, T2-1M, el manifiesto, B4, T2-0 y D1-CE)
import t2_solo_u as su  # noqa: E402  (comprueba el motor, el manifiesto, B4, D1-CE y FASE_A)
import t2_mc as mcm  # noqa: E402  (comprueba t2_calibracion, T2-0, b4_tiras y B4)

mot = riv.mot
t2 = mot.t2
b4 = mot.b4
check = b4.check

KMAX = 4
PROCS = 14
RUN_STRIPS = {1: (6, 12, 24), 2: (48,)}
SUM_STRIPS = (24, 48)
D_SUMA = 10 ** 16              # U6 (decidido por el usuario el 14/11/2026): D de la regla de suma de T2-1 (k = 4)
D_EXP = len(str(D_SUMA)) - 1
assert D_SUMA == 10 ** D_EXP
GIB = 1 << 30
SEG_C1 = 6.5 * 3600.0          # corrida 1, con el Monte Carlo
SEG_C2 = 9.0 * 3600.0
MARGEN_C1 = 900.0              # evaluacion, escritura y certificado antes del Monte Carlo
MC_CIERRE = 600.0              # R3: margen de cierre del arbol de t2_mc.py, DENTRO de las 6.5 h (no se suma)
MC_LIMPIEZA = 60.0             # R3 (v4): reserva para cerrar el arbol, DENTRO del margen de cierre (no se suma)
MEM_MAIN = 12 * GIB            # el total (24 GiB) lo vigila t2_motor.Monitor; 4 GiB por proceso, dentro de la DP
POLL = 5.0                     # segundos entre comprobaciones de techos mientras se espera a la DP
DISK_MAX = [60 * 10 ** 9]      # techo de disco: bytes en la carpeta de la corrida (trabajo/ solo se lee)
RUN_DIR = [None]
_DISK = {"t": 0.0, "bytes": 0, "max": 0}
CAP_SEGUNDOS = 15.0
N_TRABAJO = 1412
TARGETS = (F(1, 20), F(1, 10), F(1, 4))
A4_OBJ = F(1, 10)
Z95 = F(38416, 100)            # 1.96^2 / 0.1^2: 10 % al 95 %, aproximacion normal
K_PROY = 12
CERT_KEY = "%s|%s" % (mcm.RC, mcm.RM)
CERT_NAME = "TIRA_%dx6" % mcm.M_STRIP
DEADLINE = [T0_WALL + SEG_C1 - mcm.SEG_MAX - MARGEN_C1]
PLANNED = {1: ("sello", "capacidad", "trabajo", "N1_control", "N2_S2", "rival_circuito", "DP", "evaluacion",
               "llamadas"),
           2: ("sello", "previa", "capacidad", "trabajo", "rival_circuito", "DP", "evaluacion", "decisiones",
               "llamadas")}


def label(m, k, key):
    return "c%d_k%d_%s" % (m, k, key)


def sum_label(m, k):
    return "c%d_suma_k%d_D1e%d" % (m, k, D_EXP)     # el D queda en la etiqueta del trabajo


SUM_KEY = "suma_D1e%d" % D_EXP


# ------------------------------------------------------------------ regla de suma con D explicito
def odds_sum_D(cx, D):
    """t2_motor.odds_sum con D explicito: odds 1/D en todos los simbolos."""
    return {t * 8 + i: F(1, D) for t in range(len(cx["lts"])) for i in range(len(cx["CKrep"][t]))}


def precondition(nlist, k, D):
    """T = prod(D + n_b), umbral D^(NS-k) y cola = T - suma_(j<=k) e_j D^(NS-j) (enteros exactos)."""
    NS = len(nlist)
    e = t2.esym(nlist, k + 1)
    T = math.prod(D + nb for nb in nlist)
    thr = D ** (NS - k)
    return T, thr, T - sum(e[j] * D ** (NS - j) for j in range(k + 1))


def criterion_D(nlist, k, X0, D):
    """t2_calibracion.criterion con D explicito (misma formula y mismo control S >= 0)."""
    NS = len(nlist)
    T, thr, tail = precondition(nlist, k, D)
    pre = tail < thr
    S = T - X0
    check("regla de suma NS=%d k=%d: 0 <= T_cat - X0" % (NS, k), S >= 0, verbose=False)
    below = S < thr
    ver = ("COMPLETO" if below else "INCOMPLETO") if pre else "NO_DECIDIBLE"
    return {"k": k, "condicion_previa": pre, "cola_sobre_umbral": float(F(tail, thr)), "veredicto": ver,
            "margen_log10": (math.log10(thr) - math.log10(S)) if S > 0 else None,
            "peso_orden_siguiente_aprox": float(F(S, D ** (NS - k - 1)))}


def sum_interval_D(cx, k, vals, D):
    """t2_motor.sum_interval con D explicito; devuelve ademas la anchura de X0 (entero)."""
    S = 1 << mot.PREC
    Ul, Uh = vals
    DN = D ** len(cx["nclass"])
    X0_lo = (DN * Ul) >> mot.PREC
    X0_hi = -((-(DN * Uh)) >> mot.PREC)
    width = F(Uh - Ul, S)
    nl = cx["nclass"]
    return criterion_D(nl, k, X0_hi, D), criterion_D(nl, k, X0_lo, D), width, X0_hi - X0_lo


def replica_ok(cx, k, D=None):
    """La copia con D explicito = t2_motor.odds_sum, t2_motor.sum_interval y t2_calibracion.criterion con el D de
    t2_calibracion (D = None) o con el D dado (solo en la prueba de humo, cambiando t2.D y restaurandolo)."""
    D0 = t2.D
    try:
        if D is not None:
            t2.D = D
        nl = cx["nclass"]
        T, thr, tail = precondition(nl, k, t2.D)
        xs = [T - tail, T - thr // 2, T - 2 * thr, T]
        P = 1 << mot.PREC
        DN = t2.D ** len(nl)
        U = F(T - tail, DN)
        vals = [math.floor(U * P), math.ceil(U * P)]
        return (odds_sum_D(cx, t2.D) == mot.odds_sum(cx)
                and all(criterion_D(nl, k, x, t2.D) == t2.criterion(nl, k, x) for x in xs)
                and sum_interval_D(cx, k, vals, t2.D)[:3] == mot.sum_interval(cx, k, vals))
    finally:
        t2.D = D0


def preconditions_before_dp(strips, refs, D=D_SUMA):
    """Condicion previa (exacta) con D en las tiras con regla de suma, antes de lanzar la DP, y prevision de la
    anchura acumulada de U con la anchura por termino medida en T2-1M (informativa)."""
    ctl = refs["T21M"]["controles"]["C-T20"]
    per = max(ctl["circ_TIRA_%dx6_suma_k%d" % (m, k)]["anchura_U"] /
              refs["T21M"]["trabajos"]["circ_%d_suma_k%d" % (m, k)]["terminos"] for m, k in ((6, 4), (12, 4), (24, 3)))
    rows, ok = {}, True
    for m, (sc, _) in sorted(strips.items()):
        if m not in SUM_STRIPS:
            continue
        T, thr, tail = precondition(sc["cx"]["nclass"], KMAX, D)
        ok = ok and tail < thr
        est = est_terms(refs, m, KMAX)
        rows[str(m)] = {"cola_sobre_umbral": float(F(tail, thr)), "condicion_previa": tail < thr,
                        "terminos_previstos": est, "anchura_U_prevista": per * est,
                        "cota_C_anch": float(mot.SUM_W * F(1, D ** KMAX))}
    check("N3 U6: condicion previa de la regla de suma cierta con D = 1e%d y k = %d en m = %s, antes de la DP %s"
          % (len(str(D)) - 1, KMAX, ", ".join(rows), {m: "%.3g" % r["cola_sobre_umbral"] for m, r in rows.items()}),
          ok, rows)
    return {"D": str(D), "k": KMAX, "por_tira": rows, "anchura_por_termino_T2_1M_max": per,
            "nota": "prevision informativa: terminos previstos x anchura por termino de T2-1M (D = 1e14); la que "
                    "decide es la anchura medida tras la DP"}


def grids():
    ce = b4.load_all(with_b1=False)[0]
    cells = [(rc, rm) for rc in ce.GRID_P for rm in (rc, F(0))]
    grid12 = [(p, q) for p in ce.GRID_P for q in (p, p / 2, F(0))]
    return cells, grid12


def view(lo, hi):
    return [mot.dec(lo, False), mot.dec(hi, True)]


def _enc(x):
    if isinstance(x, F):
        return str(x)
    raise TypeError("no serializable: %s" % type(x).__name__)


def canon(obj):
    """Forma JSON de un objeto (la que se escribe en disco), para comparar con lo leido de un resultado."""
    return json.loads(json.dumps(obj, ensure_ascii=True, default=_enc))


def jdump(path, obj):
    with open(path, "w", encoding="ascii", newline="\n") as fo:
        fo.write(json.dumps(obj, ensure_ascii=True, indent=1, default=_enc) + "\n")
    return sha_file(path)


# ------------------------------------------------------------------ etapas, techos y paradas
class Stage:
    """Registra una etapa. Techo -> NO_EVALUADO; control fallido -> FALLA; otro error -> NO_EVALUADO."""

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
        elif issubclass(et, b4.CheckFail):
            rec.update(estado="FALLA", motivo="control fallido: %s" % ev)
        elif issubclass(et, b4.Stop):
            rec.update(estado="NO_EVALUADO", motivo="techo alcanzado: %s" % ev)
        else:
            rec.update(estado="NO_EVALUADO", motivo="error: %s: %s" % (et.__name__, ev))
        self.out["etapas"].append(rec)
        if et is not None:
            print("  etapa %s: %s (%s)" % (self.name, rec["estado"], rec.get("motivo")), flush=True)
        return False


def close_stages(out, planned):
    done = {s["etapa"] for s in out["etapas"]}
    bad = next((s for s in out["etapas"] if s["estado"] != "PASS"), None)
    for nm in planned:
        if nm not in done:
            out["etapas"].append({"etapa": nm, "estado": "NO_EVALUADO",
                                  "motivo": "no alcanzada" + (" (%s: %s)" % (bad["etapa"], bad.get("motivo"))
                                                              if bad else "")})
    st = [s["estado"] for s in out["etapas"]]
    out["estado"] = "PASS" if all(x == "PASS" for x in st) else ("FALLA" if "FALLA" in st else "NO_EVALUADO")
    if out["estado"] != "PASS":
        out["decisiones"] = {"decision": "NO_EVALUADO", "motivo": "corrida no completa: " + (
            "%s: %s" % (bad["etapa"], bad.get("motivo")) if bad else "etapas no alcanzadas")}


def disk_now(force=False):
    """Bytes en la carpeta de la corrida (lo unico que escribe T2-1: la DP en modo celdas y suma no escribe y
    trabajo/ solo se lee), medidos como mucho cada POLL s; se guarda el maximo."""
    d = RUN_DIR[0]
    if d is None:
        return 0
    now = time.time()
    if force or now - _DISK["t"] >= POLL:
        tot = 0
        for root, _, fs in os.walk(d):
            for f in fs:
                try:
                    tot += os.path.getsize(os.path.join(root, f))
                except OSError:
                    pass
        _DISK.update(t=now, bytes=tot, max=max(_DISK["max"], tot))
    return _DISK["bytes"]


def guard(where):
    """Techos de las etapas: plazo, memoria del principal, memoria total (t2_motor.Monitor marca mot.OVER) y disco.
    Lo consultan la espera de la DP, la lectura de N1 (riv.GUARD, R2) y la evaluacion."""
    if time.time() > DEADLINE[0]:
        raise b4.Stop("techo de tiempo (%s)" % where)
    m = b4.mem_now()
    if m is not None and m > MEM_MAIN:
        raise b4.Stop("techo de memoria del proceso principal (%s)" % where)
    if mot.OVER[0]:
        raise b4.Stop("techo de memoria total: %s (%s)" % (mot.OVER[0], where))
    dk = disk_now()
    if dk > DISK_MAX[0]:
        raise b4.Stop("techo de disco: %d bytes en la carpeta de la corrida (%s)" % (dk, where))


def check_seal(sello=SELLO, plan=PLAN, code_dir=UDIR):
    if not os.path.exists(sello) or not os.path.exists(plan):
        raise b4.CheckFail("plan sin sellar: faltan %s o %s" % (sello, plan))
    S = json.load(open(sello, encoding="utf-8"))
    got = {c: sha_file(os.path.join(code_dir, c)) for c in CODE}
    ok = S.get("plan") == os.path.basename(plan) and S.get("sha256") == sha_file(plan) and S.get("codigo") == got
    check("N3 sello: SHA-256 del plan y de los %d programas = SELLO_PLAN.json" % len(CODE), ok,
          {"sello": S.get("codigo"), "actual": got})
    return {"plan": S["plan"], "sha256": S["sha256"], "codigo": S["codigo"], "sellado": S.get("sealed_at"),
            "sello_sha256": sha_file(sello)}


def w_spin(sec):
    t0 = time.perf_counter()
    n = 0
    while True:
        for _ in range(20000):
            n += 1
        dt = time.perf_counter() - t0
        if dt >= sec:
            return n / dt


def capacity_probe(total=CAP_SEGUNDOS):
    """Capacidad de la maquina (sin decoder): un proceso y luego PROCS a la vez. Informativa; no cambia nada."""
    one = w_spin(total / 3)
    with mp.get_context("spawn").Pool(PROCS) as pool:
        rates = pool.map(w_spin, [2 * total / 3] * PROCS)
    return {"segundos": total, "tasa_un_proceso": one, "tasa_total": sum(rates), "procesos": PROCS,
            "nucleos_efectivos": sum(rates) / one, "nota": "bucle aritmetico; informativo"}


def w_sha(task):
    path, sha = task
    return path, sha_file(path) == sha, os.path.getsize(path)


def verify_trabajo(pool, MAN, tags=None, expect=N_TRABAJO):
    tasks = [(os.path.join(riv.WDIR, tag, e["archivo"]), e["sha256"]) for tag in sorted(MAN)
             if tags is None or tag in tags for _, e in sorted(MAN[tag].items())]
    bad, nb = [], 0
    for path, ok, size in pool.imap_unordered(w_sha, tasks, chunksize=4):
        nb += size
        if not ok:
            bad.append(path)
        guard("comprobacion de trabajo/")
    check("N3 trabajo/: %d ficheros con el SHA-256 de MANIFIESTO.json (%s...)" % (
        len(tasks), riv.SHA[riv.WDIR + "MANIFIESTO.json"][:12]), not bad and len(tasks) == expect, bad[:3])
    return {"ficheros": len(tasks), "bytes": nb, "manifiesto_sha256": riv.SHA[riv.WDIR + "MANIFIESTO.json"]}


# ------------------------------------------------------------------ N2: SOLO_U reutilizado y D-t4
def n2_s2(refs, grid12, seal_code):
    res = su.verify_output(SOLOU)
    check("N3 N2: verificador sin DP de SOLO_U (48 celdas y cortes, sin fallos)",
          res["comprobaciones"] == 48 and not res["fallos"], res)
    d = json.load(open(SOLOU, encoding="ascii"))
    check("N3 N2: SOLO_U de t2_solo_u.py sellado, sin parada y con sus controles PASS",
          d["programa_sha256"] == seal_code["t2_solo_u.py"] and "detenido_en" not in d
          and all(c["pasa"] for c in d["controles"]), d.get("detenido_en"))
    eps = riv.s2_eps(refs, grid12)
    rows, hi_max = {}, F(0)
    for p, q in grid12:
        key = "%s|%s" % (p, q)
        lo, hi = (F(x) for x in d["instancias"]["S2_34_6x8"]["celdas"][key]["3"]["eps_comp_exacto"])
        if lo <= 0:
            raise b4.CheckFail("N3 N2: eps_comp(3) de S2 sin cota inferior positiva en " + key)
        e = eps[key][3]
        R = (e / hi, e / lo)
        hi_max = max(hi_max, R[1])
        rows[key] = {"eps_ing_3": riv.show(e), "eps_comp_3": view(lo, hi), "R_3": view(*R)}
    dec = "PASA" if hi_max < 2 else "FALLA"
    return {"resultado_solo_u": {"archivo": SOLOU, "sha256": SHA[SOLOU]}, "celdas": rows,
            "D_t4": {"decision": dec, "R_3_cota_superior_max": mot.dec(hi_max, True),
                     "regla": "PASA si R_3 < 2 (cota superior) en las 12 celdas de S2",
                     "lectura": "" if dec == "PASA" else "se detiene la lectura y se revisa antes de informar"}}


# ------------------------------------------------------------------ DP
def ref_structure(refs, m, k):
    RT = refs["T20"].get("TIRA_%dx6" % m)
    if RT and str(k) in RT["estructura"]:
        e = RT["estructura"][str(k)]
        return e["polimeros"], e["estados"], e["terminos"]
    T = refs["T21M"]["trabajos"].get("circ_%d_est_k%d" % (m, k))
    if T:
        return T["polimeros"], T["estados"], T["terminos"]
    return None


def est_terms(refs, m, k):
    r = ref_structure(refs, m, k)
    if r is not None:
        return r[2]
    r = ref_structure(refs, 24, k)
    return 2 * r[2] if r else 1


def make_specs(strips, refs, cells, deadline):
    S = []
    for m, (sc, firmas) in sorted(strips.items()):
        tag = sc["tag"]
        base = {"dir": os.path.join(riv.WDIR, tag), "manifiesto": refs["MAN"][tag], "NS": sc["NS"],
                "gsize": sc["gsize"], "ngroups": sc["ngroups"], "firmas": firmas, "exacto": False,
                "limite_grupos": None, "quitar": None, "plazo": deadline}
        for k in range(1, KMAX + 1):
            est = est_terms(refs, m, k)
            for rc, rm in cells:
                key = "%s|%s" % (rc, rm)
                S.append(dict(base, etiqueta=label(m, k, key), k=k, modo="celdas",
                              celdas=[{"clave": key, "odds": mot.odds_circ(sc["cx"], rc, rm)[0]}],
                              prioridad=6.0 * est))
        if m in SUM_STRIPS:
            S.append(dict(base, etiqueta=sum_label(m, KMAX), k=KMAX, modo="suma",
                          celdas=[{"clave": SUM_KEY, "odds": odds_sum_D(sc["cx"], D_SUMA)}],
                          prioridad=2.0 * est_terms(refs, m, KMAX)))
    return sorted(S, key=lambda s: -s["prioridad"])


def collect_jobs(results, guard_fn):
    """Recoge los trabajos; un error de un trabajo o un techo detiene la etapa (NO_EVALUADO, sin parciales)."""
    R = {}
    for r in results:
        R[r["etiqueta"]] = r
        if "error" in r:
            raise b4.Stop("trabajo %s: %s" % (r["etiqueta"], r["error"]))
        print("  %s: estados=%s terminos=%d polimeros=%d (%.1fs)" % (
            r["etiqueta"], r.get("estados"), r["terminos"], r["polimeros"], r["segundos"]), flush=True)
        guard_fn("DP")
    return R


def run_dp(specs, store):
    ctx = mp.get_context("spawn")
    print("== DP: %d trabajos en %d procesos" % (len(specs), PROCS), flush=True)
    R = {}
    with ctx.Pool(PROCS, maxtasksperchild=1) as pool:   # al salir (tambien por una parada) se terminan los procesos
        mon = mot.Monitor(pool, "dp")
        try:
            R = collect_jobs(_tee(_poll(pool.imap_unordered(mot.dp_job, specs, chunksize=1), guard), store), guard)
        finally:
            mon.stop()
    check("N3 DP: %d trabajos recibidos de %d" % (len(R), len(specs)), len(R) == len(specs))
    return R


def _poll(it, guard_fn):
    """Espera los resultados comprobando los techos cada POLL segundos (no solo al terminar un trabajo)."""
    while True:
        try:
            r = it.next(timeout=POLL)
        except mp.TimeoutError:
            guard_fn("DP, en espera")
            continue
        except StopIteration:
            return
        yield r


def _tee(it, store):
    for r in it:
        store[r["etiqueta"]] = job_store(r)
        yield r


def job_store(r):
    o = mot.job_summary(r)
    o.update({x: r[x] for x in ("etiqueta", "claves", "valores", "peak") if x in r})
    return o


# ------------------------------------------------------------------ evaluacion de una tira
def combine(c_hi, c_lo):
    v = [c_hi["veredicto"], c_lo["veredicto"]]
    if v == ["COMPLETO", "COMPLETO"]:
        return "COMPLETO"
    if v == ["INCOMPLETO", "INCOMPLETO"]:
        return "INCOMPLETO"
    return "NO_DECIDIBLE"


def rival_from_evidence(sc, ev, cells):
    scan = riv.scan_from_summary(ev)
    out = {}
    for rc, rm in cells:
        od, B = mot.odds_circ(sc["cx"], rc, rm)
        terms = riv.pb_terms(riv.groups_circ(sc["cx"], rc, rm), KMAX)
        out["%s|%s" % (rc, rm)] = {"B": B, "terminos_PB": terms, "por_k": riv.rival_values(scan, od, B, terms)}
    return out


def eval_strip(m, sc, rv, jobs, refs, cells, guard_fn=None):
    """Intervalos, rival, controles y metricas de una tira. jobs: trabajos guardados de la tira (celdas y suma).
    guard_fn (R2): techos consultados por celda en la corrida (en --verificar no hay techos)."""
    tag, cx = sc["tag"], sc["cx"]
    keys = ["%s|%s" % c for c in cells]
    o = {"m": m, "estructura": {}, "regla_de_suma": {}, "H_comp": {}, "celdas": {}}
    bad_anc, st_k = 0, {}
    for k in range(1, KMAX + 1):
        st = set()
        for key in keys:
            r = jobs[label(m, k, key)]
            check("N3 %s: trabajo %s completo, modo celdas, k = %d, celda %s" % (tag, label(m, k, key), k, key),
                  r.get("completo") is True and r["modo"] == "celdas" and r["k"] == k and r["claves"] == [key]
                  and "error" not in r, verbose=False)
            st.add((r["polimeros"], r["estados"], r["terminos"]))
            bad_anc += r["anclas_malas"]
        ref = ref_structure(refs, m, k)
        check("N3 %s k=%d: polimeros, estados y terminos iguales en las %d celdas%s" % (
            tag, k, len(keys), " y = T2-0/T2-1M %s" % list(ref) if ref else " (sin referencia anterior)"),
            len(st) == 1 and (ref is None or st == {ref}), sorted(st))
        st_k[k] = next(iter(st))
        o["estructura"][str(k)] = {"polimeros_estados_terminos": list(st_k[k]),
                                   "referencia": list(ref) if ref else None}
    sum_rules = {}
    if m in SUM_STRIPS:
        k = KMAX
        r = jobs[sum_label(m, k)]
        check("N3 %s: trabajo de la regla de suma completo, modo suma, k = %d, clave %s" % (tag, k, SUM_KEY),
              r.get("completo") is True and r["modo"] == "suma" and r["k"] == k and r["claves"] == [SUM_KEY]
              and "error" not in r, verbose=False)
        bad_anc += r["anclas_malas"]
        check("N3 %s regla de suma: la copia con D explicito = t2_motor.odds_sum, sum_interval y "
              "t2_calibracion.criterion con D = 1e14 (4 valores de X0)" % tag, replica_ok(cx, k))
        Ul, Uh = r["valores"][0]
        width = F(Uh - Ul, 1 << mot.PREC)
        lim = mot.SUM_W * F(1, D_SUMA ** k)
        # anchura acumulada medida (no basta la precision de 384 bits: crece con los terminos); antes del criterio
        check("N3 C-anch %s regla de suma k=%d, D = 1e%d: anchura acumulada de U medida %.3e <= 1e-12 D^-k = %.0e "
              "(%d terminos)" % (tag, k, D_EXP, float(width), float(lim), r["terminos"]), 0 <= width <= lim,
              float(width))
        c_hi, c_lo, width2, wx = sum_interval_D(cx, k, [Ul, Uh], D_SUMA)
        T, thr, tail = precondition(cx["nclass"], k, D_SUMA)
        check("N3 %s regla de suma k=%d: misma estructura que las celdas" % (tag, k),
              (r["polimeros"], r["estados"], r["terminos"]) == st_k[k] and width2 == width)
        v = combine(c_hi, c_lo)
        sum_rules[k] = v
        o["regla_de_suma"][str(k)] = {
            "D": str(D_SUMA), "veredicto": v, "extremos": [c_hi["veredicto"], c_lo["veredicto"]],
            "condicion_previa": [c_hi["condicion_previa"], c_lo["condicion_previa"]],
            "cola_sobre_umbral": c_hi["cola_sobre_umbral"],
            "margen_log10": [c_hi["margen_log10"], c_lo["margen_log10"]],
            "peso_orden_siguiente_aprox": [c_hi["peso_orden_siguiente_aprox"], c_lo["peso_orden_siguiente_aprox"]],
            "anchura_U": float(width), "anchura_U_sobre_umbral": float(width * D_SUMA ** k),
            "anchura_X0_sobre_umbral": float(F(wx, thr)), "holgura_sobre_umbral": float(1 - F(tail, thr)),
            "anchura_U_por_termino": float(width / max(1, r["terminos"])), "terminos": r["terminos"],
            "trabajo": sum_label(m, k),
            "nota": "umbral relativo de U = D^-k; con los dos extremos de acuerdo el veredicto es certificado"}
        print("  %s regla de suma k=%d (D = 1e%d): %s; anchura de U / umbral = %.3e" % (
            tag, k, D_EXP, v, float(width * D_SUMA ** k)), flush=True)
    check("N3 C-anc %s: ningun hijo con ancla <= la de su padre en los trabajos de la tira" % tag, bad_anc == 0,
          bad_anc)
    hc = {k: riv.hcomp(tag, k, sum_rules) for k in range(1, KMAX + 1)}
    o["H_comp"] = {str(k): {"veredicto": v[0], "procedencia": v[1]} for k, v in hc.items()}
    X, maxw, inv = {}, F(0), []
    fails = {"I2": [], "I5": [], "I6": [], "I3": []}
    for (rc, rm), key in zip(cells, keys):
        if guard_fn is not None:
            guard_fn("evaluacion de %s" % tag)
        B, terms, rk = rv[key]["B"], rv[key]["terminos_PB"], rv[key]["por_k"]
        Xc, rows = {}, {}
        for k in range(1, KMAX + 1):
            r = jobs[label(m, k, key)]
            vals = r["valores"][0]
            iv = mot.cell_intervals(B, vals)
            (el, eh), (Ll, Lh), (Pl, Ph) = iv["eps_comp"], iv["L_comp"], iv["P_0"]
            if not all(iv[x][0] <= iv[x][1] for x in iv):
                inv.append([key, k])
            w = max(mot.relw(*iv[x]) for x in iv)
            maxw = max(maxw, w)
            ell, h = Ll, min(F(1), Lh + eh)
            e_comp = (h - ell) / (h + ell) if ell > 0 else None
            ei, Lc = rk[k]["eps_ing"], rk[k]["L_cat"]
            exact = hc[k][0] == "COMPLETO"
            Li = Lc if exact else None
            e_ing = ei / (2 * Li + ei) if exact and Li > 0 else None
            e_bnd = ei / (2 * Lc + ei) if not exact and Lc > 0 else None
            R = (ei / eh if eh > 0 else None, ei / el if el > 0 else None)
            if not (el >= 0 and eh <= ei):
                fails["I2"].append([key, k])
            if not ((Li <= Ll) if exact else (Lc <= Lh)):
                fails["I5"].append([key, k])
            # I6. Con L_ing = 0 el rival no da cota inferior: su error relativo es 1 y la desigualdad se cumple;
            # con L_ing > 0, I5 da ell >= L_ing > 0 y e_comp esta definido.
            if exact and Li > 0 and not (e_comp is not None and e_comp <= e_ing):
                fails["I6"].append([key, k])
            Xc[k] = {"eps_lo": el, "eps_hi": eh, "L_lo": Ll, "L_hi": Lh, "ell": ell, "h": h, "e_comp": e_comp,
                     "eps_ing": ei, "L_cat": Lc, "L_ing": Li, "e_ing": e_ing, "e_ing_cota": e_bnd, "R_lo": R[0],
                     "R_hi": R[1], "hcomp": hc[k][0]}
            row = {"trabajo": label(m, k, key), "valores_punto_fijo": vals, "precision_bits": mot.PREC,
                   "eps_comp": view(el, eh), "L_comp": view(Ll, Lh), "P_0": view(Pl, Ph),
                   "intervalo_P_sil": view(ell, h), "anchura_relativa_max": float(w),
                   "e_comp": mot.dec(e_comp, True) if e_comp is not None else None,
                   "estimacion_central": mot.dec(2 * ell * h / (ell + h), False) if ell > 0 else None,
                   "eps_ing": riv.show(ei), "H_comp": hc[k][0],
                   "L_ing": riv.show(Li) if exact else "PENDIENTE",
                   "L_cat_cota_inferior": None if exact else riv.show(Lc),
                   "e_ing": riv.show(e_ing) if e_ing is not None else ("PENDIENTE" if not exact else None),
                   "e_ing_cota_superior_con_L_cat": mot.dec(e_bnd, True) if e_bnd is not None else None,
                   "R": [mot.dec(R[0], False) if R[0] is not None else None,
                         mot.dec(R[1], True) if R[1] is not None else None]}
            if k == KMAX:
                row["intervalo_P_sil_exacto"] = [str(ell), str(h)]
            rows[str(k)] = row
        for k in range(1, KMAX):
            a, b = Xc[k], Xc[k + 1]
            if not (a["L_hi"] <= b["L_lo"] and b["eps_hi"] <= a["eps_lo"] and b["eps_ing"] <= a["eps_ing"]
                    and a["L_cat"] <= b["L_cat"]):
                fails["I3"].append([key, k])
        o["celdas"][key] = {"r_c": str(rc), "r_m": str(rm), "B": str(B), "por_k": rows, "k_estrella": kstar(Xc),
                            "muestras_10_por_ciento_95_aprox_normal": budgets(Xc, terms)}
        X[key] = Xc
    check("N3 %s: extremo inferior <= superior en los intervalos de eps_comp, L_comp y P_0 (%d filas)"
          % (tag, len(keys) * KMAX), not inv, inv[:4])
    check("N3 C-anch %s: anchura relativa de eps_comp, L_comp y P_0 <= 1e-30 en las %d filas (maxima %.3e)"
          % (tag, len(keys) * KMAX, float(maxw)), maxw <= mot.REL_W, float(maxw))
    check("N3 I2 %s: 0 <= eps_comp <= eps_ing (extremos) en las %d filas" % (tag, len(keys) * KMAX),
          not fails["I2"], fails["I2"][:4])
    check("N3 I5 %s: L_ing <= L_comp (L_ing exacto: <= L_lo; con L_cat como cota: L_cat <= L_hi)" % tag,
          not fails["I5"], fails["I5"][:4])
    check("N3 I6 %s: e_comp <= e_ing donde L_ing es exacto" % tag, not fails["I6"], fails["I6"][:4])
    check("N3 I3 %s: L no decrece y eps no crece con k en los dos rivales (extremos)" % tag, not fails["I3"],
          fails["I3"][:4])
    o["anchura_relativa_max"] = float(maxw)
    return o, X


def kstar(Xc):
    out = {}
    for t in TARGETS:
        kc = next((k for k in range(1, KMAX + 1) if Xc[k]["e_comp"] is not None and Xc[k]["e_comp"] <= t), None)
        ki, via = None, None
        for k in range(1, KMAX + 1):
            if Xc[k]["e_ing"] is not None and Xc[k]["e_ing"] <= t:
                ki, via = k, "L_ing exacto"
                break
            if Xc[k]["e_ing"] is None and Xc[k]["e_ing_cota"] is not None and Xc[k]["e_ing_cota"] <= t:
                ki, via = k, "cota superior con L_cat"
                break
        out[str(t)] = {"componentes": kc if kc is not None else "NO_DETERMINADO",
                       "ingenua": ki if ki is not None else "NO_DETERMINADO", "ingenua_via": via}
    return out


def budgets(Xc, terms):
    ell = Xc[KMAX]["ell"]
    d = {"monte_carlo_directo": math.ceil(Z95 * (1 - ell) / ell) if ell > 0 else None,
         "monte_carlo_directo_con": "L_lo de L_comp(4)"}
    L3 = Xc[3]["L_ing"] if Xc[3]["L_ing"] is not None else Xc[3]["L_cat"]
    pi = L3 / terms[3] if terms[3] > 0 else F(0)
    d["muestreo_por_subconjuntos"] = {
        "estrato_sitios": riv.F_S, "pi": float(pi), "muestras": math.ceil(Z95 * (1 - pi) / pi) if pi > 0 else None,
        "L_3": "exacto" if Xc[3]["L_ing"] is not None else "L_cat (cota inferior: muestras como cota superior)"}
    d["nota"] = "aproximacion normal; no es cobertura finita ni coste minimo; no decide"
    return d


def strip_costs(m, sc, rvres, jobs, cells):
    keys = ["%s|%s" % c for c in cells]
    out = {"certificado_por_k": {}, "rival": {}}
    for k in range(1, KMAX + 1):
        rs = [jobs[label(m, k, key)] for key in keys]
        secs = [r["segundos"] for r in rs]
        out["certificado_por_k"][str(k)] = {
            "polimeros": rs[0]["polimeros"], "estados": rs[0]["estados"], "terminos": rs[0]["terminos"],
            "segundos_max": max(secs), "segundos_medio": sum(secs) / len(secs),
            "pico_bytes_max": max(r.get("peak") or 0 for r in rs)}
    out["rival"] = {"truncacion_ingenua_decodificaciones": rvres["truncacion_ingenua"],
                    "lectura_catalogo_y_P3C": rvres["lectura"]["coste"], "leyes_por_celda": rvres["leyes"],
                    "clases": rvres["contexto"]["coste_clases"]}
    return out


def equal_time(m, sc, rvres, jobs, refs, X, cells):
    """Proyeccion a igual tiempo (informativa, no decide)."""
    tag, cx = sc["tag"], sc["cx"]
    cat = refs["T21M"]["instancias"][tag]["catalogo"]
    t_cat = cat["segundos_catalogo"] + cat["segundos_fusion"]
    t_dp = max(jobs[label(m, KMAX, "%s|%s" % c)]["segundos"] for c in cells)
    T = t_cat + t_dp
    pf = cat["por_f"][str(KMAX)]
    rate = pf["asignaciones"] / pf["segundos"]
    # leyes de una celda (favorable al rival: el certificado tambien se cuenta para una celda)
    t_law = rvres["contexto"]["coste_clases"]["segundos_clases"] + rvres["leyes"]["segundos"] / len(cells)
    e = t2.esym(cx["nclass"], K_PROY)
    kp, cum = 0, 0
    for k in range(1, K_PROY + 1):
        cum += e[k]
        if cum / rate + t_law <= T:
            kp = k
        else:
            break
    out = {"segundos_certificado_k4": T, "segundos_catalogo_T2_1M": t_cat, "segundos_dp_k4_una_celda": t_dp,
           "tasa_decoder_por_s": rate, "tasa_fuente": "catalogo de T2-1M con f = 4 (asignaciones / segundos, "
                                                         "14 procesos)", "segundos_leyes_una_celda": t_law,
           "k_prima": kp,
           "regla": "k' = mayor corte con (suma_(f<=k') e_f(n)) / tasa + clases y leyes de una celda <= tiempo "
                    "del certificado con k = 4 (catalogo de T2-1M + DP de una celda)",
           "nota": "proyeccion de coste, no el mayor corte alcanzable demostrado; informativa", "celdas": {}}
    for rc, rm in cells:
        key = "%s|%s" % (rc, rm)
        Xc = X[key]
        if kp == 0:
            v = None
        elif kp <= KMAX:
            ee = Xc[kp]["e_ing"]
            v = {"e_ing": riv.show(ee) if ee is not None else "PENDIENTE"}
        else:
            ep = riv.eps_from_terms(riv.pb_terms(riv.groups_circ(cx, rc, rm), kp), kp)
            if Xc[KMAX]["L_ing"] is None:
                v = {"eps_ing": riv.show(ep), "e_ing_cota_inferior": "PENDIENTE (sin L_ing(4) exacto)"}
            else:
                den = 2 * (Xc[KMAX]["L_ing"] + Xc[KMAX]["eps_ing"] - ep) + ep
                v = {"eps_ing": riv.show(ep), "e_ing_cota_inferior": mot.dec(ep / den, False)}
        out["celdas"][key] = {"e_comp_4": mot.dec(Xc[KMAX]["e_comp"], True) if Xc[KMAX]["e_comp"] else None,
                              "rival_en_k_prima": v}
    return out


# ------------------------------------------------------------------ decisiones (fijadas en el borrador)
def decide(X, dt4, final):
    ms = sorted(X)
    D = {"tiras": ms, "final": final}
    t1 = {}
    for m in ms:
        lo = [X[m][key][KMAX]["R_lo"] for key in X[m]]
        hi = [X[m][key][KMAX]["R_hi"] for key in X[m]]
        if any(v is None for v in lo + hi):
            t1[str(m)] = {"categoria": "NO_EVALUADO", "motivo": "R_4 sin extremos"}
            continue
        rel_y, rel_n = all(x >= 10 for x in lo), any(x < 10 for x in hi)
        nul_y, nul_n = any(x < 2 for x in hi), all(x >= 2 for x in lo)
        if not (rel_y or rel_n) or not (nul_y or nul_n):
            cat = "NO_DECIDIDA_POR_REDONDEO"
        else:
            cat = "RELEVANTE" if rel_y else ("NULA" if nul_y else "ZONA_GRIS")
        t1[str(m)] = {"categoria": cat, "celdas_R_mayor_o_igual_10": sum(x >= 10 for x in lo),
                      "celdas_R_menor_que_2": sum(x < 2 for x in hi),
                      "R_4": {key: [mot.dec(X[m][key][KMAX]["R_lo"], False), mot.dec(X[m][key][KMAX]["R_hi"], True)]
                              for key in X[m]}}
    steps, grow, nogrow, undec, slopes = [], 0, 0, 0, {}
    for a, b in zip(ms, ms[1:]):
        if b != 2 * a:
            continue
        steps.append([a, b])
        for key in X[a]:
            ra, rb = X[a][key][KMAX], X[b][key][KMAX]
            if None in (ra["R_lo"], ra["R_hi"], rb["R_lo"], rb["R_hi"]):
                undec += 1
                continue
            if rb["R_lo"] > ra["R_hi"]:
                grow += 1
            elif rb["R_hi"] <= ra["R_lo"]:
                nogrow += 1
            else:
                undec += 1
            slopes.setdefault("%d->%d" % (a, b), {})[key] = [math.log2(rb["R_lo"] / ra["R_hi"]),
                                                             math.log2(rb["R_hi"] / ra["R_lo"])]
    ncells = len(next(iter(X.values()))) if X else 0
    total = 3 * ncells
    if final and len(steps) == 3:
        c = "CRECE_CON_N" if grow == total else ("NO_DECIDIDA_POR_REDONDEO" if nogrow == 0 else "NO_CRECE")
    else:
        c = "PARCIAL (faltan tiras)"
    D["D_t1"] = {"por_tira": t1, "crecimiento": c, "pasos": steps, "crece": grow, "no_crece": nogrow,
                 "indecidibles": undec, "pendientes_log2_informativas": slopes,
                 "regla": "RELEVANTE si R_4 >= 10 en las 8 celdas; NULA si alguna < 2; ZONA_GRIS si no. "
                          "CRECE_CON_N si R_4 crece en cada paso en las 8 celdas"}
    if 12 in X:
        cnt, tot, rows = 0, 0, {}
        for m in (24, 48):
            if m not in X:
                continue
            for key in X[m]:
                a, b = X[m][key][KMAX]["e_comp"], X[12][key][KMAX]["e_comp"]
                tot += 1
                if a is None or b is None or b == 0:
                    rows["%d %s" % (m, key)] = None
                    cnt += 1
                    continue
                q = a / b
                rows["%d %s" % (m, key)] = mot.dec(q, True)
                cnt += int(q >= 2)
        full = final and 24 in X and 48 in X
        D["D_t2"] = {"decision": ("ESTABLE" if cnt == 0 else "CRECE") if full else "PARCIAL (faltan tiras)",
                     "celdas_con_cociente_mayor_o_igual_2": cnt, "celdas": tot, "cocientes": rows,
                     "regla": "ESTABLE si e_comp(4,m)/e_comp(4,12) < 2 en las 8 celdas para m = 24 y 48"}
    else:
        D["D_t2"] = {"decision": "NO_EVALUADO", "motivo": "falta m = 12"}
    t3 = {}
    for m in ms:
        if any(X[m][key][KMAX]["hcomp"] != "COMPLETO" for key in X[m]):
            t3[str(m)] = {"decision": "PENDIENTE", "motivo": "H-comp no certificado con k = 4"}
            continue
        r = {}
        for t in TARGETS:
            n = 0
            for key in X[m]:
                x = X[m][key][KMAX]
                n += int(x["e_comp"] is not None and x["e_ing"] is not None and x["e_comp"] <= t < x["e_ing"])
            r[str(t)] = n
        t3[str(m)] = {"celdas_con_ventaja": r}
    D["D_t3"] = {"por_tira": t3, "regla": "celdas con e_comp(4) <= objetivo < e_ing(4); solo con L_ing exacto"}
    D["D_t4"] = dt4
    n5 = ok5 = 0
    for m in ms:
        for key, Xc in X[m].items():
            if F(key.split("|")[0]) <= F(1, 10 ** 4):
                n5 += 1
                x = Xc[KMAX]
                ok5 += int(x["ell"] > 0 and x["e_comp"] is not None and x["e_comp"] <= A4_OBJ)
    full5 = final and ms == [6, 12, 24, 48]
    D["D_t5"] = {"decision": ("CUMPLE" if ok5 == n5 == 24 else "NO_CUMPLE") if full5 else "PARCIAL (faltan tiras)",
                 "celdas": n5, "celdas_que_cumplen": ok5,
                 "regla": "CUMPLE si en las 24 celdas con r_c <= 1e-4 ell > 0 y e_comp(4) <= 10 %"}
    return D


# ------------------------------------------------------------------ certificado de D-mc y Monte Carlo
def write_certificate(res_path, cert_path, sc6):
    with open(res_path, "rb") as fh:
        raw = fh.read()
    sha = hashlib.sha256(raw).hexdigest()
    R = json.loads(raw.decode("ascii"))
    inst = R["instancias"][CERT_NAME]
    cell = inst["celdas"][CERT_KEY]
    row = cell["por_k"][str(mcm.K_CERT)]
    lab = label(mcm.M_STRIP, mcm.K_CERT, CERT_KEY)
    job = inst["trabajos"][lab]
    B = F(cell["B"])
    vals = row["valores_punto_fijo"]
    ok = (row["trabajo"] == lab and job["k"] == mcm.K_CERT and job["modo"] == "celdas" and job["claves"] == [CERT_KEY]
          and job.get("completo") is True and "error" not in job and job["valores"][0] == vals
          and row["precision_bits"] == mot.PREC and F(cell["r_c"]) == mcm.RC and F(cell["r_m"]) == mcm.RM
          and R.get("estado") == "PASS")
    check("N3 certificado: fila %s %s k=%d del resultado escrito, con trabajo completo y metadatos correctos"
          % (CERT_NAME, CERT_KEY, mcm.K_CERT), ok)
    iv = mot.cell_intervals(B, vals)
    ell, h = iv["L_comp"][0], min(F(1), iv["L_comp"][1] + iv["eps_comp"][1])
    B2 = mot.odds_circ(sc6["cx"], mcm.RC, mcm.RM)[1]
    check("N3 certificado: extremos recalculados desde el punto fijo y B de la fila = intervalo exacto guardado; "
          "B = leyes de la celda", [str(ell), str(h)] == row["intervalo_P_sil_exacto"] and B == B2)
    org = os.path.relpath(res_path).replace("\\", "/")
    cert = {"unidad": "T2-1", "instancia": CERT_NAME, "r_c": str(mcm.RC), "r_m": str(mcm.RM), "k": mcm.K_CERT,
            "b4_tiras_sha256": mcm.SHA[mcm.B4D + "b4_tiras.py"],
            "t2_calibracion_sha256": mcm.SHA[mcm.T20 + "t2_calibracion.py"],
            "L_lo": str(ell), "H_hi": str(h), "vista": view(ell, h),
            "origen": {"archivo": org, "sha256": sha},
            "fila": ["instancias", CERT_NAME, "celdas", CERT_KEY, "por_k", str(mcm.K_CERT)], "trabajo_dp": lab,
            "punto_fijo": vals, "precision_bits": mot.PREC}
    csha = jdump(cert_path, cert)
    try:
        c = mcm.load_cert(cert_path)
    except SystemExit as ex:
        raise b4.CheckFail("N3 certificado rechazado por t2_mc.load_cert: %s" % ex)
    check("N3 certificado: t2_mc.load_cert lo acepta (L_lo, H_hi = la fila)", c["lo"] == ell and c["hi"] == h)
    return {"archivo": os.path.relpath(cert_path).replace("\\", "/"), "sha256": csha, "origen": cert["origen"],
            "intervalo": view(ell, h)}


class _BasicLimit(ctypes.Structure):
    _fields_ = [("PerProcessUserTimeLimit", ctypes.c_longlong), ("PerJobUserTimeLimit", ctypes.c_longlong),
                ("LimitFlags", ctypes.c_uint32), ("MinimumWorkingSetSize", ctypes.c_size_t),
                ("MaximumWorkingSetSize", ctypes.c_size_t), ("ActiveProcessLimit", ctypes.c_uint32),
                ("Affinity", ctypes.c_size_t), ("PriorityClass", ctypes.c_uint32), ("SchedulingClass", ctypes.c_uint32)]


class _ExtLimit(ctypes.Structure):
    _fields_ = [("BasicLimitInformation", _BasicLimit), ("IoInfo", ctypes.c_ulonglong * 6),
                ("ProcessMemoryLimit", ctypes.c_size_t), ("JobMemoryLimit", ctypes.c_size_t),
                ("PeakProcessMemoryUsed", ctypes.c_size_t), ("PeakJobMemoryUsed", ctypes.c_size_t)]


class _Accounting(ctypes.Structure):
    _fields_ = [("TotalUserTime", ctypes.c_longlong), ("TotalKernelTime", ctypes.c_longlong),
                ("ThisPeriodTotalUserTime", ctypes.c_longlong), ("ThisPeriodTotalKernelTime", ctypes.c_longlong),
                ("TotalPageFaultCount", ctypes.c_uint32), ("TotalProcesses", ctypes.c_uint32),
                ("ActiveProcesses", ctypes.c_uint32), ("TotalTerminatedProcesses", ctypes.c_uint32)]


class _PidList(ctypes.Structure):
    _fields_ = [("NumberOfAssignedProcesses", ctypes.c_uint32), ("NumberOfProcessIdsInList", ctypes.c_uint32),
                ("ProcessIdList", ctypes.c_size_t * 1024)]


class JobQueryError(OSError):
    """R3 (v4): fallo de una consulta o de una orden al Job Object, o lista de procesos incompleta. No es una
    medida: nunca se lee como cero procesos ni como memoria nula."""


class _PMC(ctypes.Structure):
    _fields_ = [("cb", ctypes.c_uint32), ("PageFaultCount", ctypes.c_uint32),
                ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]


class JobTree:
    """R3: arbol de procesos propio en un Job Object de Windows. El proceso se crea suspendido, se asigna al trabajo
    antes de reanudarse (asi no puede crear hijos fuera) y sus descendientes quedan dentro (el trabajo no permite
    salir). Se cuentan y se miden sus procesos; TerminateJobObject los termina a todos. Solo actua sobre los
    procesos de este trabajo, nunca sobre otros del equipo. Con KILL_ON_JOB_CLOSE, si el conductor termina, Windows
    termina el arbol al cerrarse el manejador.
    v4: cada consulta y cada orden comprueba su retorno. Un fallo se registra en errors (funcion y GetLastError, o
    el NTSTATUS de NtResumeProcess) y lanza JobQueryError: nunca devuelve [] ni (0, 0). Una lista de PID incompleta
    (ERROR_MORE_DATA o menos PID que procesos asignados) tambien es un fallo. peak_commit es solo diagnostico: si
    falla devuelve None (no disponible), nunca 0."""

    def __init__(self):
        from ctypes import wintypes as W
        self.W = W
        self.errors, self.h, self.seen, self.p = [], None, set(), None
        k = self.k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        k.CreateJobObjectW.restype = W.HANDLE
        k.CreateJobObjectW.argtypes = [ctypes.c_void_p, W.LPCWSTR]
        k.SetInformationJobObject.argtypes = [W.HANDLE, ctypes.c_int, ctypes.c_void_p, W.DWORD]
        k.QueryInformationJobObject.argtypes = [W.HANDLE, ctypes.c_int, ctypes.c_void_p, W.DWORD, ctypes.c_void_p]
        k.AssignProcessToJobObject.argtypes = [W.HANDLE, W.HANDLE]
        k.TerminateJobObject.argtypes = [W.HANDLE, W.UINT]
        k.IsProcessInJob.argtypes = [W.HANDLE, W.HANDLE, ctypes.POINTER(W.BOOL)]
        k.OpenProcess.restype = W.HANDLE
        k.OpenProcess.argtypes = [W.DWORD, W.BOOL, W.DWORD]
        k.K32GetProcessMemoryInfo.argtypes = [W.HANDLE, ctypes.POINTER(_PMC), W.DWORD]
        k.CloseHandle.argtypes = [W.HANDLE]
        self.nt = ctypes.WinDLL("ntdll")
        self.nt.NtResumeProcess.argtypes = [W.HANDLE]
        self.nt.NtResumeProcess.restype = ctypes.c_long
        self.h = k.CreateJobObjectW(None, None)
        if not self.h:
            self._fail("CreateJobObjectW")
        info = _ExtLimit()
        info.BasicLimitInformation.LimitFlags = 0x2000          # JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
        if not k.SetInformationJobObject(self.h, 9, ctypes.byref(info), ctypes.sizeof(info)):
            err = ctypes.get_last_error()
            self.close()
            self._fail("SetInformationJobObject(JobObjectExtendedLimitInformation)", err)

    def _fail(self, fn, err=None, detail=None, kind="GetLastError"):
        err = ctypes.get_last_error() if err is None else err
        e = {"funcion": fn, "error": err, "tipo_error": kind}
        if detail:
            e["detalle"] = detail
        self.errors.append(e)
        raise JobQueryError(err, "%s fallo (%s %d)%s" % (fn, kind, err, "; " + detail if detail else ""))

    def start(self, cmd, stdout):
        p = subprocess.Popen(cmd, stdout=stdout, stderr=subprocess.STDOUT, creationflags=0x4)   # CREATE_SUSPENDED
        try:
            hp = int(p._handle)
            if not self.k32.AssignProcessToJobObject(self.h, hp):
                self._fail("AssignProcessToJobObject")
            inj = self.W.BOOL()
            if not self.k32.IsProcessInJob(hp, self.h, ctypes.byref(inj)):
                self._fail("IsProcessInJob")
            if not inj.value:
                self._fail("IsProcessInJob", 0, "el proceso no quedo en el trabajo")
            st = self.nt.NtResumeProcess(hp)
            if st != 0:
                self._fail("NtResumeProcess", st & 0xFFFFFFFF, kind="NTSTATUS")
        except BaseException:
            p.kill()                # sigue suspendido: aun no ha ejecutado nada ni creado hijos
            p.wait()
            raise
        self.p = p
        self.seen.add(p.pid)
        return p

    def pids(self):
        """PID de los procesos activos del trabajo. La lista debe estar completa; si no, JobQueryError."""
        buf = _PidList()
        if not self.k32.QueryInformationJobObject(self.h, 3, ctypes.byref(buf), ctypes.sizeof(buf), None):
            err = ctypes.get_last_error()
            self._fail("QueryInformationJobObject(JobObjectBasicProcessIdList)", err,
                       "lista parcial: mas de 1024 procesos" if err == 234 else None)
        n_in, n_as = int(buf.NumberOfProcessIdsInList), int(buf.NumberOfAssignedProcesses)
        if n_in != n_as or n_in > 1024:
            self._fail("QueryInformationJobObject(JobObjectBasicProcessIdList)", 0,
                       "lista incompleta: %d PID de %d procesos asignados" % (n_in, n_as))
        ids = [int(buf.ProcessIdList[i]) for i in range(n_in)]
        self.seen.update(ids)
        return ids

    def counts(self):
        """(procesos activos, procesos totales) del trabajo, medidos; un fallo lanza JobQueryError."""
        acc = _Accounting()
        if not self.k32.QueryInformationJobObject(self.h, 1, ctypes.byref(acc), ctypes.sizeof(acc), None):
            self._fail("QueryInformationJobObject(JobObjectBasicAccountingInformation)")
        return int(acc.ActiveProcesses), int(acc.TotalProcesses)

    def memory(self):
        """Working set de todos los procesos activos del trabajo (bytes), medido. Un proceso de la lista que ya
        salio (OpenProcess con ERROR_INVALID_PARAMETER) cuenta 0 porque ya no ocupa memoria; cualquier otro fallo
        lanza JobQueryError."""
        tot = 0
        for pid in self.pids():
            h = self.k32.OpenProcess(0x1000 | 0x0010, False, pid)
            if not h:
                err = ctypes.get_last_error()
                if err == 87:
                    continue
                self._fail("OpenProcess", err, "pid %d" % pid)
            try:
                pmc = _PMC()
                pmc.cb = ctypes.sizeof(_PMC)
                if not self.k32.K32GetProcessMemoryInfo(h, ctypes.byref(pmc), pmc.cb):
                    self._fail("K32GetProcessMemoryInfo", None, "pid %d" % pid)
                tot += int(pmc.WorkingSetSize)
            finally:
                self.k32.CloseHandle(h)
        return tot

    def peak_commit(self):
        """Solo diagnostico: pico de memoria comprometida del trabajo, o None si no se pudo consultar (nunca 0)."""
        info = _ExtLimit()
        if not self.k32.QueryInformationJobObject(self.h, 9, ctypes.byref(info), ctypes.sizeof(info), None):
            self.errors.append({"funcion": "QueryInformationJobObject(JobObjectExtendedLimitInformation)",
                                "error": ctypes.get_last_error(), "tipo_error": "GetLastError",
                                "diagnostico": True})
            return None
        return int(info.PeakJobMemoryUsed)

    def terminate(self):
        if not self.k32.TerminateJobObject(self.h, 1):
            self._fail("TerminateJobObject")

    def wait_empty(self, wait_s):
        """True si una consulta valida da 0 procesos activos antes de wait_s s; False si siguen activos. Un fallo
        de la consulta lanza JobQueryError: no acredita un arbol vacio."""
        t = time.time()
        while True:
            if self.counts()[0] == 0:
                return True
            if time.time() - t >= wait_s:
                return False
            time.sleep(0.1)

    def close(self):
        """Cierra el manejador; con KILL_ON_JOB_CLOSE, Windows termina lo que quede. False si CloseHandle falla."""
        if not self.h:
            return True
        ok = bool(self.k32.CloseHandle(self.h))
        if not ok:
            self.errors.append({"funcion": "CloseHandle(trabajo)", "error": ctypes.get_last_error(),
                                "tipo_error": "GetLastError"})
        self.h = None
        return ok


def mc_entorno():
    """Entorno en el que corre el arbol (R3: la interfaz solo se acredita en el entorno Windows probado)."""
    w = sys.getwindowsversion()
    return {"sistema": platform.platform(), "windows": [w.major, w.minor, w.build],
            "python": platform.python_version(), "bits": platform.architecture()[0],
            "interfaz": "Job Object (kernel32) y NtResumeProcess (ntdll); acreditada solo en este entorno"}


def close_tree(job, p, tree, normal, hard, clock):
    """R3 (v4): cierre del arbol. Devuelve (confirmado, problema o None). Nunca lanza una excepcion ordinaria: cada
    consulta u orden que falla se registra como problema y se sigue, de modo que la terminacion de emergencia,
    job.close() (con KILL_ON_JOB_CLOSE termina lo que quede) y la espera de t2_mc.py siempre se ejecutan.
    normal: t2_mc.py salio por si mismo; entonces el arbol debe quedar vacio sin forzarlo. Si hay que forzarlo, o
    si no se puede comprobar, se termina igual y queda como problema. Confirmado solo si una consulta valida da 0
    procesos activos al final."""
    probs = []

    def attempt(what, fn, *a):
        try:
            return fn(*a)
        except Exception as ex:
            probs.append("%s: %s: %s" % (what, type(ex).__name__, ex))
            return None

    confirmed = p is None and job is None
    try:
        if p is not None:
            force = not normal
            if normal:
                attempt("registro de PID tras la salida", job.pids)
                empty = attempt("comprobacion del arbol tras la salida", job.wait_empty,
                                max(0.0, min(10.0, hard - clock())))
                if empty is not True:           # False (vivos) o None (consulta fallida): no acredita un arbol vacio
                    force = True
                    if empty is False:
                        n = attempt("recuento de descendientes", lambda: job.counts()[0])
                        tree["descendientes_vivos_tras_la_salida"] = n
                        attempt("registro de PID de los descendientes", job.pids)
                        probs.append("t2_mc.py salio, pero %s descendientes seguian vivos y hubo que terminar el "
                                     "arbol por fuerza; salida y recuento conservados, pendiente de inspeccion" % n)
            else:
                tree["procesos_activos_al_detener"] = attempt("recuento al detener", lambda: job.counts()[0])
                attempt("registro de PID al detener", job.pids)
            tree["terminacion_forzada"] = force
            if force:
                t_k = clock()
                attempt("TerminateJobObject", job.terminate)
                tree["segundos_terminacion"] = round(clock() - t_k, 3)
            empty = attempt("confirmacion del cierre", job.wait_empty, max(5.0, hard - clock()))
            cnt = attempt("recuento final", job.counts)
            act, total = cnt if cnt else (None, None)        # None = no medido (nunca 0 por un fallo)
            tree.update(procesos_activos_al_final=act, procesos_totales=total)
            confirmed = empty is True and act == 0
            if empty is False:
                probs.append("el arbol sigue con procesos activos tras terminarlo")
        elif job is not None:                   # el arbol se creo, pero t2_mc.py no llego a arrancar
            cnt = attempt("recuento final", job.counts)
            tree["procesos_activos_al_final"] = cnt[0] if cnt else None
            confirmed = cnt is not None and cnt[0] == 0
    finally:
        if job is not None:
            try:
                tree["pico_memoria_comprometida_arbol_bytes"] = job.peak_commit()
            except Exception:
                tree["pico_memoria_comprometida_arbol_bytes"] = None
            try:
                if not job.close():
                    probs.append("CloseHandle del trabajo fallo")
            except Exception as ex:
                probs.append("cierre del manejador: %s: %s" % (type(ex).__name__, ex))
        if p is not None:
            try:
                p.wait(timeout=30)
            except subprocess.TimeoutExpired:
                p.kill()
                p.wait()
                probs.append("t2_mc.py seguia vivo al cerrar el trabajo; terminado con kill")
    tree["cierre_confirmado"] = confirmed
    tree["problemas_de_cierre"] = probs
    return confirmed, "; ".join(probs) or None


def mc_evidence(out_path, cert_path, mc_sha, p):
    """Evidencia de la salida de t2_mc.py (archivo, SHA-256, decision propia, recuento, Clopper-Pearson y coste) y
    la decision que se aceptaria con un cierre correcto. La evidencia se registra aunque el estado operativo sea
    NO_EVALUADO (R3, v4)."""
    if p is None:
        return {}, "NO_EVALUADO", "t2_mc.py no se ejecuto"
    if not os.path.exists(out_path):
        return {}, "NO_EVALUADO", "t2_mc.py no escribio su salida"
    ev = {"salida": {"archivo": os.path.relpath(out_path).replace("\\", "/"), "sha256": sha_file(out_path)}}
    try:
        M = json.load(open(out_path, encoding="ascii"))
    except Exception as ex:
        return ev, "NO_EVALUADO", "salida de t2_mc.py ilegible: %s: %s" % (type(ex).__name__, ex)
    if not isinstance(M, dict):
        return ev, "NO_EVALUADO", "la salida de t2_mc.py no es un objeto JSON"
    dm = M.get("D_mc") if isinstance(M.get("D_mc"), dict) else {}
    ev.update(decision_de_t2_mc=dm.get("decision"), recuento=M.get("recuento_archivo"),
              clopper_pearson=M.get("clopper_pearson"), coste=M.get("coste"))
    good = (M.get("programa_sha256") == mc_sha and M.get("modo") == "Monte Carlo"
            and (M.get("certificado") or {}).get("sha256") == sha_file(cert_path))
    dmc = dm.get("decision", "NO_EVALUADO")
    if not good:
        return ev, "NO_EVALUADO", "la salida no es de t2_mc.py sellado con este certificado"
    if "detenido_en" in M:
        return ev, "NO_EVALUADO", M["detenido_en"]
    if dmc == "COMPATIBLE" and p.returncode != 0:
        return ev, "NO_EVALUADO", "codigo de retorno distinto de 0 con COMPATIBLE"
    return ev, dmc, dm.get("motivo")


def run_mc(cert_path, out_path, log_path, mc_sha, end_wall, cmd=None, work_s=None, close_s=MC_CIERRE,
           mem_total=None, job_factory=None, clock=None):
    """t2_mc.py como proceso aparte en su arbol propio (JobTree).
    Plazos (R3; v4 separa la parada del trabajo y la limpieza, sin sumar tiempo al presupuesto):
    - plazo de trabajo = techo propio de t2_mc.py (mcm.SEG_MAX, contado desde su inicio);
    - plazo duro del arbol = min(inicio + trabajo + margen de cierre close_s, end_wall);
    - plazo de aceptacion = plazo duro - reserva de limpieza (min(MC_LIMPIEZA, close_s/2)), DENTRO del margen.
    Solo se lanza si trabajo + cierre caben antes de end_wall (6.5 h de la corrida 1). Cada espera dura como mucho
    POLL s y nunca pasa del plazo de aceptacion. En cada sondeo: memoria total (principal + todos los procesos del
    arbol), memoria del principal y disco de la corrida; una medida que falla es un fallo de supervision, nunca 0.
    Un techo, el plazo de aceptacion o un fallo de supervision terminan el arbol entero.
    D-mc (de la salida de t2_mc.py) solo se acepta si, ademas de sus comprobaciones: la salida se observo dentro del
    plazo de aceptacion, el arbol quedo vacio sin forzarlo (consulta valida de 0 procesos activos), ninguna orden ni
    consulta del cierre fallo y el cierre termino antes del plazo duro. Si no, D-mc NO_EVALUADO con el motivo; la
    salida y el recuento se conservan (pendiente de inspeccion) y no se reintenta nada.
    job_factory y clock solo los sustituye la prueba de cierre (procesos sinteticos)."""
    clock = clock or time.time
    job_factory = job_factory or JobTree
    work = mcm.SEG_MAX if work_s is None else work_s
    mem_total = mot.MEM_TOTAL if mem_total is None else mem_total
    remaining = end_wall - clock()
    if remaining < work + close_s:
        return {"decision": "NO_EVALUADO", "motivo": "tiempo restante (%.0f s) menor que el plazo de trabajo de "
                "t2_mc.py (%.0f s) mas el margen de cierre (%.0f s)" % (remaining, work, close_s)}
    cmd = cmd or [sys.executable, "-I", "-B", MC_PATH, "--output", out_path, "--certificado", cert_path]
    t = clock()
    hard = min(t + work + close_s, end_wall)
    clean = min(MC_LIMPIEZA, close_s / 2.0)
    stop_at = hard - clean
    tree = {"plazo_trabajo_s": work, "margen_cierre_s": close_s, "reserva_limpieza_s": clean,
            "plazo_aceptacion_s_desde_inicio": round(stop_at - t, 3), "plazo_duro_s_desde_inicio": round(hard - t, 3),
            "entorno": mc_entorno(),
            "nota": "Job Object propio; solo sus procesos. null = no medido (una medida que falla nunca es 0)"}
    reason, peak, peak_main, p, job, t_exit = None, None, None, None, None, None
    with open(log_path, "w", encoding="utf-8", newline="\n") as lf:
        try:
            try:
                job = job_factory()
                p = job.start(cmd, lf)
            except Exception as ex:
                reason = "no se pudo iniciar el arbol de t2_mc.py (no se ejecuto): %s: %s" % (type(ex).__name__, ex)
            if p is not None:
                tree["pid_t2_mc"] = p.pid
                while True:
                    left = stop_at - clock()
                    if left <= 0:
                        tree["parada_s_desde_inicio"] = round(clock() - t, 3)
                        reason = ("plazo de aceptacion del arbol de t2_mc.py vencido (trabajo %.0f s + cierre %.0f s "
                                  "- reserva de limpieza %.0f s; el plazo duro queda dentro del techo de la corrida)"
                                  % (work, close_s, clean))
                        break
                    try:
                        p.wait(timeout=min(POLL, left))
                        t_exit = clock()
                        break
                    except subprocess.TimeoutExpired:
                        pass
                    main = b4.mem_now()
                    if main is None:
                        raise JobQueryError(0, "K32GetProcessMemoryInfo del proceso principal fallo (no medido)")
                    tot = main + job.memory()
                    peak = tot if peak is None else max(peak, tot)
                    peak_main = main if peak_main is None else max(peak_main, main)
                    dk = disk_now()
                    if tot > mem_total:
                        reason = "techo de memoria total: %d bytes (principal + arbol de t2_mc.py)" % tot
                    elif main > MEM_MAIN:
                        reason = "techo de memoria del proceso principal: %d bytes" % main
                    elif dk > DISK_MAX[0]:
                        reason = "techo de disco: %d bytes en la carpeta de la corrida" % dk
                    if reason:
                        tree["parada_s_desde_inicio"] = round(clock() - t, 3)
                        break
                if t_exit is not None:
                    tree["salida_observada_s_desde_inicio"] = round(t_exit - t, 3)
                    if t_exit > stop_at:
                        reason = ("salida de t2_mc.py observada %.3f s despues del plazo de aceptacion; no habilita "
                                  "la corrida 2" % (t_exit - stop_at))
        except Exception as ex:
            reason = "supervision del arbol de t2_mc.py interrumpida: %s: %s" % (type(ex).__name__, ex)
        finally:
            confirmed, issue = close_tree(job, p, tree, t_exit is not None, hard, clock)
    end = clock()
    if peak is not None:
        mot.MEMPEAK["monte_carlo"] = max(mot.MEMPEAK.get("monte_carlo", 0), peak)
    tree.update(pids=sorted(job.seen) if job else [], pico_memoria_total_bytes=peak,
                pico_memoria_principal_bytes=peak_main, errores_api=list(job.errors) if job else [],
                segundos=round(end - t, 3), exceso_sobre_plazo_de_trabajo_s=round(max(0.0, end - t - work), 3),
                exceso_sobre_plazo_duro_s=round(max(0.0, end - hard), 3),
                fin_menos_techo_de_la_corrida_s=round(end - end_wall, 3))
    rec = {"arbol": tree, "segundos": round(end - t, 3), "salida_consola": os.path.relpath(log_path).replace("\\", "/"),
           "codigo": None if p is None else p.returncode}
    ev, dmc, why = mc_evidence(out_path, cert_path, mc_sha, p)
    rec.update(ev)
    op = "; ".join(x for x in (reason, issue) if x) or None        # la causa de la parada y los problemas del cierre
    if op is None and not confirmed:
        op = "cierre del arbol no confirmado"
    if op is None and end > hard:
        op = "cierre del arbol terminado %.3f s despues del plazo duro" % (end - hard)
    if op:
        rec.update(decision="NO_EVALUADO", motivo=op)
        if ev.get("salida"):
            rec["pendiente_de_inspeccion"] = True
            rec["nota_recepcion"] = ("regla operativa de recepcion (R3): la salida y el recuento se conservan y se "
                                     "inspeccionan; no se declara falso el recuento; sin reintento automatico")
    else:
        rec.update(decision=dmc, motivo=why)
    return rec


# ------------------------------------------------------------------ corridas
def dependencies():
    return {**mot.SHA, **riv.SHA, **su.SHA, **mcm.SHA, **SHA}


def base_out(corrida):
    return {"unidad": "T2-1", "corrida": corrida, "programa_sha256": sha_file(os.path.abspath(__file__)),
            "dependencias": dependencies(), "precision_bits": mot.PREC,
            "D_regla_de_suma_T2_1": str(D_SUMA), "D_procedencia_T2_0_y_T2_1M": str(t2.D), "etapas": [],
            "instancias": {}}


def strip_block(m, sc, rvres, jobs, refs, cells, guard_fn=None):
    o, X = eval_strip(m, sc, rvres["valores_eval"], jobs, refs, cells, guard_fn=guard_fn)
    o["contexto"] = {x: sc[x] for x in ("tag", "NS", "gsize", "ngroups")}
    o["rival"] = {"evidencia": riv.scan_summary(rvres["lectura"]) if "lectura" in rvres else rvres["evidencia"],
                  "leyes": rvres.get("leyes"), "truncacion_ingenua": rvres.get("truncacion_ingenua")}
    o["trabajos"] = jobs
    if "lectura" in rvres:
        o["coste_igual_corte"] = strip_costs(m, sc, rvres, jobs, cells)
        o["igual_tiempo"] = equal_time(m, sc, rvres, jobs, refs, X, cells)
    return o, X


def corrida(out, n, d, previa=None):
    DEADLINE[0] = T0_WALL + (SEG_C1 - mcm.SEG_MAX - MARGEN_C1 if n == 1 else SEG_C2)
    out["techo_tiempo_etapas_s"] = DEADLINE[0] - T0_WALL
    RUN_DIR[0] = d
    _DISK.update(t=0.0, bytes=0, max=0)
    out["techo_disco"] = {"bytes": DISK_MAX[0], "medida": "bytes en la carpeta de la corrida (trabajo/ solo se lee)"}
    cells, grid12 = grids()
    ce, rect, _ = b4.load_all(with_b1=False)
    refs = riv.load_refs()
    with Stage(out, "sello"):
        seal = check_seal()
        out["sello"] = seal
    X = {}
    if n == 2:
        with Stage(out, "previa"):
            prev = load_previa(previa, seal, refs, cells, ce, rect)
            out["corrida_1"] = prev["registro"]
            X.update(prev["X"])
            out["instancias"].update(prev["instancias"])
            dt4 = prev["D_t4"]
    with Stage(out, "capacidad"):
        out["capacidad"] = capacity_probe()
    ctx = mp.get_context("spawn")
    RV = {}
    with ctx.Pool(PROCS) as pool:
        mon = mot.Monitor(pool, "lectura")
        riv.GUARD[0] = guard            # R2: la lectura de N1 y sus bucles consultan los techos de la corrida
        riv.PROGRESS.clear()
        try:
            with Stage(out, "trabajo"):
                out["trabajo"] = verify_trabajo(pool, refs["MAN"])
            if n == 1:
                with Stage(out, "N1_control"):
                    riv.prov_checks(refs)
                    out["N1_control"] = riv.fen_control(pool, ce, rect, refs, DEADLINE[0], grid12)
                with Stage(out, "N2_S2"):
                    out["N2_S2"] = n2_s2(refs, grid12, seal["codigo"])
                    dt4 = out["N2_S2"]["D_t4"]
            with Stage(out, "rival_circuito"):
                for m in RUN_STRIPS[n]:
                    RV[m] = riv.circuit_rival(pool, ce, rect, m, refs, cells, DEADLINE[0])
                    guard("rival de m = %d" % m)
        except BaseException:   # evidencia del avance de la lectura de N1 hasta la parada (sin valores)
            out["N1_lectura_hasta_la_parada"] = copy.deepcopy(riv.PROGRESS)
            raise
        finally:
            riv.GUARD[0] = None
            mon.stop()
    store = {}
    with Stage(out, "DP"):
        strips = {m: (RV[m]["contexto"], RV[m]["lectura"]["firmas"]) for m in RUN_STRIPS[n]}
        out["regla_de_suma_antes_de_la_DP"] = preconditions_before_dp(strips, refs)
        try:
            run_dp(make_specs(strips, refs, cells, DEADLINE[0]), store)
        except BaseException:   # evidencia del coste consumido; ningun valor parcial entra en filas ni decisiones
            out["trabajos_hasta_la_parada"] = {lab: {x: v for x, v in r.items() if x != "valores"}
                                               for lab, r in sorted(store.items())}
            raise
        finally:
            out["trabajos_recibidos"] = len(store)
    with Stage(out, "evaluacion"):
        for m in RUN_STRIPS[n]:
            jobs = {lab: store[lab] for lab in store if lab.startswith("c%d_" % m)}
            RV[m]["valores_eval"] = RV[m]["valores"]
            guard("evaluacion de m = %d" % m)
            o, X[m] = strip_block(m, RV[m]["contexto"], RV[m], jobs, refs, cells, guard_fn=guard)
            out["instancias"][RV[m]["contexto"]["name"]] = o
        if n == 1:
            out["decisiones"] = decide(X, dt4, final=False)
    if n == 2:
        with Stage(out, "decisiones"):
            out["decisiones"] = decide(X, dt4, final=True)
    with Stage(out, "llamadas"):
        out["llamadas_decoder_N3"] = b4.CALLS[0]
        check("N3: ninguna llamada al decoder en el conductor (S2 reutilizado de SOLO_U)", b4.CALLS[0] == 0,
              b4.CALLS[0])


def load_previa(d1, seal, refs, cells, ce, rect):
    """Corrida 2: lectura de la corrida 1 (regla de paso) y filas de m <= 24 recalculadas desde su evidencia."""
    lec_p = os.path.join(d1, "LECTURA_T2_1_C1.json")
    res_p = os.path.join(d1, "RESULTADO_T2_1_C1.json")
    L = json.load(open(lec_p, encoding="ascii"))
    rs = sha_file(res_p)
    me = sha_file(os.path.abspath(__file__))
    check("N3 previa: lectura de la corrida 1 con regla de paso LANZAR_CORRIDA_2, del mismo programa y sello",
          L.get("regla_de_paso") == "LANZAR_CORRIDA_2" and L.get("resultado", {}).get("sha256") == rs
          and L.get("programa_sha256") == me and L.get("sello_sha256") == seal["sello_sha256"],
          {"regla": L.get("regla_de_paso")})
    R = json.load(open(res_p, encoding="ascii"))
    check("N3 previa: resultado de la corrida 1 PASS, del mismo programa", R.get("estado") == "PASS"
          and R.get("programa_sha256") == me and R.get("corrida") == 1)
    issues = contract(R, seal, cells)
    check("N3 previa: entrega de la corrida 1 completa (contrato de --verificar: instancias, celdas, cortes, "
          "trabajos, decisiones, sello y programa)", not issues, issues[:6])
    X, mism = recheck(R, refs, cells, ce, rect)
    check("N3 previa: las filas de la corrida 1 se reproducen desde su evidencia (punto fijo, B y rival)",
          not mism, mism[:4])
    return {"X": X, "instancias": R["instancias"], "D_t4": R["decisiones"]["D_t4"],
            "registro": {"resultado": {"archivo": res_p.replace("\\", "/"), "sha256": rs},
                         "lectura": {"archivo": lec_p.replace("\\", "/"), "sha256": sha_file(lec_p)},
                         "D_mc": L.get("D_mc", {}).get("decision")}}


def recheck(R, refs, cells, ce, rect):
    """Recalcula cada tira de un resultado desde su evidencia guardada (punto fijo de la DP, B de las leyes y
    silenciosas y pares del rival), con todos sus controles, y la compara con lo escrito. Sin DP ni catalogos."""
    X, mism = {}, []
    for name, o in sorted(R["instancias"].items()):
        m = o["m"]
        sc = riv.strip_context(ce, rect, m, False, refs["B4"])
        rv = rival_from_evidence(sc, o["rival"]["evidencia"], cells)
        o2, X[m] = eval_strip(m, sc, rv, o["trabajos"], refs, cells)
        o2 = canon(o2)
        # R4: union de las claves, para que una celda borrada (o anadida) tambien sea una discrepancia
        mism += [name + " " + key for key in sorted(set(o["celdas"]) | set(o2["celdas"]))
                 if o["celdas"].get(key) != o2["celdas"].get(key)]
        mism += [name + " " + x for x in ("m", "estructura", "regla_de_suma", "H_comp", "anchura_relativa_max")
                 if o.get(x) != o2[x]]
    return X, mism


def expected_labels(m, cells):
    labs = {label(m, k, "%s|%s" % c) for k in range(1, KMAX + 1) for c in cells}
    if m in SUM_STRIPS:
        labs.add(sum_label(m, KMAX))
    return labs


def contract(R, seal, cells):
    """R4: lo que debe tener una entrega completa de la corrida R["corrida"]. Devuelve la lista de problemas (vacia
    si cumple). No recalcula valores (eso lo hacen recheck, N2 y decide en verify_result)."""
    if not isinstance(R, dict):
        return ["la entrega no es un objeto JSON"]
    bad = []

    def need(cond, what):
        if not cond:
            bad.append(what)

    n = R.get("corrida")
    need(R.get("unidad") == "T2-1", "unidad distinta de T2-1")
    if n not in (1, 2):
        return bad + ["corrida distinta de 1 y 2"]
    me = sha_file(os.path.abspath(__file__))
    need(R.get("programa_sha256") == me, "programa_sha256 distinto del programa que verifica")
    need(seal is not None and R.get("sello") == seal, "sello distinto del sello vigente")
    need(seal is not None and seal.get("codigo", {}).get("t2_intervalos.py") == me,
         "el programa que verifica no es el sellado")
    need(R.get("dependencias") == dependencies(), "dependencias distintas")
    need(R.get("precision_bits") == mot.PREC, "precision distinta")
    need(R.get("D_regla_de_suma_T2_1") == str(D_SUMA), "D de la regla de suma distinto de 1e16")
    need(R.get("D_procedencia_T2_0_y_T2_1M") == str(t2.D), "D de la procedencia distinto de 1e14")
    need(R.get("estado") == "PASS" and "detenido_en" not in R, "estado distinto de PASS o con parada")
    et = R.get("etapas") or []
    need([s.get("etapa") for s in et] == list(PLANNED[n]) and all(s.get("estado") == "PASS" for s in et),
         "etapas distintas de las previstas o no PASS")
    ctl = R.get("controles") or []
    need(bool(ctl) and all(c.get("pasa") is True for c in ctl), "controles ausentes o con alguno fallido")
    need(R.get("llamadas_decoder_N3") == 0, "llamadas al decoder en el conductor")
    strips = (6, 12, 24) if n == 1 else (6, 12, 24, 48)
    inst = R.get("instancias") if isinstance(R.get("instancias"), dict) else {}
    names = ["TIRA_%dx6" % m for m in strips]
    need(sorted(inst) == sorted(names), "instancias %s (se esperan %s)" % (sorted(inst), names))
    keys = {"%s|%s" % c for c in cells}
    ks = {str(k) for k in range(1, KMAX + 1)}
    rows = 0
    for m, nm in zip(strips, names):
        o = inst.get(nm)
        if not isinstance(o, dict):
            continue
        need(o.get("m") == m, nm + ": m distinto")
        cel = o.get("celdas") if isinstance(o.get("celdas"), dict) else {}
        need(set(cel) == keys, nm + ": celdas distintas de las 8 previstas")
        for key in sorted(set(cel) & keys):
            pk = set((cel[key] or {}).get("por_k") or {})
            need(pk == ks, "%s %s: cortes distintos de 1..4" % (nm, key))
            rows += len(pk & ks)
        need(set(o.get("trabajos") or {}) == expected_labels(m, cells), nm + ": trabajos distintos de los previstos")
        need(set(o.get("regla_de_suma") or {}) == ({str(KMAX)} if m in SUM_STRIPS else set()),
             nm + ": regla de suma distinta de la prevista")
        need(set(o.get("H_comp") or {}) == ks, nm + ": H_comp sin los cortes 1..4")
        need("evidencia" in (o.get("rival") or {}), nm + ": sin evidencia del rival")
    nrows = len(strips) * len(keys) * KMAX
    need(rows == nrows, "filas: %d (se esperan %d)" % (rows, nrows))
    run = RUN_STRIPS[n]
    need(R.get("trabajos_recibidos") == sum(len(expected_labels(m, cells)) for m in run),
         "trabajos recibidos distintos de las especificaciones de la corrida")
    pb = R.get("regla_de_suma_antes_de_la_DP") or {}
    pt = pb.get("por_tira") or {}
    need(pb.get("D") == str(D_SUMA) and sorted(pt) == sorted(str(m) for m in run if m in SUM_STRIPS)
         and all((v or {}).get("condicion_previa") is True for v in pt.values()),
         "condicion previa de la regla de suma antes de la DP ausente o distinta")
    D = R.get("decisiones") if isinstance(R.get("decisiones"), dict) else {}
    need(all(x in D for x in ("tiras", "final", "D_t1", "D_t2", "D_t3", "D_t4", "D_t5")), "decisiones incompletas")
    need(D.get("final") is (n == 2) and D.get("tiras") == list(strips), "decisiones: final o tiras distintas")
    if n == 1:
        need(isinstance(R.get("N2_S2"), dict), "falta N2_S2")
    else:
        c1 = (R.get("corrida_1") or {}).get("resultado") or {}
        p1 = c1.get("archivo")
        if not p1 or not os.path.exists(p1):
            bad.append("corrida_1: resultado de la corrida 1 no encontrado")
        else:
            need(sha_file(p1) == c1.get("sha256"), "corrida_1: SHA-256 del resultado de la corrida 1 distinto")
            R1 = json.load(open(p1, encoding="ascii"))
            need(all(R1.get("instancias", {}).get(x) == inst.get(x) for x in names[:3]),
                 "corrida_1: instancias distintas de las del resultado de la corrida 1")
    return bad


_VC = {}


def verify_result(path, parcial=False, sello=SELLO, plan=PLAN):
    """R4, sin DP ni lectura de catalogos. Verificacion integral (por defecto): sello vigente, contrato completo de
    la entrega (contract), recalculo de todas las filas desde la evidencia (recheck), N2 recalculado desde SOLO_U y
    decisiones recalculadas con ese D-t4 (no con el D-t4 del propio resultado). ok solo si todo cumple.
    parcial=True: inspeccion de las filas presentes; se etiqueta INSPECCION_PARCIAL y nunca es ok."""
    try:
        R = json.load(open(path, encoding="ascii"))
    except Exception as ex:
        return {"tipo": "INSPECCION_PARCIAL" if parcial else "VERIFICACION_INTEGRAL", "ok": False, "filas": 0,
                "contrato": ["no se puede leer: %s: %s" % (type(ex).__name__, ex)], "discrepancias": []}
    cells, grid12 = grids()
    if not _VC:
        _VC["ce"], _VC["rect"], _ = b4.load_all(with_b1=False)
        _VC["refs"] = riv.load_refs()
    ce, rect, refs = _VC["ce"], _VC["rect"], _VC["refs"]
    out = {"tipo": "INSPECCION_PARCIAL" if parcial else "VERIFICACION_INTEGRAL", "contrato": [], "discrepancias": []}
    seal = None
    if not parcial:
        try:
            seal = check_seal(sello, plan)
        except Exception as ex:
            out["contrato"].append("sello: %s" % ex)
        out["contrato"] += contract(R, seal, cells)
    X = None
    try:
        X, mism = recheck(R, refs, cells, ce, rect)
        out["discrepancias"] += mism
    except Exception as ex:
        out["discrepancias"].append("recalculo de las filas: %s: %s" % (type(ex).__name__, ex))
    if not parcial:
        if X is None or seal is None or not isinstance(R, dict) or R.get("corrida") not in (1, 2):
            out["discrepancias"].append("N2 y decisiones no recalculados (sin filas, sin sello o sin corrida)")
        else:
            try:
                n2 = n2_s2(refs, grid12, seal["codigo"])
                if R["corrida"] == 1 and canon(n2) != R.get("N2_S2"):
                    out["discrepancias"].append("N2_S2 distinto del recalculado desde SOLO_U")
                D2 = decide(X, n2["D_t4"], final=R["corrida"] == 2)
                if canon(D2) != R.get("decisiones"):
                    out["discrepancias"].append("decisiones distintas de las recalculadas (D-t4 desde N2)")
            except Exception as ex:
                out["discrepancias"].append("N2 o decisiones: %s: %s" % (type(ex).__name__, ex))
    inst = R.get("instancias") if isinstance(R, dict) and isinstance(R.get("instancias"), dict) else {}
    out["filas"] = sum(len((o or {}).get("celdas") or {}) * KMAX for o in inst.values())
    out["ok"] = not parcial and not out["contrato"] and not out["discrepancias"]
    return out


def verify_cli(path, parcial=False, sello=SELLO, plan=PLAN):
    """Codigo de salida de --verificar: 0 solo con verificacion integral sin problemas; 1 si hay alguno; 2 siempre
    en una inspeccion parcial."""
    r = verify_result(path, parcial, sello, plan)
    if parcial:
        print("INSPECCION PARCIAL (no es una verificacion integral): %d filas presentes; discrepancias: %s"
              % (r["filas"], r["discrepancias"] or "ninguna"), flush=True)
        return 2, r
    print("VERIFICACION INTEGRAL SIN DP: %s; %d filas; contrato: %s; discrepancias: %s" % (
        "CUMPLE" if r["ok"] else "NO CUMPLE", r["filas"], r["contrato"] or "completo",
        r["discrepancias"] or "ninguna"), flush=True)
    return (0 if r["ok"] else 1), r


# ------------------------------------------------------------------ prueba de humo (sin decoder ni DP)
def synth_scan(sc, rng):
    cx = sc["cx"]
    data = [t for t, lt in enumerate(cx["lts"]) if lt[0] == "d"]

    def sig_of(nn):
        syms = sorted(t * 8 + rng.randrange(len(cx["CKrep"][t])) for t in (rng.choice(data) for _ in range(nn)))
        s = 0
        for sy in syms:
            s = (s << 6) | (sy + 1)
        return s

    sil = {}
    for f in (3, 4):
        for _ in range(4):
            s = sig_of(f)
            sil[s] = sil.get(s, 0) + rng.randint(1, 6)
    pairs = {}
    for sg in [s for s in sil if riv.nsym(s) == 3]:
        for _ in range(2):
            kk = (sg, sig_of(1))
            pairs[kk] = pairs.get(kk, 0) + rng.randint(1, 9)
    return {"silenciosas": sil, "pares": pairs,
            "silenciosas_por_f": {str(f): sum(v for s, v in sil.items() if riv.nsym(s) == f) for f in range(1, 5)},
            "f_s": 3, "pares_compatibles": sum(pairs.values()), "singletes_que_pasan_con_clase_trivial": 0,
            "singletes_que_no_pasan": 0, "singletes_silenciosos": 0, "entradas_de_singletes_benignos": 0,
            "polimeros_silenciosos_de_3_sitios": 0, "validas_por_f": [0, 0, 0, 0], "secciones": 0, "bytes": 0,
            "coste": {"segundos_lectura": 0.0, "segundos_pares": 0.0}, "firmas": sorted(sil)}


def synth_jobs(m, rv, refs, cells, sum_mode=None, sc=None):
    S = 1 << mot.PREC
    jobs = {}
    for k in range(1, KMAX + 1):
        st = ref_structure(refs, m, k) or (100 * k, 10 * k, 1000 * k)
        for rc, rm in cells:
            key = "%s|%s" % (rc, rm)
            B, v = rv[key]["B"], rv[key]["por_k"][k]
            ec = v["eps_ing"] / (k + 1)
            Lc = v["L_cat"] * (1 + F(1, 1 << 60))
            U, O = (1 - ec) / B, Lc / B
            E = U / 2
            vals = [math.floor(U * S), math.ceil(U * S), math.floor(E * S), math.ceil(E * S), math.floor(O * S),
                    math.ceil(O * S)]
            jobs[label(m, k, key)] = {"etiqueta": label(m, k, key), "k": k, "modo": "celdas", "claves": [key],
                                      "completo": True, "polimeros": st[0], "estados": st[1], "terminos": st[2],
                                      "anclas_malas": 0, "valores": [vals], "segundos": 1.0 + k, "peak": 0}
    if sum_mode:
        st = ref_structure(refs, m, KMAX) or (100 * KMAX, 10 * KMAX, 1000 * KMAX)
        jobs[sum_label(m, KMAX)] = {"etiqueta": sum_label(m, KMAX), "k": KMAX, "modo": "suma", "claves": [SUM_KEY],
                                    "completo": True, "polimeros": st[0], "estados": st[1], "terminos": st[2],
                                    "anclas_malas": 0, "valores": [synth_sum_vals(sc["cx"], sum_mode)],
                                    "segundos": 9.0, "peak": 0}
    return jobs


def synth_sum_vals(cx, mode, D=D_SUMA, k=KMAX):
    """[U_lo, U_hi] sinteticos de la regla de suma con D: COMPLETO (X0 = T - cola: catalogo completo hasta k),
    INCOMPLETO (falta una decoracion de 3 sitios: X0 = T - cola - D^(NS-3)), BORDE (el intervalo de X0 cruza
    T - umbral: los extremos discrepan) y ANCHO (anchura de U de 1e-70, por encima de 1e-12 D^-k)."""
    S = 1 << mot.PREC
    NS = len(cx["nclass"])
    T, thr, tail = precondition(cx["nclass"], k, D)
    DN = D ** NS
    X0 = {"COMPLETO": T - tail, "INCOMPLETO": T - tail - D ** (NS - 3), "BORDE": T - thr, "ANCHO": T - tail}[mode]
    U = F(X0, DN)
    lo, hi = math.floor(U * S), math.ceil(U * S)
    if mode == "BORDE":
        lo, hi = lo - 1, hi + 1
    if mode == "ANCHO":
        lo = math.floor((U - F(1, 10 ** 70)) * S)
    return [lo, hi]


def humo(d):
    import random
    os.makedirs(d)
    rng = random.Random(20261114)
    c0 = b4.CALLS[0]
    ce, rect, _ = b4.load_all(with_b1=False)
    ce0 = ce.CALLS[0]
    refs = riv.load_refs()
    cells, grid12 = grids()
    n = [0]

    def ok(desc, cond, det=None):
        n[0] += 1
        check("humo N3: " + desc, cond, det)

    def fails(fn, exc=b4.CheckFail):
        """Texto de la excepcion esperada ('' si no se produce); quita de CHECKS los controles del negativo."""
        nc = len(b4.CHECKS)
        try:
            fn()
        except exc as ex:
            del b4.CHECKS[nc:]
            return str(ex) or type(ex).__name__
        return ""

    # 1. sello
    plan = os.path.join(d, "PLAN_SINTETICO.txt")
    with open(plan, "w", encoding="ascii", newline="\n") as fh:
        fh.write("plan sintetico de la prueba de humo\n")
    code = {c: sha_file(os.path.join(UDIR, c)) for c in CODE}
    s_ok, s_bad = os.path.join(d, "SELLO_SINTETICO.json"), os.path.join(d, "SELLO_MALO.json")
    jdump(s_ok, {"plan": os.path.basename(plan), "sha256": sha_file(plan), "codigo": code})
    jdump(s_bad, {"plan": os.path.basename(plan), "sha256": sha_file(plan), "codigo": dict(code, **{CODE[1]: "0" * 64})})
    seal = check_seal(s_ok, plan)
    ok("sello sintetico aceptado; con un SHA-256 de codigo cambiado o sin sello se rechaza",
       seal["codigo"] == code and fails(lambda: check_seal(s_bad, plan))
       and fails(lambda: check_seal(os.path.join(d, "NO_EXISTE.json"), plan)))
    # 2. capacidad
    cap = capacity_probe(1.5)
    ok("prueba de capacidad (1.5 s en la prueba de humo; 15 s en la corrida)",
       cap["tasa_un_proceso"] > 0 and cap["nucleos_efectivos"] > 0, cap)
    # 3. trabajo/ (solo fen_TIRA_6x6 en la prueba de humo) y negativo
    with mp.get_context("spawn").Pool(4) as pool:
        tv = verify_trabajo(pool, refs["MAN"], tags={"fen_TIRA_6x6"}, expect=len(refs["MAN"]["fen_TIRA_6x6"]))
        man2 = copy.deepcopy(refs["MAN"])
        g0 = sorted(man2["fen_TIRA_6x6"])[0]
        man2["fen_TIRA_6x6"][g0]["sha256"] = "0" * 64
        neg = fails(lambda: verify_trabajo(pool, man2, tags={"fen_TIRA_6x6"}, expect=len(man2["fen_TIRA_6x6"])))
    ok("trabajo/: %d ficheros de fen_TIRA_6x6 verificados; un SHA-256 cambiado se detecta; el manifiesto lista "
       "%d ficheros en total" % (tv["ficheros"], N_TRABAJO), neg
       and sum(len(v) for v in refs["MAN"].values()) == N_TRABAJO)
    # 4. tiras sinteticas: m = 6 (H-comp previo) y m = 48 (regla de suma con k = 4 y D = 1e16, sintetica
    #    INCOMPLETO: L_cat queda como cota en k = 1..4)
    X, inst = {}, {}
    for m, mode in ((6, None), (48, "INCOMPLETO")):
        sc = sc_cache(ce, rect, refs, m)
        scan = synth_scan(sc, rng)
        ev = riv.scan_summary(scan)
        rv = rival_from_evidence(sc, ev, cells)
        jobs = synth_jobs(m, rv, refs, cells, sum_mode=mode, sc=sc)
        rvres = {"valores_eval": rv, "evidencia": ev}
        o, X[m] = strip_block(m, sc, rvres, jobs, refs, cells)
        inst[sc["name"]] = o
    sc6, sc48 = sc_cache(ce, rect, refs, 6), sc_cache(ce, rect, refs, 48)
    rv6 = rv_of(inst, X, 6, sc6, cells)
    ok("m = 6: H-comp COMPLETO por procedencia; L_ing exacto, e_ing y k* calculados",
       all(inst["TIRA_6x6"]["H_comp"][str(k)]["veredicto"] == "COMPLETO" for k in range(1, 5))
       and inst["TIRA_6x6"]["celdas"][CERT_KEY]["por_k"]["4"]["e_ing"] not in (None, "PENDIENTE"))
    o48 = inst["TIRA_48x6"]
    rs48 = o48["regla_de_suma"]["4"]
    ok("m = 48, k = 4, D = 1e16: condicion previa cierta en los dos extremos (cola/umbral %.3g); regla sintetica "
       "INCOMPLETO -> L_ing PENDIENTE y L_cat como cota en k = 1..4; D y anchura sobre el umbral registrados"
       % rs48["cola_sobre_umbral"],
       rs48["veredicto"] == "INCOMPLETO" and rs48["condicion_previa"] == [True, True] and rs48["D"] == str(D_SUMA)
       and rs48["anchura_U_sobre_umbral"] < 1e-12
       and all(o48["celdas"][CERT_KEY]["por_k"][str(k)]["L_ing"] == "PENDIENTE" for k in range(1, 5)), rs48)
    # 4b. regla de suma con D = 1e16: copia frente a las funciones selladas, condicion previa, veredictos y anchura
    sc24 = sc_cache(ce, rect, refs, 24)
    rep = {"%d,%d" % (m, k): [replica_ok(sc_cache(ce, rect, refs, m)["cx"], k),
                              replica_ok(sc_cache(ce, rect, refs, m)["cx"], k, D=D_SUMA)]
           for m in (6, 24, 48) for k in (3, 4)}
    ok("la copia con D explicito = t2_motor.odds_sum, sum_interval y t2_calibracion.criterion con D = 1e14 y "
       "(cambiando t2.D solo aqui y restaurandolo) con D = 1e16, en m = 6, 24, 48 y k = 3, 4",
       all(all(v) for v in rep.values()) and t2.D == 10 ** 14, rep)
    pb = preconditions_before_dp({24: (sc24, None), 48: (sc48, None)}, refs)
    f0 = fails(lambda: preconditions_before_dp({48: (sc48, None)}, refs, D=10 ** 14))
    ok("condicion previa antes de la DP: cierta con D = 1e16 en m = 24 y 48 (%s); con D = 1e14 en m = 48 detiene"
       % {m: "%.3g" % r["cola_sobre_umbral"] for m, r in pb["por_tira"].items()},
       all(r["condicion_previa"] for r in pb["por_tira"].values()) and f0.startswith("N3 U6"), [pb, f0])
    vv = {}
    for m, sc in ((24, sc24), (48, sc48)):
        for mode in ("COMPLETO", "INCOMPLETO", "BORDE"):
            c_hi, c_lo, w, wx = sum_interval_D(sc["cx"], KMAX, synth_sum_vals(sc["cx"], mode), D_SUMA)
            vv["%d %s" % (m, mode)] = combine(c_hi, c_lo)
    ok("regla de suma sintetica con D = 1e16 en m = 24 y 48: COMPLETO e INCOMPLETO en los dos extremos; un "
       "intervalo que cruza el umbral da NO_DECIDIBLE", all(v == x.split()[1].replace("BORDE", "NO_DECIDIBLE")
                                                           for x, v in vv.items()), vv)
    rv48 = rv_of(inst, X, 48, sc48, cells)
    jc = synth_jobs(48, rv48, refs, cells, sum_mode="COMPLETO", sc=sc48)
    oc, _ = eval_strip(48, sc48, rv48, jc, refs, cells)
    ok("m = 48 con regla de suma COMPLETO (D = 1e16): H-comp COMPLETO en k = 1..4 (k <= 3 cubiertos por k = 4) y "
       "L_ing exacto", all(oc["H_comp"][str(k)]["veredicto"] == "COMPLETO" for k in range(1, 5))
       and oc["H_comp"]["3"]["procedencia"] == "regla de suma de T2-1 con k = 4"
       and all(oc["celdas"][CERT_KEY]["por_k"][str(k)]["L_ing"] != "PENDIENTE" for k in range(1, 5)), oc["H_comp"])
    ja = synth_jobs(48, rv48, refs, cells, sum_mode="ANCHO", sc=sc48)
    fa = fails(lambda: eval_strip(48, sc48, rv48, ja, refs, cells))
    ok("una anchura acumulada de U de 1e-70 (> 1e-12 D^-4 = 1e-76) detiene la evaluacion en C-anch, antes del "
       "criterio", fa.startswith("N3 C-anch circ_TIRA_48x6 regla de suma"), fa)
    k0 = label(6, 3, CERT_KEY)
    bad = copy.deepcopy(synth_jobs(6, rv6, refs, cells))
    v = bad[k0]["valores"][0]
    v[4] = v[5] = v[4] // 2                                        # L_comp(3) = L_ing(3)/2
    bad2 = copy.deepcopy(synth_jobs(6, rv6, refs, cells))
    bad2[k0]["valores"][0][5] = bad2[k0]["valores"][0][4] - 1      # O_hi < O_lo
    f1 = fails(lambda: eval_strip(6, sc6, rv6, bad, refs, cells))
    f2 = fails(lambda: eval_strip(6, sc6, rv6, bad2, refs, cells))
    ok("intervalos sinteticos que violan L_ing <= L_comp o lo <= hi detienen la evaluacion con ese control",
       f1.startswith("N3 I5") and "extremo inferior <= superior" in f2, [f1, f2])
    ng = []

    def g_eval(where):
        ng.append(where)
        if len(ng) >= 3:
            raise b4.Stop("techo sintetico (%s)" % where)

    fe = _stop_text(lambda: eval_strip(6, sc6, rv6, synth_jobs(6, rv6, refs, cells), refs, cells, guard_fn=g_eval))
    ok("R2: la evaluacion consulta los techos por celda y se detiene (%d consultas)" % len(ng),
       "techo sintetico (evaluacion de circ_TIRA_6x6)" in fe and len(ng) == 3, fe)
    # 5. decisiones parciales
    D = decide(X, {"decision": "NO_EVALUADO", "motivo": "prueba de humo"}, final=False)
    ok("decisiones parciales con m = 6 y 48: D-t2 NO_EVALUADO (falta m = 12), D-t3 PENDIENTE en m = 48, D-t5 "
       "PARCIAL", D["D_t2"]["decision"] == "NO_EVALUADO" and D["D_t3"]["por_tira"]["48"]["decision"] == "PENDIENTE"
       and D["D_t5"]["decision"].startswith("PARCIAL"))
    # 6. resultado sintetico, certificado y negativo
    out = base_out(1)
    out.update(instancias=inst, decisiones=D, estado="PASS", nota="RESULTADO SINTETICO DE LA PRUEBA DE HUMO")
    res = os.path.join(d, "RESULTADO_SINTETICO.json")
    jdump(res, out)
    cert = os.path.join(d, "CERTIFICADO_SINTETICO.json")
    ci = write_certificate(res, cert, sc6)
    ok("certificado escrito desde la fila del resultado y aceptado por t2_mc.load_cert", ci["sha256"] == sha_file(cert))
    tam = os.path.join(d, "RESULTADO_SINTETICO_ALTERADO.json")
    R2 = json.load(open(res, encoding="ascii"))
    R2["instancias"]["TIRA_6x6"]["celdas"][CERT_KEY]["por_k"]["4"]["intervalo_P_sil_exacto"][0] = "0"
    jdump(tam, R2)
    ok("una fila alterada se rechaza al escribir el certificado",
       fails(lambda: write_certificate(tam, os.path.join(d, "CERT_ALTERADO.json"), sc6)))
    with open(res, "a", encoding="ascii", newline="\n") as fh:
        fh.write("\n")
    ok("si el resultado de origen cambia despues, t2_mc.load_cert rechaza el certificado",
       fails(lambda: mcm.load_cert(cert), SystemExit))
    # 7. Monte Carlo: subproceso sustituido (no se ejecuta t2_mc.py) y plazo de aceptacion
    mo, lg = os.path.join(d, "DMC_SINTETICO.json"), os.path.join(d, "SALIDA_DMC_SINTETICA.txt")
    mc_sha = sha_file(MC_PATH)
    fake = ("import json,sys; json.dump({'programa_sha256': sys.argv[2], 'modo': 'Monte Carlo', 'certificado': "
            "{'sha256': sys.argv[3]}, 'D_mc': {'decision': 'COMPATIBLE'}}, open(sys.argv[1], 'w')); "
            "print('salida sintetica')")
    r = run_mc(cert, mo, lg, mc_sha, time.time() + 120, cmd=[sys.executable, "-I", "-c", fake, mo, mc_sha,
                                                               sha_file(cert)], work_s=60, close_s=10)
    r2 = run_mc(cert, mo + "2", lg + "2", mc_sha, time.time() + 120,
                cmd=[sys.executable, "-I", "-c", "import time; time.sleep(30)"], work_s=2, close_s=1)
    r3 = run_mc(cert, mo + "3", lg + "3", mc_sha, time.time() + mcm.SEG_MAX + MC_CIERRE - 1)
    ok("Monte Carlo sustituido: COMPATIBLE leido y comprobado; plazo de aceptacion -> NO_EVALUADO; sin tiempo para "
       "el trabajo (1.5 h) mas el margen de cierre (10 min) no se lanza", r["decision"] == "COMPATIBLE"
       and r2["decision"] == "NO_EVALUADO" and "plazo de aceptacion" in r2["motivo"] and r3["decision"] == "NO_EVALUADO"
       and "margen de cierre" in r3["motivo"] and r["arbol"]["procesos_activos_al_final"] == 0,
       [r["decision"], r2.get("motivo"), r3.get("motivo")])
    fake2 = fake.replace("sys.argv[2]", "'otro'")
    r4 = run_mc(cert, mo + "4", lg + "4", mc_sha, time.time() + 120, cmd=[
        sys.executable, "-I", "-c", fake2, mo + "4", mc_sha, sha_file(cert)], work_s=60, close_s=10)
    ok("una salida de Monte Carlo de otro programa no cuenta (NO_EVALUADO)", r4["decision"] == "NO_EVALUADO")
    # 7b. R3 (v4): cierre del arbol de t2_mc.py con procesos sinteticos, fallos de Windows inyectados y reloj de
    #     prueba (prueba_cierre_r3; la misma que --prueba-cierre)
    for desc, cond, det in prueba_cierre_r3(d, cert, mc_sha):
        ok("R3: " + desc, cond, det)
    # 7c. R2 con procesos reales: la lectura de N1 (scan_catalogue) se detiene por el Monitor del motor (memoria
    #     total), por el plazo durante una espera sin resultados y por el disco; tambien la evaluacion.
    for desc, cond, det in humo_r2(d):
        ok(desc, cond, det)
    # 8. etapas y paradas normalizadas
    so = {"etapas": []}
    for exc in (b4.Stop("techo sintetico"), b4.CheckFail("control sintetico"), ValueError("error sintetico")):
        try:
            with Stage(so, "e_%s" % type(exc).__name__):
                raise exc
        except Exception:
            pass
    close_stages(so, ("e_Stop", "e_CheckFail", "e_ValueError", "pendiente"))
    st = [s["estado"] for s in so["etapas"]]
    ok("paradas: techo -> NO_EVALUADO, control -> FALLA, error -> NO_EVALUADO, pendiente -> NO_EVALUADO; "
       "estado FALLA y decisiones NO_EVALUADO", st == ["NO_EVALUADO", "FALLA", "NO_EVALUADO", "NO_EVALUADO"]
       and so["estado"] == "FALLA" and so["decisiones"]["decision"] == "NO_EVALUADO", st)
    jobs6 = list(synth_jobs(6, rv_of(inst, X, 6, sc6, cells), refs, cells).values())
    errj = dict(jobs6[0], error="RuntimeError: techo de tiempo de la unidad")
    calls = [0]

    def g2(where):
        calls[0] += 1
        if calls[0] >= 3:
            raise b4.Stop("techo sintetico (%s)" % where)

    class SlowIt:
        """Iterador como el de imap_unordered: dos esperas sin resultado y luego los trabajos."""

        def __init__(self, xs):
            self.xs, self.waits = list(xs), 2

        def next(self, timeout=None):
            if self.waits:
                self.waits -= 1
                raise mp.TimeoutError
            if not self.xs:
                raise StopIteration
            return self.xs.pop(0)

    seen = []
    got = list(_poll(SlowIt(jobs6[:3]), seen.append))
    ok("DP: un trabajo con error o un techo en la recogida detienen la etapa (sin parciales); en espera se "
       "comprueban los techos cada %.0f s" % POLL,
       fails(lambda: collect_jobs(iter(jobs6[:2] + [errj]), lambda w: None), b4.Stop)
       and fails(lambda: collect_jobs(iter(jobs6), g2), b4.Stop) and got == jobs6[:3] and len(seen) == 2
       and fails(lambda: list(_poll(SlowIt(jobs6[:1]), g2)), b4.Stop))
    # 9. --verificar sobre el resultado sintetico (parcial: m = 6 y 48) y uno alterado. R4: la verificacion integral
    #    lo rechaza (no es una entrega completa); la inspeccion parcial se etiqueta como tal y devuelve el codigo 2
    jdump(res, out)
    nc9 = len(b4.CHECKS)
    vi = verify_result(res, sello=s_ok, plan=plan)
    ci_, v = verify_cli(res, parcial=True)
    R3 = json.load(open(res, encoding="ascii"))
    R3["instancias"]["TIRA_48x6"]["celdas"][CERT_KEY]["por_k"]["2"]["eps_comp"][0] = "0"
    alt = os.path.join(d, "RESULTADO_SINTETICO_ALTERADO_2.json")
    jdump(alt, R3)
    v2 = verify_result(alt, parcial=True)
    del b4.CHECKS[nc9:]
    ok("--verificar integral rechaza el resultado sintetico parcial (%d problemas de contrato, p. ej. %s); "
       "--parcial: INSPECCION_PARCIAL con codigo 2, 0 discrepancias en %d filas, y una vista alterada se detecta"
       % (len(vi["contrato"]), vi["contrato"][:1], v["filas"]),
       not vi["ok"] and any("instancias" in x for x in vi["contrato"]) and ci_ == 2 and not v["ok"]
       and v["tipo"] == "INSPECCION_PARCIAL" and not v["discrepancias"]
       and v2["discrepancias"] == ["TIRA_48x6 " + CERT_KEY], [vi["contrato"][:4], v2["discrepancias"]])
    # 10. corridas 1 y 2 completas del conductor con la lectura de catalogos, la DP, el control de N1 y el Monte
    #     Carlo sustituidos (sin decoder ni DP); N2 real (reproduce el D-t4 de referencia ya conocido de SOLO_U)
    for desc, cond, det in humo_corridas(d, refs, ce, rect, cells, rng, fake):
        ok(desc, cond, det)
    ok("ninguna llamada al decoder en la prueba de humo", b4.CALLS[0] == c0 and ce.CALLS[0] == ce0,
       [b4.CALLS[0] - c0, ce.CALLS[0] - ce0])
    return n[0]


class _Clock:
    """Solo pruebas: reloj real mas un desplazamiento que fija la prueba."""

    def __init__(self):
        self.shift = 0.0

    def __call__(self):
        return time.time() + self.shift


class _LateProc:
    """Solo pruebas: envuelve un Popen real; cuando su wait() observa la salida, el reloj de prueba salta jump s."""

    def __init__(self, p, clk, jump):
        self._p, self._clk, self._jump = p, clk, jump

    def wait(self, timeout=None):
        r = self._p.wait(timeout)
        self._clk.shift = self._jump
        return r

    def __getattr__(self, name):
        return getattr(self._p, name)


class _FaultK32:
    """Solo pruebas: kernel32 en el que la funcion name falla (devuelve 0 con GetLastError = err, sin llamar a
    Windows); con info_class, solo para esa clase de QueryInformationJobObject. El resto pasa sin cambios."""

    def __init__(self, real, name, err, info_class=None):
        self._real, self._name, self._err, self._cls = real, name, err, info_class
        self.fallos = 0

    def __getattr__(self, attr):
        f = getattr(self._real, attr)
        if attr != self._name:
            return f

        def g(*a):
            if self._cls is not None and a[1] != self._cls:
                return f(*a)
            self.fallos += 1
            ctypes.set_last_error(self._err)
            return 0
        return g


class _AbortoSintetico(BaseException):
    """Solo pruebas: interrupcion que no es Exception (como KeyboardInterrupt)."""


def _late_job(clk, jump):
    class _LateJob(JobTree):
        def start(self, cmd, stdout):
            return _LateProc(JobTree.start(self, cmd, stdout), clk, jump)
    return _LateJob


def _fault_job(name, err, info_class=None):
    class _FaultJob(JobTree):
        def __init__(self):
            JobTree.__init__(self)
            self.k32 = _FaultK32(self.k32, name, err, info_class)
    return _FaultJob


def _raising_job(exc):
    class _RaisingJob(JobTree):
        def memory(self):
            raise exc
    return _RaisingJob


def _no_job():
    raise JobQueryError(5, "CreateJobObjectW fallo (GetLastError 5) [sintetico, sin llamar a Windows]")


def _dead_within(pid, s=5.0):
    """True si el proceso pid ya no esta vivo (o deja de estarlo en s segundos)."""
    if not pid:
        return False
    t = time.time()
    while proc_alive(pid):
        if time.time() - t > s:
            return False
        time.sleep(0.1)
    return True


def prueba_cierre_r3(d, cert, mc_sha):
    """R3 (v4): cierre del arbol de t2_mc.py con procesos sinteticos (sin t2_mc.py, decoder, DP ni Monte Carlo) y
    sustitutos de Windows (fallos inyectados en kernel32, reloj de prueba). Devuelve (descripcion, condicion,
    detalle) por caso; cada negativo debe fallar por su causa prevista (su motivo), no por otra."""
    py, csha = sys.executable, sha_file(cert)
    fake = ("import json, sys; json.dump({'programa_sha256': sys.argv[2], 'modo': 'Monte Carlo', 'certificado': "
            "{'sha256': sys.argv[3]}, 'D_mc': {'decision': 'COMPATIBLE'}, 'recuento_archivo': "
            "{'archivo': 'recuento sintetico', 'x': 0, 'N': 1}}, open(sys.argv[1], 'w'))")
    spawn = ("import json, os, subprocess, sys, time; c = subprocess.Popen([sys.executable, '-I', '-c', "
             "'import time; time.sleep(120)']); print('hijo', c.pid, 'padre', os.getpid(), flush=True); ")
    orphan = spawn + fake.split("; ", 1)[1]           # sale bien con su salida, pero deja un hijo vivo
    sleeper = spawn + "%s time.sleep(120)"             # no sale: lo para un techo, el plazo o un fallo
    res = []

    def run(tag, script, **kw):
        mo, lg = os.path.join(d, "DMC_R3_%s.json" % tag), os.path.join(d, "SALIDA_DMC_R3_%s.txt" % tag)
        r = run_mc(cert, mo, lg, mc_sha, time.time() + 120, cmd=[py, "-I", "-c", script, mo, mc_sha, csha], **kw)
        return r, mo, lg

    def ids(lg):
        try:
            w = open(lg, encoding="utf-8").read().split()
            return int(w[1]), int(w[3])
        except Exception:
            return None, None

    def has(r, *words):
        return all(w in (r.get("motivo") or "") for w in words)

    # A. salida normal a tiempo y arbol vacio: el camino aceptado se conserva
    ra, mo, _ = run("normal", fake, work_s=60, close_s=10)
    ta = ra["arbol"]
    res.append(("salida normal a tiempo y arbol vacio: se acepta (COMPATIBLE), cierre confirmado con 0 procesos "
                "activos medidos, sin terminacion forzada ni errores de Windows y salida dentro del plazo de "
                "aceptacion", ra["decision"] == "COMPATIBLE" and ta["cierre_confirmado"] is True
                and ta["procesos_activos_al_final"] == 0 and ta["terminacion_forzada"] is False
                and ta["errores_api"] == [] and ta["problemas_de_cierre"] == [] and "pendiente_de_inspeccion" not in ra
                and ta["salida_observada_s_desde_inicio"] <= ta["plazo_aceptacion_s_desde_inicio"]
                and ta["exceso_sobre_plazo_duro_s"] == 0 and ra["salida"]["sha256"] == sha_file(mo), ta))
    # B. salida normal observada despues del plazo de aceptacion (reloj de prueba adelantado al observarla)
    clk = _Clock()
    rb, mo, _ = run("tardia", fake, work_s=60, close_s=10, clock=clk, job_factory=_late_job(clk, 120.0))
    res.append(("salida normal observada despues del plazo de aceptacion (reloj de prueba adelantado 120 s al "
                "observarla): NO_EVALUADO por la salida tardia; salida y recuento conservados, pendiente de "
                "inspeccion", rb["decision"] == "NO_EVALUADO" and has(rb, "despues del plazo de aceptacion")
                and rb.get("decision_de_t2_mc") == "COMPATIBLE" and bool(rb.get("recuento"))
                and rb.get("salida", {}).get("sha256") == sha_file(mo) and rb.get("pendiente_de_inspeccion") is True,
                [rb.get("motivo"), rb["arbol"]]))
    # C. lista de PID parcial (ERROR_MORE_DATA) durante la supervision: fallo, nunca una lista vacia ni memoria 0
    rc, _, lg = run("pids", sleeper % "", work_s=60, close_s=10,
                    job_factory=_fault_job("QueryInformationJobObject", 234, 3))
    kc, pc = ids(lg)
    tc = rc["arbol"]
    res.append(("consulta de PID fallida (ERROR_MORE_DATA, lista parcial) durante la supervision: NO_EVALUADO por "
                "fallo de supervision, registrado con la funcion y GetLastError 234; memoria total no medida (null, "
                "no 0); la terminacion de emergencia se ejecuta y cierra padre e hijo",
                rc["decision"] == "NO_EVALUADO" and has(rc, "supervision", "JobObjectBasicProcessIdList",
                                                        "lista parcial")
                and tc["pico_memoria_total_bytes"] is None
                and any(e["error"] == 234 and "ProcessIdList" in e["funcion"] for e in tc["errores_api"])
                and tc["terminacion_forzada"] is True and tc["procesos_activos_al_final"] == 0
                and _dead_within(kc) and _dead_within(pc), tc))
    # D. contadores fallidos (ERROR_INVALID_HANDLE) en el cierre: nunca 0 procesos activos medidos
    rd, mo, lg = run("contadores", orphan, work_s=60, close_s=10,
                     job_factory=_fault_job("QueryInformationJobObject", 6, 1))
    kd, _ = ids(lg)
    td = rd["arbol"]
    res.append(("contadores del trabajo fallidos (GetLastError 6) al cerrar: el arbol no se acredita vacio "
                "(procesos activos null, cierre no confirmado), NO_EVALUADO con la funcion en el motivo; aun asi se "
                "termina el arbol (el hijo muere) y la salida se conserva",
                rd["decision"] == "NO_EVALUADO" and has(rd, "JobObjectBasicAccountingInformation")
                and td["procesos_activos_al_final"] is None and td["cierre_confirmado"] is False
                and td["terminacion_forzada"] is True and _dead_within(kd)
                and rd.get("salida", {}).get("sha256") == sha_file(mo) and rd.get("pendiente_de_inspeccion") is True,
                td))
    # E. TerminateJobObject fallido con un descendiente vivo tras una salida normal (ventana de salida de 22 s: un
    #    retraso de salida de unos segundos no cambia de camino; intento 1 de la regresion R3)
    re_, _, lg = run("terminacion", orphan, work_s=20, close_s=4, job_factory=_fault_job("TerminateJobObject", 5))
    ke, _ = ids(lg)
    te = re_["arbol"]
    res.append(("TerminateJobObject fallido (GetLastError 5) con un descendiente vivo: NO_EVALUADO con la funcion "
                "en el motivo, cierre no confirmado; job.close() se ejecuta igual (KILL_ON_JOB_CLOSE) y el hijo muere",
                re_["decision"] == "NO_EVALUADO" and has(re_, "TerminateJobObject")
                and te["cierre_confirmado"] is False and te.get("descendientes_vivos_tras_la_salida") == 1
                and any(e["funcion"] == "TerminateJobObject" and e["error"] == 5 for e in te["errores_api"])
                and _dead_within(ke), te))
    # F. salida buena con un descendiente vivo: se termina por fuerza -> NO_EVALUADO, evidencia conservada
    rf, mo, lg = run("hijos", orphan, work_s=60, close_s=10)
    kf, _ = ids(lg)
    tf = rf["arbol"]
    res.append(("t2_mc.py sale bien pero deja un descendiente vivo: se termina por fuerza y D-mc queda "
                "NO_EVALUADO, pendiente de inspeccion, con la salida, la decision propia (COMPATIBLE) y el recuento "
                "conservados; el arbol queda vacio y confirmado",
                rf["decision"] == "NO_EVALUADO" and has(rf, "terminar el arbol por fuerza", "pendiente de inspeccion")
                and rf.get("pendiente_de_inspeccion") is True and rf.get("decision_de_t2_mc") == "COMPATIBLE"
                and bool(rf.get("recuento")) and rf.get("salida", {}).get("sha256") == sha_file(mo)
                and tf.get("descendientes_vivos_tras_la_salida") == 1 and tf["terminacion_forzada"] is True
                and tf["cierre_confirmado"] is True and tf["procesos_activos_al_final"] == 0 and kf in tf["pids"]
                and _dead_within(kf), tf))
    # G. excepcion ordinaria durante la supervision
    rg, _, lg = run("excepcion", sleeper % "", work_s=60, close_s=10,
                    job_factory=_raising_job(RuntimeError("fallo sintetico de supervision")))
    kg, pg = ids(lg)
    tg = rg["arbol"]
    res.append(("excepcion durante la supervision (RuntimeError sintetico): NO_EVALUADO con el motivo; la limpieza "
                "no se salta (arbol terminado, 0 procesos activos medidos, padre e hijo muertos)",
                rg["decision"] == "NO_EVALUADO" and has(rg, "supervision", "fallo sintetico de supervision")
                and tg["terminacion_forzada"] is True and tg["procesos_activos_al_final"] == 0
                and tg["cierre_confirmado"] is True and _dead_within(kg) and _dead_within(pg), tg))
    # H. interrupcion que no es Exception durante la supervision: se propaga, pero la limpieza se ejecuta
    lg_h = os.path.join(d, "SALIDA_DMC_R3_interrupcion.txt")
    propagated = False
    try:
        run("interrupcion", sleeper % "", work_s=60, close_s=10,
            job_factory=_raising_job(_AbortoSintetico("interrupcion sintetica")))
    except _AbortoSintetico:
        propagated = True
    kh, ph = ids(lg_h)
    res.append(("interrupcion durante la supervision (BaseException sintetica, como KeyboardInterrupt): se "
                "propaga y la limpieza se ejecuta igual (padre e hijo muertos)",
                propagated and _dead_within(kh) and _dead_within(ph), [propagated, kh, ph]))
    # I. plazo de aceptacion con padre e hijo; la espera no pasa del plazo; un proceso ajeno no se toca
    other = subprocess.Popen([py, "-I", "-c", "import time; time.sleep(120)"])
    try:
        ri, _, lg = run("plazo", sleeper % "", work_s=3, close_s=2)
        other_alive = other.poll() is None
    finally:
        other.kill()
        other.wait()
    ki, pi = ids(lg)
    ti = ri["arbol"]
    res.append(("plazo de aceptacion (trabajo 3 s + cierre 2 s - reserva 1 s = 4 s) con padre e hijo: la ultima "
                "espera se recorta al tiempo restante (parada a %.3f s), el arbol se termina y el cierre acaba "
                "antes del plazo duro; D-mc NO_EVALUADO; un proceso ajeno sigue vivo" % ti.get("parada_s_desde_inicio",
                                                                                              -1),
                ri["decision"] == "NO_EVALUADO" and has(ri, "plazo de aceptacion")
                and 4.0 <= ti["parada_s_desde_inicio"] < 4.0 + 1.0 and ti["exceso_sobre_plazo_duro_s"] == 0
                and ti["plazo_duro_s_desde_inicio"] <= 5.001 and ti["fin_menos_techo_de_la_corrida_s"] < 0
                and ki in ti["pids"] and pi in ti["pids"] and ti["procesos_activos_al_final"] == 0
                and ti["cierre_confirmado"] is True and _dead_within(ki) and _dead_within(pi) and other_alive
                and other.pid not in ti["pids"] and all(x in ti for x in (
                    "plazo_trabajo_s", "margen_cierre_s", "reserva_limpieza_s", "plazo_aceptacion_s_desde_inicio",
                    "exceso_sobre_plazo_de_trabajo_s", "entorno")), ti))
    # J. techo de memoria total con el arbol medido
    rj_mem = (b4.mem_now() or 0) + 150 * (1 << 20)
    rj, _, lg = run("memoria", sleeper % "b = b'x' * (400 * 2 ** 20);", work_s=60, close_s=5, mem_total=rj_mem)
    kj, pj = ids(lg)
    tj = rj["arbol"]
    res.append(("techo de memoria total con el arbol medido (padre de 400 MB, techo de principal + 150 MB): "
                "NO_EVALUADO, arbol terminado entero y pico medido por encima del techo",
                rj["decision"] == "NO_EVALUADO" and has(rj, "memoria total") and (tj["pico_memoria_total_bytes"] or 0)
                > rj_mem and tj["procesos_activos_al_final"] == 0 and _dead_within(kj) and _dead_within(pj), tj))
    # K. el arbol no se puede crear: t2_mc.py no se ejecuta
    rk, mo, _ = run("inicio", fake, work_s=60, close_s=10, job_factory=_no_job)
    res.append(("el arbol no se puede crear (fallo sintetico de CreateJobObjectW): NO_EVALUADO sin ejecutar "
                "t2_mc.py (sin salida ni codigo), con el entorno registrado",
                rk["decision"] == "NO_EVALUADO" and has(rk, "no se pudo iniciar el arbol", "CreateJobObjectW")
                and rk["codigo"] is None and "salida" not in rk and not os.path.exists(mo)
                and "entorno" in rk["arbol"], rk))
    return res


def prueba_cierre_cli(d):
    """--prueba-cierre: solo los casos de cierre de R3 (prueba_cierre_r3), sin catalogos, decoder, DP ni Monte
    Carlo. Codigo 0 si todos pasan; 1 si alguno falla (se ejecutan todos)."""
    os.makedirs(d)
    cert = os.path.join(d, "CERTIFICADO_SINTETICO_R3.txt")
    with open(cert, "w", encoding="ascii", newline="\n") as fh:
        fh.write("certificado sintetico de la prueba de cierre (run_mc solo calcula su SHA-256)\n")
    c0 = b4.CALLS[0]
    res = prueba_cierre_r3(d, cert, sha_file(MC_PATH))
    out = []
    for desc, cond, det in res:
        print("  [%s] R3 cierre: %s" % ("PASS" if cond else "FALLA", desc), flush=True)
        out.append({"control": "R3 cierre: " + desc, "pasa": bool(cond), "detalle": det})
    calls = b4.CALLS[0] - c0
    out.append({"control": "ninguna llamada al decoder", "pasa": calls == 0, "detalle": calls})
    code = 0 if all(x["pasa"] for x in out) else 1
    jdump(os.path.join(d, "CONTROLES_CIERRE_R3.json"), {"controles": out, "codigo": code})
    print("PRUEBA DE CIERRE R3: %d casos, %d PASS; %d llamadas al decoder; sin DP ni Monte Carlo; codigo %d"
          % (len(res), sum(1 for x in out[:-1] if x["pasa"]), calls, code), flush=True)
    return code


def proc_alive(pid):
    """Si el proceso pid sigue vivo (solo consulta; para las pruebas del arbol de R3)."""
    from ctypes import wintypes as W
    k = ctypes.WinDLL("kernel32", use_last_error=True)
    k.OpenProcess.restype = W.HANDLE
    k.OpenProcess.argtypes = [W.DWORD, W.BOOL, W.DWORD]
    k.GetExitCodeProcess.argtypes = [W.HANDLE, ctypes.POINTER(W.DWORD)]
    k.CloseHandle.argtypes = [W.HANDLE]
    h = k.OpenProcess(0x1000, False, pid)              # PROCESS_QUERY_LIMITED_INFORMATION
    if not h:
        return False
    try:
        code = W.DWORD()
        return bool(k.GetExitCodeProcess(h, ctypes.byref(code))) and code.value == 259     # STILL_ACTIVE
    finally:
        k.CloseHandle(h)


def w_sleep_scan(task):
    """Trabajador de la prueba de humo: no devuelve nada en 60 s (espera de N1 sin resultados)."""
    time.sleep(60)
    return {"error": "no deberia llegar", "g": task[0]}


def _stop_text(fn):
    """Texto de la parada (b4.Stop) que produce fn, o '' si no se detiene."""
    try:
        fn()
    except b4.Stop as ex:
        return str(ex)
    return ""


def humo_r2(d):
    """R2 con procesos reales y un catalogo sintetico en el formato del motor (sin leer los catalogos de T2-1M):
    la lectura de N1 se detiene por la memoria total que marca t2_motor.Monitor, por la memoria ya marcada con
    trabajadores reales, por el plazo durante una espera sin resultados; y el techo de disco detiene guard."""
    res = []
    dc = os.path.join(d, "catalogo_sintetico_R2")
    tag = "circ_SINTETICA"
    os.makedirs(os.path.join(dc, tag))
    man = {}
    for g in range(3):
        secs = [riv.section_dict(a, [(1 << a, 1 << a, [(a + 1, 1, 1, 0)])]) for a in range(g * 4, g * 4 + 4)]
        f = "g%03d.bin" % g
        man[str(g)] = {"archivo": f, "sha256": riv.write_group(os.path.join(dc, tag, f), secs)}
    MAN = {tag: man}
    sc = {"tag": tag, "NS": 12, "ngroups": 3, "gsize": 4}
    ctx = mp.get_context("spawn")
    saved = (riv.GUARD[0], mot.MEM_TOTAL, DEADLINE[0], mot.OVER[0], RUN_DIR[0], DISK_MAX[0])
    try:
        riv.GUARD[0] = guard
        DEADLINE[0] = time.time() + 3600
        RUN_DIR[0] = None
        # a. t2_motor.Monitor con un techo total sintetico de 1 byte (solo aqui; se restaura) y trabajadores que no
        #    devuelven nada: el Monitor marca OVER en su primera medida (5 s) y la espera de N1 se detiene
        mot.MEM_TOTAL, mot.OVER[0] = 1, None
        t = time.time()
        with ctx.Pool(2) as pool:
            mon = mot.Monitor(pool, "humo_R2")
            try:
                fa = _stop_text(lambda: riv.scan_catalogue(pool, sc, MAN, None, time.time() + 600, wdir=dc,
                                                           worker=w_sleep_scan))
            finally:
                mon.stop()
            t_stop = time.time() - t
        t_end = time.time() - t
        pa = copy.deepcopy(riv.PROGRESS.get(tag))
        res.append(("R2: con t2_motor.Monitor (techo total sintetico de 1 byte) y trabajadores sin resultados, la "
                    "lectura de N1 se detiene en espera por memoria total (%.1f s) y el pool se termina (%.1f s); "
                    "avance guardado (%d grupos leidos, %d esperas)" % (t_stop, t_end, pa["grupos_leidos"],
                                                                       pa["esperas"]),
                    "techo de memoria total" in fa and "en espera" in fa and t_stop < 2 * POLL + 5 and t_end < 30
                    and pa["grupos_leidos"] == 0 and pa["esperas"] >= 1 and pa["estado"] == "leyendo", [fa, pa]))
        mot.MEM_TOTAL, mot.OVER[0] = saved[1], None
        # b. memoria total ya marcada y trabajadores reales (read_group y _scan sobre el catalogo sintetico)
        mot.OVER[0] = "sintetico"
        with ctx.Pool(2) as pool:
            fb = _stop_text(lambda: riv.scan_catalogue(pool, sc, MAN, None, time.time() + 600, wdir=dc))
        pb = copy.deepcopy(riv.PROGRESS.get(tag))
        mot.OVER[0] = None
        res.append(("R2: con la memoria total ya marcada, la lectura real (read_group y _scan) se detiene en el "
                    "primer resultado o la primera espera (%d de 3 grupos leidos)" % pb["grupos_leidos"],
                    "techo de memoria total" in fb and pb["grupos_leidos"] <= 1, [fb, pb]))
        # c. plazo de la unidad que vence durante una espera sin resultados
        t = time.time()
        with ctx.Pool(2) as pool:
            fc = _stop_text(lambda: riv.scan_catalogue(pool, sc, MAN, None, time.time() + 2, wdir=dc,
                                                       worker=w_sleep_scan))
            t_c = time.time() - t
        res.append(("R2: plazo vencido durante una espera sin resultados: la lectura se detiene (%.1f s)" % t_c,
                    "techo de tiempo" in fc and t_c < 2 + POLL + 5, fc))
        # d. disco: la carpeta de la corrida por encima del techo (sintetico) detiene guard; por debajo, no
        dd = os.path.join(d, "disco_sintetico_R2")
        os.makedirs(dd)
        with open(os.path.join(dd, "x.bin"), "wb") as fh:
            fh.write(b"x" * 1000)
        RUN_DIR[0], DISK_MAX[0] = dd, 2000
        _DISK["t"] = 0.0
        fd0 = _stop_text(lambda: guard("humo, disco por debajo del techo"))
        DISK_MAX[0] = 500
        _DISK["t"] = 0.0
        fd = _stop_text(lambda: guard("humo, disco"))
        res.append(("R2: techo de disco: 1000 bytes en la carpeta de la corrida pasan con un techo de 2000 y detienen "
                    "con uno de 500", fd0 == "" and "techo de disco" in fd, [fd0, fd]))
    finally:
        riv.GUARD[0], mot.MEM_TOTAL, DEADLINE[0], mot.OVER[0], RUN_DIR[0], DISK_MAX[0] = saved
        _DISK.update(t=0.0, bytes=0, max=0)
    return res


_SC = {}


def sc_cache(ce, rect, refs, m):
    if m not in _SC:
        _SC[m] = riv.strip_context(ce, rect, m, False, refs["B4"])
    return _SC[m]


def rv_of(inst, X, m, sc, cells):
    return rival_from_evidence(sc, inst["TIRA_%dx6" % m]["rival"]["evidencia"], cells)


def humo_corridas(d, refs, ce, rect, cells, rng, fake):
    """Corridas 1 y 2 de run_main con sustitutos sinteticos de la lectura de catalogos (N1), la DP, el control
    fenomenologico, el sello, la capacidad y trabajo/. Devuelve (descripcion, condicion, detalle) para cada control."""
    G = globals()
    names = ("check_seal", "capacity_probe", "verify_trabajo", "run_dp")
    saved = {k: G[k] for k in names}
    rsaved = (riv.fen_control, riv.circuit_rival)
    sello = os.path.join(d, "SELLO_SINTETICO_CORRIDAS.json")
    plan = os.path.join(d, "PLAN_SINTETICO_CORRIDAS.txt")
    with open(plan, "w", encoding="ascii", newline="\n") as fh:
        fh.write("plan sintetico de las corridas de la prueba de humo\n")
    jdump(sello, {"nota": "sello sintetico de la prueba de humo", "plan": os.path.basename(plan),
                  "sha256": sha_file(plan), "codigo": {c: sha_file(os.path.join(UDIR, c)) for c in CODE}})
    seal = check_seal(sello, plan)
    synth, nspec, wired, stop_n1 = {}, [], [], []

    def s_rival(pool, ce_, rect_, m, refs_, cells_, deadline):
        wired.append(riv.GUARD[0] is guard)
        if stop_n1:
            riv.PROGRESS["circ_TIRA_%dx6" % m] = {"estado": "leyendo", "grupos": 256, "grupos_leidos": 17,
                                                  "nota": "avance sintetico"}
            raise b4.Stop("techo de memoria total: sintetico (N1, catalogo de circ_TIRA_%dx6, en espera)" % m)
        sc = sc_cache(ce, rect, refs, m)
        scan = synth_scan(sc, rng)
        vals = rival_from_evidence(sc, riv.scan_summary(scan), cells)
        synth[m] = (sc, vals)
        return {"contexto": sc, "leyes": riv.law_cost(sc["cx"], cells), "lectura": scan, "valores": vals,
                "truncacion_ingenua": riv.esym_counts(sc["cx"]["nclass"], KMAX)}

    def s_dp(specs, store):
        labs = [s["etiqueta"] for s in specs]
        nspec.append(len(labs))
        for m, (sc, vals) in synth.items():
            for lab, j in synth_jobs(m, vals, refs, cells, sum_mode="COMPLETO" if m in SUM_STRIPS else None,
                                     sc=sc).items():
                if lab in labs:
                    store[lab] = j
        check("humo N3: trabajos sinteticos = especificaciones de make_specs (%d)" % len(labs),
              sorted(store) == sorted(labs) and len(set(labs)) == len(labs), verbose=False)
        return store

    def mc_cmd(cert_p, mc_out):
        return [sys.executable, "-I", "-c", fake, mc_out, seal["codigo"]["t2_mc.py"], sha_file(cert_p)]

    real_seal = saved["check_seal"]

    def s_seal(sello_=sello, plan_=plan, code_dir=UDIR):
        # check_seal real sobre el sello y el plan sinteticos (intento 1 de v3: un sustituto sin argumentos
        # hacia fallar el contrato R4 por el sello en todas las verificaciones)
        return real_seal(sello_, plan_, code_dir)

    G.update(check_seal=s_seal, capacity_probe=lambda: {"nota": "sustituida en la prueba de humo"},
             verify_trabajo=lambda pool, MAN: {"nota": "sustituida en la prueba de humo"}, run_dp=s_dp)
    riv.fen_control = lambda *a: {"nota": "sustituido en la prueba de humo"}
    riv.circuit_rival = s_rival
    res = []
    try:
        d1, d2, d1b, d2b = (os.path.join(d, x) for x in ("corrida1", "corrida2", "corrida1_no_lanzar", "corrida2_b"))
        os.makedirs(d1)
        c1 = run_main(1, d1, None, mc_cmd=mc_cmd)
        L1 = json.load(open(os.path.join(d1, "LECTURA_T2_1_C1.json"), encoding="ascii"))
        R1 = json.load(open(os.path.join(d1, "RESULTADO_T2_1_C1.json"), encoding="ascii"))
        cert = mcm.load_cert(os.path.join(d1, "CERTIFICADO_DMC.json"))
        res.append(("corrida 1 sintetica: codigo 0, %d etapas PASS, %d trabajos de make_specs (96 celdas + 1 regla "
                    "de suma), certificado aceptado por t2_mc.load_cert, D-mc COMPATIBLE (sustituido), regla de "
                    "paso LANZAR_CORRIDA_2" % (len(R1["etapas"]), nspec[0]),
                    c1 == 0 and R1["estado"] == "PASS" and nspec[0] == 97 and L1["regla_de_paso"] == "LANZAR_CORRIDA_2"
                    and L1["D_mc"]["decision"] == "COMPATIBLE" and cert["lo"] > 0
                    and [s["etapa"] for s in R1["etapas"]] == list(PLANNED[1]), [c1, R1["estado"], nspec]))
        rs24 = R1["instancias"]["TIRA_24x6"]["regla_de_suma"]["4"]
        pb1 = R1["regla_de_suma_antes_de_la_DP"]["por_tira"]
        it24 = R1["instancias"]["TIRA_24x6"]["igual_tiempo"]
        res.append(("corrida 1: condicion previa comprobada antes de la DP (m = 24, cola/umbral %.3g con D = 1e16); "
                    "regla de suma sintetica COMPLETO en m = 24 -> H-comp COMPLETO en k = 1..4, D-t3 calculado y "
                    "linea de igual tiempo con L_ing(4) exacto (k' = %s); D-t4 %s (N2 real)"
                    % (pb1["24"]["cola_sobre_umbral"], it24["k_prima"], R1["N2_S2"]["D_t4"]["decision"]),
                    list(pb1) == ["24"] and pb1["24"]["condicion_previa"] is True and rs24["veredicto"] == "COMPLETO"
                    and rs24["D"] == str(D_SUMA) and R1["D_regla_de_suma_T2_1"] == str(D_SUMA)
                    and [R1["instancias"]["TIRA_24x6"]["H_comp"][str(k)]["veredicto"] for k in range(1, 5)]
                    == ["COMPLETO"] * 4
                    and "celdas_con_ventaja" in R1["decisiones"]["D_t3"]["por_tira"]["24"]
                    and all("PENDIENTE" not in json.dumps(c["rival_en_k_prima"]) for c in it24["celdas"].values()),
                    [rs24, it24["k_prima"]]))
        synth.clear()
        os.makedirs(d2)
        c2 = run_main(2, d2, d1)
        R2 = json.load(open(os.path.join(d2, "RESULTADO_T2_1.json"), encoding="ascii"))
        D = R2["decisiones"]
        res.append(("corrida 2 sintetica: codigo 0, %d trabajos (32 celdas + 1 regla de suma), filas de la corrida "
                    "1 reproducidas desde su evidencia, decisiones finales con las 4 tiras (D-t1 %s, D-t2 %s, D-t5 "
                    "%s)" % (nspec[1], D["D_t1"]["crecimiento"], D["D_t2"]["decision"], D["D_t5"]["decision"]),
                    c2 == 0 and R2["estado"] == "PASS" and nspec[1] == 33 and D["final"] is True
                    and sorted(R2["instancias"]) == ["TIRA_12x6", "TIRA_24x6", "TIRA_48x6", "TIRA_6x6"]
                    and D["D_t1"]["crecimiento"] in ("CRECE_CON_N", "NO_CRECE", "NO_DECIDIDA_POR_REDONDEO")
                    and D["D_t2"]["decision"] in ("ESTABLE", "CRECE") and D["D_t5"]["decision"] in
                    ("CUMPLE", "NO_CUMPLE") and "celdas_con_ventaja" in D["D_t3"]["por_tira"]["48"]
                    and [R2["instancias"]["TIRA_48x6"]["H_comp"][str(k)]["veredicto"] for k in range(1, 5)]
                    == ["COMPLETO"] * 4 and list(R2["regla_de_suma_antes_de_la_DP"]["por_tira"]) == ["48"],
                    [c2, R2["estado"], nspec]))
        p1, p2 = os.path.join(d1, "RESULTADO_T2_1_C1.json"), os.path.join(d2, "RESULTADO_T2_1.json")
        (k1, v1), (k2, v2) = verify_cli(p1, sello=sello, plan=plan), verify_cli(p2, sello=sello, plan=plan)
        res.append(("--verificar integral de las dos corridas sinteticas: contrato completo, %d y %d filas, 0 "
                    "discrepancias (filas, N2 recalculado desde SOLO_U y decisiones con ese D-t4), codigo 0"
                    % (v1["filas"], v2["filas"]), k1 == 0 and k2 == 0 and v1["ok"] and v2["ok"] and v1["filas"] == 96
                    and v2["filas"] == 128, [v1["contrato"], v1["discrepancias"], v2["contrato"], v2["discrepancias"]]))
        res.append(("R2: durante la corrida, la lectura de N1 tiene los techos de la corrida (riv.GUARD = guard) y "
                    "se quitan al terminar", wired == [True, True, True, True] and riv.GUARD[0] is None, wired))
        # R4: negativos por omision y por metadatos (copias alteradas; ninguna corrida ni DP)
        dn = os.path.join(d, "verificar_negativos")
        os.makedirs(dn)
        B1, B2 = (json.load(open(p, encoding="ascii")) for p in (p1, p2))
        neg = {}

        def variant(name, base, fn, want):
            # cada negativo debe fallar por SU defecto (want), no por otro comun a todos (intento 1 de v3)
            R = copy.deepcopy(base)
            fn(R)
            pth = os.path.join(dn, name + ".json")
            jdump(pth, R)
            nc = len(b4.CHECKS)
            code_, r_ = verify_cli(pth, sello=sello, plan=plan)
            del b4.CHECKS[nc:]
            probs = r_["contrato"] + r_["discrepancias"]
            neg[name] = {"codigo": code_, "esperado": want, "detectado": any(want in p for p in probs),
                         "problemas": probs}

        def alter_dt4(R):
            for x in (R["decisiones"]["D_t4"], R.get("N2_S2", {}).get("D_t4", {})):
                x["decision"] = "FALLA" if x.get("decision") == "PASA" else "PASA"

        cases1 = {
            "vacio": (lambda R: (R.clear(), R.update(instancias={})), "unidad distinta de T2-1"),
            "sin_instancias": (lambda R: R.update(instancias={}), "instancias []"),
            "sin_tira_24": (lambda R: R["instancias"].pop("TIRA_24x6"), "instancias ['TIRA_12x6', 'TIRA_6x6']"),
            "sin_celda": (lambda R: R["instancias"]["TIRA_12x6"]["celdas"].pop(CERT_KEY),
                          "TIRA_12x6: celdas distintas de las 8 previstas"),
            "sin_corte": (lambda R: R["instancias"]["TIRA_6x6"]["celdas"][CERT_KEY]["por_k"].pop("2"),
                          "cortes distintos de 1..4"),
            "sin_trabajo": (lambda R: R["instancias"]["TIRA_6x6"]["trabajos"].pop(label(6, 1, CERT_KEY)),
                            "TIRA_6x6: trabajos distintos de los previstos"),
            "sin_regla_de_suma": (lambda R: R["instancias"]["TIRA_24x6"]["regla_de_suma"].pop("4"),
                                  "TIRA_24x6: regla de suma distinta de la prevista"),
            "D_distinto": (lambda R: R.update(D_regla_de_suma_T2_1=str(10 ** 14)),
                           "D de la regla de suma distinto de 1e16"),
            "programa_distinto": (lambda R: R.update(programa_sha256="0" * 64),
                                  "programa_sha256 distinto del programa que verifica"),
            "corrida_distinta": (lambda R: R.update(corrida=2), "'TIRA_48x6']"),
            "estado_no_PASS": (lambda R: R.update(estado="NO_EVALUADO"), "estado distinto de PASS"),
            "etapa_ausente": (lambda R: R["etapas"].pop(), "etapas distintas de las previstas"),
            "D_t4_alterado_con_N2": (alter_dt4, "decisiones distintas de las recalculadas"),
            "decision_ausente": (lambda R: R["decisiones"].pop("D_t3"), "decisiones incompletas"),
            "sin_N2": (lambda R: R.pop("N2_S2"), "falta N2_S2"),
        }
        for name, (fn, want) in cases1.items():
            variant("c1_" + name, B1, fn, want)
        cases2 = {
            "sin_tira_48": (lambda R: R["instancias"].pop("TIRA_48x6"),
                            "instancias ['TIRA_12x6', 'TIRA_24x6', 'TIRA_6x6']"),
            "sin_celda_previa": (lambda R: R["instancias"]["TIRA_6x6"]["celdas"].pop(CERT_KEY),
                                 "TIRA_6x6: celdas distintas de las 8 previstas"),
            "corrida_1_sha": (lambda R: R["corrida_1"]["resultado"].update(sha256="0" * 64),
                              "SHA-256 del resultado de la corrida 1 distinto"),
            "D_t4_alterado": (alter_dt4, "decisiones distintas de las recalculadas"),
            "no_final": (lambda R: R["decisiones"].update(final=False), "decisiones: final o tiras distintas"),
        }
        for name, (fn, want) in cases2.items():
            variant("c2_" + name, B2, fn, want)
        kp, vp = verify_cli(p1, parcial=True)
        jdump(os.path.join(dn, "RESUMEN_NEGATIVOS.json"), neg)
        res.append(("R4: --verificar integral rechaza (codigo 1) las %d entregas incompletas o de otro contrato, "
                    "cada una por su defecto (las dos de partida cumplen): %s; la inspeccion parcial de una entrega "
                    "completa devuelve el codigo 2 (nunca 0)" % (len(neg), ", ".join(sorted(neg))),
                    k1 == 0 and k2 == 0 and all(v["codigo"] == 1 and v["detectado"] for v in neg.values())
                    and kp == 2 and vp["tipo"] == "INSPECCION_PARCIAL" and not vp["discrepancias"],
                    {k: v for k, v in neg.items() if not (v["codigo"] == 1 and v["detectado"])} or "todas"))
        os.makedirs(d1b)
        for f in ("RESULTADO_T2_1_C1.json", "LECTURA_T2_1_C1.json"):
            with open(os.path.join(d1, f), "rb") as fi, open(os.path.join(d1b, f), "wb") as fo:
                fo.write(fi.read())
        L1["regla_de_paso"] = "NO_LANZAR"
        jdump(os.path.join(d1b, "LECTURA_T2_1_C1.json"), L1)
        os.makedirs(d2b)
        nchk = len(b4.CHECKS)
        c3 = run_main(2, d2b, d1b)
        R3 = json.load(open(os.path.join(d2b, "RESULTADO_T2_1.json"), encoding="ascii"))
        st = {s["etapa"]: s["estado"] for s in R3["etapas"]}
        del b4.CHECKS[nchk:]
        res.append(("corrida 2 con una lectura NO_LANZAR: codigo 1, etapa previa FALLA, resto NO_EVALUADO, sin "
                    "decisiones", c3 == 1 and st["previa"] == "FALLA" and R3["estado"] == "FALLA"
                    and all(v == "NO_EVALUADO" for k, v in st.items() if k not in ("sello", "previa"))
                    and R3["decisiones"]["decision"] == "NO_EVALUADO", st))
        # R2: corrida 1 que se detiene durante la lectura de N1 (techo sintetico): etapa NO_EVALUADO, avance de la
        #     lectura guardado, sin decisiones ni Monte Carlo
        d1s = os.path.join(d, "corrida1_parada_N1")
        os.makedirs(d1s)
        stop_n1.append(True)
        nchk = len(b4.CHECKS)
        c4 = run_main(1, d1s, None, mc_cmd=mc_cmd)
        del b4.CHECKS[nchk:]
        stop_n1.clear()
        R4 = json.load(open(os.path.join(d1s, "RESULTADO_T2_1_C1.json"), encoding="ascii"))
        L4 = json.load(open(os.path.join(d1s, "LECTURA_T2_1_C1.json"), encoding="ascii"))
        st4 = {s["etapa"]: s["estado"] for s in R4["etapas"]}
        res.append(("R2: corrida 1 detenida por un techo durante la lectura de N1: codigo 1, rival_circuito "
                    "NO_EVALUADO con el motivo, avance de la lectura guardado, sin decisiones, D-mc NO_EVALUADO y "
                    "NO_LANZAR", c4 == 1 and st4["rival_circuito"] == "NO_EVALUADO"
                    and "memoria total" in next(s["motivo"] for s in R4["etapas"] if s["etapa"] == "rival_circuito")
                    and R4["N1_lectura_hasta_la_parada"]["circ_TIRA_6x6"]["grupos_leidos"] == 17
                    and R4["decisiones"]["decision"] == "NO_EVALUADO" and L4["D_mc"]["decision"] == "NO_EVALUADO"
                    and L4["regla_de_paso"] == "NO_LANZAR" and riv.GUARD[0] is None, st4))
    finally:
        G.update(saved)
        riv.fen_control, riv.circuit_rival = rsaved
    return res


# ------------------------------------------------------------------ programa
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--corrida", type=int, choices=(1, 2))
    ap.add_argument("--dir")
    ap.add_argument("--previa")
    ap.add_argument("--verificar")
    ap.add_argument("--parcial", action="store_true", help="con --verificar: inspeccion de parciales (codigo 2)")
    ap.add_argument("--humo", action="store_true")
    ap.add_argument("--prueba-cierre", action="store_true", help="solo los casos de cierre de R3 (sinteticos)")
    a = ap.parse_args()
    if a.verificar:
        raise SystemExit(verify_cli(a.verificar, a.parcial)[0])
    if not a.dir or os.path.exists(a.dir):
        raise SystemExit("--dir es obligatorio y no debe existir")
    if a.prueba_cierre:
        raise SystemExit(prueba_cierre_cli(a.dir))
    if a.humo:
        code = 1
        try:
            nchk = humo(a.dir)
            code = 0
            print("PRUEBA DE HUMO N3: %d grupos de comprobaciones, todos PASS; 0 llamadas al decoder; sin DP" % nchk,
                  flush=True)
        finally:
            jdump(os.path.join(a.dir, "CONTROLES_HUMO.json"), {"controles": b4.CHECKS, "codigo": code})
        raise SystemExit(code)
    if a.corrida is None or (a.corrida == 2 and not a.previa):
        raise SystemExit("uso: --corrida 1 --dir D | --corrida 2 --dir D --previa D1 | --verificar R [--parcial] | "
                         "--humo --dir D | --prueba-cierre --dir D")
    os.makedirs(a.dir)
    raise SystemExit(run_main(a.corrida, a.dir, a.previa))


def run_main(n, d, previa, mc_cmd=None):
    """Una corrida completa: etapas, resultado, recursos y, en la corrida 1, certificado, Monte Carlo y lectura con
    la regla de paso. mc_cmd solo lo usa la prueba de humo (sustituye t2_mc.py por un proceso sintetico)."""
    out = base_out(n)
    try:
        corrida(out, n, d, previa)
    except Exception as ex:  # la etapa ya quedo registrada (FALLA o NO_EVALUADO); nunca PASS parcial
        out["detenido_en"] = "%s: %s" % (type(ex).__name__, ex)
    close_stages(out, PLANNED[n])
    if "detenido_en" in out and out["estado"] == "PASS":   # parada fuera de una etapa: tampoco cuenta como PASS
        out["estado"] = "NO_EVALUADO"
        out["decisiones"] = {"decision": "NO_EVALUADO", "motivo": "parada fuera de una etapa: " + out["detenido_en"]}
    out["controles"] = b4.CHECKS
    out["segundos"] = round(time.time() - T0_WALL, 3)
    res = os.path.join(d, "RESULTADO_T2_1_C1.json" if n == 1 else "RESULTADO_T2_1.json")
    rsha = jdump(res, out)
    jdump(os.path.join(d, "RECURSOS_T2_1_C%d.json" % n), {
        "nota": "tiempos de reloj y memoria; no deterministas; no son resultado",
        "etapas": [{x: s[x] for x in ("etapa", "segundos") if x in s} for s in out["etapas"]],
        "capacidad": out.get("capacidad"), "pico_memoria_principal_bytes": b4.mem_now(peak=True),
        "pico_memoria_total_bytes_por_etapa": mot.MEMPEAK, "techo_memoria_total_bytes": mot.MEM_TOTAL,
        "disco_max_carpeta_corrida_bytes": max(_DISK["max"], disk_now(force=True)), "techo_disco_bytes": DISK_MAX[0],
        "segundos": out["segundos"]})
    print("CORRIDA %d: %s" % (n, out["estado"]), flush=True)
    code = 0 if out["estado"] == "PASS" else 1
    if n == 1:
        lec = {"unidad": "T2-1", "corrida": 1, "programa_sha256": out["programa_sha256"],
               "sello_sha256": out.get("sello", {}).get("sello_sha256"),
               "resultado": {"archivo": res.replace("\\", "/"), "sha256": rsha}, "estado_corrida": out["estado"],
               "etapas": []}
        dmc = {"decision": "NO_EVALUADO", "motivo": "la corrida 1 no termino con todas sus etapas PASS"}
        if out["estado"] == "PASS":
            cert_p, mc_out = os.path.join(d, "CERTIFICADO_DMC.json"), os.path.join(d, "RESULTADO_DMC.json")
            try:
                with Stage(lec, "certificado"):
                    sc6 = riv.strip_context(*b4.load_all(with_b1=False)[:2], mcm.M_STRIP, False, riv.load_refs()["B4"])
                    lec["certificado"] = write_certificate(res, cert_p, sc6)
                with Stage(lec, "monte_carlo"):
                    dmc = run_mc(cert_p, mc_out, os.path.join(d, "SALIDA_DMC.txt"), out["sello"]["codigo"]["t2_mc.py"],
                                 T0_WALL + SEG_C1, **({"cmd": mc_cmd(cert_p, mc_out), "work_s": 120, "close_s": 10}
                                                      if mc_cmd else {}))
            except Exception as ex:
                dmc = {"decision": "NO_EVALUADO", "motivo": "%s: %s" % (type(ex).__name__, ex)}
        lec["D_mc"] = dmc
        good = out["estado"] == "PASS" and all(s["estado"] == "PASS" for s in lec["etapas"])
        lec["regla_de_paso"] = "LANZAR_CORRIDA_2" if good and dmc["decision"] == "COMPATIBLE" else "NO_LANZAR"
        lec["regla"] = "la corrida 2 solo se lanza si todos los controles de la corrida 1 pasan y D-mc es COMPATIBLE"
        lec["D_t4"] = out.get("decisiones", {}).get("D_t4")
        if (lec["D_t4"] or {}).get("decision") == "FALLA":
            lec["lectura"] = "DETENIDA POR D-t4: se revisa antes de informar (no cambia la regla de paso)"
        lec["segundos_total"] = round(time.time() - T0_WALL, 3)
        jdump(os.path.join(d, "LECTURA_T2_1_C1.json"), lec)
        print("D-mc: %s; regla de paso: %s" % (dmc["decision"], lec["regla_de_paso"]), flush=True)
        code = 0 if lec["regla_de_paso"] == "LANZAR_CORRIDA_2" else 1
    return code


if __name__ == "__main__":
    main()
