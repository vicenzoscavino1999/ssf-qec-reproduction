#!/usr/bin/env python3
"""T2-2, piloto m = 6 (BORRADOR_PLAN_T2_2.txt v4): motor de las cadenas y del ancla, generico en el modelo.

Codigo preparado para revision (version 2: R1 a R3 del dictamen del codigo v1). Solo se ha ejecutado con modelos
de juguete (t2s_juguete.py) en las autopruebas y pruebas de humo; el modelo fisico (t2s_fisica.py, decoder SSF de
B4) no se ha ejecutado.

Interfaz del modelo (la cumplen t2s_juguete.Juguete y t2s_fisica.Fisico):
- NS, ncls[s] >= 1 (clases no nulas de cada sitio), resp[s][c] = (e, sindrome) enteros, resp[s][0] = (0, 0);
- laws(r): leyes exactas por sitio [l(s,0), ..., l(s,ncls[s])] (fracciones) con r_m = r_c = r;
- oracle(e, sindrome, cls): True si la configuracion esta en F (una llamada al decoder en el modelo fisico);
- anchor_sample(rng, Ba, tick): historias del bloque del ancla en orden creciente de posicion h = 1..Ba, como
  (h, pares (s, c) con c != 0 ordenados por s), solo las que tienen algun sitio con fallo.

Paso de Metropolis (2.4, algoritmo 3 de FF por clases): sitio s = randrange(NS); opcion k = randrange(ncls[s])
sobre ({0} U C_s) \\ {E_s} (b = k si k < E_s, si no k + 1); rho = l(s,b)/l(s,E_s) (float correctamente
redondeado desde la fraccion exacta); se acepta por peso si rho >= 1 o random() < rho; si la nueva
configuracion tiene menos de 3 sitios con fallo se queda sin llamar al decoder (atajo, V5); si no, se reserva
la llamada en el contador global y se mueve solo si esta en F. Cada paso da una muestra (el estado tras el
paso); E_0 no es muestra. Flujos: MT19937 de Python sembrado con el SHA-256 completo de la etiqueta.

Presupuesto unico: contadores compartidos (multiprocessing.Value) de llamadas y pasos por etapa, con reserva
previa (techo estricto: el total nunca pasa del techo); plazo comun de pared consultado antes de cada llamada,
al empezar cada bloque y cada 4096 pasos; memoria del proceso medida al empezar cada bloque (una medida que
falla detiene, nunca cuenta como 0); bytes del registro reservados antes de escribir cada bloque.
"""
import hashlib
import math
import os
import random
import sys
import time

import t2s_mat as tm
import t2s_reg as tr

TICK = 4096
GIB = 1 << 30
# Plazo compartido < 0: parada cooperativa ordenada por el conductor despues de otra parada (techo o error de
# otro proceso o del principal). Su motivo es distinto del de un techo propio (R1 del dictamen del codigo v1).
PARADA = -1.0
MOTIVO_PARADA = "parada cooperativa ordenada por el conductor (consecuencia de otra parada; no es un techo propio)"


class Agotado(Exception):
    """Techo de tiempo, llamadas, pasos o memoria: la etapa en curso queda NO_EVALUADA."""


# ------------------------------------------------------------------ flujos (C2: representacion determinista)
def stream(label):
    return random.Random(int.from_bytes(hashlib.sha256(label.encode("ascii")).digest(), "big"))


def label_anchor(master, b):
    return "T2-2|ancla|%d|%d" % (master, b)


def label_chain(master, l, i, j):
    return "T2-2|%d|%d|%d|%d" % (master, l, i, j)


def label_ctiempo(master, rep, j):
    return "T2-2|c-tiempo|%d|%d|%d" % (master, rep, j)


def label_controles(master):
    return "T2-2|controles|%d" % master


# ------------------------------------------------------------------ memoria (Windows; copia de b4.mem_now)
def mem_now():
    """Working set del proceso (bytes); None si no se puede medir."""
    try:
        import ctypes
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]
        k32 = ctypes.windll.kernel32
        k32.GetCurrentProcess.restype = wintypes.HANDLE
        k32.K32GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
        pmc = PMC()
        pmc.cb = ctypes.sizeof(PMC)
        if k32.K32GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb):
            return int(pmc.WorkingSetSize)
    except Exception:
        pass
    return None


# ------------------------------------------------------------------ presupuesto
class Budget:
    """Techos de una etapa. calls/steps: multiprocessing.Value('q') compartidos; caps estrictos."""

    def __init__(self, calls, calls_cap, steps, steps_cap, deadline, mem_cap=None, data=None, clock=None):
        self.calls, self.calls_cap = calls, calls_cap
        self.steps, self.steps_cap = steps, steps_cap
        self._dl, self.mem_cap, self.data = deadline, mem_cap, data
        self.clock = clock or time.time
        self.mem_peak = 0

    @property
    def deadline(self):
        """Plazo vigente: un multiprocessing.Value('d') compartido (el conductor lo pone a PARADA para una parada
        cooperativa) o un numero."""
        return self._dl.value if hasattr(self._dl, "value") else self._dl

    def tick(self):
        dl = self.deadline
        if self.clock() > dl:
            raise Agotado(MOTIVO_PARADA if dl < 0 else "techo de tiempo")

    def reserve_call(self):
        self.tick()
        c = self.calls
        with c.get_lock():
            if c.value >= self.calls_cap:
                raise Agotado("techo de llamadas")
            c.value += 1

    def reserve_steps(self, n):
        self.tick()
        c = self.steps
        if c is None:
            return
        with c.get_lock():
            if c.value + n > self.steps_cap:
                raise Agotado("techo de pasos")
            c.value += n

    def check_mem(self):
        if self.mem_cap is None:
            return
        m = mem_now()
        if m is None:
            raise Agotado("fallo de la medida de memoria del proceso (no medido; nunca cuenta como 0)")
        self.mem_peak = max(self.mem_peak, m)
        if m > self.mem_cap:
            raise Agotado("techo de memoria del proceso: %d bytes" % m)


# ------------------------------------------------------------------ datos por nivel
class LevelData:
    """Leyes exactas de los t + 1 niveles, tabla de claves, tablas de logaritmos y cocientes float de la
    aceptacion por nivel. Se reconstruye igual en cada proceso y se compara con los SHA-256 sellados."""

    def __init__(self, model, rates):
        self.rates = [tm.F(r) for r in rates]
        self.t = len(self.rates) - 1
        self.laws = []
        for r in self.rates:
            L = model.laws(r)
            tm.check_laws(L, model.ncls, "ley en r = %s" % tm.fstr(r))
            self.laws.append(L)
        self.kt = tm.KeyTable(self.laws)
        self.lt = tm.LogTables(self.kt)
        self.rat = []
        for L in self.laws:
            self.rat.append([[[float(row[b] / row[a]) for b in range(len(row))] for a in range(len(row))]
                             for row in L])

    def shas(self):
        return {"claves": self.kt.record()["sha256"], "logaritmos": self.lt.record()["sha256"]}


# ------------------------------------------------------------------ estado y paso de Metropolis (2.4)
class State:
    __slots__ = ("cls", "nf", "N", "sig", "e", "syn")

    def __init__(self, model, key_of, pairs):
        self.cls = [0] * model.NS
        self.N = {}
        self.e = self.syn = 0
        for s, c in pairs:
            if not 1 <= c <= model.ncls[s] or self.cls[s]:
                raise ValueError("par sitio-clase no valido: %r" % ((s, c),))
            self.cls[s] = c
            k = key_of[s][c]
            self.N[k] = self.N.get(k, 0) + 1
            e, y = model.resp[s][c]
            self.e ^= e
            self.syn ^= y
        self.nf = sum(1 for c in self.cls if c)
        self.sig = tuple(sorted(self.N.items()))

    def pairs(self):
        return [[s, c] for s, c in enumerate(self.cls) if c]


class Parcial(Agotado):
    """Agotado dentro de un bloque: lleva los contadores del bloque parcial (no se registra)."""

    def __init__(self, why, partial):
        super().__init__(why)
        self.partial = partial


def run_block(st, nsteps, rat, model, key_of, rng, budget):
    """nsteps pasos desde st (se modifica). Devuelve (histograma de firmas de las muestras, contadores)."""
    budget.reserve_steps(nsteps)
    budget.check_mem()
    NS, ncls, resp, oracle = model.NS, model.ncls, model.resp, model.oracle
    rr, rnd = rng.randrange, rng.random
    cls, N = st.cls, st.N
    sig = st.sig
    hist = {}
    acc = calls = av = trn = 0
    step = 0
    try:
        for step in range(nsteps):
            if not step % TICK:
                budget.tick()
            s = rr(NS)
            a = cls[s]
            k = rr(ncls[s])
            b = k if k < a else k + 1
            rho = rat[s][a][b]
            if rho >= 1.0 or rnd() < rho:
                acc += 1
                nf2 = st.nf + (1 if b else 0) - (1 if a else 0)
                if nf2 < 3:
                    av += 1
                else:
                    budget.reserve_call()
                    calls += 1
                    ea, ya = resp[s][a]
                    eb, yb = resp[s][b]
                    e2, y2 = st.e ^ ea ^ eb, st.syn ^ ya ^ yb
                    cls[s] = b
                    if oracle(e2, y2, cls):
                        st.e, st.syn, st.nf = e2, y2, nf2
                        if a:
                            ka = key_of[s][a]
                            v = N[ka] - 1
                            if v:
                                N[ka] = v
                            else:
                                del N[ka]
                        if b:
                            kb = key_of[s][b]
                            N[kb] = N.get(kb, 0) + 1
                        sig = tuple(sorted(N.items()))
                        trn += 1
                    else:
                        cls[s] = a
            hist[sig] = hist.get(sig, 0) + 1
    except Agotado as ex:
        st.sig = sig
        raise Parcial(str(ex), {"pasos_hechos": step, "aceptaciones": acc, "llamadas": calls, "evitadas": av,
                                "transiciones": trn})
    st.sig = sig
    return hist, {"muestras": nsteps, "aceptaciones": acc, "llamadas": calls, "evitadas": av,
                  "transiciones": trn}


# ------------------------------------------------------------------ cadena con el criterio de FF (2.8)
def run_chain(model, lvd, j, start_pairs, rng, writer, budget, lengths, block, eps):
    """Cadena del nivel j desde start_pairs. Bloques de `block` muestras; criterio de FF tras T_init y tras
    cada ampliacion (lengths). Devuelve el registro de la cadena; los bloques quedan en `writer`.
    Riesgo de memoria (no demostrado; punto 6 del dictamen del codigo v1): `blocks` conserva en memoria todos
    los histogramas de la cadena del nivel hasta su ultimo criterio. El formato compacto del disco no acredita
    que quepa en 1 GiB por proceso; si no cabe, la medida de memoria o el limite duro del Job Object detienen
    la unidad (NO_EVALUADO) y no se pide mas memoria."""
    st = State(model, lvd.kt.key_of, start_pairs)
    e0 = st.pairs()
    blocks, checks, tot = [], [], {"aceptaciones": 0, "llamadas": 0, "evitadas": 0, "transiciones": 0}
    rec = {"nivel": j, "E0": e0, "estado": None}
    t0 = time.perf_counter()
    idx = 0
    try:
        while True:
            target = lengths[idx]
            if target % block:
                raise ValueError("T no es multiplo de B")
            while len(blocks) * block < target:
                hist, cnt = run_block(st, block, lvd.rat[j], model, lvd.kt.key_of, rng, budget)
                writer.write(cnt, hist)
                blocks.append(hist)
                for k in tot:
                    tot[k] += cnt[k]
            crit = tm.ff_criterion(blocks, j, lvd.t, lvd.lt, eps)
            checks.append(crit)
            if crit["cumple"]:
                rec["estado"] = "CRITERIO_FF_CUMPLE"
                break
            idx += 1
            if idx >= len(lengths):
                rec["estado"] = "CRITERIO_FF_NO_CUMPLE"
                break
    except Parcial as ex:
        rec.update(estado="NO_EVALUADO", motivo=str(ex), bloque_parcial_no_registrado=ex.partial)
    except (Agotado, tr.DatosNoCaben) as ex:
        rec.update(estado="NO_EVALUADO", motivo=str(ex),
                   datos_no_caben=isinstance(ex, tr.DatosNoCaben))
    rec.update(T=len(blocks) * block, ampliaciones=max(0, len(checks) - 1), criterios=checks,
               contadores=dict(tot, pasos=len(blocks) * block), final=st.pairs(),
               firma_final=[list(x) for x in st.sig], segundos=round(time.perf_counter() - t0, 3))
    return rec, st


# ------------------------------------------------------------------ bloque del ancla (2.6)
def anchor_block(model, b, Ba, rng, budget):
    """Un bloque del ancla en orden logico: posiciones (1..Ba) y configuraciones de sus silenciosas."""
    pos, cfgs = [], []
    calls = av = nonempty = 0
    by_nf = {}
    cls = [0] * model.NS
    t0 = time.perf_counter()
    last = 0
    try:
        budget.tick()                  # dentro del try: un techo al empezar da un bloque interrumpido
        budget.check_mem()
        for h, pairs in model.anchor_sample(rng, Ba, budget.tick):
            if not last < h <= Ba:
                raise ValueError("posicion %d fuera de orden en el bloque %d" % (h, b))
            last = h
            nf = len(pairs)
            nonempty += 1
            by_nf[nf] = by_nf.get(nf, 0) + 1
            if nf < 3:
                av += 1
                continue
            budget.reserve_call()
            calls += 1
            e = y = 0
            for s, c in pairs:
                cls[s] = c
                ee, yy = model.resp[s][c]
                e ^= ee
                y ^= yy
            silent = model.oracle(e, y, cls)
            for s, c in pairs:
                cls[s] = 0
            if silent:
                pos.append(h)
                cfgs.append([[s, c] for s, c in pairs])
    except Agotado as ex:
        return {"b": b, "interrumpido": str(ex), "llamadas": calls, "evitadas": av,
                "segundos": round(time.perf_counter() - t0, 3), "pid": os.getpid()}
    by_nf[0] = Ba - nonempty
    return {"b": b, "historias": Ba, "posiciones": pos, "configuraciones": cfgs, "llamadas": calls,
            "evitadas": av, "por_sitios_con_fallo": {str(k): v for k, v in sorted(by_nf.items())},
            "segundos": round(time.perf_counter() - t0, 3), "pid": os.getpid()}


# ------------------------------------------------------------------ procesos de trabajo
_W = {}


def make_model(spec):
    kind = spec["tipo"]
    if kind == "juguete":
        import t2s_juguete as tj
        return tj.make(spec)
    if kind == "fisico":
        import t2s_fisica as tf
        return tf.make(spec)
    raise ValueError("modelo desconocido: %r" % kind)


def w_init(cfg, counters):
    sys.set_int_max_str_digits(0)
    model = make_model(cfg["modelo"])
    lvd = LevelData(model, cfg["escalera"]["tasas"])
    got = lvd.shas()
    if got != cfg["escalera"]["shas"]:
        raise RuntimeError("tabla de claves o de logaritmos distinta de la sellada: %r" % got)
    _W.update(cfg=cfg, model=model, lvd=lvd, counters=counters)


def _budget(stage):
    cfg, C = _W["cfg"], _W["counters"]
    caps = cfg["techos"][stage]
    data = tr.Budget(C["bytes"], cfg["techos"]["datos_bytes"]) if stage == "cadenas" else None
    return Budget(C[stage + "_llamadas"], caps["llamadas"], C.get(stage + "_pasos"), caps.get("pasos"),
                  C[stage + "_plazo"], cfg["techos"]["memoria_proceso"], data)


def anchor_task(model, cfg, b, budget):
    """Bloque b del ancla con su flujo declarado; el registro (completo o interrumpido) lleva su pico de memoria
    y la etiqueta del flujo (--verificar la compara con la semilla maestra declarada)."""
    lab = label_anchor(cfg["semilla_maestra"], b)
    r = anchor_block(model, b, cfg["ancla"]["B_a"], stream(lab), budget)
    r["memoria_pico"] = budget.mem_peak
    r["flujo"] = lab
    return r


def w_anchor(b):
    return anchor_task(_W["model"], _W["cfg"], b, _budget("ancla"))


def chain_path(run_dir, j, l, i):
    return os.path.join(run_dir, "cadenas", "n%03d_s%02d_r%d.bin" % (j, l, i))


def registry_path(run_dir, l, i):
    return os.path.join(run_dir, "cadenas", "registro_s%02d_r%d.jsonl" % (l, i))


def w_pipeline(task):
    """Repeticion i de la semilla l: niveles 0..t en orden; el nivel j parte del estado final del j - 1."""
    l, i, seed_pairs, run_dir = task
    cfg, model, lvd = _W["cfg"], _W["model"], _W["lvd"]
    budget = _budget("cadenas")
    lengths = tm.ff_lengths(cfg["cadenas"]["T_init"], cfg["cadenas"]["lambda_ini"],
                            cfg["cadenas"]["ampliaciones_max"])
    eps = tm.F(cfg["cadenas"]["eps"])
    block = cfg["cadenas"]["B"]
    start = seed_pairs
    out = {"semilla": l, "repeticion": i, "niveles": [], "pid": os.getpid()}
    for j in range(lvd.t + 1):
        path = chain_path(run_dir, j, l, i)
        w = tr.BlockWriter(path, (j, l, i), lvd.kt.nkeys, budget.data)
        rng = stream(label_chain(cfg["semilla_maestra"], l, i, j))
        try:
            rec, st = run_chain(model, lvd, j, start, rng, w, budget, lengths, block, eps)
        finally:
            fil = w.close()
        rec.update(semilla=l, repeticion=i, archivo=fil, flujo=label_chain(cfg["semilla_maestra"], l, i, j),
                   memoria_pico=budget.mem_peak, pid=os.getpid())
        tr.append_jsonl(registry_path(run_dir, l, i), rec)
        out["niveles"].append({"nivel": j, "estado": rec["estado"], "T": rec["T"]})
        if rec["estado"] == "NO_EVALUADO":
            out["interrumpido"] = rec.get("motivo")
            out["datos_no_caben"] = rec.get("datos_no_caben", False)
            break
        start = rec["final"]
    out["memoria_pico"] = budget.mem_peak
    return out


def w_ctiempo(task):
    """C-tiempo: n pasos desde la semilla 1 en el nivel j con un flujo propio; valores descartados. Mide
    segundos, llamadas y bytes que ocuparian sus bloques (se codifican y se descartan)."""
    rep, j, seed_pairs, nsteps = task
    cfg, model, lvd = _W["cfg"], _W["model"], _W["lvd"]
    budget = _budget("c_tiempo")
    block = cfg["cadenas"]["B"]
    rng = stream(label_ctiempo(cfg["semilla_maestra"], rep, j))
    st = State(model, lvd.kt.key_of, seed_pairs)
    t0 = time.perf_counter()
    nbytes = calls = done = 0
    try:
        while done < nsteps:
            hist, cnt = run_block(st, block, lvd.rat[j], model, lvd.kt.key_of, rng, budget)
            meta = dict(cnt, nivel=j, semilla=0, repeticion=rep, bloque=done // block)
            nbytes += len(tr.encode_block(meta, hist, lvd.kt.nkeys))
            calls += cnt["llamadas"]
            done += block
    except Agotado as ex:
        return {"repeticion": rep, "nivel": j, "interrumpido": str(ex), "pasos": done}
    return {"repeticion": rep, "nivel": j, "pasos": done, "llamadas": calls, "bytes": nbytes,
            "segundos": time.perf_counter() - t0, "memoria_pico": budget.mem_peak}
