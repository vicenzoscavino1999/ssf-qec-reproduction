#!/usr/bin/env python3
"""T2-2, piloto m = 6 (BORRADOR_PLAN_T2_2.txt v4, 6 y 7.1): supervisor de la corrida en un arbol de procesos.

Codigo preparado para revision (version 2: R3 del dictamen del codigo v1, disco tras la salida). Solo se ha
ejecutado con procesos sinteticos (t2s_pruebas.py, prueba del supervisor): nunca con el conductor de T2-2.

Contrato:
- JobTree se copia de t2_intervalos.py (T2-1, version revisada en DICTAMEN_N1_V2_N3_V3) y se le anaden los
  limites de memoria del trabajo: JOB_OBJECT_LIMIT_PROCESS_MEMORY (1 GiB por proceso) y JOB_OBJECT_LIMIT_JOB_MEMORY
  (8 GiB en total), ademas de KILL_ON_JOB_CLOSE. Son limites de memoria comprometida aplicados por Windows (una
  reserva que los pasa falla dentro del proceso); el sondeo mide ademas el working set de cada proceso y del
  total. El proceso se crea suspendido, entra en el trabajo antes de reanudarse y sus descendientes quedan dentro.
  Cada consulta y cada orden comprueba su retorno: un fallo se registra con GetLastError y detiene; nunca se lee
  como cero procesos ni memoria nula. Solo se actua sobre los procesos de este trabajo.
- Fase 1: t2_2.py --corrida --dir D --t0 T0, con plazo duro T0 + 8 h - 30 min (la reserva de --verificar) y plazo
  de aceptacion = plazo duro - 60 s de limpieza. Fase 2 (solo si la fase 1 termina por si misma, con el arbol
  vacio sin forzarlo y dentro de plazo): t2_2.py --verificar D, con plazo duro T0 + 8 h. Cada espera dura como
  mucho POLL = 5 s y nunca pasa del plazo de aceptacion.
- En cada sondeo: memoria total del arbol, memoria de cada proceso y bytes de la carpeta de la corrida (techo de
  disco 2e9). Un techo, el plazo de aceptacion o un fallo de supervision terminan el arbol entero. Tras una
  salida normal y el cierre del arbol se mide otra vez el disco: sobre el techo, o una medida fallida, no se
  acepta la fase.
- La corrida solo queda ACEPTADA si las dos fases salen por si mismas dentro de plazo, sus arboles quedan vacios
  sin forzarlos (consulta valida de 0 procesos activos), ninguna consulta ni orden del cierre falla y --verificar
  devuelve 0. Si no, NO_EVALUADO con el motivo, todo lo escrito conservado y pendiente de inspeccion, sin
  reintento.

Uso (desde D:\\paper ia cuantica\\qec_revision_matematica\\05_COMPROBACIONES):
  python -I -B T2_2_CODIGO_PILOTO_20261115/t2s_supervisor.py --corrida --dir <carpeta nueva>
"""
import argparse
import ctypes
import json
import os
import platform
import subprocess
import sys
import time

UDIR = "T2_2_CODIGO_PILOTO_20261115/"
sys.path.insert(0, os.path.abspath(UDIR))
import t2s_reg as tr  # noqa: E402

GIB = 1 << 30
POLL = 5.0
CLEAN = 60.0
LIM = {"global_s": 8 * 3600.0, "verificar_s": 1800.0, "memoria_proceso": GIB, "memoria_total": 8 * GIB,
       "disco_bytes": 2_000_000_000}


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


class _PMC(ctypes.Structure):
    _fields_ = [("cb", ctypes.c_uint32), ("PageFaultCount", ctypes.c_uint32),
                ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]


class JobQueryError(OSError):
    """Fallo de una consulta o de una orden al Job Object, o lista de procesos incompleta. No es una medida."""


class JobTree:
    """Copia de t2_intervalos.JobTree (T2-1) con limites de memoria por proceso y del trabajo."""

    def __init__(self, mem_proc=None, mem_job=None):
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
        flags = 0x2000                                          # KILL_ON_JOB_CLOSE
        if mem_proc:
            flags |= 0x100                                      # JOB_OBJECT_LIMIT_PROCESS_MEMORY
            info.ProcessMemoryLimit = int(mem_proc)
        if mem_job:
            flags |= 0x200                                      # JOB_OBJECT_LIMIT_JOB_MEMORY
            info.JobMemoryLimit = int(mem_job)
        info.BasicLimitInformation.LimitFlags = flags
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
            p.kill()
            p.wait()
            raise
        self.p = p
        self.seen.add(p.pid)
        return p

    def pids(self):
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
        acc = _Accounting()
        if not self.k32.QueryInformationJobObject(self.h, 1, ctypes.byref(acc), ctypes.sizeof(acc), None):
            self._fail("QueryInformationJobObject(JobObjectBasicAccountingInformation)")
        return int(acc.ActiveProcesses), int(acc.TotalProcesses)

    def memory_by_pid(self):
        """Working set de cada proceso activo del trabajo. Un proceso que ya salio (OpenProcess con
        ERROR_INVALID_PARAMETER) no ocupa memoria y se omite; cualquier otro fallo lanza JobQueryError."""
        out = {}
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
                out[pid] = int(pmc.WorkingSetSize)
            finally:
                self.k32.CloseHandle(h)
        return out

    def peak_commit(self):
        info = _ExtLimit()
        if not self.k32.QueryInformationJobObject(self.h, 9, ctypes.byref(info), ctypes.sizeof(info), None):
            self.errors.append({"funcion": "QueryInformationJobObject(JobObjectExtendedLimitInformation)",
                                "error": ctypes.get_last_error(), "tipo_error": "GetLastError",
                                "diagnostico": True})
            return None
        return {"trabajo": int(info.PeakJobMemoryUsed), "proceso": int(info.PeakProcessMemoryUsed)}

    def terminate(self):
        if not self.k32.TerminateJobObject(self.h, 1):
            self._fail("TerminateJobObject")

    def wait_empty(self, wait_s):
        t = time.time()
        while True:
            if self.counts()[0] == 0:
                return True
            if time.time() - t >= wait_s:
                return False
            time.sleep(0.1)

    def close(self):
        if not self.h:
            return True
        ok = bool(self.k32.CloseHandle(self.h))
        if not ok:
            self.errors.append({"funcion": "CloseHandle(trabajo)", "error": ctypes.get_last_error(),
                                "tipo_error": "GetLastError"})
        self.h = None
        return ok


def entorno():
    w = sys.getwindowsversion()
    return {"sistema": platform.platform(), "windows": [w.major, w.minor, w.build],
            "python": platform.python_version(), "bits": platform.architecture()[0],
            "interfaz": "Job Object (kernel32) y NtResumeProcess (ntdll); acreditada solo en este entorno"}


def close_tree(job, p, tree, normal, hard, clock):
    """Cierre del arbol (copia de t2_intervalos.close_tree, mismas reglas). Devuelve (confirmado, problema)."""
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
                if empty is not True:
                    force = True
                    if empty is False:
                        n = attempt("recuento de descendientes", lambda: job.counts()[0])
                        tree["descendientes_vivos_tras_la_salida"] = n
                        attempt("registro de PID de los descendientes", job.pids)
                        probs.append("el proceso salio, pero %s descendientes seguian vivos y hubo que terminar el "
                                     "arbol por fuerza; salida conservada, pendiente de inspeccion" % n)
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
            act, total = cnt if cnt else (None, None)
            tree.update(procesos_activos_al_final=act, procesos_totales=total)
            confirmed = empty is True and act == 0
            if empty is False:
                probs.append("el arbol sigue con procesos activos tras terminarlo")
        elif job is not None:
            cnt = attempt("recuento final", job.counts)
            tree["procesos_activos_al_final"] = cnt[0] if cnt else None
            confirmed = cnt is not None and cnt[0] == 0
    finally:
        if job is not None:
            try:
                tree["pico_memoria_comprometida"] = job.peak_commit()
            except Exception:
                tree["pico_memoria_comprometida"] = None
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
                probs.append("el proceso seguia vivo al cerrar el trabajo; terminado con kill")
    tree["cierre_confirmado"] = confirmed
    tree["problemas_de_cierre"] = probs
    return confirmed, "; ".join(probs) or None


def supervise_phase(name, cmd, log_path, hard, lim, run_dir, job_factory=None, clock=None, poll=POLL,
                    clean=CLEAN):
    """Una fase en su arbol propio. Devuelve el registro con 'aceptada' True/False y el motivo."""
    clock = clock or time.time
    job_factory = job_factory or (lambda: JobTree(lim["memoria_proceso"], lim["memoria_total"]))
    t = clock()
    stop_at = hard - min(clean, max(0.0, (hard - t) / 2.0))
    tree = {"fase": name, "plazo_duro_s_desde_inicio": round(hard - t, 3),
            "plazo_aceptacion_s_desde_inicio": round(stop_at - t, 3), "entorno": entorno(),
            "nota": "Job Object propio; null = no medido (una medida que falla nunca es 0)"}
    reason, peak, peak_proc, peak_disk, p, job, t_exit = None, None, None, None, None, None, None
    with open(log_path, "w", encoding="utf-8", newline="\n") as lf:
        try:
            try:
                job = job_factory()
                p = job.start(cmd, lf)
            except Exception as ex:
                reason = "no se pudo iniciar el arbol (no se ejecuto): %s: %s" % (type(ex).__name__, ex)
            if p is not None:
                tree["pid"] = p.pid
                while True:
                    left = stop_at - clock()
                    if left <= 0:
                        reason = "plazo de aceptacion vencido"
                        tree["parada_s_desde_inicio"] = round(clock() - t, 3)
                        break
                    try:
                        p.wait(timeout=min(poll, left))
                        t_exit = clock()
                        break
                    except subprocess.TimeoutExpired:
                        pass
                    mem = job.memory_by_pid()
                    tot = sum(mem.values())
                    mx = max(mem.values()) if mem else 0
                    dk = tr.dir_bytes(run_dir)
                    peak = tot if peak is None else max(peak, tot)
                    peak_proc = mx if peak_proc is None else max(peak_proc, mx)
                    peak_disk = dk if peak_disk is None else max(peak_disk, dk)
                    if tot > lim["memoria_total"]:
                        reason = "techo de memoria total: %d bytes en el arbol" % tot
                    elif mx > lim["memoria_proceso"]:
                        reason = "techo de memoria por proceso: %d bytes" % mx
                    elif dk > lim["disco_bytes"]:
                        reason = "techo de disco: %d bytes en la carpeta de la corrida" % dk
                    if reason:
                        tree["parada_s_desde_inicio"] = round(clock() - t, 3)
                        break
                if t_exit is not None:
                    tree["salida_observada_s_desde_inicio"] = round(t_exit - t, 3)
                    if t_exit > stop_at:
                        reason = "salida observada %.3f s despues del plazo de aceptacion" % (t_exit - stop_at)
        except Exception as ex:
            reason = "supervision interrumpida: %s: %s" % (type(ex).__name__, ex)
        finally:
            confirmed, issue = close_tree(job, p, tree, t_exit is not None and reason is None, hard, clock)
    if t_exit is not None and reason is None:
        # R3: disco tras las escrituras finales y el cierre del arbol, antes de aceptar una salida normal (el
        # sondeo periodico no ve un aumento justo antes de la salida)
        try:
            dk = tr.dir_bytes(run_dir)
            tree["disco_tras_la_salida_bytes"] = dk
            peak_disk = dk if peak_disk is None else max(peak_disk, dk)
            if dk > lim["disco_bytes"]:
                reason = "techo de disco tras la salida: %d bytes en la carpeta de la corrida" % dk
        except OSError as ex:
            tree["disco_tras_la_salida_bytes"] = None
            reason = "medida de disco tras la salida fallida (no medida; nunca cuenta como 0): %s" % ex
    end = clock()
    tree.update(pids=sorted(job.seen) if job else [], pico_memoria_total_bytes=peak,
                pico_memoria_proceso_bytes=peak_proc, pico_disco_bytes=peak_disk,
                errores_api=list(job.errors) if job else [], segundos=round(end - t, 3),
                codigo=None if p is None else p.returncode)
    op = "; ".join(x for x in (reason, issue) if x) or None
    if op is None and not confirmed:
        op = "cierre del arbol no confirmado"
    if op is None and end > hard:
        op = "cierre terminado %.3f s despues del plazo duro" % (end - hard)
    tree["aceptada"] = op is None
    tree["motivo"] = op
    return tree


def supervise_run(run_dir, t0, lim=LIM, cmd_run=None, cmd_verify=None, job_factory=None, clock=None, poll=POLL):
    clock = clock or time.time
    hard1 = t0 + lim["global_s"] - lim["verificar_s"]
    hard2 = t0 + lim["global_s"]
    cmd_run = cmd_run or [sys.executable, "-I", "-B", UDIR + "t2_2.py", "--corrida", "--dir", run_dir,
                          "--t0", repr(t0)]
    cmd_verify = cmd_verify or [sys.executable, "-I", "-B", UDIR + "t2_2.py", "--verificar", run_dir]
    out = {"unidad": "T2-2 supervision", "t0": t0, "techos": lim}
    f1 = supervise_phase("corrida", cmd_run, os.path.join(run_dir, "CONSOLA_CONDUCTOR.txt"), hard1, lim,
                         run_dir, job_factory, clock, poll)
    out["corrida"] = f1
    if f1["aceptada"]:
        f2 = supervise_phase("verificar", cmd_verify, os.path.join(run_dir, "CONSOLA_VERIFICAR.txt"), hard2, lim,
                             run_dir, job_factory, clock, poll)
        out["verificar"] = f2
    else:
        out["verificar"] = {"aceptada": False, "motivo": "no se lanzo: la fase de la corrida no fue aceptada"}
    v = out["verificar"]
    if f1["aceptada"] and v["aceptada"] and f1["codigo"] == 0 and v.get("codigo") == 0:
        out["estado"] = "ACEPTADA"
    else:
        out["estado"] = "NO_EVALUADO"
        out["motivo"] = "; ".join(x for x in (f1.get("motivo"), v.get("motivo")) if x) or (
            "codigo de retorno %s / %s" % (f1.get("codigo"), v.get("codigo")))
        out["pendiente_de_inspeccion"] = True
        out["nota"] = "todo lo escrito se conserva; sin reintento automatico"
    out["segundos"] = round(clock() - t0, 3)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--corrida", action="store_true")
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    if not a.corrida:
        raise SystemExit("uso: --corrida --dir <carpeta nueva>")
    if os.path.exists(a.dir):
        raise SystemExit("--dir ya existe")
    t0 = time.time()
    import t2_2  # noqa: E402  (solo para comprobar el sello antes de lanzar nada)
    try:
        t2_2.check_seal()
    except t2_2.CheckFail as ex:
        raise SystemExit(str(ex))
    os.makedirs(a.dir)
    out = supervise_run(a.dir, t0)
    tr.write_json(os.path.join(a.dir, "SUPERVISION_T2_2.json"), out)
    print("supervision: %s" % out["estado"], flush=True)
    raise SystemExit(0 if out["estado"] == "ACEPTADA" else 1)


if __name__ == "__main__":
    main()
