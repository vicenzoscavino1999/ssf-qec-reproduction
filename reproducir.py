"""Reproducción de Silent logical failure in small-set-flip decoding, versión limpia."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
PAPER = "Silent_logical_failure_in_small_set_flip_decoding"
CERT = ROOT / "certificados"
COMP = ROOT / "comparaciones"
SPLIT = "SPLITTING_HGP_225_VV_20261124"
DSS = "DSS_HGP_225_BERNOULLI_20261124"
STRIPS = ROOT / "tiras"
WORK = STRIPS / "workspace"
PILOT61 = "T2_2_PILOTO_M6_EJECUCION_20261115/corrida"
CODE61 = WORK / "T2_2_CODIGO_PILOTO_20261115"
ACTIONS = ("rapido","tasas","f17","sectores","piloto","path25","comparaciones",
           "splitting-225","dss-225","splitting-61","tiras-1","tiras-2","catalogos","pdf")

def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))

def scientific_equal(new, ref, keys=()):
    def normalized(value):
        if isinstance(value, dict):
            return {k:normalized(v) for k,v in value.items() if k not in keys}
        if isinstance(value, list):
            return [normalized(v) for v in value]
        return value
    if normalized(read(new)) != normalized(read(ref)):
        raise ValueError(f"Resultado diferente de la referencia: {new.name} / {ref.name}")

def integrity():
    manifest = read(ROOT/"SHA256.json")
    for name, digest in manifest.items():
        path = ROOT/name
        with path.open("rb") as source:
            if hashlib.file_digest(source,"sha256").hexdigest() != digest:
                raise ValueError(f"SHA-256 diferente: {name}")
    print(f"{len(manifest)} ficheros verificados por SHA-256.",flush=True)

def run(action, out, previous):
    logs = []
    python = [sys.executable,"-E","-s","-B","-X","utf8"]
    env = os.environ.copy()
    for name in ("PYTHONPATH", "PYTHONHOME", "PYTHONOPTIMIZE"):
        env.pop(name, None)
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1", PYTHONUTF8="1")

    def execute(title, argv, cwd):
        logfile = out / f"{len(logs)+1:02d}-{title}.log"
        logs.append(logfile.name)
        print(f"  {title}  [{logfile.name}]",flush=True)
        with logfile.open("xb") as log:
            completed = subprocess.run([str(a) for a in argv], cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT)
        if completed.returncode:
            tail = logfile.read_text(encoding="utf-8",errors="replace").splitlines()[-10:]
            raise RuntimeError("\n".join(tail)+f"\nCódigo {completed.returncode}. Registro: {logfile}")

    def py(script, *args, cwd=None):
        execute(script.stem, [*python, script, *args], cwd or script.parent)

    if action=="rapido":
        integrity()
        py(ROOT/"calculos_del_articulo.py")
        py(STRIPS/"tables/reproduce_tables.py", "--output", out/"tablas_tiras")
        py(CERT/"cero/H5A_limite_del_10_por_ciento.py")
        py(CERT/"fibras/verify_fibres.py")
        py(COMP/DSS/"comprobar_tabla.py")
        summary = read(out/logs[0])
        for size, record in summary["Bernoulli"].items():
            lo, hi = record["P_sil"]
            print(f"  {size} cúbits: P_sil en [{lo}, {hi}]",flush=True)
        strips = read(out/"tablas_tiras/POSTPROCESSING.json")
        print(f"  Tiras: {strips['circuit_rows']} filas. Figuras 4–6: "
              f"{summary['coordenadas_cotejadas_figuras_4_5_6']} valores cotejados.",flush=True)
    elif action=="tasas":
        py(CERT/"f17/verify.py")
        for name in ("verify_first_contact","verify_fibres","verificar_FASE5_circuito_propagacion",
                     "verificar_H5A_3_4_primer_apagado_de_ternas","verificar_H5A_3_4_esquina_2",
                     "verificar_H5A_3_4_esquina_3","verificar_H5A_3_4_esquina_3_circuito",
                     "verificar_FASE6_F17_circuito"):
            args = ["--global"] if name.endswith("apagado_de_ternas") else []
            py(CERT/"fibras"/(name+".py"), *args)
    elif action=="f17":
        for block, name, reference in (("h1","verificar_h1_cierre","RESULTADO_H1_CIERRE.json"),
                                       ("h3","verificar_h3","RESULTADO_H3.json")):
            folder = CERT/"supplement"/block
            py(folder/(name+".py"), "--output", out/reference)
            scientific_equal(out/reference,folder/reference)
        solver = CERT/"sector_f17/finite_sector.py"
        py(solver, "produce", out/"certificate.jsonl")
        py(solver, "verify", out/"certificate.jsonl")
        py(ROOT/"calculos_del_articulo.py")
    elif action=="sectores":
        local = Path(shutil.copytree(ROOT/"sectores",out/"codigo"))
        py(local/"reproducir.py","--out",out/"censo")
        py(local/"contrastar.py",out/"censo",ROOT/"sectores/referencia")
    elif action=="piloto":
        name = "piloto_HGP_225_esquinas_20261124.json"
        local = out/"piloto"
        local.mkdir()
        for f in COMP.glob("*.py"):
            shutil.copy2(f,local)
        py(local/"esquinas_HGP_225.py","--salida",local/name)
        scientific_equal(local/name, COMP/name,
                         ("preparacion_s","total_pared_s","pico_MiB","evaluacion_s"))
        py(local/"auditar_esquinas_HGP_225.py","--carpeta",local)
    elif action=="path25":
        local = Path(shutil.copytree(ROOT/"path25/validacion/run",out/"validacion"))
        py(local/"supervisor.py")
        scientific_equal(local/"attempt_01/RESULTS.json", ROOT/"path25/validacion/expected/RESULTS.json",
                         ("elapsed_seconds","timings_seconds"))
        py(ROOT/"path25/estimadores/reproduce_h6.py","--out",out/"estimadores")
    elif action=="comparaciones":
        py(COMP/DSS/"comprobar_tabla.py")
        py(COMP/SPLIT/"verificar_resultado.py")
        py(CODE61/"t2_2.py","--verificar",WORK/PILOT61,"--salida",out/"verificacion61.json",cwd=WORK)
    elif action in ("splitting-225","dss-225"):
        local = Path(shutil.copytree(COMP,out/"comparaciones"))
        folder=local/(SPLIT if action=="splitting-225" else DSS)
        if os.name=="nt":
            compile_cmd = ["cmd","/d","/c","compilar.cmd"]
        else:
            compile_cmd = ["g++","-O3","-std=c++17","-mbmi2","-mpopcnt","estrato_225.cpp","-o","estrato_225"]
        execute("compilar",compile_cmd,folder)
        if action=="splitting-225":
            for name in ("modelo.bin","modelo.json","ancla.txt","cadenas.tsv","histogramas.tsv",
                         "resultado_splitting.json","tiempos_splitting.json","certificado_225_tiempos.json"):
                (folder/name).unlink()
            py(folder/"ejecutar_225.py")
            py(folder/"verificar_resultado.py")
            original = read(COMP/SPLIT/"resultado_splitting.json")
            fresh = read(folder/"resultado_splitting.json")
            for key in ("estimacion","IC_nominal_95","transporte_por_semilla","diagnosticos"):
                if fresh[key]!=original[key]:
                    raise ValueError(f"Splitting 225: difiere {key}.")
            py(folder/"medir_certificado.py")
        else:
            (folder/"modelo.bin").unlink()
            execute("preparar-modelo",[*python,"-c",
                    "import runpy,sys; runpy.run_path(sys.argv[1])['build'](sys.argv[2])",
                    local/SPLIT/"modelo_225.py",folder],folder)
            result = out/"dss225.json"
            py(folder/"ejecutar.py","--output",result)
            scientific_equal(result,COMP/DSS/"resultado.json",("process_s","load_s","evaluation_s"))
    elif action=="splitting-61":
        py(CODE61/"t2s_supervisor.py","--corrida","--dir",out/"corrida",cwd=WORK)
    elif action in ("tiras-1","tiras-2","catalogos"):
        command = {"tiras-1":"run-1", "tiras-2":"run-2", "catalogos":"rebuild-engine"}[action]
        args = [command,"--out",out/"calculo","--execute-production"]
        if action=="tiras-2":
            args += ["--previous",previous/"calculo"]
        py(STRIPS/"reproduce.py",*args)
    elif action=="pdf":
        executable = shutil.which("pdflatex")
        if not executable:
            raise ValueError("Instale pdfLaTeX (MiKTeX o TeX Live) para esta acción.")
        src = Path(shutil.copy2(ROOT/"manuscrito"/(PAPER+".tex"),out))
        for i in range(3):
            execute(f"pdflatex-{i+1}",[executable,"-interaction=nonstopmode","-halt-on-error",src.name],out)
        log = src.with_suffix(".log").read_text(encoding="utf-8",errors="replace")
        if any(t in log for t in ("There were undefined references","There were undefined citations",
                                  "Rerun to get cross-references right","Label(s) may have changed")):
            raise ValueError("Persisten referencias pendientes en el log final de LaTeX.")
    return logs

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("accion",choices=ACTIONS,nargs="?",default="rapido")
    parser.add_argument("--salida",type=Path,help="Carpeta nueva para este bloque. Por defecto fuera del paquete.")
    parser.add_argument("--previa",type=Path,help="Carpeta de salida del bloque tiras-1.")
    args = parser.parse_args()
    if sys.version_info<(3,12) or sys.maxsize<=2**32 or not __debug__:
        parser.error("Se requiere Python 3.12+ de 64 bits con las aserciones activas.")
    if args.accion in ("path25","splitting-61","splitting-225","tiras-1","tiras-2","catalogos") and os.name!="nt":
        parser.error("Este productor utiliza la configuración original de Windows x64. Consulte LEEME.md.")
    if (args.previa is not None)!=(args.accion=="tiras-2"):
        parser.error("--previa se requiere exclusivamente para tiras-2.")
    if args.previa and not (args.previa/"calculo/LECTURA_T2_1_C1.json").is_file():
        parser.error("La corrida previa debe contener calculo/LECTURA_T2_1_C1.json.")
    if args.salida:
        out=args.salida.resolve()
        if out.is_relative_to(ROOT):
            parser.error("Elija una carpeta nueva fuera del paquete para las salidas.")
        if out.exists():
            parser.error("La carpeta de salida ya existe. Elija otra.")
        out.mkdir(parents=True)
    else:
        parent=ROOT.parent/"resultados_ssf"
        parent.mkdir(exist_ok=True)
        out=Path(tempfile.mkdtemp(prefix=args.accion+"-",dir=parent))
    print(f"{args.accion} -> {out}",flush=True)
    start=time.monotonic()
    result={"accion":args.accion, "python":platform.python_version(), "plataforma":platform.platform()}
    try:
        result["registros"]=run(args.accion,out,args.previa.resolve() if args.previa else None)
        result["estado"]="OK"
    except KeyboardInterrupt:
        result.update(estado="INTERRUMPIDO", detalle="Interrumpido por el usuario.")
    except Exception as error:
        result.update(estado="ERROR", detalle=str(error))
        print(error,file=sys.stderr)
    finally:
        result["segundos"]=round(time.monotonic()-start,3)
        (out/"resultado.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(f"{result.get('estado','INTERRUMPIDO')}  {out/'resultado.json'}",flush=True)
    return 0 if result.get("estado")=="OK" else 1

if __name__=="__main__":
    raise SystemExit(main())

