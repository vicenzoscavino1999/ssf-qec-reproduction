"""Corre cada fase desde cero y mide el tiempo de pared del proceso completo."""
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter

HERE=Path(__file__).resolve().parent

def run(phase, cmd, times):
    start=perf_counter()
    with (HERE/(phase+'.log')).open('w',encoding='utf-8') as log:
        with subprocess.Popen(cmd,cwd=HERE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                              text=True,errors='replace') as proc:
            for line in proc.stdout:
                print(line,end='',flush=True)
                log.write(line)
                log.flush()
        if proc.returncode:
            raise RuntimeError(f'{phase} falló con código {proc.returncode}')
    times[phase]=perf_counter()-start
    (HERE/'tiempos_splitting.json').write_text(json.dumps(times,indent=2),encoding='utf-8')

def main():
    times={}
    run('preparacion',[sys.executable,'-I','-B','-X','utf8',str(HERE/'modelo_225.py')],times)
    print('Modelo VV preparado y calendario verificado.',flush=True)
    run('ancla_motor',[str(HERE/'motor_225.exe'),'ancla'],times)
    run('cadenas_motor',[str(HERE/'motor_225.exe'),'cadenas'],times)
    run('analisis',[sys.executable,'-I','-B','-X','utf8',str(HERE/'analizar_225.py')],times)

if __name__=='__main__':
    main()
