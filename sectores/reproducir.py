"""Cold, deterministic run from the seed. Python 3.10+ and a C++17 compiler."""
import argparse
import json
import os
import platform
import subprocess
import sys
from pathlib import Path
from time import perf_counter

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import modelo
import orbitas
from intervalo import certificate


def peak_bytes():
    if os.name=='nt':
        import ctypes as c
        from ctypes import wintypes as w
        class Counters(c.Structure):
            _fields_=[('cb',w.DWORD),('faults',w.DWORD)]+[(key,c.c_size_t) for key in
                     ('peak','current','pp','p','np','n','page','peak_page')]
        out=Counters()
        out.cb=c.sizeof(out)
        handle=c.windll.kernel32.GetCurrentProcess
        handle.restype=w.HANDLE
        func=c.windll.psapi.GetProcessMemoryInfo
        func.argtypes=[w.HANDLE,c.POINTER(Counters),w.DWORD]
        if not func(handle(),c.byref(out),out.cb):
            raise c.WinError()
        return out.peak
    import resource
    v=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return v if sys.platform=='darwin' else 1024*v


def compile_if_needed():
    exe=HERE/('censo.exe' if os.name=='nt' else 'censo')
    if not exe.exists() or exe.stat().st_mtime<(HERE/'censo.cpp').stat().st_mtime:
        cmd=['cmd','/d','/c','compilar.cmd'] if os.name=='nt' else [
             'c++','-O3','-std=c++17','censo.cpp','-o',str(exe)]
        subprocess.run(cmd,cwd=HERE,check=True)
    return exe


def run(t,out,exe):
    start=perf_counter()
    model=out/f'semilla_{t}.bin'
    meta=modelo.build(t,model)
    meta['geometry_s']=perf_counter()-start
    zero,one=(out/f'orb_{t}_{k}.bin' for k in (0,1))
    ts=perf_counter()
    orbitas.prepare(t,zero,True)
    orbitas.prepare(t,one)
    meta['orbit_preparation_s']=perf_counter()-ts
    counts,marks={},{}
    total_calls=0
    comparisons=probes=0
    child_peaks=[]
    for name,args in (('pairs',[]),('triples',[zero]),('one_readout',[one])):
        result=subprocess.run([str(exe),str(model),*map(str,args)],stdout=subprocess.PIPE,
                              text=True,encoding='utf-8',check=True)
        (out/f'{meta["qubits"]}_{name}.txt').write_text(result.stdout,encoding='utf-8')
        for line in result.stdout.splitlines():
            fields=line.split()
            if line.startswith('2 '):
                ab=tuple(map(int,fields[:2]))
                counts[ab]=int(fields[4]);marks[str(ab)]=fields[9] if len(fields)>9 else ''
                total_calls+=int(fields[3]);comparisons+=int(fields[6]);probes+=int(fields[7])
            elif line.startswith('# total_s'):
                child_peaks.append(int(fields[4]))
            elif line.startswith(('triples ','orbits ')):
                population,value,representatives,comp,pro=map(int,fields[1:6])
                comparisons+=comp;probes+=pro;child_peaks.append(int(fields[7]))
                total_calls+=int(fields[8])
                if name=='triples':
                    counts[3,0]=value
                    marks['(3, 0)']=fields[9] if len(fields)>9 else ''
                    meta['no_readout_orbit_representatives']=representatives
                    # Each orbit is first tested over all permitted size/score ties.
                    total_calls+=representatives
                else:
                    assert population==meta['qubits']*(meta['qubits']-1)*(meta['qubits']-2)//6*meta['checks']
                    exceptional=value
                    meta['one_readout_orbit_representatives']=representatives
    cert=certificate(meta,counts,exceptional)
    for ab,txt in marks.items():
        import ast
        entries=[tuple(map(int,entry.split(':'))) for entry in txt.split(',') if entry]
        assert all(mark>0 and number>0 for mark,number in entries)
        assert len(set(mark for mark,_ in entries))==len(entries)
        assert sum(number for _,number in entries)==counts[ast.literal_eval(ab)]
    data={'geometry':meta,'p':'0.000001','q':'0.0001','interval':cert,
          'source_configurations_or_orbit_queries':total_calls,'winner_queries':comparisons,
          'local_table_probes':probes,'preparation_python_peak_bytes':peak_bytes(),
          'solver_cpp_peak_bytes':max(child_peaks),'cold_seconds':perf_counter()-start,
          'python':platform.python_version(),'platform':platform.platform()}
    (out/f'resultado_{meta["qubits"]}.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(meta['qubits'],cert['lower'],cert['upper'],cert['max_relative_percent'],data['cold_seconds'],flush=True)


if __name__=='__main__':
    if not __debug__:
        raise RuntimeError('Keep Python assertions enabled.')
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sheets',type=int,choices=(1,2))
    ap.add_argument('--out',type=Path,default=HERE/'run')
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    exe=compile_if_needed()
    if args.sheets:
        run(args.sheets,args.out.resolve(),exe)
    else:
        for t in (1,2):
            subprocess.run([sys.executable,'-B',str(Path(__file__).resolve()),
                            '--sheets',str(t),'--out',str(args.out.resolve())],check=True)

