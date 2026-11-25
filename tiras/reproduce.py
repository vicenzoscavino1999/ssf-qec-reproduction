"""Portable T2-1 entry point. No scientific work is run without a subcommand."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
WORK = ROOT / 'workspace'
UNIT = 'T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114'
MOTOR = 'T2_1M_MOTOR_ANCLAS_PUNTO_FIJO_20261113'
C1 = 'T2_1_CORRIDA_1_EJECUCION_20261114/corrida_1/RESULTADO_T2_1_C1.json'
C2 = 'T2_1_CORRIDA_2_EJECUCION_20261115/corrida_2/RESULTADO_T2_1.json'
SOLO = 'T2_1_PRUEBAS_ACOTADAS_EJECUCION_20261114/RESULTADO_SOLO_U.json'

def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def integrity():
    manifest = json.loads((ROOT/'MANIFEST.json').read_text(encoding='utf-8'))
    total = 0
    for i, row in enumerate(manifest['files'], 1):
        path = (ROOT/row['path']).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise ValueError('Missing or invalid path: '+row['path'])
        if path.stat().st_size != row['bytes'] or sha(path) != row['sha256']:
            raise ValueError('Integrity mismatch: '+row['path'])
        total += row['bytes']
        if i % 300 == 0: print('Checked',i,'files',flush=True)
    # Every catalogue entry must correspond to an included, hashed shard.
    items = {r['path']:r for r in manifest['files']}
    prefix = 'workspace/'+MOTOR+'/trabajo/'
    catalogues = json.loads((ROOT/prefix/'MANIFIESTO.json').read_text())
    count = 0
    for tag, groups in catalogues.items():
        for v in groups.values():
            r = items[prefix+tag+'/'+v['archivo']]
            if r['sha256'] != v['sha256'] or r['bytes'] != v['bytes']:
                raise ValueError('Catalogue manifest mismatch')
            count += 1
    if count != 1412: raise ValueError('Expected 1412 catalogue shards')
    print(f'INTEGRITY PASS: {len(items)} files, {total} bytes, {count} catalogue shards',flush=True)

def execute(script, args):
    command = [sys.executable,'-I','-B',str(WORK/script),*map(str,args)]
    return subprocess.run(command,cwd=WORK,check=False).returncode

def new_output(text):
    if not text: raise ValueError('--out is required')
    p = Path(text).resolve()
    if p == ROOT or p.is_relative_to(ROOT): raise ValueError('Output must be outside the release')
    if p.exists(): raise ValueError('Output must be new; nothing is overwritten')
    p.parent.mkdir(parents=True,exist_ok=True)
    return p

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command',choices=['integrity','tables','verify-recorded','run-1','run-2','rebuild-engine'])
    ap.add_argument('--out',help='New output folder outside the release')
    ap.add_argument('--previous',help='Completed run-1 folder, for run-2')
    ap.add_argument('--execute-production',action='store_true',help='Explicit opt-in to expensive scientific work')
    a = ap.parse_args()
    if sys.version_info < (3,12): raise ValueError('Use Python 3.12 or later')
    production = a.command in {'run-1','run-2','rebuild-engine'}
    if production and not a.execute_production:
        raise ValueError('Production was NOT started. Add --execute-production only when intending a full run')
    if a.command not in {'integrity','tables'} and os.name != 'nt':
        raise ValueError('Original production and verification environment: Windows 64-bit')
    if sys.maxsize <= 2**32: raise ValueError('64-bit Python is required')
    if not __debug__: raise ValueError('Do not disable assertions with -O')
    out = None if a.command in {'integrity','verify-recorded'} else new_output(a.out)
    previous = None
    if a.command == 'run-2':
        if not a.previous: raise ValueError('--previous is required for run-2')
        previous = Path(a.previous).resolve()
        if not (previous/'LECTURA_T2_1_C1.json').is_file(): raise ValueError('Missing run-1 reading')
    integrity()
    if a.command == 'integrity': return 0
    if a.command == 'tables':
        return subprocess.run([sys.executable,'-I','-B',str(ROOT/'tables/reproduce_tables.py'),
                               '--output',str(out)],cwd=ROOT,check=False).returncode
    if a.command == 'verify-recorded':
        jobs = [(UNIT+'/t2_solo_u.py',['--verificar',WORK/SOLO]),
                (UNIT+'/t2_intervalos.py',['--verificar',WORK/C1]),
                (UNIT+'/t2_intervalos.py',['--verificar',WORK/C2])]
        for script, args in jobs:
            rc = execute(script,args)
            if rc: return rc
        print('RECORDED RESULTS PASS; no decoder, DP or sampling rerun')
        return 0
    if a.command in {'run-1','run-2'}:
        args = ['--corrida',a.command[-1],'--dir',out]
        if previous: args += ['--previa',previous]
        return execute(UNIT+'/t2_intervalos.py',args)
    out.mkdir()
    code = execute(MOTOR+'/t2_motor.py',['--output',out/'RESULTADO_T2_1M.json',
                   '--recursos',out/'RECURSOS_T2_1M.json','--trabajo',out/'trabajo'])
    if code: return code
    result = json.loads((out/'RESULTADO_T2_1M.json').read_text())
    if result.get('status') != 'PASS' or result.get('decision_D_m1') != 'VALIDADO':
        print('Engine did not finish all required checks; evidence preserved')
        return 2
    old = json.loads((WORK/MOTOR/'trabajo/MANIFIESTO.json').read_text())
    new = json.loads((out/'trabajo/MANIFIESTO.json').read_text())
    if new != old:
        print('Fresh catalogue manifest differs; inspect before reuse, no automatic replacement')
        return 2
    print('ENGINE PASS; fresh catalogue manifest matches the supplied manifest')
    return 0

if __name__ == '__main__':
    try: raise SystemExit(main())
    except (OSError,ValueError,KeyError) as ex:
        print('STOP:',ex,file=sys.stderr)
        raise SystemExit(1)
