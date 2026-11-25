"""Reproduce the finite SSF flag and high-rate pilots against the paper's reproduction checkout."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
JOBS = [(1,3,1,2000),(1,3,2,96),(1,4,0,96),(1,4,1,96),(2,3,1,2000),(2,4,0,96)]

def invoke(args, cwd, output=None):
    if output is None:
        subprocess.run(args, cwd=cwd, check=True)
    else:
        with output.open('wb') as stream:
            subprocess.run(args, cwd=cwd, stdout=stream, check=True)

def recount(repo, work):
    sys.path.insert(0,str(repo/'sectores'))
    from modelo import build
    from orbit_tables import prepare
    shutil.copyfile(repo/'sectores/censo.cpp',work/'sector_reference.cpp')
    for name in ('censo','bounded_orbits'):
        shutil.copyfile(HERE/(name+'.cpp'),work/(name+'.cpp'))
        if os.name=='nt':
            template=(repo/'sectores/compilar.cmd').read_bytes()
            assert b'censo.exe censo.cpp' in template
            script=work/(name+'_build.cmd')
            script.write_bytes(template.replace(b'censo.exe censo.cpp',f'{name}.exe {name}.cpp'.encode('ascii')))
            invoke([os.environ.get('COMSPEC','cmd.exe'),'/d','/c',str(script.resolve())],work)
        else:
            invoke(['c++','-O3','-std=c++17',f'{name}.cpp','-o',name],work)
    suffix='.exe' if os.name=='nt' else ''
    direct=work/('censo'+suffix)
    orbits=work/('bounded_orbits'+suffix)
    for t in (1,2):
        model=work/f'model_{t}.bin'
        build(t,model)
        invoke([direct,model,'flag'],work,work/f'flag_{t}.txt')
        if t==1:
            invoke([direct,model,'four'],work,work/'four_1.txt')
    invoke([sys.executable,'-B',HERE/'audit_flag_225.py',repo,work],work)
    for t,a,b,budget in JOBS:
        name=f'{t}{a}{b}'
        table=work/f'shape_{name}.bin'
        prepare(repo,table,t,a,b)
        invoke([orbits,work/f'model_{t}.bin',table,str(budget)],work,work/f'orbit_{name}.txt')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repo',type=Path,help='path to the complete ssf-qec-reproduction checkout')
    parser.add_argument('--recount',action='store_true',help='regenerate models and all census and orbit files first')
    parser.add_argument('--output',type=Path,help='save the verified, recomputed JSON at this path')
    args=parser.parse_args()
    repo=args.repo.resolve()
    with tempfile.TemporaryDirectory(prefix='ssf_flag_rates_') as temp:
        work=Path(temp)
        if args.recount:
            recount(repo,work)
        else:
            for f in (HERE/'counts').iterdir():
                shutil.copyfile(f,work/f.name)
        invoke([sys.executable,'-B',HERE/'feasibility.py',repo,work],work,work/'check.log')
        result=json.loads((work/'pilot-results.json').read_bytes())
        ref=json.loads((HERE/'pilot-results.json').read_bytes())
        if result!=ref:
            raise RuntimeError('The recomputed interval or a census differs from the reference JSON.')
        if args.output is not None:
            shutil.copyfile(work/'pilot-results.json',args.output)
        a225=result['flag225_direct']
        a900=result['orbit_bounds']['900_at_reference']
        b=result['orbit_bounds']['225_at_high_rates']
        print('Verified both finite-code pilots.')
        print('225 and 900 conditional risk reduction factors, lower bounds:',
              a225['gain_per_kept_lower'],a900['gain_per_kept_lower'])
        print('225 high-rate silent failure interval and maximum relative error (%):',
              *b['interval'],b['relative_error_percent_upper'])

if __name__=='__main__':
    main()
