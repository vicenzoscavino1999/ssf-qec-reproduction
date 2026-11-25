"""Recompute tables and scalar certificates from frozen outputs; no decoding.

Python standard library only. This is postprocessing, not a rerun of the
catalogue producers, dynamic program or Monte Carlo.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, sys
sys.set_int_max_str_digits(100000)
ROOT=Path(__file__).resolve().parent
def load(name): return json.loads((ROOT/'results'/name).read_text(encoding='utf-8-sig'))
def up(x,d=3):
    x=F(x); n=10**d
    a=-((-x.numerator*n)//x.denominator)
    return f'{a//n}.{a%n:0{d}d}'
def close_decimal(value,exact):
    # Stored decimal summaries have 25 significant digits. All tabulated
    # bounds are recomputed from the underlying rational endpoints.
    return abs(F(str(value))-exact)<=F(1,10**23)*max(F(1),abs(exact))
def main(out):
    checks=0
    def check(ok,msg):
        nonlocal checks
        if not ok: raise AssertionError(msg)
        checks+=1
    manifest=ROOT/'MANIFEST.json'
    if manifest.exists():
        for row in json.loads(manifest.read_text()):
            check(hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path'])
    circuit=load('circuit.json'); ph=load('phenomenological.json')
    errs={};baseerrs={}; costs=[];adv=0;pass_low=0;ratios=[]
    for name,i in circuit['instancias'].items():
        m=i['m'];errs[m]={};baseerrs[m]={}
        for key,c in i['celdas'].items():
            B=F(c['B'])
            for k,row in c['por_k'].items():
                ul,uh,el,eh,ol,oh=[F(v,2**row['precision_bits']) for v in row['valores_punto_fijo']]
                check(0<=ul<=uh and 0<=el<=eh and 0<=ol<=oh,'ordered intervals')
                lo=B*ol; hi=min(F(1),B*oh+1-B*ul)
                if 'intervalo_P_sil_exacto' in row:
                    check([lo,hi]==list(map(F,row['intervalo_P_sil_exacto'])),'silent interval')
                else:
                    check(all(close_decimal(v,e) for v,e in zip(row['intervalo_P_sil'],[lo,hi])),'silent interval decimals')
                check(0<=lo<=hi<=1,'probability bounds')
                check(row['H_comp']=='COMPLETO','completeness provenance')
                e=(hi-lo)/(hi+lo) if hi+lo else F(0)
                if row['e_comp'] is None:
                    check(lo==0,'undefined relative point without lower bound')
                else:
                    check(close_decimal(row['e_comp'],e),'point error')
                if k=='4':
                    check(lo>0,'nonzero lower bound')
                    errs[m][key]=e; baseerrs[m][key]=F(row['e_ing'])
                    if e<=F(1,10) and F(row['e_ing'])>F(1,10): adv+=1
                    if F(c['r_c'])<=F(1,100000) and e<=F(1,10):pass_low+=1
                    ratios.extend(map(F,row['R']))
        cost=i['coste_igual_corte']['certificado_por_k']['4']
        costs.append({'m':m,**cost,'catalogue_seconds':i['igual_tiempo']['segundos_catalogo_T2_1M']})
    check(adv==9,'nine separating cells');check(pass_low==16,'sixteen low-rate cells')
    rows=[]
    for power in [3,4,5,6]:
        for noise in [True,False]:
            key=f'1/{10**power}|'+(f'1/{10**power}' if noise else '0')
            vals=[up(100*errs[m][key]) for m in [6,12,24,48]]
            vals.append(up(100*baseerrs[48][key]))
            rows.append('$10^{-'+str(power)+'}$ & '+('$r_c$' if noise else '$0$')+' & '+' & '.join(vals)+r'\\')
    phrows=[]
    for name,i in ph['instancias'].items():
        if not name.startswith('TIRA_'):continue
        for key,c in i['celdas'].items():
            for k,row in c['por_k'].items():
                L=F(row['L_comp']['exacto']); eps=F(row['eps_comp']['exacto'])
                hi=min(F(1),L+eps); e=(hi-L)/(hi+L) if hi+L else F(0)
                check(L==0 if row['e_comp'] is None else e==F(row['e_comp']['exacto']),'phenomenological point error')
                check(0<=L<=hi<=1,'phenomenological bounds')
                check(F(row['R']['exacto'])*eps==F(row['eps_ing']['exacto']),'tail ratio')
    for power in [3,4,5,6]:
        key=f'1/{10**power}|1/{10**power}'
        vals=[up(100*F(ph['instancias'][f'TIRA_{m}x6']['celdas'][key]['por_k']['4']['e_comp']['exacto'])) for m in [6,12,24,48]]
        vals.append(up(100*F(ph['instancias']['TIRA_48x6']['celdas'][key]['por_k']['4']['e_ing']['exacto'])))
        phrows.append('$10^{-'+str(power)+'}$ & '+' & '.join(vals)+r'\\')
    dmc=load('direct_sampling.json');rec=dmc['recuento']
    check(rec['x']==477 and rec['N']==7000000000,'sampling counts')
    check(sum(rec['por_fallos'].values())==rec['N'],'sample total')
    check(sum(rec['silenciosas_por_fallos'].values())==rec['x'],'silent total')
    cp=dmc['clopper_pearson']
    cell=circuit['instancias']['TIRA_6x6']['celdas']['1/1000|1/1000']['por_k']['4']
    lo,hi=map(F,cell['intervalo_P_sil_exacto'])
    check(max(lo,F(cp['p_inf']))<=min(hi,F(cp['p_sup'])),'sampling compatibility')
    out.mkdir(parents=True,exist_ok=False)
    (out/'circuit_rows.tex').write_text('\n'.join(rows)+'\n',encoding='ascii')
    (out/'phenomenological_rows.tex').write_text('\n'.join(phrows)+'\n',encoding='ascii')
    (out/'resources.json').write_text(json.dumps(costs,indent=2)+'\n',encoding='ascii')
    summary={'status':'PASS','checks':checks,'circuit_rows':128,'separating_cells':adv,'low_rate_cells_passing':pass_low,
             'new_decoder_calls':0,'new_DP_calls':0,'new_Monte_Carlo_histories':0,
             'scope':'Stored endpoints and scalar identities; does not reconstruct the catalogue or DP.'}
    (out/'POSTPROCESSING.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='ascii')
    print(json.dumps(summary))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True,type=Path)
    main(parser.parse_args().output)
