"""P1: exact finite reference, not a polymer estimator. Standard library only."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
from array import array
from itertools import combinations
from hashlib import sha256
import json,math,time,sys,traceback

ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parent
DESIGN=PROJECT/'DISENO_PILOTO_ESTIMADOR_20261026'
OUT=Path(sys.argv[1]) if len(sys.argv)>1 else None
BEGIN=time.perf_counter(); CASE_BEGIN=BEGIN; CURRENT='initialization'; PHASE='initialization'
PROFILE={}; COUNTS={}

def write_json(path,obj):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')
    tmp.replace(path)
def guard():
    if time.perf_counter()-BEGIN>119:raise TimeoutError('worker total limit')
    if time.perf_counter()-CASE_BEGIN>59:raise TimeoutError('worker case limit')
def phase(label):
    global PHASE
    guard();PHASE=label
    write_json(OUT/'progress.json',{'case':CURRENT,'phase':label,'case_start':CASE_BEGIN,'worker_start':BEGIN,'updated':time.perf_counter()})
def timed(name,fn,*args):
    phase(name);t=time.perf_counter();r=fn(*args);PROFILE[name]=time.perf_counter()-t;return r
def bits(v):
    while v:
        b=v&-v;yield b.bit_length()-1;v^=b
def encode(s):return sum(1<<i for i in s)
def basis(vectors):
    b={}
    for v in vectors:
        for k in sorted(b,reverse=True):
            if v>>k&1:v^=b[k]
        if v:
            k=v.bit_length()-1
            for j in list(b):
                if b[j]>>k&1:b[j]^=v
            b[k]=v
    return b
def reduce(v,b):
    for k in sorted(b,reverse=True):
        if v>>k&1:v^=b[k]
    return v
def lin(cols,v):
    r=0
    for i in bits(v):r^=cols[i]
    return r
def columns(rows,n):return [sum(1<<j for j,r in enumerate(rows) if r>>i&1) for i in range(n)]
def mismatch(kind,**data):
    write_json(OUT/'FIRST_COUNTEREXAMPLE.json',{'case':CURRENT,'kind':kind,**data})
    raise AssertionError(kind)

def build_fast(spec):
    if spec['id']=='T10':
        cols=[7,11,13,14,7]+[v<<4 for v in (7,11,13,14,7)]
        hx=[sum(1<<i for i,c in enumerate(cols) if c>>j&1) for j in range(8)]
        hz=[17|(17<<5)];n=10;actions=[1<<i for i in range(10)];beta=F(1,10)
    else:
        a=spec['seed_columns'];H=spec['seed_rows_integer_masks'];b=len(H);n=a*a+b*b
        hx=[sum(1<<(i*a+k) for k in bits(H[j]))+sum(1<<(a*a+l*b+j) for l in range(b) if H[l]>>i&1)
            for i in range(a) for j in range(b)]
        hz=[sum(1<<(k*a+j) for k in bits(H[i]))+sum(1<<(a*a+i*b+l) for l in range(b) if H[l]>>j&1)
            for i in range(b) for j in range(a)]
        cols=columns(hx,n);A=sum(r&1 for r in H);actions=set()
        for g in hz:
            f=g
            while f:
                foot=lin(cols,f)
                if foot and (not spec['regular'] or 2*foot.bit_count()>=A*f.bit_count()):actions.add(f)
                f=(f-1)&g
        actions=sorted(actions,key=lambda v:(v.bit_count(),v));beta=F(1,2)
        assert hx==spec['HX'] and hz==spec['HZ'] and actions==spec['candidate_masks_ordered']
    assert n<=25 and all((x&z).bit_count()%2==0 for x in hx for z in hz)
    return {'n':n,'hx':hx,'hz':hz,'cols':cols,'actions':actions,'foot':[lin(cols,f) for f in actions],
            'beta':beta,'checks':len(hx),'spec':spec}

def kron(A,B):
    return [[a*b for a in rowa for b in rowb] for rowa in A for rowb in B]
def eye(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def literal_syndrome(colsets,qubits):
    out=set()
    for q in qubits:out.symmetric_difference_update(colsets[q])
    return frozenset(out)

def build_literal(spec):
    # Independent matrix/candidate construction using Boolean lists and subsets.
    if spec['id']=='T10':
        colsets=[frozenset(bits(v)) for v in (7,11,13,14,7)]
        colsets+= [frozenset(j+4 for j in c) for c in colsets[:5]]
        hxs=[frozenset(i for i,c in enumerate(colsets) if j in c) for j in range(8)]
        hzs=[frozenset((0,4,5,9))];n=10;actions=[frozenset((i,)) for i in range(n)]
    else:
        a=spec['seed_columns'];h=[[int(row>>i&1) for i in range(a)] for row in spec['seed_rows_integer_masks']]
        b=len(h);ht=[list(c) for c in zip(*h)];n=a*a+b*b
        mx=[x+y for x,y in zip(kron(eye(a),h),kron(ht,eye(b)))]
        mz=[x+y for x,y in zip(kron(h,eye(a)),kron(eye(b),ht))]
        hxs=[frozenset(i for i,v in enumerate(row) if v) for row in mx]
        hzs=[frozenset(i for i,v in enumerate(row) if v) for row in mz]
        colsets=[frozenset(j for j,row in enumerate(mx) if row[i]) for i in range(n)]
        A=sum(row[0] for row in h);candidate_set=set()
        for g in hzs:
            for k in range(1,len(g)+1):
                for tup in combinations(sorted(g),k):
                    f=frozenset(tup);foot=literal_syndrome(colsets,f)
                    if foot and (not spec['regular'] or 2*len(foot)>=A*len(f)):candidate_set.add(f)
        actions=sorted(candidate_set,key=lambda f:(len(f),encode(f)))
    # Direct span for independent stabilizer membership, at most 4096 vectors.
    span={frozenset()}
    for g in hzs:span|={frozenset(x.symmetric_difference(g)) for x in list(span)}
    return {'hx':hxs,'hz':hzs,'cols':colsets,'actions':actions,
            'foot':[literal_syndrome(colsets,f) for f in actions],'span':span}

def fast_transitions(f):
    table=[];scores=0;bn=f['beta'].numerator;bd=f['beta'].denominator
    widths=[p.bit_count() for p in f['foot']]
    for s in range(1<<f['checks']):
        if s%128==0:guard()
        best=-1;num=0;den=1
        for j,(foot,d) in enumerate(zip(f['foot'],widths)):
            g=2*(s&foot).bit_count()-d;scores+=1
            if best<0 or g*den>num*d:best=j;num=g;den=d
        if best<0 or num*bd<bn*den:table.append((-1,s))
        else:
            nxt=s^f['foot'][best];assert nxt.bit_count()<s.bit_count()
            table.append((best,nxt))
    COUNTS['fast_score_evaluations']=scores
    return table

def literal_transitions(f,l):
    table=[];scores=0;bn=f['beta'].numerator;bd=f['beta'].denominator
    for s in range(1<<f['checks']):
        if s%64==0:guard()
        current=frozenset(bits(s));winner=None;gain=0;den=1
        for j,foot in enumerate(l['foot']):
            d=len(foot);g=2*len(current.intersection(foot))-d;scores+=1
            if winner is None or g*den>gain*d:winner=j;gain=g;den=d
        if winner is None or gain*bd<bn*den:table.append((-1,s))
        else:
            nxt=current.symmetric_difference(l['foot'][winner]);assert len(nxt)<len(current)
            table.append((winner,encode(nxt)))
    COUNTS['literal_score_evaluations']=scores
    return table

def resolve_and_compare(f,l,ft,lt):
    size=len(ft);cor=[0]*size;end=list(range(size));traces=[() for _ in range(size)]
    for s in range(size):
        if ft[s]!=lt[s]:mismatch('transition',syndrome=s,fast=ft[s],literal=lt[s])
    for s in sorted(range(size),key=int.bit_count):
        j,nxt=ft[s]
        if j>=0:cor[s]=f['actions'][j]^cor[nxt];end[s]=end[nxt];traces[s]=(j,)+traces[nxt]
    literal_cor=[]
    for s in range(size):
        current=s;c=set();trace=[]
        while True:
            j,nxt=lt[current]
            if j<0:break
            c.symmetric_difference_update(l['actions'][j]);trace.append(j);current=nxt
            assert len(trace)<=s.bit_count()
        if tuple(trace)!=traces[s] or encode(c)!=cor[s] or current!=end[s]:
            mismatch('trace_or_correction',syndrome=s,fast_trace=traces[s],literal_trace=trace)
        assert encode(literal_syndrome(l['cols'],c))==s^current
        literal_cor.append(frozenset(c))
    COUNTS['syndromes_compared']=size
    COUNTS['max_selected_actions']=max(map(len,traces))
    return cor,end,traces,literal_cor

def quotient(f):
    n=f['n'];b=basis(f['hz']);free=[i for i in range(n) if i not in b]
    def sig(v):
        r=reduce(v,b)
        return sum(((r>>i)&1)<<j for j,i in enumerate(free))
    qc=[sig(1<<i) for i in range(n)];bc=[f['cols'][i] for i in free]
    qr=columns(qc,len(free)) # transpose: rows of Q from its column masks
    # columns() expects row masks: this call gives d masks with n bits.
    assert len(qr)==len(free)
    assert len(basis(qr))==len(free)
    assert all(lin(qc,z)==0 for z in f['hz'])
    assert all(lin(bc,qc[i])==f['cols'][i] for i in range(n))
    assert 1<<len(free)<=8192
    return b,free,qc,bc,qr

def classify(f,l,q,resolved):
    b,free,qc,bc,qr=q;cor,end,traces,lcor=resolved
    flags=bytearray(1<<len(free));witness_detect=None;witness_silent=None
    for x in range(len(flags)):
        if x%256==0:guard()
        e=sum(1<<free[j] for j in bits(x));s=lin(bc,x)
        rc=x^lin(qc,cor[s]);terminal=lin(bc,rc)
        assert terminal==end[s]
        raw=set(bits(e));raw.symmetric_difference_update(lcor[s]);raw=frozenset(raw)
        terminal_set=literal_syndrome(l['cols'],raw)
        total=rc!=0;zero=terminal==0;silent=total and zero
        literal_total=raw not in l['span'];literal_zero=not terminal_set
        if (total,zero,silent)!=(literal_total,literal_zero,literal_total and literal_zero):
            mismatch('class_event',class_id=x,representative=e,syndrome=s)
        flags[x]=int(zero)|int(total)<<1|int(silent)<<2
        w={'error_mask':e,'error_weight':e.bit_count(),'initial_syndrome':s,
           'terminal_syndrome':terminal,'correction_mask':cor[s],
           'residual_mask':e^cor[s],'selected_action_masks':[f['actions'][j] for j in traces[s]]}
        if not zero and (witness_detect is None or w['error_weight']<witness_detect['error_weight']):witness_detect=w
        if silent and (witness_silent is None or w['error_weight']<witness_silent['error_weight']):witness_silent=w
    COUNTS['quotient_classes_compared']=len(flags)
    return flags,witness_detect,witness_silent

def exact_dp(f,q,flags):
    n=f['n'];d=len(q[1]);states=1<<d;qc=q[2]
    assert 2**n<2**64
    dp=[array('Q',[0])*states];dp[0][0]=1;scanned=0;updates=0
    for i,v in enumerate(qc):
        guard()
        nxt=[array('Q',[0])*states for _ in range(i+2)]
        for k,row in enumerate(dp):
            a=nxt[k];b=nxt[k+1]
            for x,value in enumerate(row):
                scanned+=1
                if value:a[x]+=value;b[x^v]+=value;updates+=2
        dp=nxt
    spectra={name:[0]*(n+1) for name in ('zero_syndrome','silent_failure','total_failure','detectable_failure','success')}
    for k,row in enumerate(dp):
        assert sum(row)==math.comb(n,k)
        for x,value in enumerate(row):
            mask=flags[x]
            if mask&1:spectra['zero_syndrome'][k]+=value
            if mask&2:spectra['total_failure'][k]+=value
            if mask&4:spectra['silent_failure'][k]+=value
            if not mask&1:spectra['detectable_failure'][k]+=value
            if not mask&2:spectra['success'][k]+=value
        assert spectra['zero_syndrome'][k]==spectra['success'][k]+spectra['silent_failure'][k]
        assert spectra['total_failure'][k]==spectra['detectable_failure'][k]+spectra['silent_failure'][k]
    class_size=1<<len(basis(f['hz']))
    assert all(sum(row[x] for row in dp)==class_size for x in range(states))
    digest=sha256()
    for row in dp:digest.update(row.tobytes())
    COUNTS['dp_slots_scanned']=scanned;COUNTS['dp_integer_additions']=updates
    COUNTS['dp_table_cells']=(n+1)*states
    return spectra,digest.hexdigest()

def direct_small(f,l,lcor,spectra):
    if f['n']>10:return {'performed':False,'reason':'n>10; no physical-input exhaustive run'}
    counts={key:[0]*(f['n']+1) for key in spectra}
    for k in range(f['n']+1):
        for tup in combinations(range(f['n']),k):
            e=frozenset(tup);s=encode(literal_syndrome(l['cols'],e))
            r=e.symmetric_difference(lcor[s]);zero=not literal_syndrome(l['cols'],r);total=r not in l['span']
            values={'zero_syndrome':zero,'silent_failure':zero and total,'total_failure':total,
                    'detectable_failure':not zero,'success':not total}
            for key,val in values.items():counts[key][k]+=int(val)
    if counts!=spectra:mismatch('direct_spectrum',direct=counts,dp=spectra)
    COUNTS['direct_physical_inputs']=1<<f['n']
    return {'performed':True,'inputs':1<<f['n'],'all_weight_spectra_equal':True}

def rational_certificate(f):
    if f['spec']['id']=='T10':return {'applicable':False,'reason':'historical toy regression; no MW contract asserted'}
    spec=f['spec'];kappa=F(spec['density_used']);assert kappa.numerator==1
    d=kappa.denominator;D=spec['host_max_degree'];p=F(1,10**spec['candidate_certified_noise_decade_diagnostic'])
    # e <= sum_{j=0}^8 1/j! + (1/9!)/(1-1/10) < 2.72, hence e^(3/2)<9/2.
    eu=sum((F(1,math.factorial(j)) for j in range(9)),F(0))+F(10,9*math.factorial(9))
    assert eu**3<F(9,2)**2
    alpha=F(2,9*(2*D+1));z=p/(1-p)
    # For kappa=1/d, the entropy-envelope d-th power is exact rational.
    lhs=z*F(d**d,(d-1)**(d-1));rhs=alpha**d
    assert p<=kappa and lhs<=rhs
    return {'applicable':True,'p':str(p),'compiler_margin':'1/2','kappa':str(kappa),
      'alpha':str(alpha),'e_upper':str(eu),'entropy_envelope_power_lhs':str(lhs),
      'entropy_envelope_power_rhs':str(rhs),'verified_with_integer_rational_arithmetic':True,
      'scope':'activity contract only; not an executed approximation algorithm'}

def numerical(f,spectra,cert,q):
    if f['spec']['id']=='T10':rates=[F(1,100),F(1,10),F(1,3)];witness=17
    else:rates=[F(1,100),F(1,10),F(cert['p'])];witness=f['spec']['witness_mask']
    assert lin(f['cols'],witness)==0 and lin(q[2],witness)!=0
    n=f['n'];results=[]
    def prob(counts,p):
        a=p.numerator;b=p.denominator
        return F(sum(c*a**k*(b-a)**(n-k) for k,c in enumerate(counts)),b**n)
    def formatted(x):
        with localcontext() as ctx:
            ctx.prec=16;sci=format(Decimal(x.numerator)/Decimal(x.denominator),'.12E')
        return {'exact':str(x),'decimal_display':sci}
    for p in rates:
        guard();v={name:prob(counts,p) for name,counts in spectra.items()}
        z=p/(1-p);ell=z**witness.bit_count()/(1+z**witness.bit_count())
        assert v['zero_syndrome']>0 and v['silent_failure']>=ell*v['zero_syndrome']
        assert v['success']+v['total_failure']==1
        v['conditional_silent_failure']=v['silent_failure']/v['zero_syndrome']
        results.append({'p':str(p),'probabilities':{key:formatted(x) for key,x in v.items()},
                        'conditional_witness_lower_bound':str(ell)})
    if f['spec']['id']=='T10':
        original=json.loads((PROJECT/'GAS_CONDICIONADO_20261026'/'RESULTADOS.json').read_text(encoding='utf-8'))
        for r,old in zip(results,original['toy_models']):
            assert F(r['p'])==F(old['p'])
            for new,key in [('zero_syndrome','P0'),('silent_failure','P_silent'),('conditional_silent_failure','conditional_failure')]:
                assert F(r['probabilities'][new]['exact'])==F(old[key])
    return results

def run_case(spec):
    global CURRENT,CASE_BEGIN,PROFILE,COUNTS
    CURRENT=spec['id'];CASE_BEGIN=time.perf_counter();PROFILE={};COUNTS={}
    f=timed('construct_fast',build_fast,spec)
    l=timed('construct_literal',build_literal,spec)
    assert [encode(x) for x in l['hx']]==f['hx'] and [encode(x) for x in l['hz']]==f['hz']
    assert [encode(x) for x in l['actions']]==f['actions']
    assert [encode(x) for x in l['foot']]==f['foot']
    ft=timed('fast_transitions',fast_transitions,f)
    lt=timed('literal_transitions',literal_transitions,f,l)
    resolved=timed('all_transcripts',resolve_and_compare,f,l,ft,lt)
    q=timed('quotient_linear_algebra',quotient,f)
    flags,detect,silent=timed('all_class_events',classify,f,l,q,resolved)
    spectra,tablehash=timed('exact_weight_dp',exact_dp,f,q,flags)
    direct=timed('direct_small_reference',direct_small,f,l,resolved[3],spectra)
    cert=timed('rational_activity_certificate',rational_certificate,f)
    probs=timed('exact_probability_evaluation',numerical,f,spectra,cert,q)
    if CURRENT in ('T10','MOD25'):
        assert detect is not None and all(F(r['probabilities']['zero_syndrome']['exact'])<1 for r in probs)
    assert silent is not None
    profiles=dict(PROFILE)
    result={'id':CURRENT,'status':'PASS_EXACT_REFERENCE','n':f['n'],'checks':f['checks'],
      'candidate_count':len(f['actions']),'beta':str(f['beta']),'quotient_dimension':len(q[1]),
      'syndrome_count':len(ft),'class_count':len(flags),'counts':dict(COUNTS),'timings_seconds':profiles,
      'case_elapsed_seconds':time.perf_counter()-CASE_BEGIN,'weight_spectra':spectra,
      'Q_S_rows':q[4],'B_columns':q[3],'representative_coordinate_positions':q[1],
      'normalization_by_weight':True,'each_class_cardinality':1<<len(basis(f['hz'])),
      'dp_table_sha256_native_unsigned64':tablehash,'classification_sha256':sha256(flags).hexdigest(),
      'direct_check':direct,'detectable_input_witness':detect,'silent_input_witness':silent,
      'activity_certificate':cert,'probabilities':probs,
      'scope':'exact finite reference and decoder verification, not an estimator or performance advantage'}
    write_json(OUT/(CURRENT+'.json'),result)
    return {'id':CURRENT,'status':result['status'],'elapsed_seconds':result['case_elapsed_seconds'],'counts':result['counts']}

def main():
    assert OUT is not None and OUT.is_dir()
    phase('input_verification')
    tex=Path('D:/paper ia cuantica/qec/BORRADOR_INTRODUCCION_Y_RESULTADOS_EN_v2.tex')
    expected='ef4d0be1ff2e7b66916d6f0b38fe1fa84b69fa9b0cfbc733291859ae078de5b7'
    assert sha256(tex.read_bytes()).hexdigest()==expected
    preflight=json.loads((DESIGN/'PREFLIGHT.json').read_text(encoding='utf-8'))
    specs={s['id']:s for s in preflight['fixtures']}
    jobs=[specs['FGL5'],{'id':'T10'},specs['FGL10'],specs['MOD25']]
    results=[]
    for spec in jobs:results.append(run_case(spec))
    assert sha256(tex.read_bytes()).hexdigest()==expected
    report={'status':'PASS_P1_ALL_FOUR_CASES','worker_elapsed_seconds':time.perf_counter()-BEGIN,
       'cases':results,'source_sha256':expected,'manuscript_unchanged':True,'agents':0,
       'P2_started':False,'P3_started':False,'scope':'exact references only'}
    write_json(OUT/'RESULTS.json',report)
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    try:main()
    except BaseException as ex:
        if OUT is not None:
            write_json(OUT/'FAILURE.json',{'case':CURRENT,'phase':PHASE,'error':repr(ex),
                       'traceback':traceback.format_exc(),'elapsed_seconds':time.perf_counter()-BEGIN})
        raise
