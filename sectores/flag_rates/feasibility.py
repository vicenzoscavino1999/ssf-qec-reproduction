"""Finite-code feasibility calculation. Arguments: reproduction checkout, directory containing native tallies."""
import json,sys
from pathlib import Path
from fractions import Fraction as Q
from math import comb
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from functools import lru_cache
repo,work=map(Path,sys.argv[1:3])
sys.path.insert(0,str(repo/'sectores'))
from mayorante import false_perfect,no_perfect

def fmt(x,upper=False):
    with localcontext() as ctx:
        ctx.prec=13
        ctx.rounding=ROUND_CEILING if upper else ROUND_FLOOR
        return str(Decimal(x.numerator)/Decimal(x.denominator))

@lru_cache(None)
def reference(n):
    d=json.loads((repo/'sectores/referencia'/f'resultado_{n}.json').read_bytes())
    geo,ref=d['geometry'],d['interval']
    m=geo['checks'];marks=geo['marked_pairs']
    counts={tuple(map(int,k.split(','))):v for k,v in ref['positive_census'].items()}
    cap=dict(counts)
    def sizes(a,b):return comb(n,a)*comb(m,b)
    cap[3,1]=min(sizes(3,1),ref['first_step_exceptional_triples']+(n-2)*counts[2,1])
    for a,bs in ((3,(2,3)),(4,(0,1))):
        for b in bs:
            danger=false_perfect(n,m,a,b)+no_perfect(n,m,a,b)*sizes(a,b)
            smaller=comb(n-2,a-2)*cap.get((2,b),marks*comb(m,b))
            if a==4: smaller+=(n-3)*cap[3,b]
            cap[a,b]=min(sizes(a,b),danger+smaller)
    return m,marks,counts,cap,ref

def weights(n,p,q):
    m,marks,counts,cap,ref=reference(n)
    def weight(a,b):return p**a*(1-p)**(n-a)*q**b*(1-q)**(m-b)
    return weight

def evaluate(n,p,q,c40=None,affine=False,refined_caps=None):
    m,marks,known,caps,ref=reference(n)
    c,cap=known.copy(),caps.copy()
    if refined_caps is not None:
        cap.update(refined_caps)
    u=weights(n,p,q)
    if c40 is not None:
        assert n==225 and 0<=c40<=cap[4,0]
        c[4,0]=cap[4,0]=c40
    A,B=(12,9) if n==225 else (24,18)
    pair3=marks*p**2*(1-p)**(n-2)*q**3*(B*comb(m-1,2)+3*A*comb(12,2)*(m-2))
    if affine:
        assert n==225
        # Necessary condition in the independent 225-qubit pilot: one specified bad readout
        # and eight correct ones in the VV row, with the other 99 arbitrary.
        pair_all=known[2,1]*p**2*(1-p)**(n-2)*q*(1-q)**8
        pair_exact=known[2,1]*u(2,1)+known[2,2]*u(2,2)
        assert pair_all>=pair_exact
        pair3=min(pair3,pair_all-pair_exact)
        cap[3,1]=min(cap[3,1],101088)
    tails={'2,>=3':pair3,
           '3,>=4':comb(n,3)*p**3*comb(m,4)*q**4,
           '4,>=2':comb(n,4)*p**4*comb(m,2)*q**2,
           '>=5,*':comb(n,5)*p**5}
    lower=sum(k*u(a,b) for (a,b),k in c.items())
    gaps={f'{a},{b}':(v-c.get((a,b),0))*u(a,b) for (a,b),v in cap.items() if v!=c.get((a,b),0)}
    gaps.update(tails)
    upper=lower+sum(gaps.values())
    err=(upper-lower)/(upper+lower)
    return lower, upper, err, gaps

def binom_tail(n,p,k):
    return 1-sum(comb(n,a)*p**a*(1-p)**(n-a) for a in range(k))
def coordinates(n,p,q,c40=None,upper40=None,upper50=None):
    m,_,known,cap,ref=reference(n)
    u=weights(n,p,q)
    t=1 if n==225 else 2
    rows=[tuple(map(int,line.split())) for line in (work/f'flag_{t}.txt').read_text().splitlines()]
    if n==225:
        assert rows==list(map(tuple,json.loads((work/'audit_225.json').read_bytes())))
    assert len(rows)==6
    for e,b,total,flag,si,silent in rows:
        assert total==comb(n,e)*comb(m,b) and flag<=si<=total and silent==0
    tail=1-sum(count*u(e,b) for e,b,count,flag,si,silent in rows)
    abort_low=sum(flag*u(e,b) for e,b,_,flag,_,_ in rows)
    abort_hi=abort_low+tail
    reject0=sum(si*u(e,b) for e,b,_,flag,si,_ in rows)
    static=binom_tail(n,p,6)*(1-(1-q)**m)
    quad=cap[4,0] if c40 is None else c40
    if upper40 is not None:
        assert c40 is None and upper40<=quad
        quad=upper40
    omitted=(1-q)**m*binom_tail(n,p,5)
    if upper50 is not None:
        assert 0<=upper50<=comb(n,5)
        omitted=upper50*u(5,0)+(1-q)**m*binom_tail(n,p,6)
    hidden_hi=known[3,0]*u(3,0)+quad*u(4,0)+omitted+static
    hidden_lo=known[3,0]*u(3,0)+(c40*u(4,0) if c40 is not None else 0)
    # The gain uses only exhibited silent inputs. The optional affine upper bound
    # is included in silent_all as a comparison and does not enter either gain.
    lo=sum(count*u(a,b) for (a,b),count in known.items())
    if c40 is not None:
        lo+=c40*u(4,0)
    interval_lo,hi,_,_=evaluate(n,p,q,c40,affine=(n==225))
    assert lo==interval_lo
    return {'aborted':[fmt(abort_low),fmt(abort_hi,1)],'input_reject':[fmt(reject0),fmt(reject0+tail,1)],
            'silent_all':[fmt(lo),fmt(hi,1)],'silent_unflagged':[fmt(hidden_lo),fmt(hidden_hi,1)],
            'silent_per_accepted_upper':fmt(hidden_hi/(1-abort_hi),1),
            'gain_per_attempt_lower':fmt(lo/hidden_hi),
            'gain_per_kept_lower':fmt(lo*(1-abort_hi)/hidden_hi),
            'discard_ratio_lower':fmt(reject0/abort_hi),
            'discard_ratio_upper':fmt((reject0+tail)/abort_low,1),
            'static_nonzero_readout_upper':fmt(static,1),
            'postselection_silent_joint_upper':fmt(binom_tail(n,p,6)*(1-q)**m,1),
            'missed_e5plus':fmt(omitted,1),'c40_used':str(quad)}

def rectangle225(c40):
    """Bound the harmonic error over the entire rate rectangle via positive-coefficient odds polynomials."""
    n=225
    m,_,known,caps,_=reference(n)
    assert known[2,2]==(m-9)*known[2,1]
    p0,q0=Q('1e-5'),Q('3e-4')
    x,y=p0/(1-p0),q0/(1-q0)
    # With the new exact C40, uncertainty splits into terms divisible by T21 or by T30.
    # The common factor (1-p)^n (1-q)^m has been removed. Each normalized
    # polynomial increases with both odds, so the upper corner suffices also on the axes.
    assert c40 is not None and 0<=c40<=caps[4,0]
    caps=caps.copy()
    caps[3,1]=min(caps[3,1],101088)
    c21,c30=known[2,1],known[3,0]
    parts={'two_data_extra_readouts':(1+y)**(m-9)-1-(m-9)*y}
    for a,b in [(3,1),(3,2),(3,3),(4,1)]:
        parts[f'exact_weight_{a}_{b}']=Q(caps[a,b]-known.get((a,b),0),c21)*x**(a-2)*y**(b-1)
    for a,b in [(3,4),(4,2)]:
        parts[f'tail_{a}_{b}']=Q(comb(n,a)*comb(m,b),c21)*x**(a-2)*y**(b-1)*(1+x)**(n-a)*(1+y)**(m-b)
    k21=sum(parts.values())
    k30=Q(comb(n,5),c30)*x**2*(1+x)**(n-5)*(1+y)**m
    k=max(k21,k30)
    bound=k/(2+k)
    assert bound<Q('0.10')
    for ratio_p,ratio_q in [(1,1),(100,1),(1,10**3),(5*10**2,10**2),(1,None)]:
        p=p0/ratio_p
        q=q0/ratio_q if ratio_q else Q(0)
        L,U,err,gaps=evaluate(n,p,q,c40,True)
        w=weights(n,p,q)
        assert U-L<=k21*c21*w(2,1)+k30*c30*w(3,0) and err<=bound
    return {'p_max':str(p0),'q_max':str(q0),'K21_upper':fmt(k21,True),'K30_upper':fmt(k30,True),
            'harmonic_error_percent_upper':fmt(100*bound,True),'exact_zero_if_p_zero':True}

def remaining_high_rate(c40):
    lo,hi,_,gap=evaluate(225,Q('1e-4'),Q('1e-3'),c40,True)
    # Both census values can be arbitrarily small or large inside their upper counts.
    # After closing them the true new interval width is the sum of the other gaps
    # and its lower endpoint can only increase.
    remaining=hi-lo-gap['3,2']-gap['4,1']
    bound=remaining/(2*lo+remaining)
    return {'p':'1e-4','q':'1e-3','sectors_to_resolve_exactly':['3,2','4,1'],
            'upper_remaining_gap_over_current_lower':fmt(remaining/lo,True),
            'harmonic_error_percent_upper_after_both_exact':fmt(100*bound,True)}

def flag_region(n,c40=None,upper40=None,upper50=None):
    """At the paper's odds cone, bound unflagged silent probability divided by T30."""
    p0,q0=Q('1e-6'),Q('1e-4')
    x0,y0=p0/(1-p0),q0/(1-q0)
    m,_,known,caps,_=reference(n)
    c30,c21=known[3,0],known[2,1]
    quad=caps[4,0] if c40 is None else c40
    if upper40 is not None:
        assert c40 is None and upper40<=quad
        quad=upper40
    # Zero readout failure can arise from 3 data, 4 data, or at least 5 data.
    # Positive readout and no flip requires a logical E, hence at least 6 data.
    k0=1+Q(quad,c30)*x0+Q(comb(n,5),c30)*x0**2*(1+x0)**(n-5)
    if upper50 is not None:
        k0=1+Q(quad,c30)*x0+Q(upper50,c30)*x0**2+Q(comb(n,6),c30)*x0**3*(1+x0)**(n-6)
    kstatic=Q(comb(n,6),c30)*x0**3*(1+x0)**(n-6)*((1+y0)**m-1)
    # Drop clean-site factors in the low-weight abort terms. The >=3 tail is monotone.
    t=1 if n==225 else 2
    rows=[tuple(map(int,line.split())) for line in (work/f'flag_{t}.txt').read_text().splitlines()]
    w=weights(n,p0,q0)
    tail=1-sum(count*w(e,b) for e,b,count,flag,si,silent in rows)
    abort_upper=tail+sum(flag*p0**e*q0**b for e,b,_,flag,_,_ in rows)
    assert abort_upper<1
    # For 0 < x/y <= x0/y0 and y <= y0 the unfiltered lower bound is T30 + T21.
    # The numerator below is the weakest value of 1+(C21/C30)*y/x on this region.
    gain=(1-abort_upper)*(1+Q(c21,c30)*y0/x0)/(k0+kstatic)
    return {'C21_over_C30':fmt(Q(c21,c30)),'odds_ratio_max':str(x0/y0),'q_max':str(q0),
            'abort_upper_on_cone':fmt(abort_upper,True),'unflagged_over_T30_upper':fmt(k0+kstatic,True),
            'gain_per_kept_lower_on_cone':fmt(gain)}



def orbit_upper_counts(n):
    m,_,known,caps,ref=reference(n)
    a_b = [(3,1),(3,2),(4,0),(4,1)] if n==225 else [(3,1),(4,0)]
    caps=caps.copy()
    t=1 if n==225 else 2
    record={}
    for a,b in a_b:
        i,j,total,exception,reps,budget=map(int,(work/f'orbit_{t}{a}{b}.txt').read_text().split())
        assert (a,b)==(i,j) and total==comb(n,a)*comb(m,b) and 0<=exception<=total
        if (a,b)==(3,1):
            assert exception==ref['first_step_exceptional_triples']
        # At d=6 zero or one residual data fault cannot fail for any readout mask.
        # After a genuine first move leaving r>=2 faults, the remaining execution is the
        # prescribed physical decoding of that residual with the same D.
        charge=sum(comb(n-r,a-r)*caps[r,b] for r in range(2,a))
        bound=min(caps[a,b],exception+charge)
        if n==225 and (a,b)==(3,1):
            bound=min(bound,101088)  # archived independent stratified proof
        caps[a,b]=bound
        record[f'{a},{b}']={'physical_configurations':total,'orbits_checked':reps,'exceptional_upper':exception,
                              'lower_sector_charge':charge,'count_upper':bound}
    return caps,record


def spurious_five_unweighted(n,m):
    """Union over all generators and distinct nonempty FGL footprints at zero readout.
    Each external qubit meets at most one cell in the twelve-cell rectangle of these seeds.
    Count additions to the seven internal data qubits by parity and degree up to five."""
    def footprint(sub):
        return sum(1<<(3*i+j) for i in range(4) for j in range(3)
                   if ((sub>>i)^(sub>>(4+j))) & 1)
    shapes={}
    for sub in range(1,128):
        h=footprint(sub)
        if h and 2*h.bit_count()>=3*sub.bit_count():
            shapes[h]=sub
    @lru_cache(None)
    def completions(d,odd,left):
        coeff=[1]+[0]*left
        for cell in range(d):
            choices=range(1 if cell<odd else 0,6,2)
            coeff=[sum(comb(5,i)*coeff[k-i] for i in choices if i<=k)
                   for k in range(left+1)]
        return sum(v*comb(n-7-5*d,left-k) for k,v in enumerate(coeff))
    count=0
    for h,representative in shapes.items():
        for pattern in range(128):
            size=pattern.bit_count()
            if size>5 or pattern & representative == representative or pattern & (127^representative)==(127^representative):
                continue
            missing=(h&~footprint(pattern)).bit_count()
            count+=completions(h.bit_count(),missing,5-size)
    return m*count

def tighten_five900(caps):
    n=900
    m,_,known,_,_=reference(n)
    spurious=spurious_five_unweighted(n,m)
    old_spurious=false_perfect(n,m,5,0)
    assert old_spurious<=spurious<=4*old_spurious
    # If none of the five originals is fully lit, a distinguished original overlaps
    # at least one of the other four. A qubit has at most 24 possible neighbours.
    blocked=min(Q(1),Q(4*24,n-1))*comb(n,5)
    extended=comb(n-4,1)*caps[4,0]+comb(n-3,2)*known[3,0]
    count=min(comb(n,5),spurious+blocked+extended)
    assert count.denominator==1
    return int(count),{'spurious_perfect_upper':spurious,'no_perfect_original_upper':int(blocked),
                        'smaller_residual_charge':extended,'count_upper':int(count),'all_five':comb(n,5)}


def check_geometry(t):
    from modelo import seed, bits, mask
    columns=seed(t)
    nline,npoint=12*t,9*t
    lines_at_point=[[j for j,col in enumerate(columns) if (col>>point)&1] for point in range(npoint)]
    data_syndromes=[col<<(npoint*left) for left in range(nline) for col in columns]
    data_syndromes.extend(mask(left*npoint+p for left in touching) for touching in lines_at_point for p in range(npoint))
    assert max(sum(bool(s&u) for u in data_syndromes)-1 for s in data_syndromes)<=24
    for r,touching in enumerate(lines_at_point):
        for j,col in enumerate(columns):
            cells=[i*npoint+p for i in touching for p in bits(col)]
            support=set(i*nline+j for i in touching) | set(nline*nline+r*npoint+p for p in bits(col))
            rectangle=mask(cells)
            assert len(cells)==12 and len(support)==7
            counts={cell:0 for cell in cells}
            for site,syn in enumerate(data_syndromes):
                if site not in support:
                    inside=syn&rectangle
                    assert inside.bit_count()<=1
                    if inside: counts[inside.bit_length()-1]+=1
            assert set(counts.values())=={5}

def main():
    p0,q0=Q('1e-6'),Q('1e-4')
    archived_affine=json.loads((repo/'comparaciones/piloto_HGP_225_esquinas_20261124.json').read_bytes())
    assert archived_affine['trios_Bernoulli_una_lectura_max']==101088
    out={}
    for n in (225,900):
        m,_,known,cap,ref=reference(n)
        L,U,err,_=evaluate(n,p0,q0)
        assert Q(ref['lower'])<=L<=U<=Q(ref['upper'])
        assert abs(100*err-Q(ref['max_relative_percent']))<Q('1e-11')
        out[f'flag{n}_archived']=coordinates(n,p0,q0)
    f=work/'four_1.txt';newcount=None
    if f.is_file() and f.stat().st_size:
        total,newcount=map(int,f.read_text().split())
        assert total==comb(225,4)
        out['C40_direct225']={'total':total,'silent':newcount}
        out['flag225_direct']=coordinates(225,p0,q0,newcount)
    scan={}
    for pp,qq in [('1e-6','1e-4'),('1e-5','0'),('1e-5','1e-5'),('1e-5','1e-4'),('1e-5','3e-4'),('1e-4','1e-3')]:
        p,q=Q(pp),Q(qq);res={}
        for name,c40,aff in [('archive',None,False),('combined_affine',None,True)]+([('combined_with_C40',newcount,True)] if newcount is not None else []):
            L,U,e,g=evaluate(225,p,q,c40,aff)
            res[name]={'error_percent':fmt(100*e,1),'interval':[fmt(L),fmt(U,1)]}
            if pp=='1e-4':
                res[name]['excess_over_lower']={i:fmt(x/L,1) for i,x in sorted(g.items(),key=lambda t:-t[1])}
        scan[pp+','+qq]=res
    out['scan225']=scan
    if newcount is not None:
        out['finite_region_225']=rectangle225(newcount)
        out['conditional_plan_for_high_rates']=remaining_high_rate(newcount)
        out['flag_regional']={'225':flag_region(225,newcount),'900':flag_region(900)}
    if all((work/f'orbit_{key}.txt').exists() for key in ['131','132','140','141','231','240']):
        hr={}
        for n in (225,900):
            check_geometry(1 if n==225 else 2)
            caps,record=orbit_upper_counts(n)
            hr[str(n)]=record
            if n==225:
                assert newcount<=caps[4,0]
                p,q=Q('1e-4'),Q('1e-3')
                L,U,e,gaps=evaluate(n,p,q,newcount,True,caps)
                hr['225_at_high_rates']={'p':str(p),'q':str(q),'interval':[fmt(L),fmt(U,True)],
                        'estimate':fmt(2*L*U/(L+U)),'relative_error_percent_upper':fmt(100*e,True),
                        'excess_over_lower':{k:fmt(x/L,True) for k,x in gaps.items()}}
            else:
                upper50,record50=tighten_five900(caps)
                hr['900']['5,0']=record50
                hr['900_at_reference']=coordinates(n,p0,q0,upper40=caps[4,0],upper50=upper50)
                hr['900_gain_on_cone']=flag_region(n,upper40=caps[4,0],upper50=upper50)
        out['orbit_bounds']=hr
    (work/'pilot-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    for n in (225,900):
        print('A',n,out[f'flag{n}_archived'])
    if newcount is not None: print('A exact C40',newcount,out['flag225_direct'])
    for key,recs in scan.items(): print('B',key,[(x,v['error_percent']) for x,v in recs.items()])
    print('B high-rate gap',scan['1e-4,1e-3'])
if __name__=='__main__': main()
