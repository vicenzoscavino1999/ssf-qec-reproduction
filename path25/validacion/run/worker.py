"""PATH25: one bounded finite reference/local diagnostic. No global estimator.
Reuses pinned P1 Boolean-set and bitmask implementations; local oracle adapted
from P2, removing its clique assumption. No calls into original run entrypoints.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from collections import Counter,defaultdict
import hashlib,importlib.util,json,math,sys,time,traceback
D=Path(__file__).resolve().parent;OUT=Path(sys.argv[1]);START=time.perf_counter()
COUNTS=Counter();PHASE='init';TIMES={};PHASE_START=START
def save(name,obj):
 data=json.dumps(obj,ensure_ascii=False,indent=2)+'\n'
 assert len(data.encode('utf-8'))<16*1024**2
 (OUT/name).write_text(data,encoding='utf-8')
def guard():
 if time.perf_counter()-START>58:raise TimeoutError('internal 58s stop before external 60s cap')
 assert COUNTS['oracle_queries']<=19714
def phase(name):
 global PHASE,PHASE_START
 guard();now=time.perf_counter();TIMES[PHASE]=now-PHASE_START;PHASE=name;PHASE_START=now
 save('progress.json',{'phase':name,'elapsed_seconds':now-START,'counts':dict(COUNTS)})
def require(ok,kind,**data):
 if not ok:
  save('FIRST_COUNTEREXAMPLE.json',{'kind':kind,'phase':PHASE,**data})
  raise AssertionError(kind)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
P1=D/'inputs'/'P1_WORKER.py'
assert sha(P1)=='57ea8e5af6cb01f92be0944586615438e49dc50a65dfaceccfb6d1251601e231'
specmod=importlib.util.spec_from_file_location('p1_frozen',P1)
p1=importlib.util.module_from_spec(specmod);specmod.loader.exec_module(p1)
p1.guard=guard;p1.OUT=OUT;p1.CURRENT='PATH25';p1.COUNTS={}
bits=p1.bits;lin=p1.lin

def connected_parts(Y,adj):
 out=[];remaining=Y
 while remaining:
  seen=0;todo=remaining&-remaining
  while todo:
   bit=todo&-todo;todo^=bit;i=bit.bit_length()-1
   seen|=bit;todo|=adj[i]&remaining&~seen
  remaining&=~seen;out.append(seen)
 return out
def construct_host(f):
 n=f['n'];psi=f['cols']+[1<<j for j in range(f['checks'])];adj=[0]*len(psi)
 inverse=[[] for _ in range(f['checks'])]
 for j,(action,foot) in enumerate(zip(f['actions'],f['foot'])):
  for c in bits(foot):inverse[c].append(j)
  readers=sum(1<<i for i,x in enumerate(psi) if x&foot)
  for i in bits(readers):adj[i]|=readers|action
  for i in bits(action):adj[i]|=readers
 for i in range(len(adj)):adj[i]&=~(1<<i)
 assert all(adj[j]>>i&1 for i in range(len(adj)) for j in bits(adj[i]))
 return adj,inverse

def main():
 phase('construction_and_reference')
 pf=json.loads((D/'inputs'/'PREFLIGHT.json').read_text())
 sp={'id':'PATH25','regular':False,'seed_columns':4,'seed_rows_integer_masks':[3,6,12],
     'HX':pf['HX'],'HZ':pf['HZ']}
 literal=p1.build_literal(sp)
 sp['candidate_masks_ordered']=[p1.encode(x) for x in literal['actions']]
 f=p1.build_fast(sp);adj,inverse=construct_host(f)
 assert adj==pf['host_adjacency'] and len(f['actions'])==103
 assert [p1.encode(x) for x in literal['hx']]==f['hx']
 assert [p1.encode(x) for x in literal['hz']]==f['hz']
 assert [p1.encode(x) for x in literal['foot']]==f['foot']
 ft=p1.fast_transitions(f);lt=p1.literal_transitions(f,literal)
 resolved=p1.resolve_and_compare(f,literal,ft,lt)
 cor,end,traces,lcor=resolved
 touched=[]
 for tr in traces:
  mask=0
  for j in tr:mask|=f['actions'][j]
  touched.append(mask)
 q=p1.quotient(f);qc=q[2]
 flags,detect,silent=p1.classify(f,literal,q,resolved)
 phase('quotient_weight_reference')
 spectra,dphash=p1.exact_dp(f,q,flags)
 assert p1.COUNTS['fast_score_evaluations']==p1.COUNTS['literal_score_evaluations']==421888
 assert p1.COUNTS['dp_integer_additions']<=5324800
 save('REFERENCE.json',{'fixture':'PATH25','counts':dict(p1.COUNTS),
    'spectra':spectra,'quotient_rows':q[4],'representative_positions':q[1],
    'stabilizer_rank':len(q[0]),'class_count':len(flags),
    'classification_sha256':hashlib.sha256(flags).hexdigest(),'dp_native_uint64_sha256':dphash,
    'detectable_witness':detect,'silent_witness':silent,
    'DP_is_one_exact_method_not_two_independent_DP_implementations':True})
 save('REFERENCE_TABLE.json',{'actions':f['actions'],'footprints':f['foot'],
    'transitions':ft,'correction':cor,'terminal':end,'trace':traces,'touched':touched})
 cache={}
 def query(e):
  COUNTS['reference_table_lookups']+=1
  if e not in cache:
   s=lin(f['cols'],e)
   cache[e]={'E':e,'f':e.bit_count(),'s0':s,'Y':e|touched[s],
      'R':e^cor[s],'syndrome':end[s],'class':lin(qc,e^cor[s]),'trace':list(traces[s])}
  return cache[e]
 def local(Y,E,readers):
  COUNTS['oracle_queries']+=1;guard()
  s=lin(f['cols'],E);seen=E;c=0;tr=[]
  while True:
   best=-1;gn=0;gd=1
   for j in readers:
    foot=f['foot'][j];d=foot.bit_count();g=2*(s&foot).bit_count()-d
    COUNTS['oracle_score_evaluations']+=1
    if best<0 or g*gd>gn*d:best=j;gn=g;gd=d
   if best<0 or 2*gn<gd:return ('accept' if seen==Y else 'underfilled',tr,c,s,seen)
   a=f['actions'][best];tr.append(best)
   if a&~Y:return ('exit',tr,c,s,seen)
   old=s.bit_count();c^=a;seen|=a;s^=f['foot'][best]
   assert s.bit_count()<old
 def expected_prefix(Y,r):
  s=r['s0'];seen=r['E'];c=0;tr=[]
  for j in r['trace']:
   tr.append(j);a=f['actions'][j]
   if a&~Y:return ('exit',tr,c,s,seen)
   c^=a;s^=f['foot'][j];seen|=a
  return ('accept' if seen==Y else 'underfilled',tr,c,s,seen)

 phase('bounded_local_oracle')
 supports=[]
 for t in range(1,4):
  for tup in combinations(range(25),t):
   Y=sum(1<<i for i in tup)
   if len(connected_parts(Y,adj))==1:supports.append(Y)
 supports += [15<<(4*i) for i in range(4)]
 assert len(set(supports))==len(supports)
 cat=[];outcomes=Counter();examples={};support_sizes=Counter();partial_by_Y={}
 for Y in supports:
  support_sizes[Y.bit_count()]+=1
  sy=0
  for i in bits(Y):sy|=f['cols'][i]
  readers=sorted({j for c in bits(sy) for j in inverse[c]})
  E=Y;accepted=[]
  while True:
   got=local(Y,E,readers);actual=query(E);expected=expected_prefix(Y,actual)
   require(got==expected,'local_oracle_mismatch',Y=Y,E=E,got=got,expected=expected)
   outcomes[got[0]]+=1
   if got[0] not in examples:examples[got[0]]={'Y':Y,'E':E,'outcome':got,'full_execution':actual}
   if got[0]=='accept':
    assert E and actual['Y']==Y
    cat.append(actual);accepted.append(actual)
   if not E:break
   E=(E-1)&Y
  partial_by_Y[Y]=accepted
 assert len({(r['Y'],r['E']) for r in cat})==len(cat)
 save('CATALOGUE.json',cat)
 save('ORACLE.json',{'queried_supports_by_size':dict(support_sizes),'queries':COUNTS['oracle_queries'],
    'outcomes':dict(outcomes),'examples':examples,'catalogue_decorations':len(cat),
    'complete_only_for':'all connected data supports of size <=3, plus four fixed VV row supports of size 4',
    'whole_instance_catalogue_complete':False})

 phase('restriction_and_canonical_recomposition')
 low=[0]+[sum(1<<i for i in c) for t in (1,2) for c in combinations(range(25),t)]
 assert len(low)==326
 fixed_pairs=pf['compatible_zero_syndrome_stabilizer_pairs']+pf['compatible_zero_syndrome_logical_pairs']
 selected=set(low)|{r['E'] for r in cat}|{x['Y']|x['Z'] for x in fixed_pairs}
 low_rows=[query(e) for e in low];save('LOW_WEIGHT_SOURCES.json',low_rows)
 recomposed=[];multi_adaptive=[]
 for E in sorted(selected,key=lambda e:(e.bit_count(),e)):
  guard();r=query(E);parts=connected_parts(r['Y'],adj);local_rows=[]
  xor=0
  for Y in parts:
   rr=query(E&Y);local_rows.append(rr);xor^=rr['R']
   expected=[j for j in r['trace'] if f['actions'][j]&Y]
   require(rr['Y']==Y and rr['trace']==expected,'restriction',E=E,Y=Y,local=rr,global_trace=r['trace'])
   require(all(not (f['actions'][j]&~Y) for j in expected),'action_crosses_components',E=E,Y=Y)
  require(xor==r['R'],'residual_xor',E=E)
  # Merge only local current heads by exact scores and the original global rank.
  ss=[rr['s0'] for rr in local_rows];merged=[]
  while True:
   heads=[]
   for h,s in enumerate(ss):
    j,_=ft[s]
    if j>=0:
     foot=f['foot'][j];heads.append((F(2*(s&foot).bit_count()-foot.bit_count(),foot.bit_count()),-j,h,j))
   if not heads:break
   _,_,h,j=max(heads);merged.append(j);ss[h]^=f['foot'][j]
   require(len(merged)<=12,'merge_termination',E=E)
  require(merged==r['trace'],'canonical_merge',E=E,merged=merged,reference=r['trace'])
  row={'E':E,'Y':r['Y'],'components':parts,'local_sources':[rr['E'] for rr in local_rows],
       'local_residuals':[rr['R'] for rr in local_rows],'trace':r['trace'],
       'local_traces':[rr['trace'] for rr in local_rows],'initial_syndrome':r['s0']}
  if len(parts)>1:
   recomposed.append(row)
   if r['trace']:multi_adaptive.append(row)
 COUNTS['recomposition_sources']=len(selected)
 save('RECOMPOSITION.json',{'selected_sources':len(selected),'multicomponent_cases':len(recomposed),
    'adaptive_multicomponent_cases':len(multi_adaptive),'multicomponent_records':recomposed,
    'first_adaptive_multicomponent_witness':multi_adaptive[0] if multi_adaptive else None,
    'fixed_pair_count':len(fixed_pairs),'scope':'selected sources only; no exhaustive family enumeration'})
 low_spectra={key:[0,0,0] for key in spectra}
 for r in low_rows:
  k=r['f'];total=r['class']!=0;zero=r['syndrome']==0
  for key,v in {'zero_syndrome':zero,'silent_failure':zero and total,'total_failure':total,
                'detectable_failure':not zero,'success':not total}.items():low_spectra[key][k]+=v
 require(all(v==spectra[k][:3] for k,v in low_spectra.items()),'low_weight_vs_dp',direct=low_spectra)

 phase('activity_and_candidate_diagnostic')
 # This grouping is analysis only; neither reference engine used reduced candidates.
 leaders={};groups=defaultdict(list)
 for j,foot in enumerate(f['foot']):leaders.setdefault(foot,j);groups[foot].append(j)
 require(all(j<0 or leaders[f['foot'][j]]==j for j,_ in ft),'dominated_winner')
 def min_gain(d):
  g=(d+1)//2
  if (g-d)%2:g+=1
  return g
 efficiencies=[F(f['actions'][j].bit_count(),min_gain(foot.bit_count())) for foot,j in leaders.items()]
 eta=max(efficiencies);dx=max(x.bit_count() for x in f['cols']);wz=max(x.bit_count() for x in f['hz'])
 k_generic=F(1,1+wz*max(1,dx));k_leaders=1/(1+eta*max(1,dx))
 z=F(1,99);activity=[]
 for Y,rs in partial_by_Y.items():
  t=Y.bit_count();hist=Counter(r['f'] for r in rs)
  ht=Counter(r['f'] for r in rs if not r['syndrome'])
  hb=Counter(r['f'] for r in rs if not r['syndrome'] and r['class'])
  w=sum((v*z**f0 for f0,v in hist.items()),F(0))
  wt=sum((v*z**f0 for f0,v in ht.items()),F(0));ws=sum((v*z**f0 for f0,v in hb.items()),F(0))
  m=min(hist) if hist else None
  if m is not None:assert F(m,t)>=k_leaders
  bm=sum((math.comb(t,i)*z**i for i in range(m,t+1)),F(0)) if m else None
  m0=math.ceil(k_generic*t);bk=sum((math.comb(t,i)*z**i for i in range(m0,t+1)),F(0))
  activity.append({'Y':Y,'t':t,'N_f':dict(hist),'N_terminal_f':dict(ht),'N_silent_f':dict(hb),
    'minimum_f':m,'W0':str(w),'WT':str(wt),'Ws':str(ws),'p':'1/100',
    'binomial_generic':str(bk),'binomial_actual_minimum':str(bm) if bm is not None else None})
 save('ACTIVITIES.json',activity)
 low_nonempty=[r for r in low_rows if r['E']]
 worst=min(low_nonempty,key=lambda r:F(r['f'],r['Y'].bit_count()))
 diagnostic={'generic_density':str(k_generic),'fixed_order_leader_density':str(k_leaders),
    'd_X':dx,'w_Z':wz,'eta_max_action_over_minimum_integer_decrease':str(eta),
    'candidate_count':103,'distinct_footprints':len(leaders),'never_winning_candidates':103-len(leaders),
    'groups':[{'footprint':foot,'leader':j,'leader_action':f['actions'][j],
               'members':groups[foot],'minimum_integer_decrease':min_gain(foot.bit_count())}
              for foot,j in leaders.items()],
    'no_decoder_modification':True,'low_weight_min_f_over_Y':str(F(worst['f'],worst['Y'].bit_count())),
    'low_weight_worst_witness':worst,'catalogue_min_f_over_Y':str(min(F(r['f'],r['Y'].bit_count()) for r in cat)),
    'catalogue_incomplete':True,'claimed_global_density_basis':'candidate dominance plus syndrome descent; requires separate proof review'}
 save('DIAGNOSTIC.json',diagnostic)
 phase('complete')
 report={'status':'PASS_BOUNDED_PATH25','fixture':'PATH25','elapsed_seconds':time.perf_counter()-START,
    'timings_seconds':TIMES,'reference_counts':dict(p1.COUNTS),'counts':dict(COUNTS),
    'queried_supports':len(supports),'catalogue_decorations':len(cat),
    'adaptive_multicomponent_cases':len(multi_adaptive),'distinct_low_sources':326,
    'new_global_estimator':False,'full_catalogue':False,'agents':0,
    'minimum_total_failure_weight':next(i for i,c in enumerate(spectra['total_failure']) if c),
    'minimum_silent_failure_weight':next(i for i,c in enumerate(spectra['silent_failure']) if c),
    'manuscript_edited':False}
 save('RESULTS.json',report);print(json.dumps(report),flush=True)

if __name__=='__main__':
 try:main()
 except BaseException as ex:
  save('FAILURE.json',{'phase':PHASE,'error':repr(ex),'traceback':traceback.format_exc(),
       'elapsed_seconds':time.perf_counter()-START,'counts':dict(COUNTS),'reference_counts':dict(p1.COUNTS)})
  raise
