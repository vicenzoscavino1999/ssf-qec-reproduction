from pathlib import Path
import itertools,runpy,sys,json
repo,stage=map(Path,sys.argv[1:3])
aud=runpy.run_path(str(repo/'comparaciones/auditar_esquinas_HGP_225.py'))
_, checks, gens=aud['incidence'](); tags, run, _ = aud['decoder'](checks,gens)
res=[]
for e,b in [(0,0),(0,1),(0,2),(1,0),(1,1),(2,0)]:
    tot=flag=initial=fail=0
    for E in itertools.combinations(range(len(checks)),e):
        he=aud['bit_sum'](checks[v] for v in E); te=aud['bit_sum'](tags[v] for v in E)
        for D in itertools.combinations(range(108),b):
            d=sum(1<<v for v in D);s=he^d; end,cor=run(s)
            tot+=1; flag+=bool(end and end!=s); initial+=bool(s); fail+=bool(end==d and te^cor)
    res.append([e,b,tot,flag,initial,fail])
(stage/'audit_225.json').write_text(json.dumps(res))
print(res)
# Static logical input and one bad readout leaves the positive-score decoder idle.
E=range(6)
syndrome=aud['bit_sum'](checks[i] for i in E)
logical_tag=aud['bit_sum'](tags[i] for i in E)
assert syndrome==0 and logical_tag!=0 and run(1)==(1,0)

# Check the 324 reduced patterns against the independent global decoder with its physical tie order.
row_sectors = runpy.run_path(str(repo/'comparaciones/piloto_HGP_225.py'))['row_sectors']
patterns = []
assert row_sectors(checks,tags,run,patterns=patterns) == 324
expected = set()
for row in range(12):
    for direction in (1,2,3):
        parallel = set(range(3*direction,3*direction+3))
        for left,right in itertools.combinations(sorted(parallel),2):
            third, = parallel - {left,right}
            for point in range(9):
                if checks[third] & (1 << point):
                    expected.add((row,12*row+left,12*row+right,1 << (9*row+point)))
assert set(patterns) == expected
print('Affine row lemma: all 7,920 reduced inputs checked with the physical decoder.')

# The high-rate bound uses this earlier, tighter census of the sector with three
# data faults and one wrong readout, including every maximum-score tie.
H = [sum(1 << j for j in range(12) if checks[j] & (1 << point))
     for point in range(9)]
single_readout_triples = runpy.run_path(
    str(repo/'comparaciones/esquinas_HGP_225.py'))['single_readout_triples']
assert single_readout_triples(H, checks, tags, run) == 101088
print('225-qubit single-readout census: C31 <= 101088.')
