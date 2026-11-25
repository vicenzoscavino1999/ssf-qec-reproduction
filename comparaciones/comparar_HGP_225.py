"""Comparación de los dos certificados deterministas para el piloto afín."""
from argparse import ArgumentParser
from decimal import Decimal, localcontext, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb, factorial
from pathlib import Path
import json
import platform
from runpy import run_path
from statistics import median
import subprocess
import sys
from time import perf_counter


pilot = run_path(str(Path(__file__).with_name("piloto_HGP_225.py")))
bits, mask = pilot["bits"], pilot["mask"]
P = Q(375, 10**9)
Q_MEAS = Q(1, 1000)
TARGET = Q(1, 10)


def affine_maps(lines):
    """Permutaciones de puntos y rectas inducidas por AGL(2,3)."""
    line_index = {line: j for j, line in enumerate(lines)}
    maps = []
    for a,b,c,d in product(range(3), repeat=4):
        if (a*d-b*c) % 3 == 0:
            continue
        for x0,y0 in product(range(3), repeat=2):
            points = tuple(3*((a*x+b*y+x0)%3)+(c*x+d*y+y0)%3
                           for x,y in product(range(3), repeat=2))
            mapped = tuple(line_index[mask(points[t] for t in bits(line))] for line in lines)
            maps.append((mapped,points))
    assert len(maps)==len(set(maps))==432
    return maps


def axis_orbits(nvv, maps):
    """Órbitas de tres coordenadas ordenadas, primero rectas y luego puntos."""
    domains = [range(12)]*nvv + [range(9)]*(3-nvv)
    seen, reps = set(), []
    for labels in product(*domains):
        if labels in seen:
            continue
        orbit = {tuple(m[int(i>=nvv)][v] for i,v in enumerate(labels)) for m in maps}
        assert not (seen & orbit)
        seen.update(orbit)
        reps.append((labels,len(orbit)))
    assert len(seen)==12**nvv*9**(3-nvv)
    return reps


def physical_actions(H, checks, tags):
    """Todas las huellas y clases de corrección, antes de desempatar."""
    classes = {}
    for r,row in enumerate(H):
        for j,line in enumerate(checks[:12]):
            qubits = [12*i+j for i in bits(row)] + [144+9*r+t for t in bits(line)]
            for sub in range(1,128):
                sites = [qubits[i] for i in bits(sub)]
                phi = logical = 0
                for site in sites:
                    phi ^= checks[site]
                    logical ^= tags[site]
                if phi and 2*phi.bit_count() >= 3*len(sites):
                    classes.setdefault(phi,set()).add(logical)
    assert len(classes)==6273
    # Cualquier elección de representante de una misma huella da la misma
    # clase en esta semilla. Se comprueba antes de usar la simetría.
    assert all(len(choices)==1 for choices in classes.values())
    return [(phi,phi.bit_count(),next(iter(choices))) for phi,choices in classes.items()]


def possible_corrections(actions):
    """Clases alcanzables al terminar en cero para algún desempate de máximos."""
    @lru_cache(maxsize=None)
    def terminal_zero(s):
        if s==0:
            return frozenset([0])
        choices = []
        best_lit, best_size = 1,2
        for phi,size,logical in actions:
            lit = (s & phi).bit_count()
            if 100*(2*lit-size)<size:
                continue
            difference = lit*best_size-best_lit*size
            if difference>0:
                best_lit,best_size=lit,size
                choices=[(phi,logical)]
            elif difference==0:
                choices.append((phi,logical))
        return frozenset(logical^after for phi,logical in choices for after in terminal_zero(s^phi))
    return terminal_zero


def safe_measurements(sites, checks, tags, actions):
    """Máscara de bits exteriores con limpieza perfecta para todo desempate."""
    phis = [checks[i] for i in sites]
    if any(a & b for a,b in combinations(phis,2)):
        return 0
    original = phis[0] | phis[1] | phis[2]
    explained = {0:0}
    for i in sites:
        explained.update({s^checks[i]:t^tags[i] for s,t in explained.items()})
    forbidden = original
    for phi,_,logical in actions:
        extra = phi & ~original
        if extra==0:
            if explained.get(phi)!=logical:
                return 0
        elif extra & (extra-1)==0:
            forbidden |= extra
    return ((1 << 108)-1) ^ forbidden


def triple_orbit_bounds(H, checks, tags, contrast=None):
    maps = affine_maps(checks[:12])
    actions = physical_actions(H,checks,tags)
    corrections = possible_corrections(actions)
    upper_zero = safe_one = count = orbit_pairs = 0
    for nvv in range(4):
        reps = axis_orbits(nvv,maps)
        ordered_upper = ordered_safe = ordered_total = 0
        for (left,m1),(right,m2) in product(reps,repeat=2):
            sites = tuple(12*r+c if i<nvv else 144+9*r+c for i,(r,c) in enumerate(zip(left,right)))
            if len(set(sites))!=3:
                continue
            orbit_pairs += 1
            mult = m1*m2
            ordered_total += mult
            syn = checks[sites[0]] ^ checks[sites[1]] ^ checks[sites[2]]
            logical = tags[sites[0]] ^ tags[sites[1]] ^ tags[sites[2]]
            terminal_classes = corrections(syn)
            if any(t!=logical for t in terminal_classes):
                ordered_upper += mult
            safe_bits = safe_measurements(sites,checks,tags,actions)
            ordered_safe += mult*safe_bits.bit_count()
            if contrast is not None:
                run,counts = contrast
                end,correction = run(syn)
                assert end!=0 or correction in terminal_classes
                counts["triples_sin_medida"] += 1
                for j in bits(safe_bits):
                    d = 1 << j
                    assert run(syn ^ d)==(d,logical)
                    counts["triples_y_medida_exterior"] += 1
        permutations = factorial(nvv)*factorial(3-nvv)
        assert all(x%permutations==0 for x in (ordered_total,ordered_upper,ordered_safe))
        count += ordered_total//permutations
        upper_zero += ordered_upper//permutations
        safe_one += ordered_safe//permutations
    assert count==comb(225,3)
    return upper_zero,safe_one,orbit_pairs,corrections.cache_info().currsize


def verify_exterior_bit():
    """Un único bit exterior conserva toda ejecución restringida a una fila."""
    stripes = [7 << (3*i) for i in range(4)]
    stripes += [mask(3*i+j for i in range(4)) for j in range(3)]
    footprints = set()
    for subset in range(1,128):
        phi = pilot["xor_at"](stripes,subset)
        if phi and 2*phi.bit_count() >= 3*subset.bit_count():
            footprints.add(phi)
    assert len(footprints)==63
    for row in stripes[:4]:
        for phi in footprints:
            if phi & ~row:
                n = phi.bit_count()
                for r in range(4):
                    h = min(r,(row & phi).bit_count())+1
                    assert 2*h<=n if r<=1 else 3*h<r*n


def weights(k):
    return ([P**a*(1-P)**(225-a) for a in range(k+1)],
            [Q_MEAS**b*(1-Q_MEAS)**(108-b) for b in range(k+1)])


def retained_noise(k,w,m):
    return sum(comb(225,a)*w[a] * sum(comb(108,b)*m[b] for b in range(k-a+1))
               for a in range(k+1))


def interval(method,pairs,upper_zero,safe_one):
    triples = comb(225,3)
    k = {"sectores":3,"global8":8,"global12":12}[method]
    w,m = weights(k)
    low = pairs*w[2]*(m[1]+99*m[2])
    if method=="sectores":
        pairs_high = pairs*w[2]*Q_MEAS*(1-Q_MEAS)**8
        triples_high = w[3]*(upper_zero*m[0]+triples*(1-m[0])-safe_one*m[1])
        higher = 1-sum(comb(225,a)*w[a] for a in range(4))
        return low,pairs_high+triples_high+higher,{}
    pairs_high = pairs*w[2]*sum(comb(99,b-1)*m[b] for b in range(1,k-1))
    triples_high = w[3]*(upper_zero*m[0]+(108*triples-safe_one)*m[1]
                         +triples*sum(comb(108,b)*m[b] for b in range(2,k-2)))
    higher = sum(comb(225,a)*w[a]*sum(comb(108,b)*m[b] for b in range(k-a+1))
                 for a in range(4,k+1))
    retained_high = pairs_high+triples_high+higher
    tail = 1-retained_noise(k,w,m)
    tail7 = 1-retained_noise(7,w,m)
    # A_7 <= P_sil. El techo por sectores permite descartar incluso el A_7 exacto.
    _,sector_high,_ = interval("sectores",pairs,upper_zero,safe_one)
    assert tail7 > Q(2,9)*sector_high
    assert 0 < tail and tail <= Q(2,9)*pairs*w[2]*m[1]
    assert low<=retained_high<=sector_high<=retained_high+tail
    return low,retained_high+tail,{
        "k":k, "configuraciones_retenidas":sum(comb(333,j) for j in range(k+1)),
        "cola_incondicional_sup":decimal_text(tail,ROUND_CEILING),
        "cola_k7_inf":decimal_text(tail7,ROUND_FLOOR),
        "error_minimo_k7_pct_inf":decimal_text(100*tail7/(2*sector_high+tail7),ROUND_FLOOR),
        "masa_retenida_silenciosa_sup":decimal_text(retained_high,ROUND_CEILING)}


def literal_cut(checks,tags,run):
    """Censo por ubicación física, sin clasificación de filas ni simetrías."""
    w,m = weights(3)
    low = Q(0)
    counts = {}
    total = 0
    for a in range(4):
        for b in range(4-a):
            measurements = [mask(c) for c in combinations(range(108),b)]
            count = hits = 0
            for data in combinations(range(225),a):
                syn = logical = 0
                for i in data:
                    syn ^= checks[i]
                    logical ^= tags[i]
                for d in measurements:
                    end,correction = run(syn ^ d)
                    hits += end==d and correction!=logical
                    count += 1
            assert count==comb(225,a)*comb(108,b)
            assert hits=={(2,1):324,(3,0):36}.get((a,b),0)
            counts[f"{a},{b}"] = {"entradas":count,"silenciosas":hits}
            total += count
            low += hits*w[a]*m[b]
    tail = 1-retained_noise(3,w,m)
    assert total==sum(comb(333,j) for j in range(4))
    return low,low+tail,{"k":3,"configuraciones_evaluadas":total,"sectores":counts,
                         "cola_incondicional_sup":decimal_text(tail,ROUND_CEILING)}


def decimal_text(value,rounding):
    with localcontext() as ctx:
        ctx.prec,ctx.rounding = 14,rounding
        return format(Decimal(value.numerator)/Decimal(value.denominator),".13E")


def worker(method):
    start = perf_counter()
    H,checks,tags,_,_,run = pilot["build"]()
    if method=="literal3":
        prepared = perf_counter()
        low,high,details = literal_cut(checks,tags,run)
    else:
        pairs = pilot["row_sectors"](checks,tags,run)
        verify_exterior_bit()
        upper_zero,safe_one,reps,states = triple_orbit_bounds(H,checks,tags)
        assert (pairs,upper_zero,safe_one,reps)==(324,48,134073792,292)
        prepared = perf_counter()
        low,high,details = interval(method,pairs,upper_zero,safe_one)
        details.update(patrones_de_fila=pairs,fallos_dos_medidas=pairs*99,
                       orbitas_triples_ordenados_por_tipo=reps,estados_empates=states,
                       tres_datos_sin_medida_silenciosos_max=upper_zero,
                       tres_datos_una_medida_sin_fallo_min=safe_one)
    estimate = 2*low*high/(low+high)
    estimate_text = decimal_text(estimate,ROUND_HALF_EVEN)
    point = Q(estimate_text)
    error = max(point/low-1,1-point/high)
    assert (error <= TARGET)==(method!="literal3")
    cert = {"L_inf":decimal_text(low,ROUND_FLOOR),"U_sup":decimal_text(high,ROUND_CEILING),
            "estimacion":estimate_text,"error_relativo_pct_sup":decimal_text(100*error,ROUND_CEILING),**details}
    done = perf_counter()
    return {"metodo":method,"certificado":cert,"preparacion_s":prepared-start,
            "evaluacion_s":done-prepared,"pico_bytes":pilot["peak_memory"]()}


def timed_worker(method):
    cmd = [sys.executable,"-I","-B","-X","utf8",str(Path(__file__).resolve()),"--metodo",method]
    start = perf_counter()
    result = subprocess.run(cmd,capture_output=True,encoding="utf-8")
    elapsed = perf_counter()-start
    if result.returncode:
        raise RuntimeError(result.stderr)
    data = json.loads(result.stdout)
    data["pared_s"] = elapsed
    return data


def benchmark(repetitions,literal):
    methods = ("sectores","global8","global12")
    records = {method:[] for method in methods}
    for i in range(repetitions):
        order = methods[i%3:]+methods[:i%3]
        if i//3%2:
            order = tuple(reversed(order))
        for method in order:
            records[method].append(timed_worker(method))
    if literal:
        records["literal3"] = [timed_worker("literal3")]
    results = {}
    for method,runs in records.items():
        certificate = runs[0]["certificado"]
        assert all(r["certificado"]==certificate for r in runs)
        costs = [{k:r[k] for k in ("pared_s","pico_bytes","preparacion_s","evaluacion_s")} for r in runs]
        results[method] = {"certificado":certificate,"corridas":costs,
                           "mediana_pared_s":median(r["pared_s"] for r in costs),
                           "mediana_pico_MiB":median(r["pico_bytes"] for r in costs)/2**20}
    return {"p":"3.75e-7","q":"1e-3","n_datos":225,"n_medidas":108,
            "H_hex":["249","292","324","4a1","50a","454","911","862","88c"],
            "decodificador":"SSF FGL, beta=1/100, maximo global normalizado, desempate por mayor |F| y menor mascara fisica",
            "objetivo_error_relativo":"0.10","evento":"sindrome fisico cero y clase X no trivial",
            "python":sys.version,"sistema":platform.platform(),"procesador":platform.processor(),
            "costes":"Procesos nuevos secuenciales. Pared exterior, pico del proceso trabajador, incluida la preparacion.",
            "resultados":results,"contraste_geometrico_posterior":check_geometry()}


def check_geometry():
    """Contrastar las reducciones sobre representantes con el SSF de orden físico fijo."""
    H,checks,tags,_,_,run = pilot["build"]()
    patterns = []
    pilot["row_sectors"](checks,tags,run,patterns=patterns)
    verify_exterior_bit()
    counts = {"parejas_y_segunda_medida":0,"triples_sin_medida":0,"triples_y_medida_exterior":0}
    for i,u,v,d0 in patterns:
        for j in range(108):
            if j//9!=i:
                d = d0 | (1 << j)
                end,correction = run(checks[u] ^ checks[v] ^ d)
                assert end==d and correction!=(tags[u] ^ tags[v])
                counts["parejas_y_segunda_medida"] += 1
    assert counts["parejas_y_segunda_medida"]==32076
    upper,safe,reps,_ = triple_orbit_bounds(H,checks,tags,(run,counts))
    assert (upper,safe,reps)==(48,134073792,292)
    return counts


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para activar las comprobaciones.")
    parser = ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--metodo",choices=["sectores","global8","global12","literal3"],
                        help="Evaluar un único certificado y emitirlo en JSON.")
    mode.add_argument("--contraste-geometrico",action="store_true",
                      help="Validar las clases certificadas con el desempate físico.")
    parser.add_argument("--repeticiones",type=int,default=6,help="Procesos por cada método agrupado.")
    parser.add_argument("--sin-literal",action="store_true",help="Omitir el control exhaustivo hasta k=3.")
    parser.add_argument("--salida",type=Path,help="Escribir los resultados en este archivo JSON.")
    args = parser.parse_args()
    if args.repeticiones<1:
        parser.error("Indicar al menos una repetición.")
    if args.contraste_geometrico:
        data = check_geometry()
    elif args.metodo:
        data = worker(args.metodo)
    else:
        data = benchmark(args.repeticiones,not args.sin_literal)
    output = json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+"\n"
    if args.salida:
        args.salida.write_text(output,encoding="utf-8")
    else:
        print(output,end="")
