"""HGP binario (3,4). Episodios en checks compartidos. Python 3.10+, biblioteca estándar."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import factorial
from pathlib import Path
from runpy import run_path

v11 = run_path(str(Path(__file__).with_name("verificar_H5A_3_4_ramificaciones_dos_ejes.py")))
prop, old = v11["prop"], v11["prior"]
effective, upper = v11["effective"], v11["upper"]
ROWS, COLS = prop["ROWS"], prop["COLS"]
THETA_MAX = Q(".0012")
PAIR_PROB = (Q(".152"), Q(".167"), Q(".528"))

def submasks(mask):
    chosen = mask
    while chosen:
        yield chosen
        chosen = (chosen-1) & mask
    yield 0


def check_local():
    _, stable, reps = v11["previous"]("contacto_de_checks")["geometry"]()
    winners = prop["winners"](reps)
    cuts = {phi for lit,phi in winners if lit==phi}
    assert len(cuts)==63
    classes = Counter()
    for lit,phi in winners:
        births = phi ^ lit
        delta = lit.bit_count()-births.bit_count()
        for lines in (ROWS,COLS):
            export = sum(max(0,(births&s).bit_count()-(lit&s).bit_count()) for s in lines)
            assert export<=delta
        if delta!=1:
            continue
        if phi in ROWS:
            classes["fila"]+=1
        else:
            crosses=[(row,col) for row,col in product(ROWS,COLS) if row^col==phi]
            assert len(crosses)==1
            row,col=crosses[0]
            assert ((lit&row).bit_count(),(lit&col).bit_count(),
                    (births&row).bit_count(),(births&col).bit_count())==(1,2,1,1)
            classes["cruz"]+=1
    assert len(winners)==771 and classes=={"fila":12,"cruz":72}
    # En una unión de una fila y una columna sólo cabe como huella mixta
    # la cruz formada por ambas. Las demás huellas posibles son franjas.
    assert all(phi in (row,col,row^col) for row,col in product(ROWS,COLS)
               for phi in cuts if not phi&~(row|col))
    assert all(phi in (a,b,a|b) for lines in (ROWS,COLS) for a,b in product(lines,repeat=2)
               for phi in cuts if not phi&~(a|b))
    observed=[]
    for i,j in ((0,0),(0,1),(1,1)):
        left=(ROWS,COLS)[i][0]
        right=(ROWS,COLS)[j][int(i==j)]
        peak=Q(0)
        for a,b in product(submasks(left),submasks(right)):
            if i!=j and a and b and ((a&right) or (b&left)):
                continue
            inside=a|b
            events={phi&~inside for phi in cuts if phi&inside and phi&~inside and stable[phi&~inside]}
            minimal=[e for e in events if not any(t!=e and (t&e)==t for t in events)]
            peak=max(peak,sum(THETA_MAX**e.bit_count() for e in minimal))
        prob=(6+i+j)*THETA_MAX+((30,36)[i]+(30,36)[j])*peak
        assert prob<=PAIR_PROB[i+j]
        observed.append(prob)
    print("Ganadores de descenso uno: 12 filas y 72 cruces. Techos locales de incidencias:",
          *(f"{float(s):.9g}" for s in observed),flush=True)


def prices(Y,xi,theta,ball,discounts,hist):
    assert theta<=THETA_MAX and 0<Y<xi<=1
    regular,physical=v11["polynomials"](Y,ball,discounts,hist)
    def clean_pair(i,j,k):
        return Y**2/xi**(k+2)+PAIR_PROB[i+j]*Y**(-k)
    same=((576,54,9),(384,72,20))
    for i in (0,1):
        key=(1-i,i)
        rate=Y**2/xi+(7,12)[i]*THETA_MAX*Y
        regular[i][key]=min(regular[i][key],discounts[key]*rate)
        key=(2-2*i,2*i)
        zero,one,two=same[i]
        rate=zero*clean_pair(i,i,0)+two/Y**2
        rate+=one*(clean_pair(i,i,1) if i else 1/Y)
        if i:
            rate+=48+24/Y
        regular[i][key]=min(regular[i][key],discounts[key]*rate)
        aligned=old["anchor"](ball,i,2,fiber=True)
        physical[i][key]=min(physical[i][key],discounts[key]*two*aligned)
        j=1-i
        strong,overlap,paths=(117,9,324) if i==0 else (208,16,576)
        rate=(strong+24)/Y+overlap+paths*(clean_pair(i,j,0)+clean_pair(i,j,1)/2)
        regular[i][1,1]=min(regular[i][1,1],discounts[1,1]*rate)
    return regular,physical


def low_degree_tables(tables,discounts,Y,xi,theta):
    """Representación positiva de los dos precios para reutilizar las filas exactas de v11."""
    ball,halos,series,geometry,hist,trees=tables
    regular,physical=prices(Y,xi,theta,ball,discounts,hist)
    moments=[row.copy() for row in hist]
    for i in (0,1):
        for key in hist[i]:
            if sum(key)>2:
                continue
            factor=discounts[key]
            two=physical[i][key]/(factor*old["anchor"](ball,i,2))
            zero=regular[i][key]/factor-two/Y**2
            assert min(zero,two)>=0
            moments[i][key]={0:zero,2:two}
    assert v11["polynomials"](Y,ball,discounts,moments)==(regular,physical)
    return ball,halos,series,geometry,moments,trees


def projected(witness,tables,discounts):
    a,xi=Q(witness["a"]),Q(witness["xi"])
    theta=Q(witness["outer"][0])
    inputs=low_degree_tables(tables,discounts,a**6,xi,theta)
    return v11["fit_row"](witness,inputs,discounts)


def rectangle(witness,tables,discounts):
    ball,halos,series,_,hist,_=tables
    a,s=Q(witness["a"]),Q(witness["s"])
    Y=a**6
    theta=Q(witness["outer"][0])
    reg,phys=prices(Y,Q(1),theta,ball,discounts,hist)
    # Grupo único, ancla en generador, desde once datos (v11 §2).
    m=max(map(Q,witness["m"]))
    rate=upper(Q(1))*1297*m*Q(12,11)**9
    assert rate<1
    phys[0][0,0]=2*2000*11**9*Q(11**9,factorial(10))*m*(1297*m)**10/(1-rate)
    d,e=[tuple(map(Q,witness[key])) for key in ("d","e")]
    kernel=prop["kernels"](ball,Y,s)
    effective["certify"](Y,s,halos,reg,phys,kernel,
                         lambda x,i:effective["full_tail"](series,a,x,i,d,e),witness)

RECTANGLE = dict(a=".555", s=".825", name="Rectángulo por episodios",
    rates=("3.4e-7",".00115",".004",".00061","1.01"),
    m=("6.626303568224e-7","7.801234943708e-7"), u=("4.559063931966e-7","4.905216663034e-7"),
    h=("2.988992071134e-8","2.012296976029e-7"), physical=("1.828158752374e-9","1.131887995484e-8"),
    outer=(".001156586821741","3.511203304195e-7","8.210412876932e-6"),
    d=("5.004323975292e-5",".0001235454832456"),e=("5.004323975292e-5",".0003570503180436"))

SPARSE = dict(a=".525", xi="1", b=".00115", scale=".13", group=Q(215,13), refund=0,
    m=("2.482575860034e-7","2.90117098157e-7"), u=("2.85398066865e-6","3.0138253057e-6"),
    h=("5.813057526991e-8","3.495802219194e-7"), physical=("2.415124652791e-9","1.17937578785e-8"),
    outer=(".001149738751494","4.877103090167e-8","1.2663895096e-6"),
    d=("2.466878082344e-5","6.024478013391e-5"),e=("2.466878082344e-5",".0002022380176613"))

AXES = dict(a=".64", xi=".225", b=".001153", scale="1", group=Q(33,20), refund=1,
    m=("2.020336858678e-6","2.399648710067e-6"), u=("2.260806311085e-6","2.406909963778e-6"),
    h=("1.431097469476e-7","3.554932113185e-7"), physical=("1.212752630153e-8","2.681120885208e-8"),
    outer=(".001164556705849","3.355787400666e-7","1.140359719069e-5"),
    d=("5.282325393427e-5",".0001065450308925"),e=("5.282325393427e-5",".0002453572210939"))


def main():
    check_local()
    tables=v11["data_tables"]()
    discounts=v11["geometry_discounts"](tables[3])
    rectangle(RECTANGLE,tables,discounts)
    sparse=projected(SPARSE,tables,discounts)
    angular=projected(AXES,tables,discounts)
    tau,xi=AXES["group"],Q(AXES["xi"])
    # Para m>=3 y 3H<=m+2, los casos H=0, H=1 y H>=2.
    minimum=tau**4*xi**2
    assert tau**3*xi>1 and 0<minimum<min(tau**3,tau**3*xi,1)
    assert Q("2.15")**8*Q(".13")**3>1
    extra=sparse+angular/minimum
    total=Q(".002903264")+extra
    assert extra<Q(".000505949") and total<Q(".003409213")
    assert total+Q(".00059078")<Q(".004")
    print(f"G3, dos clases del resto v10: {float(extra):.12g}; subtotal calculado {float(total):.12g}.")
    print("Subtotal < 0.003409213. Presupuesto suficiente para el complemento: 0.00059078. G3 abierto.")


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para comprobar las desigualdades.")
    main()
