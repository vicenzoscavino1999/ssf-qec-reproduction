"""Verify rational outer bounds for the specified finite-sector cost recipe."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext,ROUND_FLOOR,ROUND_CEILING
from math import comb
import json,sys
if hasattr(sys,'set_int_max_str_digits'):
    sys.set_int_max_str_digits(100000)
ROOT=Path(__file__).resolve().parent
def f(o): return F(int(o['numerator']),int(o['denominator']))
def bounds(a):
    result=[]
    for mode in (ROUND_FLOOR,ROUND_CEILING):
        with localcontext() as c:
            c.prec=90;c.rounding=mode
            result.append(F(str(D(a.numerator)/D(a.denominator))))
    assert result[0]<=a<=result[1]
    return result
def display(a):
    with localcontext() as c:
        c.prec=20
        return format(D(a.numerator)/D(a.denominator),'.12E')
def polynomial(x,y):
    return 1+2*x+21*x*x+100*x**3+66*x**4+220*x**5+y*(3*x+33*x*x+78*x**3+286*x**4)

def main():
    data=json.loads((ROOT/'rational_inputs.json').read_text(encoding='utf-8'))
    rc,rm=f(data['rc']),f(data['rm'])
    a0=F(2,3)*rc;b0=F(8,15)*rc;a=a0/(1-a0);b=b0/(1-b0)
    theta=(1-(1-2*a0)**19*((15-16*rc)/(15-8*rc))**35*(1-2*rm))/2
    tl,tu=bounds(theta)
    rl=bounds((17*a*tl+b/2)/(17*a+b))[0];ru=bounds((17*a*tu+b/2)/(17*a+b))[1]
    sl=bounds((2*a*tl+b/2)/(2*a+b))[0];su=bounds((2*a*tu+b/2)/(2*a+b))[1]
    ql=(1-ru)**14*(1-su);qu=(1-rl)**14*(1-sl)
    pl=polynomial(rl/(1-rl),sl/(1-sl));pu=polynomial(ru/(1-ru),su/(1-su))
    # Bernoulli lower bounds avoid forming enormous powers exactly.
    omega=6471288*a0+3095190*b0;z0=(17*a+b)**14*(2*a+b);m=88434
    assert 0<omega<1 and (m-18)*tu<1
    fl=306*(1-omega)*z0*tl**3*(1-(m-18)*tu);fu=306*z0*tu**3
    lower=fl*ql*pl;upper=fu*qu*pu
    lamv=F(506,15)*rc;lamc=F(504,15)*rc;check=22*rc+rm*(1-4*rc/3)
    def et(t):return sum(F(comb(93636,i)*comb(83521,t-i))*lamv**i*lamc**(t-i) for i in range(t+1))
    zbar=sum(et(t)*comb(m,d)*check**d for t,d in ((5,0),(4,5),(3,6)))
    arg=2000*zbar/lower;budget=(arg.numerator+arg.denominator-1)//arg.denominator
    assert lower==f(data['L_lower']) and upper==f(data['L_upper'])
    assert zbar==f(data['Zbar']) and budget==int(data['sufficient_attempts'])
    assert lower>F('7.98339561055e-1258')
    assert F('5.38398e944')<budget<F('5.38400e944')
    print(json.dumps({'status':'PASS','L_lower':display(lower),'L_upper':display(upper),
                      'Zbar':display(zbar),'sufficient_attempts':display(F(budget)),
                      'relative_error':'1/10','confidence':'19/20','global_probability_computed':False}))
if __name__=='__main__': main()
