"""Exact-rational oracle tests; no floating point is used to judge correctness."""
from fractions import Fraction as F
from random import Random
from decimal_interval import Box,ZERO,ONE,ldlt,solve_spd,matrix_minor_pivots
from second_order_ad import Jet,N,PAIRS

def self_tests():
    rng=Random(16092026);checks=0
    for _ in range(500):
        aa=sorted([F(rng.randrange(-1000,1001),rng.randrange(1,101)) for _ in range(2)])
        bb=sorted([F(rng.randrange(-1000,1001),rng.randrange(1,101)) for _ in range(2)])
        a=Box(*aa);b=Box(*bb)
        ops=[(a+b,lambda u,v:u+v),(a-b,lambda u,v:u-v),(a*b,lambda u,v:u*v)]
        if not bb[0]<=0<=bb[1]:ops.append((a/b,lambda u,v:u/v))
        for ans,op in ops:
            for u in [aa[0],sum(aa)/2,aa[1]]:
                for v in [bb[0],sum(bb)/2,bb[1]]:
                    if not ans.contains(op(u,v)):raise AssertionError('interval containment failed')
                    checks+=1
    oracle_checks=checks
    jetchecks=0
    for x,y,z in [(F(2,3),F(-4,5),F(7,11)),(F(9,7),F(3,8),F(-1,3)),(F(3,2),F(-7,9),F(2,5))]:
        xx=Jet.variable(Box(x),0);yy=Jet.variable(Box(y),1);zz=Jet.variable(Box(z),2)
        r=(xx*yy+1)/(xx+2)
        grad=[F(0)]*N;h={}
        grad[0]=(2*y-1)/(x+2)**2;grad[1]=x/(x+2)
        h[0,0]=-2*(2*y-1)/(x+2)**3;h[0,1]=F(2)/(x+2)**2
        if not r.value.contains((x*y+1)/(x+2)):raise AssertionError('AD value oracle failed')
        for i in range(N):
            if not r.gradient[i].contains(grad[i]):raise AssertionError('AD gradient oracle failed')
        for i,j in PAIRS:
            if not r.h(i,j).contains(h.get((i,j),F(0))):raise AssertionError('AD Hessian oracle failed')
        jetchecks+=1+N+len(PAIRS)
        # A scalar Schur complement, tested by an independently written formula.
        r=zz-yy*(xx.reciprocal())*yy
        grad=[F(0)]*N;grad[0]=y*y/(x*x);grad[1]=-2*y/x;grad[2]=1
        h={(0,0):-2*y*y/(x**3),(0,1):2*y/(x*x),(1,1):-2/x}
        if not r.value.contains(z-y*y/x):raise AssertionError('Schur value oracle failed')
        for i in range(N):
            if not r.gradient[i].contains(grad[i]):raise AssertionError('Schur gradient oracle failed')
        for i,j in PAIRS:
            if not r.h(i,j).contains(h.get((i,j),F(0))):raise AssertionError('Schur Hessian oracle failed')
        jetchecks+=1+N+len(PAIRS)
    sol,piv=solve_spd([[Box(2),Box(1)],[Box(1),Box(2)]],[[Box(1)],[Box(0)]])
    if not sol[0][0].contains(F(2,3)) or not sol[1][0].contains(F(-1,3)):raise AssertionError('LDL solve test')
    if not piv[0].contains(2) or not piv[1].contains(F(3,2)):raise AssertionError('LDL pivots test')
    return {'interval_scalar_oracle_containments':oracle_checks,'AD_value_gradient_hessian_oracle_checks':jetchecks,'LDL_oracle_checks':4,'seed':16092026,'all_passed':True}
