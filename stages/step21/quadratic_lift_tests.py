"""Exact checks for the new second-order Schur lifting construction."""
from fractions import Fraction as F
import random
from exact_algebra import *


def block(A,E,D):
    return [a+e for a,e in zip(A,E)]+[e+d for e,d in zip(transpose(E),D)]


def run():
    # Exact curved model g(x,z)=(z,x-z-x^2).
    J=matrix([[0,1],[1,-1]]);R=inverse(J)
    h0=[F(1,5),F(4,5)];hi=[F(1,2),F(1,2)]
    d0=mv(R,h0);g2=[F(0),F(-2)];d2=mv(R,[h-v/2 for h,v in zip(hi,g2)])
    if d0!=[F(1),F(1,5)] or d2!=[F(2),F(1,2)]:raise AssertionError('quadratic lift directions')
    xp=[F(0),F(1),F(2)];zp=[F(0),F(1,5),F(1,2)]
    h1=zp;h2=padd(padd(xp,pscale(zp,-1)),pscale(pmul(xp,xp),-1))
    if h2!=[F(0),F(4,5),F(1,2),F(-4),F(-4)]:raise AssertionError('Schur path coefficients')
    # t=5/x+10. Multiply each determinant by x>0 to avoid rational functions.
    fac=padd([F(5)],pscale(xp,10))
    num1=padd(pmul(fac,h1),pscale(xp,-1))
    num2=padd(pmul(fac,h2),pscale(xp,-4))
    if num1!=[F(0),F(0),F(5,2),F(9),F(10)]:raise AssertionError('border 1 identity')
    if num2!=[F(0),F(0),F(5,2),F(1),F(-50),F(-120),F(-80)]:raise AssertionError('border 2 identity')
    r0=F(1,100)
    margin1=F(4,5)-4*r0*r0-4*r0**3
    margin2=F(5,2)-50*r0*r0-120*r0**3-80*r0**4
    if margin1<=0 or margin2<=0:raise AssertionError('continuous polynomial positivity')
    samples=[]
    for r in [r0,F(1,10**4),F(1,10**8),F(1,10**12)]:
        x=peval(xp,r);z=peval(zp,r);t=5/x+10
        certs=[psd([[z,1],[1,t]]),psd([[x-z-x*x,2],[2,t]])]
        samples.append({'r':r,'delta':x,'t':t,'certificates':certs})
    # An independent power-series derivation of g'' compared to fixed-kernel formula.
    rng=random.Random(210908);count=0;entries=0
    for case in range(24):
        n=1+case%3;k=1+(case//3)%3
        M=matrix([[rng.randint(-2,2) for j in range(n)] for i in range(n)])
        A0=add(mul(transpose(M),M),eye(n));Ai=inverse(A0)
        def symmetric(m):
            a=matrix([[rng.randint(-2,2) for j in range(m)] for i in range(m)])
            return add(a,transpose(a))
        A1=symmetric(n);A2=symmetric(n)
        E0=matrix([[rng.randint(-2,2) for j in range(k)] for i in range(n)])
        E1=matrix([[rng.randint(-2,2) for j in range(k)] for i in range(n)])
        E2=matrix([[rng.randint(-2,2) for j in range(k)] for i in range(n)])
        D0=mul(mul(transpose(E0),Ai),E0);D1=symmetric(k);D2=symmetric(k)
        X=mul(Ai,E0);K=scale(X,-1)+eye(k)
        M1=block(A1,E1,D1);M2=block(A2,E2,D2)
        T=mul(M1,K)[:n]
        direct=add(scale(mul(mul(transpose(K),M2),K),2),scale(mul(mul(transpose(T),Ai),T),-2))
        Ai1=scale(mul(mul(Ai,A1),Ai),-1)
        Ai2=add(mul(mul(mul(mul(Ai,A1),Ai),A1),Ai),scale(mul(mul(Ai,A2),Ai),-1))
        AA=[Ai,Ai1,Ai2];EE=[E0,E1,E2];S2=D2
        for i in range(3):
            for j in range(3):
                kk=2-i-j
                if 0<=kk<3:S2=add(S2,scale(mul(mul(transpose(EE[i]),AA[j]),EE[kk]),-1))
        if direct!=scale(S2,2):raise AssertionError('second Schur derivative differs from series')
        count+=1;entries+=k*k
    return {'curved_model':{'J':J,'d0':d0,'d2':d2,'x_polynomial':xp,'z_polynomial':zp,
              'r_interval_upper':r0,'parameter_delta_upper':peval(xp,r0),
              'all_positive_r_in_interval_complete':True,'chosen_t':'5/delta+10',
              'cleared_determinant_polynomials':[num1,num2],
              'uniform_positive_lower_coefficients':[margin1,margin2],'samples':samples},
            'second_Schur_derivative_reference_cases':count,'matrix_entries_compared':entries,
            'identity':'Sdd=K^T Mdd K-2*(Md K)_range^T A^-1*(Md K)_range',
            'scope':'The interval positivity conclusion follows from exact coefficient bounds, not just the four displayed samples.'}
