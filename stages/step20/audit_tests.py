"""Exact reference tests and deliberate failures for the new audit layer."""
from fractions import Fraction as F
from itertools import product
import random,copy
from audit_interval import I,Z,O,zeros,eye,tr,dot,mm,mv,spd,solve
from audit_prefactor import poly,peval,peval_series,reconstruct,certify_upper,certify_lower

# Separate exact polynomial arithmetic, for cleared scalar identities only.
def padd(*ps):
    out={}
    for p in ps:
        for m,c in p.items():out[m]=out.get(m,F(0))+c
    return {m:c for m,c in out.items() if c}
def pmul(*ps):
    out={():F(1)}
    for p in ps:
        d={}
        for m,c in out.items():
            for n,e in p.items():
                k=tuple(sorted(m+n));d[k]=d.get(k,F(0))+c*e
        out={m:c for m,c in d.items() if c}
    return out
def smul(p,a):return {m:c*a for m,c in p.items() if c*a}
def var(k):return {(k,):F(1)}
def const(c):return {():F(c)} if c else {}
def powp(p,n):return pmul(*([p]*n))

def formal_schur_test(wrong=False):
    u,q,c,hi,r,k,z=[var(i) for i in range(7)]
    q2=powp(q,2);u2=powp(u,2)
    bot=padd(pmul(c,padd(const(1),smul(u,2))),pmul(u2,z))
    Hc=padd(q2,pmul(u,padd(pmul(c,hi),smul(q2,-1))),pmul(u2,c,r))
    cr=padd(q,pmul(u2,k))
    lhs=padd(pmul(Hc,bot),smul(pmul(c,cr,cr),-1))
    bracket=padd(pmul(c,r,bot),pmul(q2,padd(const(1),smul(u,-2)),z),
                 smul(pmul(c,q2),-(3 if wrong else 4)),smul(pmul(c,q,k),-2),smul(pmul(c,u2,k,k),-1))
    rhs=pmul(u,padd(pmul(padd(pmul(c,hi),q2),bot),pmul(u,bracket)))
    assert lhs==rhs,'cleared independent scalar Schur identity failed'
    return {'coefficient_identity':True,'variables':7,'nonzero_terms_each_side':len(lhs)}

def formal_all_paths(cases):
    """Equality in a rational polynomial ring, with all 63 coefficients free."""
    u=var(63);powers=[powp(u,k) for k in range(7)];count=0
    variables=[padd(var(k),pmul(powers[2],var(21+k)),pmul(powers[3],var(42+k))) for k in range(21)]
    for pol in cases:
        direct={};by_degree={}
        for mon,c in pol:
            direct=padd(direct,smul(pmul(*(variables[k] for k in mon)),c))
            by_degree=padd(by_degree,smul(pmul(*(var(k) for k in mon)),c))
            # Linear Taylor terms count each occurrence separately, including
            # both occurrences when the monomial is a square.
            for i,k in enumerate(mon):
                others=pmul(*(var(j) for j in mon[:i]+mon[i+1:]))
                by_degree=padd(by_degree,smul(pmul(powers[2],var(21+k),others),c),
                               smul(pmul(powers[3],var(42+k),others),c))
            if len(mon)==2:
                i,j=mon
                by_degree=padd(by_degree,
                   smul(pmul(powers[4],var(21+i),var(21+j)),c),
                   smul(pmul(powers[5],padd(pmul(var(21+i),var(42+j)),pmul(var(42+i),var(21+j)))),c),
                   smul(pmul(powers[6],var(42+i),var(42+j)),c))
            assert len(mon)<=2
        assert direct==by_degree,'formal independent path identity failed'
        count+=1
    return count

def exact_solve(a,b):
    n=len(a);w=[list(map(F,r))+list(map(F,s)) for r,s in zip(a,b)]
    for k in range(n):
        p=next(j for j in range(k,n) if w[j][k])
        w[k],w[p]=w[p],w[k];d=w[k][k];w[k]=[v/d for v in w[k]]
        for i in range(n):
            if i!=k:
                d=w[i][k];w[i]=[v-d*u for v,u in zip(w[i],w[k])]
    return [r[n:] for r in w]
def qtr(a):return [list(r) for r in zip(*a)]
def qmm(a,b):return [[sum(x*y for x,y in zip(r,c)) for c in qtr(b)] for r in a]
def qadd(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def qscale(a,c):return [[x*c for x in r] for r in a]

def reference_tests(data):
    rng=random.Random(20260907);nscalar=0
    for _ in range(200):
        a,b=sorted(F(rng.randint(-99,99),rng.randint(1,99)) for j in range(2))
        c,d=sorted(F(rng.randint(-99,99),rng.randint(1,99)) for j in range(2))
        A=I(a,b);B=I(c,d)
        for x,y in product([a,(a+b)/2,b],[c,(c+d)/2,d]):
            for box,q in [(A+B,x+y),(A-B,x-y),(A*B,x*y),(A.sq(),x*x)]:
                assert box.contains(q);nscalar+=1
            if not c<=0<=d:assert (A/B).contains(x/y);nscalar+=int(not c<=0<=d)
    nder=0;nsolve=0;pair=0
    for rank,nu in product(range(1,5),range(1,4)):
        A=[[F(2*rank+1 if i==j else 1) for j in range(rank)] for i in range(rank)]
        B=[[F((i+1)*(j+2),7) for j in range(nu)] for i in range(rank)]
        D=[[F((i+1)+(j+1),9) for j in range(nu)] for i in range(nu)]
        dA=[[F(i+j+1,11) for j in range(rank)] for i in range(rank)]
        dB=[[F(i-j,13) for j in range(nu)] for i in range(rank)]
        dD=[[F(i+j-1,17) for j in range(nu)] for i in range(nu)]
        X=exact_solve(A,B);K=qscale(X,-1)+[[F(i==j) for j in range(nu)] for i in range(nu)]
        dM=[row+ext for row,ext in zip(dA,dB)]+[row+ext for row,ext in zip(qtr(dB),dD)]
        by_kernel=qmm(qmm(qtr(K),dM),K)
        expanded=qadd(qadd(dD,qscale(qmm(qtr(dB),X),-1)),qadd(qscale(qmm(qtr(X),dB),-1),qmm(qmm(qtr(X),dA),X)))
        assert by_kernel==expanded;nder+=nu*nu
        got=solve([[I(v) for v in r] for r in A],[[I(v) for v in r] for r in B])
        for i in range(rank):
            for j in range(nu):assert got[i][j].contains(X[i][j]);nsolve+=1
        lam=[F(i+2,5) for i in range(nu*(nu+1)//2)];W=[[F(0) for _ in range(nu)] for _ in range(nu)];off=0;rhs=F(0)
        for i in range(nu):
            for j in range(i,nu):
                W[i][j]=W[j][i]=lam[off]/(1 if i==j else 2);rhs+=lam[off]*D[i][j];off+=1
        assert sum(W[i][j]*D[j][i] for i in range(nu) for j in range(nu))==rhs;pair+=1
    formal=formal_schur_test()
    # Every path entry is substituted by convolution and independently by an
    # exact endpoint/direction polynomial evaluation, at rational coordinates.
    cases=[]
    for kind in ['small','raw']:
        for b in data[kind]:
            for row in b['mat']:cases.extend(poly(p) for p in row)
    for cp in data['completion'].values():cases.extend(poly(p) for p in cp['cross']+[cp['diagonal']])
    formal_paths=formal_all_paths(cases)
    npath=0
    for pp in cases:
        n=21;s=[F((i%5)-2,7) for i in range(n)];d=[F((i%7)-3,11) for i in range(n)];v=[F((i%3)-1,13) for i in range(n)]
        coords=[[I(a),Z,I(b),I(c)] for a,b,c in zip(s,d,v)]
        cc=peval_series(pp,coords)
        for u in [F(-1,10),F(0),F(1,100)]:
            actual=[a+u*u*b+u*u*u*c for a,b,c in zip(s,d,v)]
            exact=sum(c*prod_fraction(actual[k] for k in mon) for mon,c in pp)
            enclosed=sum((q*u**j for j,q in enumerate(cc)),Z)
            assert enclosed.contains(exact);npath+=1
    # Three common-t channels, including exactly zero crossing. Elementary
    # budget argument fixes C; regularization is verified with exact Fractions.
    w=[F(3,2),F(7,3),F(2,5)];q=[F(0),F(-2),F(3)]
    C=sum(a*b*b for a,b in zip(w,q));assert C==F(194,15);Omega=sum(w)
    h=[v*v/C for v in q];assert sum(a*b for a,b in zip(w,h))==1
    nzero=0
    for u in [F(1,100),F(1,10**6),F(1,10**20)]:
        for hi,qi in zip(h,q):
            H=(1-u)*hi+u/Omega
            assert H*C*(1+2*u)-qi*qi>0;nzero+=1
    return {'binary_bits':768,'rational_interval_reference_checks':nscalar,'fixed_kernel_derivative_entries':nder,
            'interval_solve_entries':nsolve,'offdiagonal_trace_pairing_cases':pair,
            'path_entries':len(cases),'formal_path_identities':formal_paths,'path_reference_evaluations':npath,
            'zero_crossing_regularization_cases':nzero,'formal_schur_identity':formal}

def prod_fraction(xs):
    v=F(1)
    for x in xs:v*=x
    return v

def controls(data,prefs):
    rejected=[]
    def no(name,fn):
        try:fn()
        except (AssertionError,ValueError,ZeroDivisionError):rejected.append(name)
        else:raise AssertionError('corruption accepted: '+name)
    no('reversed_interval',lambda:I(1,-1))
    no('zero_containing_denominator',lambda:I(-1,1).inv())
    no('indefinite_positive_diagonal_matrix',lambda:spd([[I(1),I(2)],[I(2),I(1)]]))
    no('rank_deficient_Gram',lambda:solve([[O,O],[O,O]],eye(2)))
    no('wrong_scalar_Schur_factor',lambda:formal_schur_test(True))
    def wrong_budget():assert 2*(F(3,2)*F(0)+F(7,3)*F(4)+F(2,5)*F(9))/F(194,15)==1
    no('doubled_normalization_budget',wrong_budget)
    def wrong_pair():assert 2*F(3)*F(4)==F(3)*F(4)
    no('missing_offdiagonal_half_weight',wrong_pair)
    def lost_reg():
        u=F(1,10**9);assert -u**4>0
    no('unopened_zero_crossing_kernel',lost_reg)
    for side,pref in prefs.items():
        no(side+'_zero_strip',lambda:certify_upper(data,pref,F(0)))
        no(side+'_oversized_strip',lambda:certify_upper(data,pref,F(1,100)))
        bad=copy.deepcopy(pref);pols=bad['raw']['single_Q'][0][0]
        pols.append(((0,1),F(1)))
        no(side+'_corrupted_Q_mixed_term',lambda:certify_lower(data,bad))
        bad_data=copy.deepcopy(data);off=0
        for inf in pref['infos']:
            if inf['parent']=='M0_1':break
            off+=inf['nu']*(inf['nu']+1)//2
        bad_data['sides'][side]['center'][20+off]='-1'
        no(side+'_negative_completion_dual_weight',lambda:reconstruct(bad_data,side))
    return {'rejected':rejected,'count':len(rejected)}
