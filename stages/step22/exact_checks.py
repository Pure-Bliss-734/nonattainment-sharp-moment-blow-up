"""Independent small Fraction/polynomial checks of the new identities.

The continuous theorem comes from the written identities plus interval bounds;
finite reference tests do not establish a universal mathematical claim.
"""
from fractions import Fraction as F
import random

def padd(*ps):
    out={}
    for p in ps:
        for m,c in p.items():out[m]=out.get(m,F(0))+c
    return {m:c for m,c in out.items() if c}
def pscale(p,a):return {m:c*a for m,c in p.items() if c*a}
def pmul(p,q):
    out={}
    for m,c in p.items():
        for n,d in q.items():
            key=tuple(sorted(m+n));out[key]=out.get(key,F(0))+c*d
    return {m:c for m,c in out.items() if c}
def pconst(a):return {():F(a)} if a else {}
def var(i):return {(i,):F(1)}
def parse(a):
    p={}
    for m,c in a:p=padd(p,{tuple(m):F(c)})
    return p
def compose(p,coords):
    out={}
    for mon,c in p.items():
        term=pconst(c)
        for k in mon:term=pmul(term,coords[k])
        out=padd(out,term)
    return out
def value(p,coords):
    out=F(0)
    for mon,c in p.items():
        for k in mon:c*=coords[k]
        out+=c
    return out

def symmetry_and_path(data):
    # Constant rational parity charts, verified over formal endpoint variables.
    raw_to=[parse(p) for p in data['to_raw']];nchart=entries=0
    rawsym={}
    for b in data['raw']:
        P=[[compose(parse(q),raw_to) for q in row] for row in b['mat']];rawsym[b['name']]=P
        children=[v for v in data['small'] if v['parent']==b['name']]
        cols=[];expected=[]
        for child in children:cols += list(zip(*[[F(v) for v in row] for row in child['T']]))
        n=len(P)
        if len(cols)!=n:raise AssertionError('missing parity column')
        C=[[F(0) for _ in range(n)] for _ in range(n)];offset=0
        exp=[[{} for _ in range(n)] for _ in range(n)]
        for child in children:
            for i,row in enumerate(child['mat']):
                for j,p in enumerate(row):exp[offset+i][offset+j]=parse(p)
            offset+=len(child['mat'])
        # Verify invertibility over Fraction, independently of interval pivoting.
        inverse([list(row) for row in zip(*cols)]);nchart+=1
        for i in range(n):
            for j in range(n):
                Q={}
                for k in range(n):
                    if not cols[i][k]:continue
                    for l in range(n):
                        if cols[j][l]:Q=padd(Q,pscale(P[k][l],cols[i][k]*cols[j][l]))
                if Q!=exp[i][j]:raise AssertionError('parity polynomial congruence')
                entries+=1
    pols=[]
    for b in data['small']:
        for row in b['mat']:pols += [parse(p) for p in row]
    for rows in rawsym.values():
        for row in rows:pols+=row
    for cp in data['completion'].values():
        pols += [compose(parse(p),raw_to) for p in cp['cross']]
        pols += [compose(parse(cp['diagonal']),raw_to)]
    # r uses index 60, centers 0..19, d0 20..39, d2 40..59.
    coords=[padd(var(j),pmul(var(60),var(20+j)),pmul(pmul(var(60),var(60)),var(40+j))) for j in range(20)]
    from bounded_model import ap
    vals=[F((j*7)%19-9,j+2) for j in range(60)];references=0
    for p in pols:
        direct=compose(p,coords);by_degree=[{} for _ in range(5)]
        for mon,c in p.items():
            if len(mon)==0:by_degree[0]=padd(by_degree[0],pconst(c))
            elif len(mon)==1:
                i=mon[0]
                for k in range(3):by_degree[k]=padd(by_degree[k],pscale(var(i+20*k),c))
            elif len(mon)==2:
                i,j=mon
                for a in range(3):
                    for b in range(3):by_degree[a+b]=padd(by_degree[a+b],pscale(pmul(var(i+20*a),var(j+20*b)),c))
            else:raise AssertionError('nonquadratic pencil')
        reassembled={}
        for k,q in enumerate(by_degree):reassembled=padd(reassembled,pmul(q,{(60,)*k:F(1)}))
        if direct!=reassembled:raise AssertionError('formal quadratic-path expansion')
        z=[[ap.I(vals[j]),ap.I(vals[j+20]),ap.I(vals[j+40])] for j in range(20)]
        intervals=ap.peval_series(list(p.items()),z)
        for k in range(7):
            ref=value(by_degree[k],vals) if k<5 else F(0)
            if not intervals[k].contains(ref):raise AssertionError('path interval excluded exact coefficient')
            references+=1
    return {'formal_parity_charts':nchart,'polynomial_congruence_entries':entries,
            'formal_quadratic_path_entries':len(pols),'exact_coefficient_inclusion_checks':references}

def inverse(A):
    n=len(A);w=[list(map(F,row))+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for k in range(n):
        found=next((j for j in range(k,n) if w[j][k]),None)
        if found is None:raise ValueError('singular')
        w[k],w[found]=w[found],w[k];p=w[k][k];w[k]=[v/p for v in w[k]]
        for i in range(n):
            if i==k:continue
            v=w[i][k];w[i]=[a-v*b for a,b in zip(w[i],w[k])]
    return [r[n:] for r in w]
def mt(A):return list(map(list,zip(*A)))
def mm(A,B):return [[sum((a*b for a,b in zip(row,col)),F(0)) for col in zip(*B)] for row in A]
def plus(A,B):return [[a+b for a,b in zip(r,s)] for r,s in zip(A,B)]
def times(A,c):return [[a*c for a in r] for r in A]
def minus(A,B):return plus(A,times(B,-1))
def eye(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def evalmat(ms,r):
    out=times(ms[0],0)
    for a in reversed(ms):out=plus(a,times(out,r))
    return out

def divided_difference_tests():
    rng=random.Random(22022);entries=0
    for case in range(24):
        rank=1+case%3;nu=1+case%2
        rand=lambda n,m:[[F(rng.randrange(-5,6),7) for _ in range(m)] for _ in range(n)]
        V=rand(rank,rank);A0=plus(eye(rank),mm(mt(V),V));As=[A0]
        for _ in range(4):
            a=rand(rank,rank);As.append(times(plus(a,mt(a)),F(1,2)))
        Es=[rand(rank,nu) for _ in range(4)]
        Ds=[]
        for _ in range(4):
            d=rand(nu,nu);Ds.append(times(plus(d,mt(d)),F(1,2)))
        r=F(1,100+case);A=evalmat(As,r);Ai=inverse(A);A0i=inverse(A0)
        Et=evalmat(Es,r);Eh=evalmat(Es[1:],r);Ah=evalmat(As[1:],r);G0=mm(A0i,Es[0]);G=mm(Ai,Et)
        Dg=mm(Ai,minus(Eh,mm(Ah,G0)))
        if minus(G,G0)!=times(Dg,r):raise AssertionError('inverse divided difference')
        R=minus(evalmat(Ds[2:],r),plus(mm(mt(Eh),G),mm(mt(Es[0]),Dg)))
        H1=Ds[0];H2=minus(Ds[1],mm(mt(Es[0]),G0))
        actual=minus(times(evalmat(Ds,r),r),times(mm(mt(Et),G),r*r))
        expected=plus(times(H1,r),plus(times(H2,r*r),times(R,r**3)))
        if actual!=expected:raise AssertionError('Schur cancellation identity')
        entries+=nu*nu
    # Fully formal scalar scaled Schur identity; C*h0=q^2 is substituted.
    r,C,q,k,hI,R,m,Z=map(var,range(8));one=pconst(1)
    B=padd(C,pmul(r,padd(m,Z)))
    top=padd(pmul(q,q),pmul(C,pmul(r,hI)),pmul(C,pmul(pmul(r,r),R)))
    left=padd(pmul(B,top),pscale(pmul(C,pmul(padd(q,pmul(r,k)),padd(q,pmul(r,k)))),-1))
    inn=padd(pmul(pmul(q,q),padd(m,Z)),pscale(pmul(C,pmul(q,k)),-2),pscale(pmul(C,pmul(r,pmul(k,k))),-1))
    right=padd(pmul(pmul(C,pmul(r,padd(hI,pmul(r,R)))),B),pmul(r,inn))
    if left!=right:raise AssertionError('formal scalar Schur identity')
    # Exact conversion identity: (C/r+m)-(C/[r(1+ar)]+m+Ca/(1+ar))=0.
    a=var(8);factor=padd(one,pmul(a,r))
    conv=padd(pmul(C,factor),pscale(C,-1),pscale(pmul(pmul(C,a),r),-1))
    if conv:raise AssertionError('parameter conversion identity')
    return {'rational_matrix_cases':24,'Schur_identity_entries':entries,
            'formal_scaled_Schur_identity':True,'formal_parameter_conversion_identity':True}

def run(data):return {'structure':symmetry_and_path(data),'divided_differences':divided_difference_tests()}
