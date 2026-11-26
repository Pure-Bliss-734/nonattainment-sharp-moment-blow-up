"""Small exact reference checks for the NEW identities; not a formal prover."""
from fractions import Fraction as F
from random import Random


def pa(*ps):
    out={}
    for p in ps:
        for m,c in p.items():out[m]=out.get(m,F(0))+c
    return {m:c for m,c in out.items() if c}
def pm(p,q):
    out={}
    for a,c in p.items():
        for b,d in q.items():
            m=tuple(sorted(a+b));out[m]=out.get(m,F(0))+c*d
    return {m:c for m,c in out.items() if c}
def ps(p,c):return {m:v*c for m,v in p.items() if v*c}
def var(k):return {(k,):F(1)}
def const(c):return {():F(c)} if c else {}

def matmul(A,B):return [[sum((x*y for x,y in zip(row,col)),F(0)) for col in zip(*B)] for row in A]
def transpose(A):return [list(r) for r in zip(*A)]
def trace(A):return sum((A[i][i] for i in range(len(A))),F(0))
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def eye(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def inverse(A):
    n=len(A);M=[list(row)+e for row,e in zip(A,eye(n))]
    for k in range(n):
        pivot=next((i for i in range(k,n) if M[i][k]),None)
        if pivot is None:raise ValueError('singular reference matrix')
        M[k],M[pivot]=M[pivot],M[k];d=M[k][k];M[k]=[v/d for v in M[k]]
        for i in range(n):
            if i!=k:
                d=M[i][k];M[i]=[a-d*b for a,b in zip(M[i],M[k])]
    return [r[n:] for r in M]

def formal_square_identity(half=F(1,2)):
    C,delta,D,r,k,t=[var(i) for i in range(6)]
    A=pa(ps(pm(C,pa(const(1),D)),half),pm(delta,r))
    T=pm(delta,D);B=pa(pm(C,t),k)
    lhs=pa(pm(pm(delta,delta),B),ps(pm(pm(C,delta),A),-2),pm(pm(C,C),T))
    rhs=pa(pm(pm(delta,delta),B),ps(pm(pm(C,C),delta),-1),ps(pm(pm(pm(C,delta),delta),r),-2))
    if pa(lhs,ps(rhs,-1)):raise AssertionError('formal completed-square identity failed')
    return True

def run():
    formal_square_identity();counts={'formal_completed_square_identity':1}
    # Root factorization checked with xi,p1,p2,z kept formal.
    X,xi,z,a,b,c,d=[var(i) for i in range(7)]
    p1=pa(a,pm(b,z));p2=pa(c,pm(d,z))
    p=pa(pm(pa(X,ps(xi,-1)),p1),pm(pa(pm(X,X),ps(pm(xi,xi),-1)),p2))
    q=pm(pa(X,ps(xi,-1)),pa(p1,pm(pa(X,xi),p2)))
    if pa(p,ps(q,-1)):raise AssertionError('formal root quotient')
    counts['formal_anchored_quotient_identity']=1
    rng=Random(20260907);nproj=ntrace=nsolve=0
    for case in range(24):
        n=3+(case%2);r=1+case%2
        Fmat=[[F(rng.randint(-4,4),rng.randint(1,5)) for j in range(n)] for i in range(n)]
        G=matmul(transpose(Fmat),Fmat)
        for i in range(n):G[i][i]+=1
        Rmat=[[F(rng.randint(-3,3),rng.randint(1,4)) for j in range(r)] for i in range(r)]
        W=matmul(transpose(Rmat),Rmat)
        for i in range(r):W[i][i]+=1
        K=[[F(rng.randint(-5,5),7) for j in range(r)] for i in range(n)]
        V=[[F(rng.randint(-5,5),11) for j in range(r)] for i in range(n)]
        lam=F(case+1,13)
        T=trace(matmul(W,matmul(transpose(K),matmul(G,K))))
        A=trace(matmul(W,matmul(transpose(K),matmul(G,V))))
        B=trace(matmul(W,matmul(transpose(V),matmul(G,V))))
        Q=[[v-lam*k for v,k in zip(vrow,krow)] for vrow,krow in zip(V,K)]
        sq=trace(matmul(W,matmul(transpose(Q),matmul(G,Q))))
        if sq!=B-2*lam*A+lam*lam*T or sq<0:raise AssertionError('exact weighted Gram identity')
        ntrace+=1
        E=matmul(V,W);V2=matmul(E,inverse(W))
        if V2!=V:raise AssertionError('V W=E orientation')
        pairing=trace(matmul(transpose(K),matmul(G,E)))
        if pairing!=A:raise AssertionError('off-diagonal trace pairing')
        ntrace+=1
        # H has a nonsingular leading square block; P is strictly positive.
        H=[[F(1),F(0)],[F(0),F(1)]]+[[F(rng.randint(-4,4),3),F(rng.randint(-4,4),5)] for _ in range(3)]
        P=[[F(0) for _ in range(5)] for _ in range(5)]
        for i in range(5):P[i][i]=F(i+1,case+1)
        HP=matmul(transpose(H),P);N=matmul(HP,H)
        g=[[F(case+2,17)],[F(3-case,19)]]
        al=matmul(matmul(P,H),matmul(inverse(N),g))
        if matmul(transpose(H),al)!=g:raise AssertionError('weighted inverse equality')
        nproj+=1
        if matmul(N,inverse(N))!=eye(2):raise AssertionError('reference inverse')
        nsolve+=4
    counts.update({'exact_weighted_trace_checks':ntrace,'exact_right_inverse_checks':nproj,'reference_inverse_entries':nsolve})
    # Contrasting signed D: the proof uses a square, not division by D.
    for Dv in [F(0),F(1,3),F(1),F(2),F(7,2)]:
        C=F(3,7);delta=F(1,1000);rr=F(-1,13);kk=F(5,11);tt=F(1234)
        Av=C*(1+Dv)/2+delta*rr
        Sv=C*tt+kk-2*C/delta*Av+C*C/delta*Dv
        if Sv!=C*(tt-C/delta+kk/C-2*rr):raise AssertionError('zero-D branch identity')
    counts['D_zero_and_positive_scalar_checks']=5
    bad=[]
    for h in [F(0),F(1,4),F(3,4),F(1)]:
        try:formal_square_identity(h)
        except AssertionError:bad.append(str(h))
        else:raise AssertionError('wrong midpoint coefficient accepted')
    counts['incorrect_formal_midpoint_coefficients_rejected']=bad
    return counts
