"""Independent rational reference tests for the shifted partner formula."""
from fractions import Fraction as F

def solve(A,b):
    n=len(A);M=[[F(x) for x in row]+[F(b[i])] for i,row in enumerate(A)]
    for j in range(n):
        k=next((k for k in range(j,n) if M[k][j]),None)
        if k is None:raise ValueError('singular')
        M[j],M[k]=M[k],M[j];p=M[j][j];M[j]=[x/p for x in M[j]]
        for k in range(n):
            if k!=j:
                q=M[k][j];M[k]=[x-q*y for x,y in zip(M[k],M[j])]
    return [row[-1] for row in M]

def run():
    constraints=gradients=variations=0
    for k in range(1,13):
        q=[F(k+1),F(k+2,3),F(k+4,5)];P=[1/v for v in q]
        H=[[F(1),F(2)],[F(2),F(1)],[F(k,7),F(-k,11)]]
        l=[F(k,13),F(-1,3),F(2,5)];g=[F(1,7),F(k,17)]
        PH=[[P[i]*v for v in row] for i,row in enumerate(H)]
        Gram=[[sum(H[i][a]*PH[i][b] for i in range(3)) for b in range(2)] for a in range(2)]
        rhs=[g[a]+sum(H[i][a]*P[i]*l[i] for i in range(3)) for a in range(2)]
        v=solve(Gram,rhs)
        alpha=[-P[i]*l[i]+sum(PH[i][a]*v[a] for a in range(2)) for i in range(3)]
        for a in range(2):
            assert sum(H[i][a]*alpha[i] for i in range(3))==g[a];constraints+=1
        for i in range(3):
            assert q[i]*alpha[i]+l[i]==sum(H[i][a]*v[a] for a in range(2));gradients+=1
        n=solve([[H[i][a] for i in range(2)] for a in range(2)],[-H[2][a] for a in range(2)])+[F(1)]
        def obj(z):return sum(q[i]*z[i]**2+2*l[i]*z[i] for i in range(3))
        for s in (F(-2),F(1,7),F(3)):
            z=[a+s*b for a,b in zip(alpha,n)]
            assert obj(z)-obj(alpha)==s*s*sum(q[i]*n[i]**2 for i in range(3))>0;variations+=1
    return {'rational_constraint_identities':constraints,'rational_stationarity_identities':gradients,
            'positive_quadratic_variation_checks':variations,'role':'Reference tests, not quantifier elimination or a formal theorem-prover proof.'}
