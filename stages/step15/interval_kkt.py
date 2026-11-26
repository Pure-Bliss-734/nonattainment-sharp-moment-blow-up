"""Rational interval evaluation of the complete symmetry-reduced KKT system.

Differentiation formulas are Schur-complement identities; no equations are
selected or discarded numerically in this verifier.
"""
from fractions import Fraction as F
from dyadic_intervals import I, zmat,eye,tr,mm,mv,dot,sub,ip,pd_pivots,inverse
from symmetric_model import small_blocks

class KKTIntervals:
    def __init__(self,side,layout,nc):
        self.side=side;self.layout=layout;self.nc=nc;self.blocks=small_blocks();self.nv=20
        self.polys=[]
        for bl in self.blocks:
            n=len(bl['mat']);P=[[[[F(0) for _ in range(n)] for __ in range(n)] for k in range(20)] for d in range(3)]
            for i,row in enumerate(bl['mat']):
                for j,form in enumerate(row):
                    for k,p in enumerate(form):
                        if len(p)>3:raise ValueError('unexpected polynomial degree')
                        for r,c in enumerate(p):P[r][k][i][j]=c
            self.polys.append(P)
        assert sum(it['nu']*(it['nu']+1)//2 for it in layout)==nc
    def matrices(self,index,y):
        P=self.polys[index];x=y[0];n=len(P[0][0]);z=[I(1)]+y[1:20]
        coeff=[];dercoeff=[]
        for k in range(20):
            A=zmat(n);D=zmat(n)
            for i in range(n):
                for j in range(i,n):
                    A[i][j]=A[j][i]=I(P[0][k][i][j])+x*P[1][k][i][j]+x*x*P[2][k][i][j]
                    D[i][j]=D[j][i]=I(P[1][k][i][j])+2*x*P[2][k][i][j]
            coeff.append(A);dercoeff.append(D)
        M=zmat(n);Dx=zmat(n);Hxx=zmat(n)
        for i in range(n):
            for j in range(i,n):
                M[i][j]=M[j][i]=sum((z[k]*coeff[k][i][j] for k in range(20)),I(0))
                Dx[i][j]=Dx[j][i]=sum((z[k]*dercoeff[k][i][j] for k in range(20)),I(0))
                Hxx[i][j]=Hxx[j][i]=sum((z[k]*2*P[2][k][i][j] for k in range(20)),I(0))
        return M,[Dx]+coeff[1:],[Hxx]+dercoeff[1:]
    def evaluate(self,v):
        y=v[:20];lam=v[20:];J=[];fs=[];HW=zmat(20);infos=[];off=0
        for it in self.layout:
            M,D,HX=self.matrices(it['index'],y);perm=it['perm'];r=it['rank'];nu=it['nu'];n=len(M)
            assert r+nu==n and sorted(perm)==list(range(n))
            def re(A):return [[A[i][j] for j in perm] for i in perm]
            PM=re(M);PD=[re(a) for a in D];PH=[re(a) for a in HX]
            A=[row[:r] for row in PM[:r]]
            if r:
                Ai=inverse(A);W=mm(Ai,[row[r:] for row in PM[:r]])
                E=[[-v for v in row] for row in W]+eye(nu)
                Sch=sub([row[r:] for row in PM[r:]],mm([row[:r] for row in PM[r:]],W))
            else:E=eye(nu);Sch=PM
            ET=tr(E);grad=[mm(mm(ET,a),E) for a in PD]
            Z=zmat(nu);kk=0
            for i in range(nu):
                for j in range(i,nu):
                    fs.append(Sch[i][j]);J.append([gg[i][j] for gg in grad]);Z[i][j]=Z[j][i]=lam[off+kk]/(1 if i==j else 2);kk+=1
            hrow=[ip(Z,mm(mm(ET,a),E)) for a in PH]
            for k in range(20):
                HW[0][k]=HW[0][k]+hrow[k]
                if k:HW[k][0]=HW[k][0]+hrow[k]
            if r:
                V=[sub([row[r:] for row in a[:r]],mm([row[:r] for row in a[:r]],W)) for a in PD]
                VZ=[mm(a,Z) for a in V];AiV=[mm(Ai,a) for a in V]
                for i in range(20):
                    for j in range(i,20):
                        h=2*ip(VZ[i],AiV[j])
                        HW[i][j]=HW[i][j]-h
                        if j!=i:HW[j][i]=HW[j][i]-h
            K=zmat(n,nu)
            for i,ii in enumerate(perm):K[ii]=E[i]
            infos.append({'item':it,'matrix':M,'A':A,'K':K,'multiplier':Z,'dual':mm(mm(K,Z),tr(K)),'Schur':Sch})
            off+=kk
        assert len(fs)==self.nc and off==self.nc
        sigma=1 if self.side=='lower' else -1
        station=[I(sigma if j==0 else 0)-dot([row[j] for row in J],lam) for j in range(20)]
        Jac=[row+[I(0)]*self.nc for row in J]
        Jac += [[-v for v in HW[j]]+[-J[k][j] for k in range(self.nc)] for j in range(20)]
        return fs+station,Jac,infos
