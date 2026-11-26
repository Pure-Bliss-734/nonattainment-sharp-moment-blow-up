"""Exact finite-completion model. Proof path: Python standard library only.
The representative has z20=0; t >= 0 is the eighth-coordinate moment.
"""
from __future__ import annotations
from functools import lru_cache
from fractions import Fraction as F
from upstream_loader import ROOT,UP
from step12_model import build,evaluate_matrix,evaluate_form,peval,linear,zero,wa,md,check
from step12_layout import projected
from verify_step10 import exact_psd_rank
RDIR=[F(0)]*14+[F(1)]*7

def dot(a,b):return sum((v*w for v,w in zip(a,b)),F(0))
def mv(A,v):return [dot(row,v) for row in A]
def transpose(A):return [list(row) for row in zip(*A)]
def mm(A,B):return [[dot(row,col) for col in transpose(B)] for row in A]
def solve(A,b):
    """Rational solve of A*x=b. Free coordinates are set to zero.
    Raise ValueError on an inconsistent system; no small pivot tolerances.
    """
    if len(A)!=len(b):raise ValueError('RHS dimension mismatch')
    n=len(A[0]) if A else 0
    C=[list(map(F,row))+[F(v)] for row,v in zip(A,b)];r=0;pivots=[]
    for j in range(n):
        k=next((k for k in range(r,len(C)) if C[k][j]),None)
        if k is None:continue
        C[r],C[k]=C[k],C[r];p=C[r][j];C[r]=[v/p for v in C[r]]
        for i in range(len(C)):
            if i==r:continue
            p=C[i][j]
            if p:C[i]=[a-p*v for a,v in zip(C[i],C[r])]
        pivots.append(j);r+=1
    for row in C:
        if not any(row[:n]) and row[n]:raise ValueError('RHS not in matrix range')
    x=[F(0)]*n
    for i,j in enumerate(pivots):x[j]=C[i][n]
    check(mv(A,x)==list(b),'Exact solve failed')
    return x

def operator_linear(coefs,ops):
    out={}
    for c,op in zip(coefs,ops):
        for w,v in op.items():wa.add_term(out,w,c*v)
    return wa.clean(out)

def operator_json(op):
    return [{'word':list(w),'coefficient':str(v)} for w,v in sorted(op.items())]

@lru_cache(None)
def model():
    b=build();P={v['name']:v for v in projected()};out={}
    for name in ('M0_1','M2_None'):
        raw=next(bl for bl in b['blocks'] if bl['name']==name)
        W=next(op for op in raw['basis'] if all(len(w)==4 and all(i<2 for i in w) for w in op) and sum(op.values())==1)
        ops=P[name]['basis']
        out[name]={'raw':raw,'projected':P[name],'W':W,'cross':[b['old_form'](wa.inner(op,W)) for op in ops],
                   'diagonal':b['old_form'](wa.inner(W,W))}
    return b,P,out

def augmented(name,x,z,t):
    _,_,C=model();c=C[name]
    B=evaluate_matrix(c['projected'],x,z)
    v=[evaluate_form(f,x,z) for f in c['cross']]
    a=evaluate_form(c['diagonal'],x,z)+t
    return [row+[v[i]] for i,row in enumerate(B)]+[v+[a]]

def thresholds(x,z):
    check(len(z)==21 and z[20]==0,'Use representative z20=0')
    b,P,C=model();out={}
    for bl in P.values():exact_psd_rank(evaluate_matrix(bl,x,z))
    for name,c in C.items():
        B=evaluate_matrix(c['projected'],x,z);v=[evaluate_form(f,x,z) for f in c['cross']]
        y=solve(B,v)
        value=dot(v,y)-evaluate_form(c['diagonal'],x,z)
        out[name]={'value':value,'solution':y,'cross':v}
    out['minimum']=max(F(0),*(out[n]['value'] for n in C))
    return out

def full_psd(x,z,t):
    b,_,_=model();full=[v+t*d for v,d in zip(z,RDIR)]
    blocks=[]
    for bl in b['blocks']:
        pp=exact_psd_rank(evaluate_matrix(bl,x,full))
        blocks.append({'name':bl['name'],'rank':len(pp),'positive_pivots':list(map(str,pp))})
    return blocks
