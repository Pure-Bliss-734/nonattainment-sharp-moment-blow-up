"""Exact reflection reduction w=z2+a/2 (discovery and verification model)."""
from functools import lru_cache
from fractions import Fraction as F
import json
import load_frozen
from step12_model import build
from step12_layout import projected

# Form representation: rows are constant, then the 19 retained coordinates;
# each coefficient is a polynomial in x of degree at most two.
KEEP=[i for i in range(20) if i!=2]
def padd(p,q):
    a=[F(0)]*max(len(p),len(q))
    for i,v in enumerate(p):a[i]+=v
    for i,v in enumerate(q):a[i]+=v
    while len(a)>1 and a[-1]==0:a.pop()
    return tuple(a)
def pscale(p,c):return tuple(c*a for a in p)
def pmul(p,q):
    a=[F(0)]*(len(p)+len(q)-1)
    for i,v in enumerate(p):
        for j,w in enumerate(q):a[i+j]+=v*w
    while len(a)>1 and not a[-1]:a.pop()
    return tuple(a)
def fscale(f,c):return tuple(pscale(a,c) for a in f)
def fadd(f,g):return tuple(padd(a,b) for a,b in zip(f,g))
def zero():return tuple((F(0),) for _ in range(20))
def symmetrize(f,a):
    # a is independent of z2. Substitute z2=-a/2, and omit z20=0.
    assert len(f)==22 and not any(a[3])
    ff=[padd(f[i],pscale(pmul(f[3],a[i]),F(-1,2))) for i in range(21)]
    assert not any(ff[3]) or ff[3]==f[3] # removed below
    return tuple([ff[0]]+[ff[i+1] for i in KEEP])
def congruence(M,T):
    n=len(M);m=len(T[0]);out=[[zero() for _ in range(m)] for _ in range(m)]
    for i in range(m):
        for j in range(i,m):
            v=zero()
            for k in range(n):
                if not T[k][i]:continue
                for l in range(n):
                    if T[l][j]:v=fadd(v,fscale(M[k][l],T[k][i]*T[l][j]))
            out[i][j]=out[j][i]=v
    return out
@lru_cache(None)
def small_blocks():
    b=build();a=b['obs']['a'];sym=json.loads((load_frozen.ROOT/'inputs/symmetry.json').read_text());out=[]
    for bl in projected():
        M=[[symmetrize(f,a) for f in row] for row in bl['mat']]
        cfg=sym[bl['name']];T=[list(map(F,row)) for row in cfg['T']]
        C=congruence(M,T);r=cfg['plus'];n=len(C)
        for i in range(r):
            for j in range(r,n):assert C[i][j]==zero()
        for suffix,idx in [('plus',list(range(r))),('minus',list(range(r,n)))]:
            if not idx:continue
            out.append({'name':bl['name']+'_'+suffix,'parent':bl['name'],'T':[[row[j] for j in idx] for row in T], 'mat':[[C[i][j] for j in idx] for i in idx]})
    return out

def encode_blocks():
 return [{'name':b['name'],'parent':b['parent'],'T':[[str(x) for x in r] for r in b['T']], 'mat':[[[[str(q) for q in p] for p in f] for f in r] for r in b['mat']]} for b in small_blocks()]
