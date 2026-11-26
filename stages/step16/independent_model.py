"""Sparse-polynomial reconstruction of Step-15 symmetry and block pencils.

Monomials are sorted tuples of variable indices. Polynomial multiplication and
substitution, matrix congruence and exact rank are implemented here from scratch.
The input is the earlier word quotient, deliberately a shared premise.
"""
from fractions import Fraction as F
from pathlib import Path
import json
from decimal_interval import Box,ZERO,ONE
from second_order_ad import Jet
ROOT=Path(__file__).resolve().parent
KEEP=tuple(i for i in range(20) if i!=2)

def pc(c):return {} if not c else {():F(c)}
def pv(i):return {(i,):F(1)}
def plus(a,b):
    c=a.copy()
    for k,v in b.items():
        c[k]=c.get(k,F(0))+v
        if not c[k]:del c[k]
    return c
def times(a,b):
    c={}
    for p,v in a.items():
        for q,w in b.items():
            k=tuple(sorted(p+q));c[k]=c.get(k,F(0))+v*w
    return {k:v for k,v in c.items() if v}
def scaling(a,s):return {k:v*F(s) for k,v in a.items() if v*F(s)}
def from_form(f,nvariables=21):
    # Raw form index k>0 corresponds to z_(k-1), at variable index k.
    # z20 is set to zero in the visible section and thus omitted.
    p={}
    for k,poly in enumerate(f):
        if k>=nvariables:continue
        for d,c in enumerate(poly):
            c=F(c)
            if not c:continue
            m=(0,)*d+((k,) if k else ())
            p=plus(p,{tuple(sorted(m)):c})
    return p

def substitute(p,variables):
    out={}
    for mon,c in p.items():
        term=pc(c)
        for i in mon:term=times(term,variables[i])
        out=plus(out,term)
    return out

def derivative(p,i):
    o={}
    for mon,c in p.items():
        m=list(mon);n=m.count(i)
        if n:
            m.remove(i);o=plus(o,{tuple(m):n*c})
    return o

def evaluate(p,coords):
    # Works with intervals, jets or Fractions; no derivative formulas here.
    out=coords[0]*0
    for mon,c in p.items():
        if not mon:out=out+c;continue
        term=coords[mon[0]]*c
        for i in mon[1:]:term=term*coords[i]
        out=out+term
    return out

def congruence(m,t):
    n=len(m);r=len(t[0]);out=[[{} for _ in range(r)] for __ in range(r)]
    for i in range(r):
        for j in range(r):
            q={}
            for k in range(n):
                if not t[k][i]:continue
                for l in range(n):
                    if t[l][j]:q=plus(q,scaling(m[k][l],t[k][i]*t[l][j]))
            out[i][j]=q
    return out

def rational_rank(a):
    a=[list(map(F,r)) for r in a];rank=0
    if not a:return 0
    for j in range(len(a[0])):
        k=next((i for i in range(rank,len(a)) if a[i][j]),None)
        if k is None:continue
        a[rank],a[k]=a[k],a[rank];q=a[rank][j];a[rank]=[v/q for v in a[rank]]
        for i in range(rank+1,len(a)):
            q=a[i][j]
            if q:a[i]=[v-q*w for v,w in zip(a[i],a[rank])]
        rank+=1
        if rank==len(a):break
    return rank

def encode_poly(p):return [{'monomial':list(k),'coefficient':str(v)} for k,v in sorted(p.items())]
class Model:
    def __init__(self,data=None,cfg=None):
        self.data=data if data is not None else json.loads((ROOT/'inputs/pencils.json').read_text())
        cfg=cfg if cfg is not None else json.loads((ROOT/'inputs/symmetry.json').read_text())
        self.raw=[];self.small=[];self.symmetry=[]
        a=from_form(self.data['observables']['a'])
        if any(3 in mon for mon in a):raise AssertionError('a depends on reflected coordinate')
        # New coordinates: x, 19 retained z's, w. z2=w-a/2.
        replacements={0:pv(0)}
        for j,k in enumerate(KEEP):replacements[k+1]=pv(j+1)
        replacements[3]=pc(0)
        newa=substitute(a,replacements)
        replacements[3]=plus(pv(20),scaling(newa,F(-1,2)))
        self.to_symmetric=replacements
        entries=terms=0
        for b in self.data['raw_blocks']:
            name=b['name'];m=[[from_form(f) for f in row] for row in b['mat']];n=len(m)
            if any(m[i][j]!=m[j][i] for i in range(n) for j in range(n)):raise AssertionError('raw symmetry lost')
            self.raw.append({'name':name,'mat':m})
            mnew=[[substitute(f,replacements) for f in row] for row in m]
            sc=cfg[name];s=[list(map(F,row)) for row in sc['S']];t=[list(map(F,row)) for row in sc['T']]
            if any(sum(s[i][k]*s[k][j] for k in range(n))!=int(i==j) for i in range(n) for j in range(n)):raise AssertionError('involution failed')
            if rational_rank(t)!=n:raise AssertionError('eigenbasis singular')
            r=sc['plus']
            if r+sc['minus']!=n:raise AssertionError('split size')
            if any(sum(s[i][k]*t[k][j] for k in range(n))!=(1 if j<r else -1)*t[i][j] for i in range(n) for j in range(n)):raise AssertionError('wrong eigenspaces')
            reverse=[ [{mon:(-c if mon.count(20)%2 else c) for mon,c in f.items()} for f in row] for row in mnew]
            if congruence(mnew,s)!=reverse:raise AssertionError('full polynomial reflection identity failed')
            even=[[{mon:c for mon,c in f.items() if 20 not in mon} for f in row] for row in mnew]
            cm=congruence(even,t)
            if any(cm[i][j] for i in range(r) for j in range(r,n)):raise AssertionError('nonzero mixed parity block')
            # The odd derivative has only mixed-parity entries, hence every
            # lifted block-diagonal dual has exactly zero w-stationarity.
            odd=congruence([[derivative(f,20) for f in row] for row in mnew],t)
            if any(odd[i][j] for inds in (range(r),range(r,n)) for i in inds for j in inds):raise AssertionError('missing full-coordinate stationarity')
            for suffix,inds in [('plus',list(range(r))),('minus',list(range(r,n)))]:
                if not inds:continue
                sm=[[cm[i][j] for j in inds] for i in inds]
                if any(len(mon)>2 for row in sm for f in row for mon in f):raise AssertionError('not quadratic in joint parameters')
                self.small.append({'name':name+'_'+suffix,'parent':name,'T':[[row[j] for j in inds] for row in t],'mat':sm})
            entries+=n*n;terms+=sum(len(f) for row in mnew for f in row)
            self.symmetry.append({'name':name,'dimension':n,'plus':r,'minus':n-r,'full_polynomial_involution':True,'odd_dual_stationarity_exact':True})
        old=self.data['legacy_small_blocks_for_comparison']
        if len(old)!=len(self.small):raise AssertionError('small block count differs')
        compared=0
        for oldb,b in zip(old,self.small):
            if oldb['name']!=b['name']:raise AssertionError('small block order differs')
            for i,row in enumerate(b['mat']):
                for j,p in enumerate(row):
                    # In old reduced forms: k>0 is y_k, not raw z.
                    if from_form(oldb['mat'][i][j],20)!=p:raise AssertionError('independent reduction disagrees')
                    compared+=1
        self.report={'original_blocks':len(self.raw),'small_blocks':len(self.small),'raw_entry_identities':entries,'new_polynomial_terms':terms,'reduced_matrix_entries_compared':compared,'reduction_by_sparse_polynomial_substitution':True,'checks':self.symmetry}
    def matrix(self,i,y):return [[evaluate(p,y) for p in row] for row in self.small[i]['mat']]
    def raw_coordinates(self,y):
        return [evaluate(self.to_symmetric[k],list(y)+[ZERO]) for k in range(21)]
    def observable(self,name,y):return evaluate(from_form(self.data['observables'][name]),self.raw_coordinates(y))
    def raw_matrix(self,i,coordinates):return [[evaluate(p,coordinates) for p in row] for row in self.raw[i]['mat']]
    def derivatives(self,i,y):
        b=self.small[i]
        return [[[evaluate(derivative(p,k),y) for p in row] for row in b['mat']] for k in range(20)]
