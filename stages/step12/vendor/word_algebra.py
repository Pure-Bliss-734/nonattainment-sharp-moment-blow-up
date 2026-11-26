"""Independent exact planar word algebra. No imports from previous bootstrap code.

Letters are integers: z=0, zb=1, pi=2, pib=3. pi=-i P_phys.
Equations are sparse polynomials in products of traces, with Fraction coefficients.
The Hamiltonian derivation is built from a cyclic gradient, not a stored EOM list.
"""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from typing import Dict, Tuple

Word = Tuple[int, ...]
WEIGHT = (1,1,2,2)
CHARGE = (1,-1,1,-1)
NAMES = ('Z','Zb','P','Pb')

def level(word): return sum(WEIGHT[l] for l in word)
def charge(word): return sum(CHARGE[l] for l in word)
def ref(word): return tuple(l ^ 1 for l in word)
def dagger(word): return (-1)**sum(l>=2 for l in word),ref(word[::-1])

@lru_cache(None)
def words(L, q=None):
    # Enumerate by length, filter weighted degree. Deliberately not the old BFS.
    return tuple(sorted((w for n in range(L+1) for w in product(range(4),repeat=n)
                         if level(w)<=L and (q is None or charge(w)==q)),
                        key=lambda w:(level(w),len(w),w)))

def clean(p): return {k:Q(v) for k,v in p.items() if v}
def trmono(*ws): return tuple(sorted(w for w in ws if w))
def add_term(out,key,value): out[key]=out.get(key,Q(0))+value

def linear_combine(*terms):
    out={}
    for a,p in terms:
        for k,v in p.items(): add_term(out,k,Q(a)*v)
    return clean(out)

def cyclic_gradient(H,letter):
    out={}
    for w,c in H.items():
        for i,l in enumerate(w):
            if l==letter: add_term(out,w[i+1:]+w[:i],c)
    return clean(out)

def hamiltonian(mass2=Q(1),coupling=Q(1)):
    # H/N² = tr[-pi pib + m² z zb + g(z zb z zb - z z zb zb)].
    return {(2,3):Q(-1),(0,1):Q(mass2),
            (0,1,0,1):Q(coupling),(0,0,1,1):-Q(coupling)}

def h_derivation(mass2=Q(1),coupling=Q(1)):
    H=hamiltonian(mass2,coupling)
    # [H,z]=-partial_{pib}H; [H,pi]=partial_{zb}H, etc.
    partners=(3,2,1,0)
    return tuple({w:((-1) if a<2 else 1)*c
                  for w,c in cyclic_gradient(H,partners[a]).items()}
                 for a in range(4))

def commutator(word,mass2=Q(1),coupling=Q(1)):
    ds=h_derivation(mass2,coupling);out={}
    for i,l in enumerate(word):
        for w,c in ds[l].items():add_term(out,word[:i]+w+word[i+1:],c)
    return clean(out)

def rotation_relation(word):
    if len(word)<2:return {}
    head,tail=word[0],word[1:];out={}
    add_term(out,trmono(word),Q(1));add_term(out,trmono(tail+(head,)),Q(-1))
    # [z,pib]=[zb,pi]=1/N, reversed pairs carry a minus sign.
    partners=(3,2,1,0);sign=1 if head<2 else -1
    for i,l in enumerate(tail):
        u,v=tail[:i],tail[i+1:]
        if l==partners[head] and charge(u)==0 and charge(v)==0:
            add_term(out,trmono(u,v),Q(-sign))
    return clean(out)

def gauge_relation(w):
    return clean({**{trmono(w+(0,3)):Q(1),trmono(w+(3,0)):Q(-1),
                       trmono(w+(1,2)):Q(1),trmono(w+(2,1)):Q(-1)},
                  trmono(w):Q(-2)})

def equations(L,mass2=Q(1),coupling=Q(1)):
    out=[]
    for w in words(L,0):
        out.append(('cyclic',w,rotation_relation(w)))
        out.append(('reflection',w,linear_combine((1,{trmono(ref(w)):Q(1)}),(-1,{trmono(w):Q(1)}))))
        sign=(-1)**sum(l>=2 for l in w)
        out.append(('reverse',w,linear_combine((sign,{trmono(w[::-1]):Q(1)}),(-1,{trmono(w):Q(1)}))))
    if L>=3:
        out.extend(('gauge',w,gauge_relation(w)) for w in words(L-3,0))
    for w in words(L-1,0):
        out.append(('eom',w,{trmono(v):c for v,c in commutator(w,mass2,coupling).items()}))
    return tuple((kind,w,p) for kind,w,p in out if p)

def sparse_rref(rows,nunknown):
    """Fraction-only sparse elimination with explicit affine RHS columns.

    Pivots are restricted to [0,nunknown). RHS-only nonzero rows are inconsistent.
    Uses row dictionaries, not SymPy/DomainMatrix's elimination routines.
    """
    piv={}
    for source in rows:
        row=clean(source.copy())
        while True:
            candidates=[i for i in row if 0<=i<nunknown]
            if not candidates:
                if row: raise ValueError(f'Inconsistent affine equations: {row}')
                break
            p=min(candidates)
            if p not in piv:
                fac=row[p];row={i:v/fac for i,v in row.items()}
                piv[p]=row;break
            fac=row[p]
            for i,v in piv[p].items():
                new=row.get(i,Q(0))-fac*v
                if new:row[i]=new
                elif i in row:del row[i]
    for p in sorted(piv,reverse=True):
        for k in sorted(i for i in piv if i<p):
            if p not in piv[k]:continue
            fac=piv[k][p]
            for i,v in piv[p].items():
                new=piv[k].get(i,Q(0))-fac*v
                if new:piv[k][i]=new
                elif i in piv[k]:del piv[k][i]
    return piv

def low_moments():
    """Derive degree<=3 neutral expectations from exact equations and x=<tr z zb>."""
    ws=[w for w in words(3,0) if w];idx={w:i for i,w in enumerate(ws)};n=len(ws)
    rows=[]
    for _,_,poly in equations(3):
        row={}
        for mon,c in poly.items():
            if len(mon)>1:raise ValueError('unexpected low nonlinear term')
            add_term(row,idx[mon[0]] if mon else n,c)
        rows.append(row)
    rows.append({idx[(0,1)]:Q(1),n+1:Q(-1)})
    piv=sparse_rref(rows,n)
    if len(piv)!=n:raise ValueError('low moments not determined by x')
    return {w:{i-n:-v for i,v in piv[j].items() if i>=n} for w,j in idx.items()}

def linearized_system(L=8,mass2=Q(1),coupling=Q(1)):
    ws=[w for w in words(L,0) if w];idx={w:i for i,w in enumerate(ws)};n=len(ws);low=low_moments()
    rows=[]
    for kind,w,poly in equations(L,mass2,coupling):
        row={}
        for mon,c in poly.items():
            if not mon:add_term(row,n,c)
            elif len(mon)==1:add_term(row,idx[mon[0]],c)
            elif len(mon)==2:
                u,v=mon
                for i,ci in low[u].items():
                    for j,cj in low[v].items():add_term(row,n+i+j,c*ci*cj)
            else:raise ValueError('unexpected product of >2 traces')
        rows.append(clean(row))
    rows.append({idx[(0,1)]:Q(1),n+1:Q(-1)})
    return ws,rows

def quotient(L=8,mass2=Q(1),coupling=Q(1)):
    ws,rows=linearized_system(L,mass2,coupling);n=len(ws);piv=sparse_rref(rows,n)
    free=[i for i in range(n) if i not in piv]
    forms={():{n:Q(1)}}
    for i,w in enumerate(ws):
        forms[w]=({j:-v for j,v in piv[i].items() if j!=i} if i in piv else {i:Q(1)})
    return {'words':ws,'rows':rows,'rank':len(piv),'free':free,'forms':forms,'n':n}

def operator_basis(L,q,parity=None):
    ws=words(L,q)
    if q!=0:return [{w:Q(1)} for w in ws]
    used=set();out=[]
    for w in ws:
        if w in used:continue
        w1=ref(w);used.update([w,w1])
        if w==w1:
            if parity==1:out.append({w:Q(1)})
        else:out.append({w:Q(1,2),w1:Q(parity,2)})
    return out

def inner(op1,op2,ground=False,mass2=Q(1),coupling=Q(1)):
    out={}
    for u,cu in op1.items():
        sign,du=dagger(u)
        for v,cv in op2.items():
            tail=commutator(v,mass2,coupling) if ground else {v:Q(1)}
            for t,c in tail.items():add_term(out,du+t,sign*cu*cv*c)
    return clean(out)

def evaluate_trace_polynomial(poly,forms):
    out={}
    for w,c in poly.items():
        for k,v in forms[w].items():add_term(out,k,c*v)
    return clean(out)

def mixed_word(w):return tuple(NAMES[l] for l in w)
def encode_old(w):return tuple(NAMES.index(l) for l in w)
