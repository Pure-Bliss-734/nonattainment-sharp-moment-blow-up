"""Exact parameterized model for the single-Q extension; standard library only.
An affine form is a length-22 list of univariate polynomials in x (low to high).
Its value is form[0] + sum(form[i+1]*z[i]). All coefficients are Fraction.
"""
from __future__ import annotations
import sys
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
sys.path.insert(0,str(Path(__file__).resolve().parent/'vendor'))
import model as md
from model import wa,check
from verify_step11 import QOP

def padd(a,b):
 out=[F(0)]*max(len(a),len(b))
 for i,c in enumerate(a):out[i]+=c
 for i,c in enumerate(b):out[i]+=c
 while len(out)>1 and out[-1]==0:out.pop()
 return tuple(out)
def pscale(a,c):return tuple(c*v for v in a)
def pmul(a,b):
 out=[F(0)]*(len(a)+len(b)-1)
 for i,c in enumerate(a):
  for j,d in enumerate(b):out[i+j]+=c*d
 return padd(out,(F(0),))
def peval(a,x):
 out=F(0)
 for c in reversed(a):out=out*x+c
 return out
def zero(n):return tuple((F(0),) for _ in range(n))
def linear(terms,n):
 out=zero(n)
 for c,v in terms:out=tuple(padd(a,pscale(b,c)) for a,b in zip(out,v))
 return out
def times_x(v):return tuple(pmul(p,(F(0),F(1))) for p in v)

def exact_right_inverse(rows):
 """For a full-row-rank Fraction matrix R, return T with R*T=I and ker R."""
 nr=len(rows);nc=len(rows[0]);eq=[]
 for i,r in enumerate(rows):
  eq.append({**{j:v for j,v in enumerate(r) if v},nc+i:F(-1)})
 rr=wa.sparse_rref(eq,nc);check(len(rr)==nr,'Bridge is not full row rank')
 T=[[F(0)]*nr for _ in range(nc)]
 for j,row in rr.items():
  for i in range(nr):T[j][i]=-row.get(nc+i,F(0))
 free=[i for i in range(nc) if i not in rr];ker=[]
 for k in free:
  v=[F(0)]*nc;v[k]=F(1)
  for j,row in rr.items():v[j]=-row.get(k,F(0))
  ker.append(v)
 for i,r in enumerate(rows):
  for j in range(nr):check(sum((r[k]*T[k][j] for k in range(nc)),F(0))==int(i==j),'Right inverse mismatch')
 for v in ker:
  check(all(sum((a*b for a,b in zip(r,v)),F(0))==0 for r in rows),'Bridge null vector mismatch')
 return T,ker

@lru_cache(None)
def build():
 m8=wa.quotient(8);free=m8['free'];nvar=len(free);nform=nvar+1;check(nvar==21,'Wrong level8 dimension')
 idx={v:i+1 for i,v in enumerate(free)}
 def old_form(poly):
  raw=wa.evaluate_trace_polynomial(poly,m8['forms']);out=list(zero(nform))
  for k,c in raw.items():
   if k<m8['n']:out[idx[k]]=padd(out[idx[k]],(c,))
   else:
    p=[F(0)]*(k-m8['n']+1);p[-1]=c;out[0]=padd(out[0],p)
  return tuple(out)
 obs={name:old_form(p) for name,p in {
  'one':{():F(1)},'x':{(0,1):F(1)},'p':{(2,3):F(-1)},'a':{(0,0,1,1):F(1)},
  'd':wa.inner(md.D_OPERATOR,md.D_OPERATOR),
  'j':{w:c/2 for w,c in wa.inner(md.C_OPERATOR,md.C_OPERATOR).items()}}.items()}
 g=md.build_global9();go=md.observables(g)
 oldwords=[m8['words'][j] for j in free]
 bridge_forms=[md.affine({w:F(1)},g) for w in oldwords]+[go[k] for k in ('x','q','u','v')]
 targets=[old_form({w:F(1)}) for w in oldwords]+[obs['x'],times_x(obs['x']),times_x(obs['p']),times_x(obs['a'])]
 R=[list(v[1:]) for v in bridge_forms];T,ker=exact_right_inverse(R)
 differences=[linear([(1,t),(-v[0],obs['one'])],nform) for t,v in zip(targets,bridge_forms)]
 gcoords=[linear(list(zip(r,differences)),nform) for r in T]
 def lifted(v):return linear([(v[0],obs['one'])]+list(zip(v[1:],gcoords)),nform)
 # Every old word agrees, not only chosen coordinates.
 for w in m8['forms']:check(lifted(md.affine({w:F(1)},g))==old_form({w:F(1)}),'Old moment bridge failed')
 for name,target in zip(('x','q','u','v'),targets[-4:]):check(lifted(go[name])==target,'Physical lift identity failed')
 nqg=md.affine(wa.inner(QOP,QOP,True),g)
 check(all(sum((a*b for a,b in zip(nqg[1:],v)),F(0))==0 for v in ker),'Q depends on an invisible coordinate')
 nq=lifted(nqg);obs['NQ']=nq
 blocks=[]
 for ground,L in [(False,4),(True,3)]:
  for charge,parity in [(0,-1),(0,1)]+[(q,None) for q in range(1,L+1)]:
   basis=wa.operator_basis(L,charge,parity)
   mat=[[old_form(wa.inner(u,v,ground)) for v in basis] for u in basis]
   check(all(mat[i][j]==mat[j][i] for i in range(len(basis)) for j in range(len(basis))),'Gram not symmetric')
   blocks.append(dict(name=('N' if ground else 'M')+str(charge)+'_'+str(parity),ground=ground,L=L,charge=charge,parity=parity,basis=basis,mat=mat))
 blocks.append(dict(name='single_Q',ground=True,L=4,charge=1,parity=None,basis=[QOP],mat=[[nq]]))
 return dict(m8=m8,g=g,old_form=old_form,obs=obs,blocks=blocks,bridge=R,right_inverse=T,lineality=ker,gcoords=gcoords,oldwords=oldwords)

def evaluate_form(form,x,z=None):
 a=[peval(p,x) for p in form]
 return a if z is None else a[0]+sum((c*v for c,v in zip(a[1:],z)),F(0))
def evaluate_matrix(block,x,z):return [[evaluate_form(v,x,z) for v in row] for row in block['mat']]

if __name__=='__main__':
 b=build();print('Exact bridge:',len(b['g']['free']),'global lifted coordinates ->',len(b['bridge']),'visible combinations +',len(b['lineality']),'invisible directions')
 print('Q polynomial affine form:')
 for i,p in enumerate(b['obs']['NQ']):
  if any(p):print(i,[str(v) for v in p])
