#!/usr/bin/env python3
"""Reconstruct the finite word system and its Gram matrices."""
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import argparse, hashlib, json, time
WGT=(1,1,2,2); CH=(1,-1,1,-1); PART=(3,2,1,0)
def check(ok,msg):
    if not ok: raise RuntimeError(msg)
def clean(p):return {k:v for k,v in p.items() if v}
def plus(out,k,c):
    v=out.get(k,F(0))+c
    if v:out[k]=v
    else:out.pop(k,None)
def lin(*terms):
    p={}
    for c,q in terms:
        for k,v in q.items():plus(p,k,c*v)
    return p
def wt(w):return sum(WGT[i] for i in w)
def charge(w):return sum(CH[i] for i in w)
def reflect(w):return tuple(i^1 for i in w)
def dagger(w):return (-1)**sum(i>=2 for i in w),reflect(w[::-1])
def mono(*ws):return tuple(sorted(w for w in ws if w))
@lru_cache(None)
def words(L,q=0):
    out=[]; stack=[()]
    while stack:
        w=stack.pop()
        if q is None or charge(w)==q:out.append(w)
        for a in range(4):
            if wt(w)+WGT[a]<=L:stack.append(w+(a,))
    return tuple(sorted(out,key=lambda w:(wt(w),len(w),w)))
H={(2,3):F(-1),(0,1):F(1),(0,1,0,1):F(1),(0,0,1,1):F(-1)}
def gradient(a):
    out={}
    for w,c in H.items():
        for i,l in enumerate(w):
            if a==l:plus(out,w[i+1:]+w[:i],c)
    return out
DER=tuple({w:(-1 if a<2 else 1)*c for w,c in gradient(PART[a]).items()} for a in range(4))
def deriv(w):
    out={}
    for i,a in enumerate(w):
        for v,c in DER[a].items():plus(out,w[:i]+v+w[i+1:],c)
    return out
def equations(L):
    out=[]
    for w in words(L):
        if len(w)>1:
            a,v=w[0],w[1:];p=lin((F(1),{mono(w):F(1)}),(F(-1),{mono(v+(a,)):F(1)}))
            for j,b in enumerate(v):
                if b==PART[a] and charge(v[:j])==charge(v[j+1:])==0:
                    plus(p,mono(v[:j],v[j+1:]),F(-1 if a<2 else 1))
            if p:out.append(p)
        p=lin((F(1),{mono(reflect(w)):F(1)}),(F(-1),{mono(w):F(1)}))
        if p:out.append(p)
        p=lin((F((-1)**sum(a>=2 for a in w)),{mono(w[::-1]):F(1)}),(F(-1),{mono(w):F(1)}))
        if p:out.append(p)
    if L>=3:
        for w in words(L-3):
            p={}
            for c,v in [(1,(0,3)),(-1,(3,0)),(1,(1,2)),(-1,(2,1)),(-2,())]:plus(p,mono(w+v),F(c))
            if p:out.append(p)
    for w in words(L-1):
        p={mono(v):c for v,c in deriv(w).items()}
        if p:out.append(p)
    return out

def rref(rows,n):
    piv={}
    for source in rows:
        row=dict(source)
        while row:
            choices=[j for j in row if j<n]
            if not choices:
                check(not row,'Inconsistent residual equation');break
            p=min(choices)
            if p not in piv:
                a=row[p];piv[p]={j:c/a for j,c in row.items()};break
            a=row[p]
            for j,c in piv[p].items():plus(row,j,-a*c)
    for p in sorted(piv,reverse=True):
        for k in sorted(k for k in piv if k<p):
            a=piv[k].get(p,F(0))
            if a:
                for j,c in piv[p].items():plus(piv[k],j,-a*c)
    return piv

def low_forms():
    ws=[w for w in words(3) if w];n=len(ws);idx={w:i for i,w in enumerate(ws)}
    rows=[]
    for p in equations(3):
        r={}
        for m,c in p.items():
            check(len(m)<=1,'Low nonlinear relation');plus(r,n if not m else idx[m[0]],c)
        rows.append(r)
    rows.append({idx[(0,1)]:F(1),n+1:F(-1)})
    rr=rref(rows,n);check(len(rr)==n,'Low moments not fixed')
    return {w:{j-n:-c for j,c in rr[i].items() if j>=n} for i,w in enumerate(ws)}
LOW=low_forms()
def level8():
    ws=[w for w in words(8) if w];n=len(ws);idx={w:i for i,w in enumerate(ws)};rows=[]
    for p in equations(8):
        r={}
        for m,c in p.items():
            if not m:plus(r,n,c)
            elif len(m)==1:plus(r,idx[m[0]],c)
            else:
                check(len(m)==2 and m[0] in LOW and m[1] in LOW,'Unexpected level8 product')
                for a,ca in LOW[m[0]].items():
                    for b,cb in LOW[m[1]].items():plus(r,n+a+b,c*ca*cb)
        rows.append(r)
    rows.append({idx[(0,1)]:F(1),n+1:F(-1)})
    rr=rref(rows,n);free=[j for j in range(n) if j not in rr];fidx={j:i for i,j in enumerate(free)}
    forms={():{(0,-1):F(1)}}
    for j,w in enumerate(ws):
        vals={k:-c for k,c in rr[j].items() if k!=j} if j in rr else {j:F(1)}
        forms[w]={(k-n,-1) if k>=n else (0,fidx[k]):c for k,c in vals.items()}
    return ws,rows,free,forms

def basis(L,q,par=None):
    ws=words(L,q)
    if q:return [{w:F(1)} for w in ws]
    used=set();out=[]
    for w in ws:
        if w in used:continue
        v=reflect(w);used.update([w,v])
        if w==v:
            if par==1:out.append({w:F(1)})
        else:out.append({w:F(1,2),v:F(par,2)})
    return out

def inner(u,v,ground=False):
    out={}
    for w,c in u.items():
        s,dw=dagger(w)
        for t,e in v.items():
            for a,b in (deriv(t) if ground else {t:F(1)}).items():plus(out,dw+a,s*c*e*b)
    return out

def reduce(p,forms):
    out={}
    for w,c in p.items():
        for k,v in forms[w].items():plus(out,k,c*v)
    return out

def global9():
    ws=[w for w in words(9) if w];nraw=len(ws);n=nraw+3;idx={w:i for i,w in enumerate(ws)};xx,xp,xa=range(nraw,n);rows=[]
    for p in equations(9):
        row={}
        for m,c in p.items():
            if not m:plus(row,n,c)
            elif len(m)==1:plus(row,idx[m[0]],c)
            else:
                check(len(m)==2,'More than two trace factors')
                u,v=sorted(m,key=wt);check(u in LOW,'Unknown low factor')
                lo=LOW[u]
                if 0 in lo:
                    check(len(lo)==1,'Mixed low constant');plus(row,idx[v],c*lo[0])
                else:
                    check(lo=={1:F(1)},'Unexpected low factor')
                    if v in LOW:
                        for power,val in LOW[v].items():plus(row,idx[(0,1)] if power==0 else xx,c*val)
                    elif len(v)==2:
                        check(set(v)=={2,3},'Expected kinetic factor');plus(row,xp,-c)
                    else:
                        check(wt(v)==4 and all(a<2 for a in v),'Unexpected quartic factor')
                        canon=min(v[i:]+v[:i] for i in range(4));check(canon in [(0,0,1,1),(0,1,0,1)],'Unknown fourth moment')
                        plus(row,xa,c)
                        if canon==(0,1,0,1):plus(row,xp,c/2);plus(row,xx,-c/2)
        rows.append(row)
    rr=rref(rows,n);free=[i for i in range(n) if i not in rr]
    forms={():{n:F(1)}}
    vals={j:({k:-c for k,c in rr[j].items() if k!=j} if j in rr else {j:F(1)}) for j in range(n)}
    forms.update({w:vals[j] for j,w in enumerate(ws)})
    # Recheck every row as an exact identity.
    for row in rows:
        p={}
        for j,c in row.items():
            for k,v in ({n:F(1)} if j==n else vals[j]).items():plus(p,k,c*v)
        check(not p,'Failed generated level9 identity')
    return ws,n,rows,free,forms,vals

QOP={(0,1,0):F(3),(0,1,2):F(-1),(2,1,0):F(-1)}
GAMMA={(0,3):F(1),(3,0):F(-1),(1,2):F(1),(2,1):F(-1),():F(-2)}
def op(s):return {tuple(map(int,s)):F(1)}
def pm(s,sign):
    w=tuple(map(int,s));v=reflect(w);return lin((F(1,2),{w:F(1)}),(F(sign,2),{v:F(1)}))
def emm(s):return pm(s,-1)
def epp(s):return pm(s,1)
def mulop(u,v):
    out={}
    for w,c in u.items():
        for t,d in v.items():plus(out,w+t,c*d)
    return out

def main(outdir):
    t0=time.monotonic(); outdir.mkdir(parents=True,exist_ok=True)
    ws,rows,free,forms=level8(); print('level8',len(ws),len(free),flush=True)
    expected=['3210','3300','311000','3232','3322','312100','313000','320011','320101','321001','321010','321100','330100','331000','10101010','11001100','11010010','11010100','11100100','11101000','11110000']
    check([''.join(map(str,ws[i])) for i in free]==expected,'Free-word table mismatch')
    check(len(ws)==614 and len(free)==21,'Level8 counts')
    ws9,n9,rows9,free9,forms9,vals9=global9(); print('level9',len(ws9),len(rows9),len(free9),flush=True)
    check((len(ws9),n9,len(rows9),len(free9))==(1574,1577,5417,35),'Level9 counts')
    def R(p):return reduce(p,forms)
    obs={k:R(p) for k,p in {'p':lin((F(-1),op('23'))),'a':op('0011'),'b':op('0101')}.items()}
    check(not lin((F(1),obs['p']),(F(-1),{(1,-1):F(1)}),(F(-2),obs['b']),(F(2),obs['a'])),'Virial coefficient')
    ptab=[0,0,0,2,-2,60,4,4,28,-36,32,-88,-12,8,8,32,70,-116,36,-28]
    atab=[0,0,0,0,1,-10,-2,-2,-2,2,0,12,6,-4,0,-4,-9,14,-4,2]
    for k,constant,coeffs in [('p',9,ptab),('a',-4,atab)]:
        check(obs[k].get((1,-1))==constant,'Wrong '+k+' constant')
        check(all(obs[k].get((0,j),F(0))==c for j,c in enumerate(coeffs)),'Wrong '+k+' table')
        check(sum(obs[k].get((0,j),F(0)) for j in range(14,21))==0,k+' ray invariance')
    # Common null columns, tested from the full word definition.
    null_checks=0
    for L,q,par,ground,nulls in [(4,0,1,False,[GAMMA]),(4,1,None,False,[mulop(op('0'),GAMMA),lin((F(1),mulop(GAMMA,op('0'))),(F(1),op('0')))]),(3,0,1,True,[{():F(1)},GAMMA])]:
        for u in basis(L,q,par):
            for v in nulls:
                check(not R(inner(u,v,ground)),'Nonzero common-null column');null_checks+=1
    dims=[];symmetric=0
    for ground,L in [(False,4),(True,3)]:
        for q,par in [(0,-1),(0,1)]+[(q,None) for q in range(1,L+1)]:
            bs=basis(L,q,par);dims.append(len(bs))
            A=[[R(inner(u,v,ground)) for v in bs] for u in bs]
            for i in range(len(A)):
                for j in range(len(A)):
                    check(A[i][j]==A[j][i],'Non-Hermitian Gram');symmetric+=1
    check(dims==[7,8,14,8,4,1,3,4,5,3,1],'Raw block layout')
    # Recovery bridge rank, with all 35 lifted free coordinates.
    targets=[forms9[w] for w in [ws[j] for j in free]]+[forms9[(0,1)],vals9[n9-3],vals9[n9-2],vals9[n9-1]]
    idxf={j:i for i,j in enumerate(free9)}
    bridge=[{idxf[j]:c for j,c in row.items() if j!=n9} for row in targets]
    # Homogeneous rank extraction: no affine RHS columns.
    rank=len(rref(bridge,35));check(rank==25,'Recovery bridge rank')
    # Q coefficient row factors through this bridge: solve augmented span by rank.
    nq=reduce(inner(QOP,QOP,True),forms9)
    nqrow={idxf[j]:c for j,c in nq.items() if j!=n9}
    check(len(rref(bridge+[nqrow],35))==rank,'Q depends on invisible coordinates')
    # Reconstruct global free coordinates as rational forms in old phi,x,x²,xp,xa.
    nr=len(bridge);br=[]
    for i,row in enumerate(bridge):br.append({**row,35+i:F(-1)})
    brr=rref(br,35)
    target_forms=[dict(forms[w]) for w in [ws[j] for j in free]]+[{(1,-1):F(1)},{(2,-1):F(1)}, {(p+1,j):c for (p,j),c in obs['p'].items()}, {(p+1,j):c for (p,j),c in obs['a'].items()}]
    for i,tgt in enumerate(targets):
        plus(target_forms[i],(0,-1),-tgt.get(n9,F(0)))
    rec=[{} for _ in range(35)]
    for j,row in brr.items():
        rec[j]=lin(*[(-c,target_forms[k-35]) for k,c in row.items() if k>=35])
    def lift(form):return lin((form.get(n9,F(0)),{(0,-1):F(1)}),*[(c,rec[idxf[j]]) for j,c in form.items() if j!=n9])
    for w in forms:check(lift(forms9[w])==forms[w],'Old-word recovery inconsistency')
    qform=lift(nq)
    Lform=lin((F(7,2),obs['p']),(F(3),obs['a']))
    etab=[1,0,0,12,0,254,0,0,151,-194,198,-388,10,0,36,F(291,2),F(653,2),-527,170,-151]
    Qexpected=lin((F(12),{(2,-1):F(1)}),(F(-13),{(1,-1):F(1)}),(F(1),{(0,j):F(c) for j,c in enumerate(etab) if c}),(F(1),{(p+1,j):c for (p,j),c in Lform.items()}))
    # Move from phi to z,t; t coefficient must vanish.
    qray=lin(*[(F(1),{(p,-1):c}) for (p,j),c in qform.items() if j>=14])
    check(not qray,'Q varies along the eighth ray')
    
    difference=lin((F(1),{k:c for k,c in qform.items() if k[1]!=20}),(F(-1),{k:c for k,c in Qexpected.items() if k[1]!=20}))
    if difference:print('Q difference',difference,flush=True)
    check(not difference,'Q coordinate table mismatch')
    def projected_form(p):return {k:c for k,c in R(p).items() if k[1]!=20}
    W0=epp('0011');W2=op('0001')
    for wi,j,coef in [(W0,15,F(1,2)),(W2,19,F(1))]:
        p=R(inner(wi,wi));check(sum(c for (power,k),c in p.items() if k>=14)==1,'Completion t coefficient')
        check(projected_form(inner(wi,wi))=={(0,j):coef},'Completion shift mismatch')
    # Exact compactness identity with the eighteen printed operator tests.
    atoms=[
      (F(1,5),False,emm('0110')),
      (F(3,5),False,lin((F(1),emm('03')),(F(-1),emm('23')),(F(1),emm('0110')))),
      (F(671,1020),False,lin((F(1),op('001')),(F(-1),op('100')))),
      (F(4,17),False,lin((F(1),op('0010')),(F(-1),op('0100')))),
      (F(122,255),True,lin((F(1),epp('01')),(F(-1),epp('03')))),
      (F(273,680),True,lin((F(1),op('0')),(F(-1),op('2')))),
      (F(631,680),True,lin((F(1),op('00')),(F(-1),op('02')))),
      (F(193,2040),True,lin((F(1),op('00')),(F(-1),op('20')))),
      (F(1,5),False,lin((F(1),emm('03')),(F(-1),emm('23')),(F(1),emm('0011')),(F(-1),emm('0101')))),
      (F(3,5),False,lin((F(1),epp('01')),(F(-1),epp('03')),(F(-1),epp('0011')),(F(1),epp('0101')))),
      (F(11,30),False,lin((F(1),op('0')),(F(-1),op('003')),(F(-1),op('012')),(F(1),op('201')))),
      (F(3,5),False,lin((F(1),op('2')),(F(-1),op('003')),(F(-1),op('012')),(F(1),op('021')))),
      (F(1,5),False,lin((F(1),op('001')),(F(-1),op('100')),(F(-1),op('003')),(F(1),op('102')))),
      (F(1,30),False,lin((F(1),op('010')),(F(-1),op('003')),(F(-1),op('012')),(F(1),op('201')))),
      (F(1,5),False,lin((F(1),op('012')),(F(-1),op('021')),(F(-1),op('120')),(F(1),op('210')))),
      (F(14,85),False,lin((F(1),op('02')),(F(-1),op('20')),(F(-1),op('0010')),(F(1),op('0100')))),
      (F(587,2040),True,lin((F(1),op('0')),(F(-1),op('2')),(F(1),op('010')))),
      (F(1,6),True,QOP)]
    pp={}
    for c,g,u in atoms:
        p=lift(reduce(inner(u,u,g),forms9))
        pp=lin((F(1),pp),(c,p))
    x={(1,-1):F(1)};one={(0,-1):F(1)}
    scalars=[(F(324883,54400),lin((F(2,5),one),(F(-1),x))),
      (F(3853,1360),{(0,-1):F(-3,20),(1,-1):F(31,40),(2,-1):F(-1)})]
    for c,p in scalars:pp=lin((F(1),pp),(c,p))
    for name,c in [('p',F(7,12)),('a',F(1,2))]:
        p=obs[name];xp1={(a+1,j):v for (a,j),v in p.items()}
        pp=lin((F(1),pp),(c,lin((F(2,5),p),(F(-1),xp1))))
    D=lin((F(1),op('210')),(F(1),op('300')))
    expected=lin((F(7751,6375),one),(F(-1),R(inner(D,D))))
    check(pp==expected,'Printed completed compactness identity is not exact')
    # The projected bases printed in Appendix A.
    pbases=[
      [emm(v) for v in ['01','03','21','23','0011','0101','0110']],
      [{():F(1)},epp('01'),epp('03'),epp('23'),lin((F(1),epp('0101')),(F(-1),epp('0011'))),lin((F(1),epp('0110')),(F(-1),epp('0011')))],
      [op(v) for v in ['0','2','001','010','100','003','012','021','102','120','201','210']],
      [op(v) for v in ['00','02','20','22']]+[lin((F(1),op(v)),(F(-1),op('0001'))) for v in ['0010','0100','1000']],
      [op(v) for v in ['000','002','020','200']],
      [emm(v) for v in ['01','03','21']],
      [epp(v) for v in ['01','03']],
      [op(v) for v in ['0','2','001','010','100']],
      [op(v) for v in ['00','02','20']],
      [op('000')], [QOP]]
    names=['M0_-1','M0_1','M1_None','M2_None','M3_None','N0_-1','N0_1','N1_None','N2_None','N3_None','single_Q']
    matrices=[]
    for ib,bs in enumerate(pbases):
        A=[]
        for u in bs:
            ar=[]
            for v in bs:
                form=lift(reduce(inner(u,v,ib>=5),forms9))
                for power in (0,1,2):check(sum(form.get((power,k),F(0)) for k in range(14,21))==0,'Projected block retains ray dependence')
                ar.append({k:c for k,c in form.items() if k[1]!=20})
            A.append(ar)
        for i in range(len(A)):
            for j in range(len(A)):check(A[i][j]==A[j][i],'Projected symmetry mismatch')
        matrices.append(A)
    check([len(A) for A in matrices]==[7,6,12,7,4,3,2,5,3,1,1],'Projected layout')
    def tangent(form):
        check(all(power==0 for power,k in form if k>=0),'Unexpected x-dependent annihilation row')
        return {k:c for (power,k),c in form.items() if k>=0}
    ann=[tangent(obs['p']),tangent(obs['a']),tangent(R(inner(D,D)))]
    parents=[0,0,2,3,6,7,8,8,0,1,2,2,2,2,2,3,7]
    for (c,g,u),parent in zip(atoms[:-1],parents):
        for v in pbases[parent]:ann.append(tangent(R(inner(v,u,g))))
    for v in pbases[2]:ann.append(tangent(R(inner(v,D))))
    arr=rref(ann,21)
    check(len(ann)==146 and len(arr)==20,'Projected recession rank')
    free_r=[j for j in range(21) if j not in arr];check(len(free_r)==1,'Unexpected ray dimension')
    ray=[F(0)]*21;ray[free_r[0]]=1
    for j,row in arr.items():ray[j]=-row.get(free_r[0],F(0))
    check(ray==[F(0)]*14+[F(1)]*7,'Wrong recession ray')
    # Export the exact reconstructed projected matrices, plus complete columns.
    exported=[]
    def encode_form(f):return [[[power,k],str(c)] for (power,k),c in sorted(f.items())]
    for name,A in zip(names,matrices):exported.append({'name':name,'matrix':[[encode_form(f) for f in row] for row in A]})
    cps={}
    for name,idx,wi in [('M0_1',1,W0),('M2_None',3,W2)]:
        cross=[{k:c for k,c in R(inner(u,wi)).items() if k[1]!=20} for u in pbases[idx]]
        shift=projected_form(inner(wi,wi))
        cps[name]={'cross':[encode_form(f) for f in cross],'shift':encode_form(shift)}
    (outdir/'reconstructed_model.json').write_text(json.dumps({'representation':'sum coefficient * x^power * (1 if coordinate=-1 else z[coordinate])','projected':exported,'completion':cps},indent=2)+'\n')

    # Scalar reference cases for the asymptotic identities (finite checks only).
    scalar_count=0
    for r in [F(1,10),F(1,100),F(1,1000)]:
        for C,q,P,e in [(F(2),F(3),F(5),F(7)),(F(7,3),F(2,5),F(11,2),F(-3,7))]:
            Q0=q*q/P;Q1=q*e/P;Q2=e*e/P
            lhs=(q+r*e)**2/(r*q*q/C+r*r*P)
            rhs=C/r+(2*C*Q1+Q0*Q2-Q1*Q1-C*C+C*r*Q2)/(Q0+C*r)
            check(lhs==rhs,'Rank-one inverse identity');scalar_count+=1
    for A,B,R0 in [(211,223,F(1,100)),(754,1072,F(1,40))]:
        d=F(1,10**9);check(R0-B*d>0 and (R0-B*d)**2-A*A*d>0,'Neighborhood strict inequality')
    report={
        'status': 'passed',
        'Q_annihilates_invisible_directions': True,
        'affine_coordinate_table_exact': True,
        'bridge_rank': rank,
        'common_null_column_entries_checked': null_checks,
        'compactness_identity_exact': True,
        'completion_shifts_exact': True,
        'invisible_dimension': 35-rank,
        'level8_free_coordinates': len(free),
        'level8_nonempty_words': len(ws),
        'level9_free_lifted_coordinates': len(free9),
        'level9_generated_relations_verified': len(rows9),
        'level9_lifted_variables': n9,
        'level9_nonempty_words': len(ws9),
        'neighborhood_rational_comparisons': 2,
        'printed_compactness_atom_count': len(atoms),
        'projected_matrix_entries_checked': sum(len(A)**2 for A in matrices),
        'rank_one_inverse_reference_checks': scalar_count,
        'raw_gram_dimensions_without_Q': dims,
        'raw_gram_symmetry_entries_checked': symmetric,
        'recession_annihilation_rows': len(ann),
        'recession_rank': len(arr),
    }
    (outdir/'exact_structure_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,default=Path('audit_results'));a=ap.parse_args();main(a.output_dir)
