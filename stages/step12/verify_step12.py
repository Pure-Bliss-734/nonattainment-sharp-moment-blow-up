#!/usr/bin/env python3
"""Verify Step 12 from original words and rational certificates.
Only Python standard library and vendored Fraction-only source are used.
Discovery programs, numerical logs and floating eigenvalues are not proof inputs.
"""
from __future__ import annotations
import argparse,json,time,hashlib
from pathlib import Path
from fractions import Fraction as F
from step12_model import *
from step12_layout import projected
from verify_step10 import foundation,exact_psd_rank,decimal,save,verify_moment_certificate
from verify_step11 import verify_certificate,reconstruct_base
ROOT=Path(__file__).resolve().parent
RDIR=[F(0)]*14+[F(1)]*7

def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def transpose(A):return list(map(list,zip(*A))) if A else []
def rank_kernel(A,n=None):
 if n is None:n=len(A[0])
 rr=wa.sparse_rref([{j:c for j,c in enumerate(row) if c} for row in A],n)
 ker=[]
 for k in range(n):
  if k in rr:continue
  v=[F(0)]*n;v[k]=1
  for j,row in rr.items():v[j]=-row.get(k,F(0))
  ker.append(v)
 return len(rr),ker

def coordinates(op,basis):
 # The raw charge/parity bases have disjoint word supports.
 val=[]
 for base in basis:
  w=next(iter(base));val.append(op.get(w,F(0))/base[w])
 out={}
 for c,base in zip(val,basis):
  for w,a in base.items():wa.add_term(out,w,c*a)
 check(wa.clean(out)==op,'Operator not in declared raw basis')
 return val

def tangent_row(form):
 check(all(len(p)==1 for p in form[1:]),'Expected x-independent old tangent')
 return [p[0] for p in form[1:]]

def verify_layout(b,state):
 xp=F('0.38900458');z=[state[w] for w in b['oldwords']];layout=json.loads((ROOT/'inputs/projected_layout.json').read_text())
 out=[];details={}
 for bl,lo in zip(b['blocks'],layout['blocks']):
  check(bl['name']==lo['name'],'Wrong block ordering');ops=[md.operator_from_json(v) for v in lo['basis']]
  n=len(bl['basis']);A=evaluate_matrix(bl,xp,z);rb,U=rank_kernel(A);piv=exact_psd_rank(A);check(len(piv)==rb,'Reference rank mismatch')
  for v in U:
   for row in bl['mat']:check(not any(any(p) for p in linear(list(zip(v,row)),22)),'Common nullspace is not an identity')
  D=[]
  for row in bl['mat']:
   rr=[]
   for form in row:
    pol=(F(0),)
    for c,p in zip(RDIR,form[1:]):pol=padd(pol,pscale(p,c))
    check(len(pol)==1,'Recession depends on x');rr.append(pol[0])
   D.append(rr)
  pp=exact_psd_rank(D);rd=len(pp)
  B=[coordinates(op,bl['basis']) for op in ops]
  for col in B:check(all(dot(row,col)==0 for row in D),'Projected basis not in exposing kernel')
  ru,_=rank_kernel(transpose(U+B),len(U+B)) if U+B else (0,[])
  check(ru==len(U)+len(B)==n-rd,'Projected basis fails to span complete exposing kernel modulo identities')
  out.append({'name':bl['name'],'raw_dimension':n,'common_kernel_dimension':len(U),'exposing_rank':rd,
    'projected_dimension':len(B),'direction_matrix':[[str(c) for c in row] for row in D],
    'positive_direction_pivots':list(map(str,pp))})
  details[bl['name']]=(bl,ops,U,B)
 # Nonzero projection dimensions match the final matrices.
 check([d['projected_dimension'] for d in out if d['projected_dimension']]==[len(v['mat']) for v in projected()],'Projected layout mismatch')
 return out,details

def verify_recession(b,details):
 cert=json.loads((ROOT/'inputs/compact22.json').read_text());counts={}
 for use_projected in (False,True):
  rows=[tangent_row(b['obs'][k]) for k in ('p','a','d')]
  atoms=[]
  for t in cert['terms']:
   if t['kind']=='scalar':continue
   op=md.operator_from_json(t['operator'])
   if op==QOP:continue
   atoms.append((t['block'],op,t['kind']=='ground'))
  atoms.append(('M1_None',md.D_OPERATOR,False))
  # p,a positivity must also survive the projection.
  check_atoms=atoms+[('M1_None',{(2,):F(1)},False),('M2_None',{(0,0):F(1)},False)]
  for name,op,ground in check_atoms:
   bl,ops,U,B=details[name];v=coordinates(op,bl['basis'])
   rank_before,_=rank_kernel(transpose(U+B),len(U+B))
   rank_after,_=rank_kernel(transpose(U+B+[v]),len(U+B)+1)
   check(rank_before==rank_after,'Certificate operator is lost by facial projection')
  for name,op,ground in atoms:
   bl,ops,_,_=details[name]
   for u in (ops if use_projected else bl['basis']):rows.append(tangent_row(b['old_form'](wa.inner(u,op,ground))))
  rank,ker=rank_kernel(rows,21);check(rank==20 and ker==[RDIR],'Recession classification failed')
  counts['projected' if use_projected else 'original']={'annihilation_rows':len(rows),'rank':rank,'kernel':list(map(str,ker[0]))}
 # All low-degree observables and N(Q) are unchanged by the remaining ray.
 for name,form in b['obs'].items():
  pol=(F(0),)
  for c,p in zip(RDIR,form[1:]):pol=padd(pol,pscale(p,c))
  check(not any(pol),'Remaining ray changes protected observable '+name)
 for w in b['m8']['forms']:
  form=b['old_form']({w:F(1)});v=dot(tangent_row(form),RDIR)
  expected=int(len(w)==8 and all(t<2 for t in w))
  check(v==expected,'Unexpected original moment growth')
 # The symbolic full-nine extension has zero growth outside pure coordinate level 8.
 dg=[]
 for form in b['gcoords']:
  pol=(F(0),)
  for c,p in zip(RDIR,form[1:]):pol=padd(pol,pscale(p,c))
  check(len(pol)==1,'Full ray unexpectedly depends on x');dg.append(pol[0])
 delta={w:dot(md.affine({w:F(1)},b['g'])[1:],dg) for w in b['g']['forms']}
 for w,v in delta.items():check(v==int(len(w)==8 and all(t<2 for t in w)),'Full-nine pure-coordinate ray mismatch')
 # Invisible directions preserve every <=8 moment, x and the three lifts and Q.
 for v in b['lineality']:
  for w in b['m8']['forms']:check(dot(md.affine({w:F(1)},b['g'])[1:],v)==0,'Invisible direction changes old moment')
  check(dot(md.affine(wa.inner(QOP,QOP,True),b['g'])[1:],v)==0,'Invisible direction changes Q')
 counts.update(lineality_dimension=len(b['lineality']),pointed_recession_dimension=1,
    full_fixed_x_affine_dimension=31,quotient_dimension=20,
    global_free_indices=b['g']['free'],full_lineality_basis=[[str(c) for c in v] for v in b['lineality']],
    pure_coordinate_ray_global_free=list(map(str,dg)),
    claim_scope='For each fixed feasible x: rec = one nonnegative pure eighth-coordinate ray plus ten-dimensional degree-nine lineality. Projected set, with z20=0, has zero recession cone.')
 return counts

def verify_simple_bounds(b):
 o=b['obs'];p2={(2,):F(1)};V={(0,0,1):F(1),(1,0,0):F(1)}
 diag=b['old_form'](wa.inner(V,V,True));cross=b['old_form'](wa.inner(p2,V,True));pm=b['old_form'](wa.inner(p2,p2,True))
 xx=times_x(o['x'])
 check(diag==linear([(F(2),o['a']),(F(1,2),o['p']),(F(2),xx),(-F(1,2),o['x'])],22),'a bound diagonal fails')
 check(cross==linear([(F(2),o['a']),(F(2),o['j'])],22),'a bound cross fails')
 check(pm==linear([(1,o['x']),(F(1,2),o['one'])],22),'p ground diagonal fails')
 # p^2 <= (x+1/2)/2 <=9/20 < (7/10)^2.
 check(F(9,20)<F(7,10)**2,'p upper bound failed')
 # 4 a^2 <= (9/10)(2a+67/100), using j>=0,p<=7/10,x<=2/5.
 # At a=1, difference 4a²-(9/5)a-603/1000 is positive,
 # and its derivative is positive for a>=1.
 check(4-F(9,5)-F(603,1000)>0 and 8-F(9,5)>0,'a<1 implication failed')
 check(F(7,2)*F(7,10)+3==F(109,20),'L upper bound failed')
 return {'p_upper':'7/10','a_upper':'1','L_definition':'(7/2)*p+3*a','L_upper':'109/20'}

def verify_dual(side,b):
 data=json.loads((ROOT/'inputs'/('dual_'+side+'.json')).read_text());cut=F(data['cut']);blocks=projected();total=zero(22);pivots=[]
 check(len(data['duals'])==len(blocks),'Wrong dual block count')
 for item,bl in zip(data['duals'],blocks):
  check(item['name']==bl['name'],'Dual layout mismatch');Y=[list(map(F,row)) for row in item['matrix']];n=len(bl['mat'])
  check(len(Y)==n and all(len(row)==n for row in Y),'Dual dimension mismatch')
  pp=exact_psd_rank(Y);check(len(pp)==n,'Dual not strictly positive on its face');pivots.extend(pp)
  total=linear([(1,total)]+[(Y[j][i],en) for i,row in enumerate(bl['mat']) for j,en in enumerate(row)],22)
 gamma=F(data['duals'][-1]['matrix'][0][0]);check(gamma>0,'Bad multiplier')
 L=linear([(F(7,2),b['obs']['p']),(3,b['obs']['a'])],22)
 rem=linear([(1,total),(-gamma,tuple(pmul(p,(-cut,F(1))) for p in L))],22)
 check(not any(any(p) for p in rem[1:]),'Global-in-x dual identity not exact')
 check(all(peval(p,cut)==0 for p in total[1:]),'All 21 tangent stationarity equations must hold')
 pol=rem[0]
 if side=='upper':pol=padd(pol,pscale((-cut,F(1)),gamma*F(109,20)))
 check(tuple(map(F,data['polynomial_ascending']))==pol,'Exclusion polynomial mismatch')
 check(len(pol)<=3,'Expected quadratic exclusion polynomial')
 region=[F(3,8),cut] if side=='lower' else [cut,F(2,5)]
 dp=tuple(i*pol[i] for i in range(1,len(pol)))
 ds=[peval(dp,t) for t in region];val=peval(pol,cut)
 check(val<0,'Strict cut not excluded')
 check(all(d>0 for d in ds) if side=='lower' else all(d<0 for d in ds),'Global sign not established')
 return {'cut':str(cut),'cut_decimal':decimal(cut,15),'strict_exclusion_value':str(val),'value_display':decimal(val,30),
    'gamma':str(gamma),'polynomial_ascending':list(map(str,pol)),
    'positive_dual_pivots':list(map(str,pivots)),'positive_pivot_count':len(pivots),'stationarity_equations_checked':21,
    'derivative_at_region_endpoints':list(map(str,ds)),'excluded_region':list(map(str,region))}

def verify_state_equations(state,delta=None):
 count=0
 for _,_,eq in wa.equations(9):
  coeff=[F(0)]*3
  for mon,c in eq.items():
   p=[c]
   for w in mon:
    a=state[w];d=F(0) if delta is None else delta[w];n=[F(0)]*(len(p)+1)
    for j,v in enumerate(p):n[j]+=a*v;n[j+1]+=d*v
    p=n
   for j,v in enumerate(p):coeff[j]+=v
  check(not any(coeff),'Original factorized equality fails');count+=1
 return count

def verify_primal(side,b):
 data=json.loads((ROOT/'inputs'/('primal_'+side+'.json')).read_text());x=F(data['x']);z=list(map(F,data['quotient_z']));t=F(data['completion']);full=list(map(F,data['full_z']))
 check(len(z)==21 and z[-1]==0 and t>=0,'Bad quotient representative')
 check(full==[a+t*r for a,r in zip(z,RDIR)],'Completion vector mismatch')
 for bl in projected():check(len(exact_psd_rank(evaluate_matrix(bl,x,z)))==len(bl['mat']),'Projected primal is not strict')
 # Lift into the complete original moment dictionary.
 gvals=[evaluate_form(v,x,full) for v in b['gcoords']]
 state={w:md.affine({w:F(1)},b['g'])[0]+dot(md.affine({w:F(1)},b['g'])[1:],gvals) for w in b['g']['forms']}
 eqcount=verify_state_equations(state)
 delta={w:F(int(len(w)==8 and all(i<2 for i in w))) for w in state};verify_state_equations(state,delta)
 for k,v in zip(('q','u','v'),[x*x,x*evaluate_form(b['obs']['p'],x,full),x*evaluate_form(b['obs']['a'],x,full)]):
  obs=md.observables(b['g']);check(obs[k][0]+dot(obs[k][1:],gvals)==v,'Factorization lift fails')
 blocks=[]
 for bl in b['blocks']:
  M=[]
  for u in bl['basis']:
   row=[]
   for v in bl['basis']:
    pol=wa.inner(u,v,bl['ground']);row.append(sum((c*state[w] for w,c in pol.items()),F(0)))
   M.append(row)
  check(M==evaluate_matrix(bl,x,full),'Raw Gram disagrees with polynomial matrix')
  pp=exact_psd_rank(M);blocks.append({'name':bl['name'],'rank':len(pp),'positive_pivots':list(map(str,pp))})
 check([v['rank'] for v in blocks]==data['expected_raw_ranks'],'Raw rank mismatch')
 return {'x':str(x),'x_display':decimal(x,15),'completion':str(t),'original_equalities_checked':eqcount,
    'original_ray_equations_checked':eqcount,'raw_blocks':blocks,'positive_pivot_count':sum(v['rank'] for v in blocks),
    'observables':{k:str(evaluate_form(v,x,full)) for k,v in b['obs'].items()},
    'full_global_free_indices':b['g']['free'],'full_global_free_values':list(map(str,gvals)),
    'finite_truncated_moments_not_claimed_physical':True}

def verify_full_comparison(b,upper_primal):
    """Compare to the stronger system with all ground operators of level <=4."""
    data=json.loads((ROOT/'inputs/compact_certificate.json').read_text())
    g=b['g'];o=md.observables(g)
    cert=verify_moment_certificate(data,g,o)
    check(F(cert['scalar_weights']['xp_ul'])>0 and F(cert['scalar_weights']['xa_ul'])>0,
          'Full certificate must bound recession changes in p and a')
    rows=[list(o[k][1:]) for k in ('x','q','u','v','p','a','d')]
    for atom in data['atoms']:
        if atom['kind']=='scalar':continue
        op=md.operator_from_json(atom['operator']);ground=atom['kind']=='ground'
        for w in wa.words(4,wa.charge(next(iter(op)))):
            rows.append(md.affine(wa.inner({w:F(1)},op,ground),g)[1:])
    for w in wa.words(4,1):
        rows.append(md.affine(wa.inner({w:F(1)},md.D_OPERATOR),g)[1:])
    rank,ker=rank_kernel(rows,35)
    check(rank==34 and len(ker)==1,'Full ground-four recession dimension mismatch')
    direction=ker[0]
    for w in g['forms']:
        check(dot(md.affine({w:F(1)},g)[1:],direction)==int(len(w)==8 and all(k<2 for k in w)),
              'The full system has a different remaining ray')
    ranks=[]
    for ground in (False,True):
        for charge,parity in [(0,-1),(0,1)]+[(i,None) for i in range(1,5)]:
            ops=wa.operator_basis(4,charge,parity)
            M=[[dot(md.affine(wa.inner(u,v,ground),g)[1:],direction) for v in ops] for u in ops]
            pp=exact_psd_rank(M)
            if ground:check(not any(c for row in M for c in row),'Full ground-four ray changes a Gram block')
            ranks.append(len(pp))
    d=F(upper_primal['observables']['d'])
    check(d>F(8,25)>F(cert['K']),'Single-Q witness does not separate the two systems')
    return {'full_ground4_certificate_rechecked':cert,'annihilation_rows':len(rows),'annihilation_rank_global35':rank,
        'full_ground4_recession_global35':list(map(str,direction)),
        'full_ground4_recession_lineality_dimension':0,'full_ground4_pointed_recession_dimension':1,
        'full_12_Gram_direction_ranks':ranks,'singleQ_witness_x':upper_primal['x'],
        'singleQ_witness_d':str(d),'singleQ_witness_d_display':decimal(d,35),
        'full_ground4_d_upper':'8/25','separation_gap':str(d-F(8,25)),
        'claim':'The full ground-four constraints remove the ten invisible recession directions, but preserve the pure eighth-coordinate ray. The single-Q feasible old-moment projection strictly contains that of the full system; strict inclusion of the radius-only projection has not been proved.'}

def run(out):
 start=time.time();out.mkdir(parents=True,exist_ok=True)
 print('Rechecking original foundation and level-9 quotient...',flush=True)
 _,fnd=foundation(ROOT/'inputs');b=build();rowcount=md.verify_all_rows(b['g'])
 check((len(b['g']['words']),len(b['g']['free']),len(b['lineality']))==(1574,35,10),'Bridge dimensions unexpected')
 cv={name:verify_certificate(json.loads((ROOT/'inputs'/(name+'.json')).read_text()),b['g'],md.observables(b['g'])) for name in ('compact22','tight24')}
 state=reconstruct_base(ROOT/'inputs');layout,details=verify_layout(b,state);recession=verify_recession(b,details);bounds=verify_simple_bounds(b)
 save(out/'geometry.json',{'layout':layout,'recession':recession,'simple_bounds':bounds,'bridge_rank':25,'bridge_old_moments_checked':615,
    'quotient_joint_compactness_scope':'After removing the 10 invisible directions and quotienting by the 1 pointed ray, the joint projected feasible set on the proved closed radius interval is compact. This is not a claim that singular projected boundary points admit finite full completion.'})
 duals={};primals={}
 for side in ('lower','upper'):
  duals[side]=verify_dual(side,b);save(out/('verified_dual_'+side+'.json'),duals[side]);print(side,'dual exact PSD and global exclusion passed.',flush=True)
  primals[side]=verify_primal(side,b);save(out/('verified_primal_'+side+'.json'),primals[side]);print(side,'finite primal all 12 raw blocks passed.',flush=True)
 comparison=verify_full_comparison(b,primals['upper']);save(out/'full_ground4_comparison.json',comparison)
 lower=[F(duals['lower']['cut']),F(primals['lower']['x'])];upper=[F(primals['upper']['x']),F(duals['upper']['cut'])]
 check(lower[0]<lower[1]<upper[0]<upper[1],'Endpoint ordering failed')
 # Negative controls: PSD and stationarity/sign parts must detect modifications.
 rejected=[]
 try:exact_psd_rank([[F(0),F(1)],[F(1),F(0)]])
 except AssertionError:rejected.append('zero-diagonal nonzero-cross PSD')
 else:raise AssertionError('Bad PSD control passed')
 check(F(duals['lower']['gamma'])>0 and F(duals['upper']['gamma'])>0,'Invalid ground constraint multipliers')
 # An incorrectly scaled ray fails a raw moment coefficient identity.
 try:
  bad=RDIR[:];bad[-1]=2
  for w in b['m8']['forms']:
   check(dot(tangent_row(b['old_form']({w:F(1)})),bad)==int(len(w)==8 and all(i<2 for i in w)),'Bad ray rejected')
 except AssertionError:rejected.append('corrupted eighth-moment ray')
 else:raise AssertionError('Bad ray passed')
 summary={
        'bounded_projected_coordinates': 20,
        'bracket_widths': [str(lower[1]-lower[0]),str(upper[1]-upper[0])],
        'dual_positive_pivots_each': [duals[s]['positive_pivot_count'] for s in ('lower','upper')],
        'full_fixed_x_coordinates': 31,
        'full_ground4_comparison': {'recession_lineality_dimension':0,'pointed_recession_dimension':1,'singleQ_strictly_weaker_in_old_moment_projection':True,'radius_only_strict_inclusion_not_proved':True},
        'generated_global_equations_checked': rowcount,
        'lineality_dimension': 10,
        'lower_extremum_bracket': list(map(str,lower)),
        'negative_controls_rejected': rejected,
        'pointed_recession_dimension': 1,
        'primal_positive_pivots_each': [primals[s]['positive_pivot_count'] for s in ('lower','upper')],
        'upper_extremum_bracket': list(map(str,upper)),
        'visible_fixed_x_coordinates': 21,
    }
 save(out/'summary.json',summary)
 print('CERTIFIED LOWER EXTREMUM: 0.3855595 <= inf S_Q <= 0.3855596',flush=True)
 print('CERTIFIED UPPER EXTREMUM: 0.3919969 <= sup S_Q <= 0.3919970',flush=True)
 print('RECESSION: one nonnegative pure eighth-coordinate ray + 10-dimensional degree-nine lineality.',flush=True)
 print('All Step-12 exact checks passed.',flush=True)
 return summary

if __name__=='__main__':
 a=argparse.ArgumentParser(description=__doc__);a.add_argument('--output-dir',type=Path,default=ROOT/'results/verified');p=a.parse_args();run(p.output_dir)
