#!/usr/bin/env python3
"""Step 15 exact verification (Python standard library only).

Reconstructs the pinned word model, verifies the reflection reduction, certifies
both square KKT roots, checks all primal/dual ranks, global exclusions, inner
approximation, uniqueness of each full projected optimal face, and failure of
finite completion. Numerical KKT residuals are never proof criteria.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import argparse, json, subprocess, sys, time
import load_frozen
from symmetric_model import small_blocks, symmetrize,congruence,zero,fscale
from step12_model import build
from step12_layout import projected
from completion_model import solve
from dyadic_intervals import I,eye,mm,mv,sub,dec,record,pd_pivots,nonsingular_pivots,self_test
from interval_kkt import KKTIntervals
from endpoint_analysis import certify_consequences

ROOT=load_frozen.ROOT

def check(ok,msg):
    if not ok:raise AssertionError(msg)

def rational_rank(A):
    R=[list(map(F,row)) for row in A];r=0
    if not R:return 0
    for j in range(len(R[0])):
        k=next((k for k in range(r,len(R)) if R[k][j]),None)
        if k is None:continue
        R[r],R[k]=R[k],R[r];p=R[r][j];R[r]=[v/p for v in R[r]]
        for k in range(r+1,len(R)):
            c=R[k][j]
            if c:R[k]=[u-c*v for u,v in zip(R[k],R[r])]
        r+=1
        if r==len(R):break
    return r

def verify_symmetry(cfg=None):
    if cfg is None:cfg=json.loads((ROOT/'inputs/symmetry.json').read_text())
    b=build();a=b['obs']['a'];obs=b['obs'];count=0;info=[]
    # p,a are invariant under changing only z2; L0 is (39/2)*x.
    check(not any(a[3]) and not any(obs['p'][3]),'p,a depend on reflected coordinate')
    Lform=[tuple(F(7,2)*q+3*r for q,r in zip(padded(obs['p'][i]),padded(a[i]))) for i in range(22)]
    check(Lform[0]==(0,F(39,2),0),'incorrect constant part of L')
    for ob in projected():
        n=len(ob['mat']);cc=cfg[ob['name']];S=[list(map(F,row)) for row in cc['S']];T=[list(map(F,row)) for row in cc['T']]
        check(len(S)==n and all(len(row)==n for row in S),'S dimension')
        for i in range(n):
            for j in range(n):check(sum(S[i][k]*S[k][j] for k in range(n))==int(i==j),'S is not an involution')
        check(rational_rank(T)==n,'T singular')
        r=cc['plus'];check(r+cc['minus']==n,'split dimension')
        for i in range(n):
            for j in range(n):check(sum(S[i][k]*T[k][j] for k in range(n))==(1 if j<r else -1)*T[i][j],'incorrect eigenbasis')
        even=[[symmetrize(f,a) for f in row] for row in ob['mat']]
        odd=[[tuple([f[3]]+[(F(0),)]*19) for f in row] for row in ob['mat']]
        check(congruence(even,S)==even,'even form not preserved')
        check(congruence(odd,S)==[[fscale(f,-1) for f in row] for row in odd],'odd form not reversed')
        for row in ob['mat']:
            for f in row:
                for k in range(1,22):
                    ff=padded(f[k]);ll=Lform[k]
                    if ob['name']=='single_Q':
                        check(ff[1]==ll[0] and ff[2]==0,'Q mixed term structure')
                    else:check(ff[1]==ff[2]==0,'unexpected x-dependent tangent')
        count+=n*n
        info.append({'block':ob['name'],'dimensions':[r,n-r],'involution_verified':True,'eigenbasis_rank':n})
    sm=small_blocks();check(len(sm)==18,'wrong total split block count')
    return {'reflected_coordinate':'w=z2+a/2','transformation':'z2 -> -a-z2; all other visible coordinates fixed','original_projected_blocks':11,'small_blocks':18,'original_visible_coordinates':20,'symmetric_visible_coordinates':19,'entry_checks_each_identity':count,'blocks':info,'exact_congruences_and_parameter_structure_passed':True,'preserves_radius_projection_of_closed_relaxation':True,'does_not_by_itself_assume_completion_at_singular_points':True}

def padded(p):return tuple(list(p)+[F(0)]*(3-len(p)))

def certify_root(side,data):
    c=list(map(F,data['center']));r=F(data['radius']);n=len(c)
    check(n==20+data['nc'],'KKT dimensions')
    C=[[I(F(q)) for q in row] for row in data['preconditioner']]
    check(len(C)==n and all(len(row)==n for row in C),'preconditioner dimension')
    P=KKTIntervals(side,data['layout'],data['nc'])
    fc,_,_=P.evaluate(list(map(I,c)));box=[I(q-r,q+r) for q in c]
    _,JB,infos=P.evaluate(box)
    eta=max(q.absupper() for q in mv(C,fc));D=sub(eye(n),mm(C,JB))
    rho=max(sum(q.absupper() for q in row) for row in D)
    check(rho<1,'contraction bound failed')
    check(eta+rho*r<r,'self-map inclusion failed')
    pivot_info=[];np=nd=0;seen=set()
    for d in infos:
        pp=pd_pivots(d['A']) if d['A'] else [];dd=pd_pivots(d['multiplier']);np+=len(pp);nd+=len(dd)
        seen.add(d['item']['index'])
        pivot_info.append({'block':d['item']['name'],'active':True,'primal_rank':len(pp),'kernel_dimension':len(dd),'primal_leading_pivots':[record(q) for q in pp],'dual_kernel_pivots':[record(q) for q in dd]})
    for ib,bl in enumerate(P.blocks):
        if ib in seen:continue
        M,_,_=P.matrices(ib,box[:20]);pp=pd_pivots(M);np+=len(pp)
        pivot_info.append({'block':bl['name'],'active':False,'primal_rank':len(pp),'primal_leading_pivots':[record(q) for q in pp]})
    check(np+nd==51,'rank partition wrong')
    report={'side':side,'kkt_variables':n,'all_active_schur_equations':data['nc'],'root_interval':{'lo':str(c[0]-r),'hi':str(c[0]+r),'center_decimal':data['center'][0],'radius':str(r)},'eta_upper':str(eta),'rho_upper':str(rho),'eta_decimal':dec(eta,40),'rho_decimal':dec(rho,40),'selfmap_ratio_upper':str(eta/r+rho),'unique_root_in_full_kkt_box':True,'all_active_equations_retained':True,'root_jacobian_nonsingular':True,'primal_positive_pivots':np,'dual_positive_pivots':nd,'blocks':pivot_info,'interpretation':'The actual root, not the rational center, defines the endpoint.'}
    return report,certify_consequences(side,data,P,box,infos)

def negative_controls():
    out=[]
    bad=json.loads((ROOT/'inputs/symmetry.json').read_text());bad['M0_-1']['S'][0][0]='2'
    try:verify_symmetry(bad)
    except AssertionError:out.append('wrong_involution_rejected')
    else:raise AssertionError('negative control not rejected')
    for side in ['lower','upper']:
        d=json.loads((ROOT/f'inputs/kkt_{side}.json').read_text());d['center'][0]=str(F(d['center'][0])+F(1,10**8))
        try:certify_root(side,d)
        except (AssertionError,ValueError,ZeroDivisionError):out.append(side+'_displaced_root_box_rejected')
        else:raise AssertionError('displaced box accepted')
    try:pd_pivots([[I(-1)]])
    except ValueError:out.append('negative_dual_multiplier_rejected')
    else:raise AssertionError('negative multiplier accepted')
    try:nonsingular_pivots([[I(1),I(2)],[I(1),I(2)]])
    except ValueError:out.append('dependent_face_rows_rejected')
    else:raise AssertionError('singular face minor accepted')
    return {'controls':out,'all_rejected':True,'count':len(out)}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    a = parser.parse_args()
    dest=a.output_dir;dest.mkdir(parents=True,exist_ok=True)
    tests=self_test();sym=verify_symmetry();(dest/'symmetry.json').write_text(json.dumps(sym,indent=2)+'\n')
    print('EXACT SYMMETRY: 11 projected blocks -> 18 smaller blocks; 20 -> 19 visible coordinates.',flush=True)
    summaries={}
    for side in ['lower','upper']:
        tic=time.time();data=json.loads((ROOT/f'inputs/kkt_{side}.json').read_text());root,cons=certify_root(side,data)
        (dest/f'{side}_root.json').write_text(json.dumps(root,indent=2)+'\n')
        (dest/f'{side}_endpoint.json').write_text(json.dumps(cons,indent=2)+'\n')
        summaries[side]={'root_interval':root['root_interval'],'kkt_dimension':root['kkt_variables'],'contraction_rho_upper':root['rho_upper'],'primal_rank':root['primal_positive_pivots'],'dual_rank':root['dual_positive_pivots'],'whole_optimal_face_rank':cons['optimal_face']['rank'],'cross_interval':cons['nonattainment']['cross'],'infimum_or_supremum_certified':True,'not_attained_by_finite_complete_moments':True}
        print(side.upper(),'KKT ROOT CERTIFIED; rho <=',root['rho_decimal'],flush=True)
        print('GLOBAL EXCLUSION + INNER APPROXIMATION PASSED; UNIQUE OPTIMAL PROJECTED POINT.',flush=True)
        print('CROSS > 1/250; NO FINITE COMPLETION. Seconds:',round(time.time()-tic,2),flush=True)
    neg=negative_controls();(dest/'negative_controls.json').write_text(json.dumps(neg,indent=2)+'\n')
    result={
        'both_exact_extrema_nonattained': True,
        'dyadic_precision_bits': 720,
        'endpoints': summaries,
        'interval_self_checks': tests,
        'negative_controls': neg['count'],
    }
    (dest/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Both true radius extrema are not attained by finite complete moment sequences.')
    print('All Step-15 exact checks passed.')
if __name__=='__main__':main()
