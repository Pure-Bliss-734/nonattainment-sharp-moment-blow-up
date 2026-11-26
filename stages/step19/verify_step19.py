#!/usr/bin/env python3
"""Certify the globally optimal leading prefactor of the eighth-moment cost.

No numerical optimizer, numerical derivative, or floating PSD tolerance is used.
Roots, global consequences, Stage-17 estimates, tangent rank, and continuous
Puiseux completions are regenerated. The original word quotient and completion
interpretation remain shared premises documented in the proof in the cited paper.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,sys,copy
from frozen_source import prepare,ROOT as HERE,SHA18
P18,P17,P16,HASHES=prepare()
from independent_model import Model,ROOT as SRC
from verify_step16 import prove_root
from independent_consequences import audit_consequences
from neighborhood_model import bound_data,rowwise_bounds
from verify_step17 import check_parameters
from prefactor_model import formal_path_checks,build_prefactor,continuous_family,quantitative_lower,full_sample
from decimal_interval import Box,display,ldlt,solve_spd,mul,transpose,dot,ONE,ZERO

def encode(v):
    if isinstance(v,F):return str(v)
    if isinstance(v,Box):return v.record()
    if isinstance(v,dict):return {str(k):encode(w) for k,w in v.items()}
    if isinstance(v,(list,tuple)):return [encode(w) for w in v]
    return v

def dump(path,v):path.write_text(json.dumps(encode(v),indent=2,ensure_ascii=False)+'\n')

def scalar_tests():
    # A two-block exact reference with a common completion t.
    w=[F(2),F(5)];q=[F(3),F(4)];C=sum(a*b*b for a,b in zip(w,q));Omega=sum(w)
    assert C==98
    slopes=[v*v/C for v in q]
    assert sum(a*b for a,b in zip(w,slopes))==1
    checks=0
    for delta in [F(1,10),F(1,1000),F(1,10**12)]:
        t=C/delta
        for a,b in zip(slopes,q):
            assert (delta*a)*t-b*b==0;checks+=1
    for u in [F(1,100),F(1,10**3),F(1,10**9)]:
        for a,b in zip(slopes,q):
            H=(1-u)*a+u/Omega
            det=H*C*(1+2*u)-b*b
            exact=u*(C/Omega+b*b)+u*u*(2*C*(F(1,Omega)-a))
            assert det==exact and det>0;checks+=1
    # A rectangular, exactly known Schur Jacobian right inverse.
    J=[[Box(1),Box(0),Box(1)],[Box(0),Box(2),Box(0)]]
    gram=mul(J,transpose(J));gi,_=solve_spd(gram,[[ONE,ZERO],[ZERO,ONE]])
    R=mul(transpose(J),gi);ident=mul(J,R)
    for i in range(2):
        for j in range(2):assert ident[i][j].contains(int(i==j));checks+=1
    return {'exact_reference_checks':checks,'reference_C':'98','reference_cases':'two common-t completions, positive regularization, rectangular right inverse'}

def finite_lower_sample(pref,lo,delta):
    # These samples use delta=10^(-even number) so sqrt(delta) is rational.
    from math import isqrt
    a=isqrt(delta.numerator);b=isqrt(delta.denominator)
    if a*a!=delta.numerator or b*b!=delta.denominator:raise AssertionError('sample sqrt is not rational')
    u=F(a,b)
    E=sum(c*u**i for i,c in enumerate(lo['E_coefficients_in_u']))
    den=sum(c*u**i for i,c in enumerate(lo['T_over_delta_upper_coefficients']))
    lower=(pref['C'].lo()-E)/den-lo['common_diagonal_upper']*delta
    upper=pref['C'].hi()*(1+2*u)
    return {'delta':delta,'scaled_infimum_lower':lower,'scaled_infimum_upper':upper,
            'relative_gap_upper':(upper-lower)/pref['C'].lo(),
            'lower_decimal':display(lower,45),'upper_decimal':display(upper,45)}

def main():
    p=argparse.ArgumentParser();p.add_argument('--side',choices=['lower','upper','both'],default='both')
    p.add_argument('--output-dir',default='results/verified');args=p.parse_args()
    out=Path(args.output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    model=Model();formals=formal_path_checks(model);dump(out/'formal_identities.json',formals)
    tests=scalar_tests();negative=[]
    def reject(name,fn):
        try:fn()
        except (AssertionError,ValueError,ZeroDivisionError):negative.append(name)
        else:raise AssertionError('false proof input accepted: '+name)
    reject('nonpositive_Gram',lambda:ldlt([[Box(0),Box(1)],[Box(1),Box(0)]]))
    reject('rank_deficient_right_inverse',lambda:solve_spd([[Box(1),Box(1)],[Box(1),Box(1)]],[[ONE,ZERO],[ZERO,ONE]]))
    oldparams=json.loads((P17/'inputs/neighborhood_parameters.json').read_text())
    params=json.loads((HERE/'inputs/parameters.json').read_text())
    reports={}
    for side in (['lower','upper'] if args.side=='both' else [args.side]):
        data=json.loads((SRC/f'inputs/kkt_{side}.json').read_text())
        root,box,infos=prove_root(model,data)
        consequences=audit_consequences(model,data,box,infos)
        bd=bound_data(model,data,box,infos);rr=rowwise_bounds(model,bd,F(oldparams[side]['R']))
        lower17=check_parameters(bd,rr,oldparams,side)
        print(side.upper()+': root, global uniqueness/nonattainment, and Stage-17 uniform estimates re-certified.',flush=True)
        pref=build_prefactor(model,data,box,infos)
        family=continuous_family(model,data,box,infos,pref,F(params[side]['u_limit']))
        lower=quantitative_lower(model,box,pref,consequences,oldparams,side)
        actual_samples=[full_sample(model,box,pref,F(d)) for d in params['full_sample_distances']]
        dump(out/f'{side}_full_samples.json',actual_samples)
        if family['distance_limit']>lower['radius']:raise AssertionError('new upper strip outside lower strip')
        samples=[finite_lower_sample(pref,lower,F(d)) for d in params['sample_distances']]
        if any(s['delta']>family['distance_limit'] for s in samples):raise AssertionError('sample upper outside strip')
        for s in samples:
            if s['scaled_infimum_lower']>=s['scaled_infimum_upper']:raise AssertionError('lower/upper contradiction')
        previous=json.loads((P18/f'results/verified/{side}_family.json').read_text())
        old=previous['minimal_completion_along_this_path_poles']['M0_1']['distance_times_minimal_completion_limit']
        ratio=Box(F(old['lower']),F(old['upper']))/pref['C']
        # Errors in interval inputs must actually fail, not merely be logged.
        reject(side+'_zero_u_strip',lambda:continuous_family(model,data,box,infos,pref,F(0)))
        reject(side+'_excessive_u_strip',lambda:continuous_family(model,data,box,infos,pref,F(1,100)))
        bad=copy.deepcopy(infos)
        first=next(a for a in bad if model.small[a['item']['index']]['parent']=='M0_1')
        first['multiplier'][0][0]=Box(-1)
        reject(side+'_negative_completion_weight',lambda:build_prefactor(model,data,box,bad))
        # Deliberately breaking normalization changes d_x away from sigma.
        if not (pref['d0'][0]*2-Box(pref['sigma'])).contains(0):negative.append(side+'_doubled_budget_rejected')
        else:raise AssertionError('normalization bad control not detected')
        dump(out/f'{side}_root_and_baseline.json',{'root':root,'consequences':consequences,'neighborhood17':lower17})
        dump(out/f'{side}_prefactor.json',pref)
        dump(out/f'{side}_continuous_family.json',family)
        dump(out/f'{side}_universal_lower.json',lower)
        dump(out/f'{side}_finite_brackets.json',samples)
        reports[side]={'C':pref['C'],'exact_definition':'sum over i=0,2 of omega_i*(nu_i^T*b_i_star)^2, normalized by KKT x derivative sigma',
          'tangent_schur_rank':pref['tangent_rank'],'u_limit':family['u_limit'],'distance_limit':family['distance_limit'],
          'old_fixed_path_prefactor_divided_by_optimal_C':ratio,
          'limit':'lim_{delta->0+} delta*tau_side(delta)=C',
          'explicit_attaining_family':'y=y_star+delta*d0+delta^(3/2)*(dI-d0); t=C*(1+2*sqrt(delta))/delta',
          'universal_lower_proved_without_differentiable_path_assumption':True,'samples':samples}
        print('  CERTIFIED C in ['+display(pref['C'].lo(),60)+', '+display(pref['C'].hi(),60)+']',flush=True)
        print('  ALL 0<delta<= '+str(family['distance_limit'])+': full feasible t=C*(1+2*sqrt(delta))/delta.',flush=True)
    tests['negative_controls']=negative;tests['negative_control_count']=len(negative);dump(out/'tests.json',tests)
    forbidden=[k for k in ['numpy','scipy','sympy','mpmath'] if k in sys.modules]
    if forbidden:raise AssertionError('third-party library imported: '+str(forbidden))
    report={
        'endpoints': reports,
        'exact_formal_checks': formals,
        'main_result': 'Both globally minimal scaled eighth-moment costs have limits, equal to the certified dual budget C.',
        'tests': tests,
    }
    dump(out/'summary.json',report)
    print('Global prefactor equality proved by a universal lower bound and an explicit continuous attaining family.')
    print('All Step-19 exact checks passed.')
if __name__=='__main__':main()
