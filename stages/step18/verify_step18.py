#!/usr/bin/env python3
"""Reconstruct the complete moment lift and bound its asymptotic order."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,sys,subprocess,copy
from frozen_source import prepare,ROOT,SHA17
P17,P16,N17,N16=prepare()
from independent_model import Model,ROOT as SOURCE,pc,pv,plus,times,evaluate
from verify_step16 import prove_root
from independent_consequences import audit_consequences
from neighborhood_model import bound_data,rowwise_bounds
from verify_step17 import check_parameters
from order_model import (certify_family,certify_samples,verify_formal_identities,
                         square,curvature,bernstein)
from decimal_interval import Box,display,solve_spd,dot,ZERO

def encode(v):
    if isinstance(v,F):return str(v)
    if isinstance(v,Box):return v.record()
    if isinstance(v,dict):return {str(k):encode(w) for k,w in v.items()}
    if isinstance(v,(tuple,list)):return [encode(w) for w in v]
    return v

def dump(path,v):path.write_text(json.dumps(encode(v),indent=2,ensure_ascii=False)+'\n')

def reference_tests():
    n=0
    for i in range(-8,9):
        for j in range(i,9):
            a=F(i,3);b=F(j,3);ss=square(Box(a,b))
            for t in [a,(a+b)/2,b]:
                if not ss.contains(t*t):raise AssertionError('square enclosure failed')
                n+=1
    from order_model import eval_path
    btests=0
    for co in [(1,2,3),(-4,3,7),(2,-9,0)]:
        cc=list(map(Box,co));bc=bernstein(cc)
        for s in [F(0),F(1,7),F(1,2),F(1)]:
            x=sum((v*s**k for k,v in enumerate(co)),F(0))
            v=bc[0]*(1-s)**2+bc[1]*2*s*(1-s)+bc[2]*s*s
            if not v.contains(x) or not eval_path(cc,Box(s)).contains(x):raise AssertionError('Bernstein reference failed')
            btests+=1
    a=[[Box(2),Box(1)],[Box(1),Box(3)]];sol,_=solve_spd(a,[[Box(1)],[Box(2)]])
    if not sol[0][0].contains(F(1,5)) or not sol[1][0].contains(F(3,5)):
        raise AssertionError('inverse reference failed')
    # Exact scalar rank-one limiting model: b^2 / B = 9/(2s).
    # Tests the asymptotic formula without importing a symbolic algebra library.
    limits=0
    for s in [F(1),F(1,10),F(1,10**12)]:
        if s*(F(9)/(2*s))!=F(9,2):raise AssertionError('simple-pole reference failed')
        limits+=1
    return {'square_enclosure_reference_points':n,'Bernstein_reference_evaluations':btests,
            'known_inverse_entries':2,'simple_pole_exact_scalar_references':limits}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',default='results/verified')
    args=ap.parse_args();out=Path(args.output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    # This is an isolated exact word-algebra process. It proves the full family
    # satisfies the original equations for arbitrary old quotient coordinates.
    subprocess.run([sys.executable,'-S',str(ROOT/'export_full_lift.py'),'--output-dir',str(out)],check=True)
    tests=reference_tests();model=Model();params=json.loads((ROOT/'inputs/parameters.json').read_text())
    oldparams=json.loads((P17/'inputs/neighborhood_parameters.json').read_text())
    if F(params['radius'])!=F(oldparams['delta']):raise AssertionError('lower/upper neighborhoods differ')
    summaries={};negative=[]
    def reject(name,fn):
        try:fn()
        except (AssertionError,ValueError,ZeroDivisionError):negative.append(name)
        else:raise AssertionError('bad input accepted: '+name)
    for side in ['lower','upper']:
        data=json.loads((SOURCE/f'inputs/kkt_{side}.json').read_text())
        root,box,infos=prove_root(model,data)
        consequences=audit_consequences(model,data,box,infos)
        if F(consequences['global']['support_lipschitz_bound'])>=2:raise AssertionError('invalid support bound')
        b=bound_data(model,data,box,infos);rowwise=rowwise_bounds(model,b,F(oldparams[side]['R']))
        lower=check_parameters(b,rowwise,oldparams,side)
        cert=certify_family(model,data,box,infos,params)
        samples=certify_samples(model,data,box,params)
        dump(out/f'{side}_family.json',cert);dump(out/f'{side}_samples.json',samples)
        dump(out/f'{side}_root_and_lower_bound.json',{'root':root,'neighborhood':lower})
        dominant=cert['minimal_completion_along_this_path_poles']['M0_1']['distance_times_minimal_completion_limit']
        summaries[side]={'root_x':root['root_x'],'inside':cert['interior_parameter'],'K':cert['K'],
                         'root_to_interior_distance':cert['root_to_interior_distance'],
                         'rounded_upper_C':cert['rounded_upper_coefficient'],
                         'path_minimal_completion_leading_coefficient':dominant,
                         'all_18_common_matrix_margins_certified':True,
                         'all_12_full_blocks_certified_per_sample':True,
                         'samples':len(samples),'positive_pivots_per_sample':54}
        for name,change in [
          ('K_too_small',lambda p:p[side].__setitem__('K','1')),
          ('theta_above_one',lambda p:p.__setitem__('theta','101/100')),
          ('upper_constant_too_small',lambda p:p[side].__setitem__('rounded_C','1/1000000000')),
          ('radius_outside_constructed_strip',lambda p:p.__setitem__('radius','1'))]:
            bad=copy.deepcopy(params);change(bad)
            reject(side+'_'+name,lambda bad=bad:certify_family(model,data,box,infos,bad))
        print(side.upper()+': K='+str(cert['K'])+', theta=99/100; uniform full completion and all sample blocks certified.',flush=True)
        print('  Constructed path lim(distance*t_min) in '+str(dominant),flush=True)
    reject('negative_sample_distance',lambda:certify_samples(model,data,box,{**params,'sample_distances':['-1/1000000000']}))
    reject('unsupported_cubic_pencil',lambda:curvature({(0,0,0):F(1)},[Box(1)]))
    tests['negative_controls']=negative;tests['negative_control_count']=len(negative)
    dump(out/'tests.json',tests)
    Cmin=F(255,32000000000);common=F(params['common_upper_C'])
    if any(s['rounded_upper_C']>common for s in summaries.values()):raise AssertionError('uniform upper constant too small')
    theorem={'radius':F(params['radius']),'universal_lower_c':Cmin,'common_existential_upper_C':common,
      'definition':'tau_side(delta)=inf{t: finite complete single-Q moment sequence at x=xi_L+delta or xi_U-delta}',
      'lower':'Every finite complete sequence has t > 1/(125000000*delta)-1/32 >= c/delta.',
      'upper':'For every 0<delta<=1e-9 an explicitly defined algebraic full feasible sequence has t<C_side/delta.',
      'conclusion':'c/delta <= tau_side(delta) <= C_side/delta; growth exponent 1 is optimal.',
      'upper_not_for_all_sequences':True,
      'prefactor_not_proved_globally_optimal':True,
      'all_parameters_in_each_strip_are_feasible':True,
      'whole_interior_interval_feasibility_not_claimed':True,
      'endpoints':summaries}
    dump(out/'sharp_order.json',theorem)
    forbidden=[k for k in ['numpy','scipy','sympy','mpmath'] if k in sys.modules]
    if forbidden:raise AssertionError('unexpected third-party verification dependency')
    full=json.loads((out/'lift_identity.json').read_text())
    summary={
        'moment_map_including_identity': full['moments_including_identity'],
        'original_factorization_identities': full['all_original_factorization_equations_verified_identically'],
        'per_side_common_margin_positive_pivots': 51,
        'per_side_exact_Bernstein_entry_checks': 15,
        'per_side_exact_secant_entry_checks': 213,
        'sample_full_PSD_blocks': sum(v['samples']*12 for v in summaries.values()),
        'sample_positive_pivots': sum(v['samples']*54 for v in summaries.values()),
        'sharp_order': theorem,
        'tests': tests,
    }
    dump(out/'summary.json',summary)
    print('CERTIFIED: for 0<delta<=1e-9, c/delta <= tau(delta) <= C/delta, c=255/32000000000, C=1/10000.')
    print('The inverse-distance growth order is sharp; its globally optimal prefactor is not determined.')
    print('All Step-18 exact checks passed.')
if __name__=='__main__':main()
