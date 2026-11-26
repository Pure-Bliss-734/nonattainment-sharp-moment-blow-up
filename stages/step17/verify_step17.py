#!/usr/bin/env python3
"""Verify the full-fiber neighborhood estimates."""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import argparse,copy,json,sys,subprocess
from prepare_source import ROOT,prepare,SHA
SRC, SOURCE_COUNT=prepare()
from neighborhood_model import (Model,Box,ZERO,ONE,inverse,eye,mul,display,
    derivative,from_form,bound_data,rowwise_bounds,sqrt_up)
from verify_step16 import prove_root
from independent_consequences import audit_consequences

def encode(v):
    if isinstance(v,F):return str(v)
    if isinstance(v,Box):return v.record()
    if isinstance(v,Path):return str(v)
    if isinstance(v,list):return [encode(x) for x in v]
    if isinstance(v,tuple):return [encode(x) for x in v]
    if isinstance(v,dict):return {str(k):encode(x) for k,x in v.items()}
    return v

def dump(path,obj):
    path.write_text(json.dumps(encode(obj),indent=2,ensure_ascii=False)+'\n')

def check_parameters(b,q,params,side):
    p=params[side]
    d=F(params['delta']);R=F(p['R']);smin=F(params['cross_floor'])
    cmax=F(params['c0_upper']);lambda_min=F(params['kernel_weight_lower'])
    if not (F(0)<d<=F(1,100)):
        raise AssertionError('invalid neighborhood radius')
    if q['R']!=R:raise AssertionError('coordinate radius mismatch')
    vals={k:F(p[k]) for k in ['coordinate_A','coordinate_B','cross_A','cross_B','root_cross_lower']}
    for key,new in [('coordinate_A','coord_A'),('coordinate_B','coord_B'),('cross_A','cross_A'),('cross_B','cross_B')]:
        if vals[key]<q[new]:raise AssertionError('underestimated '+key)
    if b['cross_at_root'].lo()<=vals['root_cross_lower']:
        raise AssertionError('uncertified lower bound on root crossing')
    if b['weight'].lo()<=lambda_min:raise AssertionError('weight lower bound fails')
    # For r>=R, all row-bound ratios are decreasing in r. This strict inequality
    # contradicts r<=||H^-1 h||, so r<R for every feasible point in the x strip.
    A=vals['coordinate_A'];B=vals['coordinate_B']
    positive_side=R-B*d
    coordinate_margin=positive_side**2-A*A*d
    if positive_side<=0 or coordinate_margin<=0:
        raise AssertionError('coordinate contraction inequality failed')
    gap=vals['root_cross_lower']-smin-vals['cross_B']*d
    cross_margin=gap**2-vals['cross_A']**2*d
    if gap<=0 or cross_margin<=0:raise AssertionError('crossing gap failed')
    c0root=b['raw_root'][16]/2
    if c0root.hi()+R/2>=cmax:raise AssertionError('c0 bound failed')
    # 0<=support<=2*distance and lambda>1/250 imply rho<500*distance.
    if F(params['support_lipschitz_upper'])!=2:raise AssertionError('support normalization changed')
    denominator=F(2)/lambda_min/(smin*smin)
    if denominator!=125000000:raise AssertionError('wrong pole coefficient')
    pure_coefficient=F(1,denominator)-cmax*d
    if pure_coefficient!=F(255,32000000000):raise AssertionError('wrong unshifted coefficient')
    return {'side':side,'delta':d,'R':R,'rounded_constants':vals,
            'squared_coordinate_margin':coordinate_margin,'squared_cross_margin':cross_margin,
            'c0_root_interval':c0root,'c0_uniform_upper':cmax,
            'kernel_norm_factor_upper':F(2)/lambda_min,
            'shifted_pole_denominator':denominator,'unshifted_pole_coefficient':pure_coefficient,
            'coordinate_norm_strictly_below_R':True,'cross_strictly_above_floor':True,
            'all_feasible_points_in_strip_covered':True}

def tests():
    from math import isqrt
    count=0
    for n in range(80):
        for den in [1,3,7,1001]:
            q=F(n,den);u=sqrt_up(q,1000)
            assert u*u>=q and (u==0 or (u-F(1,1000))**2<q)
            count+=1
    matrices=[([[2,1],[1,3]],[[F(3,5),F(-1,5)],[F(-1,5),F(2,5)]]),
              ([[0,2],[3,4]],[[F(-2,3),F(1,3)],[F(1,2),F(0)]]),
              ([[5]],[[F(1,5)]])]
    cells=0
    for a,ref in matrices:
        inv,_=inverse([[Box(x) for x in r] for r in a])
        for i,row in enumerate(ref):
            for j,v in enumerate(row):
                assert inv[i][j].contains(v);cells+=1
    rejected=[]
    try:inverse([[Box(1),Box(1)],[Box(1),Box(1)]])
    except ValueError:rejected.append('singular_inverse')
    else:raise AssertionError('singular matrix accepted')
    try:sqrt_up(F(-1))
    except ValueError:rejected.append('negative_square_root')
    else:raise AssertionError('negative radicand accepted')
    return {'rational_sqrt_checks':count,'known_inverse_entries_checked':cells,'negative_controls':rejected}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output-dir',default='results/verified')
    args=ap.parse_args();out=Path(args.output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    test=tests();model=Model();params=json.loads((ROOT/'inputs/neighborhood_parameters.json').read_text())
    count=0
    # All raw blocks must be affine in z, with only x^2 and x*z parameter terms.
    for b in model.raw:
        for row in b['mat']:
            for poly in row:
                if any(sum(i>0 for i in mon)>1 or len(mon)>2 for mon in poly):
                    raise AssertionError('raw matrix outside quadratic-affine class')
                count+=1
    if from_form(model.data['completion']['M0_1']['diagonal'])!={(16,):F(1,2)}:
        raise AssertionError('c0 != z15/2')
    reports={};negative=list(test['negative_controls'])
    for side in ['lower','upper']:
        data=json.loads((SRC/f'inputs/kkt_{side}.json').read_text())
        root,values,infos=prove_root(model,data)
        old=audit_consequences(model,data,values,infos)
        if F(old['global']['support_lipschitz_bound'])>=2:raise AssertionError('support > 2')
        b=bound_data(model,data,values,infos)
        q=rowwise_bounds(model,b,F(params[side]['R']))
        report=check_parameters(b,q,params,side)
        report['endpoint_box']=root['root_x']
        report['root_proof']={'dimension':root['KKT_dimension'],'eta':root['eta_upper'],'rho':root['rho_upper']}
        report['support_lipschitz_actual_bound']=old['global']['support_lipschitz_bound']
        report['computed_constants']={k:q[k] for k in ['coord_A','coord_B','cross_A','cross_B']}
        dump(out/f'{side}_neighborhood.json',report)
        dump(out/f'{side}_coefficient_certificate.json',{'bounds':b,'rowwise':q})
        # Negative controls exercise exactly the gates used in the main theorem.
        for label,change in [
            ('underestimated_cross_A',lambda p:p[side].__setitem__('cross_A','0')),
            ('enlarged_delta',lambda p:p.__setitem__('delta','1/10000000')),
            ('false_root_cross_bound',lambda p:p[side].__setitem__('root_cross_lower','1/10')),
            ('false_c0_bound',lambda p:p.__setitem__('c0_upper','1/1000000'))]:
            bad=copy.deepcopy(params);change(bad)
            try:check_parameters(b,q,bad,side)
            except AssertionError:negative.append(side+'_'+label)
            else:raise AssertionError('bad parameter set accepted: '+label)
        reports[side]=report
        print(side.upper()+': root re-proved; coefficient inequalities and delta=1e-9 certified.',flush=True)
    delta=F(params['delta']);coefficient=F(255,32000000000)
    examples=[]
    for d in [F(1,10**9),F(1,10**12),F(1,10**16)]:
        bound=F(1,125000000)/d-F(1,32)
        examples.append({'endpoint_distance':d,'strict_shifted_free_lower_bound':bound,'display':display(bound)})
    global_cap={'domain':'any finite single-Q feasible moment sequence; cap t<=T is an additional conditional hypothesis',
                'necessary_endpoint_distance':'min(1/1000000000, 1/(125000000*(T+1/32)))',
                'valid_for':'T>=0 and either certified endpoint; no physical upper cap is assumed'}
    dump(out/'moment_bounds.json',{'examples':examples,'conditional_moment_cap':global_cap,
         'uniform_eighth_moment_bound':'t > 1/(125000000*d)-1/32 >= 255/(32000000000*d), 0<d<=1e-9'})
    test['negative_controls']=negative;test['all_bad_inputs_rejected']=True
    dump(out/'tests.json',test)
    forbidden=[s for s in ['numpy','scipy','sympy','mpmath'] if s in sys.modules]
    if forbidden:raise AssertionError('unexpected third-party proof dependency')
    summary={
        'c0_upper': F(1,32),
        'coordinate_count': 20,
        'cross_floor': F(1,500),
        'delta': delta,
        'endpoints': reports,
        'raw_polynomial_entries_checked': count,
        'tests': test,
        'unshifted_pole_coefficient': coefficient,
    }
    dump(out/'summary.json',summary)
    print('CERTIFIED: 0 < abs(x-xi) <= 1/1000000000 implies s > 1/500 and c0 < 1/32.')
    print('CERTIFIED: t > 1/(125000000*abs(x-xi)) - 1/32 >= 255/(32000000000*abs(x-xi)).')
    print('All Step-17 exact checks passed.')
if __name__=='__main__':main()
