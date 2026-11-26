#!/usr/bin/env python3
"""Verify the full-coordinate lower remainder bound."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,subprocess,sys,shutil
from stage23_sources import ROOT,SOURCE_SHA,PENCILS_SHA,prepare
from lower_model import reconstruct,bound_deficit,coefficient_report,bm,ap,I,Z,O,zeros,tr,mm,scale,sub,mv,val,plus,poly,sc,DATA,INTEGRITY,zero_diagnostic
from exact_tests import run as reference_tests

CORE=['exact_checks.json','negative_controls.json','summary.json']+[f'{s}_{k}.json' for s in ('lower','upper') for k in ('coefficients','universal_lower','full_samples')]

def dump(path,obj):path.write_text(json.dumps(bm.encode(obj),indent=2,ensure_ascii=False)+'\n')

def sample_certificate(model,curve,upper_certificate,delta):
    point=bm.sample_quadratic(DATA,curve,upper_certificate,F(delta))
    coords=point['raw_coordinates'];t=point['t'];C=model['C'];d=F(delta)
    square=Z
    for part in model['partners']:
        B=ap.matval(model['pref']['raw'][part['parent']],coords)
        K=part['K'];V=part['V']
        if part['completion_bottom'] is not None:
            cp=DATA['completion'][part['parent']]
            b=[ap.peval(ap.poly(p),coords) for p in cp['cross']]
            c=ap.peval(ap.poly(cp['diagonal']),coords)
            B=[row+[v] for row,v in zip(B,b)]+[b+[t+c]]
            K=K+[[Z]];V=V+[part['completion_bottom']]
        shifted=sub(V,scale(K,C/d))
        value=mm(part['W'],mm(tr(shifted),mm(B,shifted)))
        square=square+sum((value[i][i] for i in range(len(value))),Z)
    S=val(model['deficit'],coords)
    # Stable excess comes from the exact parameter conversion, not pole subtraction.
    ex=point['excess_over_C_div_delta'];rhs=C*(ex+S)
    zero_diagnostic(square-rhs,'full completed-square sample')
    if square.lo()<=0 or (ex+S).lo()<=0:raise AssertionError('sample square positivity not certified')
    return {'delta':d,'x':point['x'],'actual_family_excess':ex,'deficit_polynomial':S,
            'sum_of_positive_Gram_forms':square,'C_times_excess_plus_deficit':rhs,
            'full_blocks':point['blocks'],'block_count':len(point['blocks']),
            'positive_pivots':sum(len(b['pivots']) for b in point['blocks']),
            'role':'Regression on algebraic feasible points, not the proof for a continuous neighborhood.'}


def controls(models,pars):
    rejected=[]
    def fail(label,call):
        try:call()
        except (AssertionError,ValueError,ZeroDivisionError):rejected.append(label)
        else:raise AssertionError('bad input accepted: '+label)
    for side in ('lower','upper'):
        m=models[side];p=pars[side]
        fail(side+': zero deficit bound',lambda:bound_deficit(m,F(p['distance_limit']),0))
        fail(side+': outside inherited neighborhood',lambda:bound_deficit(m,F(1,10**8),10**10))
        fail(side+': wrong one-half cancellation coefficient',lambda:reconstruct(side,F(p['metric_shift']),half=F(3,4)))
    fail('nonpositive auxiliary metric',lambda:reconstruct('lower',0))
    # Dependent annihilator columns must fail, regardless of how small residuals look.
    from audit_interval import spd
    fail('singular weighted normal matrix',lambda:spd([[I(1),I(1)],[I(1),I(1)]]))
    return {'count':len(rejected),'rejected':rejected,
            'scope':'Rejection of these deliberate errors is not formal proof-assistant certification.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    a = parser.parse_args()
    if sys.flags.optimize:raise RuntimeError('Do not disable inherited assertions with -O')
    out=Path(a.output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    pars=json.loads((ROOT/'inputs/parameters.json').read_text());upper_par=json.loads((prepare()[0]/'inputs/parameters.json').read_text())
    refs=reference_tests();dump(out/'exact_checks.json',refs)
    models={};bounds={};counts={}
    for side in ('lower','upper'):
        p=pars[side];model=reconstruct(side,F(p['metric_shift']));models[side]=model
        bound=bound_deficit(model,F(p['distance_limit']),F(p['lower_remainder_bound']));bounds[side]=bound
        dump(out/f'{side}_coefficients.json',coefficient_report(model));dump(out/f'{side}_universal_lower.json',bound)
        # Recompute, not just read, the Stage22 continuous upper certificate.
        curve=bm.build_quadratic(DATA,side);up=upper_par[side]
        est=bm.estimate_strip(curve,F(up['r_limit']));uc=bm.certify_constant(curve,est,F(up['completion_m']),F(up['M']))
        if uc['distance_limit']<F(p['distance_limit']) or uc['M']>F(p['upper_remainder_bound']):raise AssertionError('upper/lower domains do not match')
        samples=[sample_certificate(model,curve,uc,F(d)) for d in up['sample_deltas']]
        dump(out/f'{side}_full_samples.json',samples)
        counts[side]={'face_rows':len(model['rows']),'face_rank':len(model['gram_pivots']),
                       'weight_positive_pivots':sum(len(i['weight_pivots']) for i in model['pref']['infos']),
                       'metric_positive_pivots':sum(len(i['pivots']) for i in model['metric_pivots']),
                       'samples':len(samples),'sample_blocks':sum(s['block_count'] for s in samples),
                       'sample_positive_pivots':sum(s['positive_pivots'] for s in samples)}
        print('CERTIFIED '+side+': C/delta - '+p['lower_remainder_bound']+' < t for EVERY feasible sequence, 0<delta<='+p['distance_limit'],flush=True)
        print('COMBINED '+side+': C/delta - '+p['lower_remainder_bound']+' <= tau(delta) <= C/delta + '+p['upper_remainder_bound'],flush=True)
    bad=controls(models,pars);dump(out/'negative_controls.json',bad)
    forbidden=[name for name in ('numpy','scipy','sympy','mpmath') if name in sys.modules]
    if forbidden:raise AssertionError('third-party final proof dependency')
    report={
        'counts': counts,
        'model_negative_controls': bad['count'],
        'quantifiers': 'Lower bound: EVERY finite complete feasible truncated moment sequence. Upper bound: existence by the rechecked Stage22 family. Neither assumes attainment of the infimum.',
        'reference_checks': refs,
        'sides': {s:{'C':models[s]['C'],
            'distance_limit':F(pars[s]['distance_limit']),
            'lower_remainder_bound':F(pars[s]['lower_remainder_bound']),
            'upper_remainder_bound':F(pars[s]['upper_remainder_bound']),
            'computed_deficit_upper':bounds[s]['computed_deficit_upper'],
            'absolute_remainder_bound':max(F(pars[s]['lower_remainder_bound']),F(pars[s]['upper_remainder_bound']))} for s in models},
    }
    dump(out/'summary.json',report)
    print('Both-sided bounded remainders proved. All Step-23 exact/interval checks passed.',flush=True)
if __name__=='__main__':main()
