#!/usr/bin/env python3
"""Verify the finite constant-term enclosures and polynomial structure."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,sys
from stage24_sources import ROOT
import partner_model as pm
from semialgebraic_model import export_model,formulas,examples
from stage24_exact_tests import run as exact_tests
bm=pm.bm
CORE=['model_structure.json','semialgebraic_formulas.json','examples.json','exact_tests.json','negative_controls.json','summary.json']+[f'{s}_{k}.json' for s in ('lower','upper') for k in ('partner_coefficients','finite_lower','finite_upper')]

def dump(p,obj):p.write_text(json.dumps(bm.encode(obj),indent=2,ensure_ascii=False)+'\n')

def controls(models,pars):
    rejected=[]
    def reject(label,fn):
        try:fn()
        except (AssertionError,ValueError,ZeroDivisionError):rejected.append(label)
        else:raise AssertionError('bad input accepted: '+label)
    for side in ('lower','upper'):
        m=models[side];par=pars[side]
        reject(side+': impossible lower-constant bound',lambda:pm.bound_deficit(m,F(par['lower_distance_limit']),-100))
        reject(side+': outside inherited concentration range',lambda:pm.bound_deficit(m,F(1,10**8),100))
    reject('zero auxiliary regularization',lambda:pm.reconstruct('lower',0))
    reject('wrong half coefficient',lambda:pm.reconstruct('upper',F(1,10**8),half=F(3,4)))
    bad=json.loads(json.dumps(pm.DATA));bad['raw'][0]['mat'][0][0].append([[0,0,0],'1'])
    reject('nonquadratic finite pencil',lambda:export_model(bad))
    bad2=json.loads(json.dumps(pm.DATA));bad2['raw'][0]['mat'][0][1].append([[],'1'])
    reject('nonsymmetric finite pencil',lambda:export_model(bad2))
    return {'count':len(rejected),'rejected':rejected}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    v = parser.parse_args()
    if sys.flags.optimize:raise RuntimeError('Do not run with -O')
    out=Path(v.output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    pars=json.loads((ROOT/'inputs/parameters.json').read_text())
    structure=export_model(pm.DATA);dump(out/'model_structure.json',structure)
    dump(out/'semialgebraic_formulas.json',formulas());ex=examples();dump(out/'examples.json',ex)
    ref=exact_tests();dump(out/'exact_tests.json',ref)
    models={};reports={}
    for side in ('lower','upper'):
        p=pars[side];m=pm.reconstruct(side,F(p['metric_shift']));models[side]=m
        lo=pm.bound_deficit(m,F(p['lower_distance_limit']),-F(p['D_lower']))
        curve=bm.build_quadratic(pm.DATA,side)
        est=bm.estimate_strip(curve,F(p['upper_distance_limit']))
        up=bm.certify_constant(curve,est,F(p['completion_m']),F(p['D_upper']))
        dump(out/f'{side}_partner_coefficients.json',pm.coefficient_report(m))
        dump(out/f'{side}_finite_lower.json',lo);dump(out/f'{side}_finite_upper.json',up)
        Bstar=pm.val(m['deficit'],m['pref']['rawstar'])
        reports[side]={'C':m['C'],'D_lower':F(p['D_lower']),'D_upper':F(p['D_upper']),
                      'finite_joint_distance_limit':min(lo['distance_limit'],up['distance_limit']),
                      'endpoint_deficit_interval':Bstar,'endpoint_lower_interval':-Bstar,
                      'lower_certificate_margin':lo['strict_margin'],
                      'upper_certificate_margin':F(p['D_upper'])-up['unrounded_M_upper'],
                      'face_rows':len(m['rows']),'weighted_Gram_positive_pivots':len(m['gram_pivots']),
                      'regularized_choice_not_global_subleading_optimizer':True,
                      'new_lower_identity':'sum tr W(U-C/delta*K)^T G(U-C/delta*K)=C*(t-C/delta+deficit)>=0'}
        print(f"CERTIFIED {side}: {p['D_lower']} <= D_{side} <= {p['D_upper']}",flush=True)
        print(f"Finite-distance joint bounds valid for every 0<delta<={reports[side]['finite_joint_distance_limit']}",flush=True)
    neg=controls(models,pars);dump(out/'negative_controls.json',neg)
    forbidden=[m for m in ('numpy','scipy','sympy','mpmath') if m in sys.modules]
    if forbidden:raise AssertionError('third-party proof dependency')
    summary={
        'all_nonempty_principal_minors_count': structure['all_principal_minors_count'],
        'completed_entries': structure['complete_entries'],
        'examples': ex,
        'full_blocks': len(structure['blocks']),
        'negative_controls': neg,
        'raw_entries': structure['raw_entries'],
        'reference_tests': ref,
        'sides': reports,
    }
    dump(out/'summary.json',summary)
    print('All Stage-24 computational checks passed. Limit existence uses the written semialgebraic proof.',flush=True)
if __name__=='__main__':main()
