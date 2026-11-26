#!/usr/bin/env python3
"""Verify the feasible quadratic family and the upper remainder."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import argparse,json,subprocess,sys,copy
from source_tools import ROOT,prepare,SOURCE_SHA,PENCILS_SHA
from bounded_model import (I,encode,load_data,build_quadratic,estimate_strip,certify_constant,
                           sample_quadratic,sqrt_interval,scale)
from exact_checks import run as run_checks

def dump(path,obj):path.write_text(json.dumps(encode(obj),indent=2,ensure_ascii=False)+'\n')
def controls(data,curves,ests,params):
    rejected=[]
    def fail(label,fn):
        try:fn()
        except (AssertionError,ValueError,ZeroDivisionError):rejected.append(label)
        else:raise AssertionError('negative control accepted: '+label)
    for side in ['lower','upper']:
        c=curves[side];e=ests[side];p=params[side]
        fail(side+': m=0 insufficient for proposed scaled Schur certificate',lambda:certify_constant(c,e,0))
        fail(side+': M=m ignores positive distance conversion',lambda:certify_constant(c,e,F(p['completion_m']),F(p['completion_m'])))
        broken=copy.deepcopy(c);broken['pref']['C']=broken['pref']['C']*F(99,100)
        fail(side+': corrupted pole coefficient',lambda:certify_constant(broken,e,F(p['completion_m']),F(p['M'])))
        broken=copy.deepcopy(c);first=next(v for v in broken['projected'] if v['active']);first['HI']=scale(first['HI'],2)
        fail(side+': wrong second-order target',lambda:estimate_strip(broken,F(p['r_limit'])))
    fail('negative parameter strip',lambda:estimate_strip(curves['lower'],-1))
    fail('negative completion shift',lambda:certify_constant(curves['lower'],ests['lower'],-1))
    fail('negative square-root argument',lambda:sqrt_interval(I(-1)))
    sqrtchecks=0
    for i in range(65):
        v=F(i,19);q=sqrt_interval(I(v*v))
        if not q.contains(v):raise AssertionError('sqrt inclusion')
        sqrtchecks+=1
    return {'rejected_count':len(rejected),'rejected':rejected,'square_root_reference_checks':sqrtchecks,
            'scope':'A failed sufficient certificate is not a proof that a wider interval/family is infeasible.'}

def main():
    ar=argparse.ArgumentParser();ar.add_argument('--output-dir',default='results/verified')
    args=ar.parse_args();out=Path(args.output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    old,_,integrity=prepare();data=load_data();params=json.loads((ROOT/'inputs/parameters.json').read_text())
    checks=run_checks(data);dump(out/'exact_checks.json',checks)
    print('Exact parity, path-coefficient and Schur divided-difference checks passed.',flush=True)
    curves={};ests={};sides={}
    for side in ('lower','upper'):
        c=build_quadratic(data,side);curves[side]=c;p=params[side]
        prev=json.loads((old/'results/verified'/f'{side}_quadratic_direction.json').read_text())
        if encode(c['d2'])!=prev['d2']:raise AssertionError('recomputed Stage21 d2 differs')
        est=estimate_strip(c,F(p['r_limit']));ests[side]=est
        cert=certify_constant(c,est,F(p['completion_m']),F(p['M']))
        direction={'side':side,'C':c['pref']['C'],'Omega':c['pref']['Omega'],'d0':c['pref']['d0'],
                   'd2':c['d2'],'a':c['a'],'second_derivatives':c['curvatures'],
                   'matches_stage21_direction_bytes':True,'parameter':'delta=r+a*r^2, a>0; x=xi+sigma*delta',
                   'exact_direction_definition':'d0=J^T(JJ^T)^-1 H0; d2=J^T(JJ^T)^-1[HI-(1/2)D^2g(d0,d0)]',
                   'identity_basis':'Certified root equations and matrix inverse identities; zero-containing intervals are diagnostics only.'}
        samples=[sample_quadratic(data,c,cert,F(d)) for d in p['sample_deltas']]
        dump(out/f'{side}_direction.json',direction);dump(out/f'{side}_continuous_upper.json',cert)
        dump(out/f'{side}_full_samples.json',samples)
        rec={'distance_limit':cert['distance_limit'],'M':cert['M'],'completion_m':cert['completion_m'],
             'a':c['a'],'C':c['pref']['C'],'family_upper_excess_limit':cert['completion_m']+c['pref']['C']*c['a'],
             'unrounded_M_upper':cert['unrounded_M_upper'],
             'number_of_projected_blocks':len(cert['projected']),
             'active_blocks':sum(v['type']=='active' for v in cert['projected']),
             'minimum_active_second_order_margin':min(v['HI_minus_r_normR_lower'] for v in cert['projected'] if v['type']=='active'),
             'complete_blocks_per_sample':12,'sample_count':len(samples),
             'positive_sample_pivots':sum(len(b['pivots']) for s in samples for b in s['blocks']),
             'theorem':'tau(delta)<=C/delta+M for every positive delta up to the stated distance_limit'}
        sides[side]=rec
        print('CERTIFIED '+side+': tau(delta) <= C/delta + '+str(cert['M'])+' for all 0<delta<='+str(cert['distance_limit']),flush=True)
    bad=controls(data,curves,ests,params);dump(out/'negative_controls.json',bad)
    forbidden=[m for m in ['numpy','scipy','sympy','mpmath'] if m in sys.modules]
    if forbidden:raise AssertionError('third-party proof dependency: '+str(forbidden))
    summary={
        'sides': sides,
        'exact_checks': checks,
        'negative_controls': bad['rejected_count'],
    }
    dump(out/'summary.json',summary)
    print('All Step-22 exact checks passed.',flush=True)
if __name__=='__main__':main()
