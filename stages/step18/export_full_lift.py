#!/usr/bin/env python3
"""Rebuild every neutral level-9 moment and verify the original equations.

This process deliberately shares the historical word algebra, independently of
new interval bounds. It exports a complete algebraic construction, not a flag
asserting that any particular point is feasible.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,subprocess,argparse
from frozen_source import prepare,ROOT
P17,P16,_,_=prepare()
old=ROOT.parent/'step15'
sys.path.insert(0,str(old))
import load_frozen
from step12_model import build,linear
from completion_model import RDIR
from step12_model import md,wa

def add(a,b):
    o=a.copy()
    for k,v in b.items():
        o[k]=o.get(k,F(0))+v
        if not o[k]:del o[k]
    return o

def scale(a,c):return {k:v*c for k,v in a.items() if v*c}

def mul(a,b):
    o={}
    for m,c in a.items():
        for n,d in b.items():
            k=tuple(sorted(m+n));o[k]=o.get(k,F(0))+c*d
    return {k:v for k,v in o.items() if v}

def poly(form):
    p={}
    for i,row in enumerate(form):
        for power,c in enumerate(row):
            if not c:continue
            mon=(0,)*power+((i,) if i else ())
            p=add(p,{tuple(sorted(mon)):c})
    return p

def encode(p):return [{'variables':list(k),'coefficient':str(v)} for k,v in sorted(p.items())]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',default=str(ROOT/'results/verified'))
    arg=ap.parse_args();out=Path(arg.output_dir);out.mkdir(parents=True,exist_ok=True)
    b=build();g=b['g'];md.verify_all_rows(g)
    moments={}
    for w in [()]+g['words']:
        av=md.affine({w:F(1)},g)
        lifted=linear([(av[0],b['obs']['one'])]+list(zip(av[1:],b['gcoords'])),22)
        moments[w]=poly(lifted)
    count=0
    for _,_,equation in wa.equations(9):
        value={}
        for mon,c in equation.items():
            term={():c}
            for w in mon:term=mul(term,moments[w])
            value=add(value,term)
        if value:raise AssertionError('unlifted original factorization equation failed')
        count+=1
    if count!=5417 or len(moments)!=1575:raise AssertionError('unexpected model size')
    if moments[(1,1,1,1,0,0,0,0)]!={(21,):F(1)}:
        raise AssertionError('eighth moment is not z20')
    r=list(RDIR)
    for name in ('x','p','a','j','d','NQ'):
        p=poly(b['obs'][name]);der={}
        for mon,c in p.items():
            for k in range(1,22):
                if k not in mon or not r[k-1]:continue
                n=list(mon);n.remove(k)
                der=add(der,{tuple(n):c*r[k-1]*mon.count(k)})
        if der:raise AssertionError('completion direction changed '+name)
    report={
      'coordinate_order':['x']+[f'z{k}' for k in range(21)],
      'word_alphabet':{'0':'Z','1':'Zbar','2':'Pi','3':'Pibar'},
      'construction':'Given symmetric family y, recover z2=-a/2, set z20=0, then add t*ray; evaluate these polynomials. The 10 unused level-9 directions are set to zero by the exact right inverse.',
      'ray':list(map(str,r)),
      'moments':[{'word':list(w),'polynomial':encode(p)} for w,p in moments.items()],
      'all_original_factorization_equations_verified_identically':count,
      'identity_domain':'all real x,z0,...,z20; no positivity assumed in equality verification',
      'invisible_directions':len(b['lineality']),
      'moment_t_is_exactly_z20':True,
      'observables_invariant_under_completion_ray':['x','p','a','j','d','NQ'],
      'shared_word_algebra_not_an_independent_reimplementation':True}
    (out/'full_moment_map.json').write_text(json.dumps(report,indent=2)+'\n')
    brief={k:v for k,v in report.items() if k!='moments'};brief['moments_including_identity']=len(moments)
    (out/'lift_identity.json').write_text(json.dumps(brief,indent=2)+'\n')
    print(f'FULL LIFT: {len(moments)} moments; {count} original equations hold as exact multivariate polynomial identities.')
if __name__=='__main__':main()
