#!/usr/bin/env python3
"""Verify the endpoint certificates with independent interval and AD routines."""
from pathlib import Path
from fractions import Fraction as F
import argparse,copy,json,sys,subprocess
from decimal_interval import Box,ZERO,ONE,display,dot,mul,sub,eye,ldlt,matrix_minor_pivots,DIGITS
from independent_model import Model,ROOT
from ad_kkt import evaluate_kkt,primal_dual_data
from independent_consequences import audit_consequences
from test_independent import self_tests
from second_order_ad import Jet

def dump(path,obj):path.write_text(json.dumps(obj,indent=2)+'\n')
def prove_root(model,data):
    center=list(map(F,data['center']));r=F(data['radius']);dim=20+data['nc']
    if len(center)!=dim or r<=0:raise AssertionError('invalid root box')
    box=[Box(q-r,q+r) for q in center]
    fc,_=evaluate_kkt(model,data,[Box(q) for q in center])
    _,jac=evaluate_kkt(model,data,box)
    c=[[Box(F(q)) for q in row] for row in data['preconditioner']]
    if len(c)!=dim or any(len(row)!=dim for row in c):raise AssertionError('invalid preconditioner shape')
    eta=max(dot(row,fc).norm() for row in c)
    residual=sub(eye(dim),mul(c,jac))
    rho=max(sum(q.norm() for q in row) for row in residual)
    if rho>=1:raise AssertionError('interval contraction failed')
    if eta+rho*r>=r:raise AssertionError('interval selfmap failed')
    infos,ranks=primal_dual_data(model,data,box)
    info={'method':'scalar-elimination second-order forward AD + fixed decimal interval Banach proof',
          'decimal_interval_digits':DIGITS,'KKT_dimension':dim,'primal_coordinates':20,'active_equations':data['nc'],
          'root_x':{'lower':str(center[0]-r),'upper':str(center[0]+r),'lower_decimal':display(center[0]-r,80),'upper_decimal':display(center[0]+r,80)},
          'box_radius':str(r),'eta_upper':str(eta),'rho_upper':str(rho),'eta_decimal':display(eta),'rho_decimal':display(rho),'selfmap_ratio':str(eta/r+rho),
          'unique_zero_of_all_KKT_equations_in_box':True,'rank_positivity':ranks,
          'function_and_jacobian_entries':dim+dim*dim}
    return info,box,infos

def rejected_controls(model):
    out=[]
    def reject(name,fn):
        try:fn()
        except (AssertionError,ValueError,ZeroDivisionError):out.append(name)
        else:raise AssertionError('bad input accepted: '+name)
    bad=json.loads((ROOT/'inputs/symmetry.json').read_text());bad['M0_-1']['S'][0][0]='2'
    reject('false_reflection',lambda:Model(cfg=bad))
    d=copy.deepcopy(model.data);d['raw_blocks'][-1]['mat'][0][0][1][0]=str(F(d['raw_blocks'][-1]['mat'][0][0][1][0])+1)
    reject('changed_Q_pencil_against_frozen_reduction',lambda:Model(data=d))
    for side in ['lower','upper']:
        b=json.loads((ROOT/f'inputs/kkt_{side}.json').read_text());b['center'][0]=str(F(b['center'][0])+F(1,10**8))
        reject(side+'_displaced_root',lambda b=b:prove_root(model,b))
    b=json.loads((ROOT/'inputs/kkt_lower.json').read_text());b['preconditioner']=[['0']*b['dimension'] for _ in range(b['dimension'])]
    reject('zero_preconditioner',lambda:prove_root(model,b))
    reject('indefinite_PSD_block',lambda:ldlt([[Box(0),Box(1)],[Box(1),Box(3)]]))
    reject('negative_dual_weight',lambda:ldlt([[Box(-1)]]))
    reject('dependent_optimal_face_rows',lambda:matrix_minor_pivots([[Box(1),Box(2)],[Box(1),Box(2)]]))
    reject('zero_crossing_denominator',lambda:Box(-1,1).reciprocal())
    # Deliberately corrupt the second-derivative chain rule, while keeping
    # values and first derivatives. The root-Jacobian proof must reject it.
    original=Jet.reciprocal
    def broken(self):
        j=original(self);return Jet(j.value,j.gradient,[ZERO]*len(j.hessian))
    Jet.reciprocal=broken
    try:
        good=json.loads((ROOT/'inputs/kkt_lower.json').read_text())
        reject('broken_reciprocal_second_derivative',lambda:prove_root(model,good))
    finally:Jet.reciprocal=original
    return {'count':len(out),'rejected':out,'all_rejected':True}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    arg = parser.parse_args()
    dest=arg.output_dir;dest.mkdir(parents=True,exist_ok=True)
    tests=self_tests();dump(dest/'oracle_tests.json',tests)
    model=Model();dump(dest/'symmetry.json',model.report)
    print('INDEPENDENT SYMMETRY: 11 -> 18 blocks; 213 reduced entries agree exactly.',flush=True)
    reports={}
    for side in ['lower','upper']:
        data=json.loads((ROOT/f'inputs/kkt_{side}.json').read_text())
        root,box,infos=prove_root(model,data)
        cons=audit_consequences(model,data,box,infos)
        dump(dest/f'{side}_root.json',root);dump(dest/f'{side}_consequences.json',cons)
        reports[side]={'root_x':root['root_x'],'rho_upper':root['rho_upper'],'rho_decimal':root['rho_decimal'],'dimension':root['KKT_dimension'],'primal_rank':root['rank_positivity']['primal_rank'],'dual_rank':root['rank_positivity']['dual_rank'],'global_exclusion_passed':True,'all_corrected_inner_blocks_positive':True,'full_optimal_face_rank':cons['optimal_face']['proved_rank'],'completion_cross':cons['nonattainment']['cross'],'nonattainment_recertified':True}
        print(side.upper(),'independent root certified; rho <=',root['rho_decimal'],flush=True)
        print('GLOBAL SIGNS, CURVATURE, FULL FACE RANK AND CROSS > 1/250 VERIFIED.',flush=True)
    controls=rejected_controls(model)
    dump(dest/'negative_controls.json',controls)
    forbidden=['dyadic_intervals','interval_kkt','endpoint_analysis','symmetric_model','sympy','numpy','scipy','mpmath']
    loaded=[s for s in forbidden if s in sys.modules]
    if loaded:raise AssertionError('independent proof imported forbidden old/backend modules: '+str(loaded))
    report={
        'both_endpoint_and_nonattainment_claims_recertified': True,
        'endpoints': reports,
        'negative_controls': controls,
        'oracle_tests': tests,
        'forbidden_modules_loaded': loaded,
    }
    dump(dest/'summary.json',report)
    print('Both endpoint and nonattainment claims independently re-certified.')
    print('All Step-16 independent exact checks passed.')
if __name__=='__main__':main()
