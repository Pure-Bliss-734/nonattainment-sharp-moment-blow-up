#!/usr/bin/env python3
"""Reconstruct the sharp prefactor and continuous interval inequalities."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import argparse,json,subprocess,sys,time
from audit_interval import I,overlap
from audit_prefactor import reconstruct,certify_upper,certify_lower,sample_full
from audit_tests import reference_tests,controls
ROOT=Path(__file__).resolve().parent

def enc(x):
    if isinstance(x,I):return x.record()
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [enc(v) for v in x]
    return x

def dump(path,obj):path.write_text(json.dumps(enc(obj),indent=2,ensure_ascii=False)+'\n')
def dec(x,n=45):
    with localcontext() as c:
        c.prec=n;return str(Decimal(x.numerator)/Decimal(x.denominator))

def compare_to_previous(root,pref,side):
    # Interval overlap is a diagnostic, not an equality proof. Algebraic
    # equality follows from the independently derived Schur derivative identity.
    q=json.loads((root/f'results/verified/{side}_prefactor.json').read_text());counts={}
    def old(v):return I(F(v['lower']),F(v['upper']))
    assert overlap(pref['C'],old(q['C']));counts['C']=1
    for key in ['J','R']:
        prev=q['J' if key=='J' else 'right_inverse'];current=pref[key]
        assert len(current)==len(prev) and len(current[0])==len(prev[0]);n=0
        for row,pr in zip(current,prev):
            for a,b in zip(row,pr):assert overlap(a,old(b)),key+' enclosure mismatch';n+=1
        counts[key]=n
    for key in ['d0','dI']:
        assert all(overlap(a,old(b)) for a,b in zip(pref[key],q[key]));counts[key]=len(pref[key])
    for name,cp in pref['completion'].items():
        for key in ['cross','weight','contribution']:assert overlap(cp[key],old(q['completion'][name][key]))
    counts['completion_scalars']=6
    return counts

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output-dir',default='results/verified');args=ap.parse_args()
    out=Path(args.output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    from hashlib import sha256
    payload=(ROOT/'inputs/shared_pencils.json').read_bytes()
    assert sha256(payload).hexdigest()=='35f905894a2c61f6488bb0480c6c6551b229f63b5800a5bf92833dde277ba39e','shared input fingerprint mismatch'
    data=json.loads(payload)
    tests=reference_tests(data);dump(out/'reference_tests.json',tests)
    prefs={};summary={}
    old=ROOT.parent/'step19'
    for side,u0 in [('lower',F(1,10**9)),('upper',F(1,10**6))]:
        pref=reconstruct(data,side);prefs[side]=pref
        upper=certify_upper(data,pref,u0);lower=certify_lower(data,pref)
        compare=compare_to_previous(old,pref,side)
        samples=[sample_full(data,pref,F(1,10**k)) for k in [20,24,28]]
        newlo=[]
        for exponent in [20,24,28]:
            u=F(1,10**(exponent//2));delta=u*u
            numerator=pref['C'].lo()-lower['e_coefficients'][1]*u-lower['e_coefficients'][2]*delta
            denom=1+lower['denominator_coefficients'][1]*u+lower['denominator_coefficients'][2]*delta
            low=numerator/denom-delta;high=pref['C'].hi()*(1+2*u)
            assert low<=high
            # On the negative-numerator branch this lower bound is trivial
            # from t>=0, not obtained by dividing a signed inequality.
            newlo.append({'delta':delta,'lower':low,'upper':high,
                          'lower_display':dec(low),'upper_display':dec(high),
                          'positive_numerator':numerator>0})
        # A readable finite error rate follows without claiming a next coefficient.
        e1,e2=lower['e_coefficients'][1:];b1,b2=lower['denominator_coefficients'][1:]
        # C-δτ <= E+C*(den-1)+δ provided C-E>=0; if negative, same bound
        # follows from δτ>=0 since E>=C. The denominator is at least one.
        m=e1+pref['C'].hi()*b1+(e2+pref['C'].hi()*b2+1)*u0
        rate=max(m,2*pref['C'].hi());rate_int=(rate.numerator+rate.denominator-1)//rate.denominator
        dump(out/f'{side}_independent_prefactor.json',{
            'C':pref['C'],'Omega':pref['Omega'],'completion':pref['completion'],'J':pref['J'],
            'right_inverse':pref['R'],'gram_pivots':pref['gram_pivots'],'d0':pref['d0'],'dI':pref['dI'],
            'different_derivative_method':'fixed-kernel compression K^T(dM)K, no automatic differentiation',
            'baseline_overlap_diagnostics_not_equality_proofs':compare})
        dump(out/f'{side}_continuous_upper.json',upper);dump(out/f'{side}_universal_lower.json',lower)
        dump(out/f'{side}_full_samples.json',samples);dump(out/f'{side}_finite_brackets.json',newlo)
        summary[side]={'C':pref['C'],'C_display':dec(pref['C'].lo(),60),'u_limit':u0,'delta_limit':u0*u0,
                       'Schur_jacobian_rank':len(pref['J']),
                       'all_path_lower_and_continuous_upper_reestablished':True,
                       'completion_Schur_margin_lowers':{k:v['Schur_over_u_lower'] for k,v in upper['completion'].items()},
                       'new_error_numerator_degree_in_sqrt_delta':2,'common_diagonal_upper':1,
                       'scaled_error_bound_integer':rate_int,
                       'scaled_error_bound':'|delta*tau(delta)-C| <= M*sqrt(delta) on the stated constructed strip; M is conservative',
                       'full_sample_blocks':sum(len(v['blocks']) for v in samples),
                       'sample_positive_pivots':sum(len(b['pivots']) for v in samples for b in v['blocks'])}
        print(side.upper()+': C='+summary[side]['C_display']+'; new continuous Schur checks pass.',flush=True)
        print('  all 0<delta<= '+str(u0*u0)+'; quantitative modulus M='+str(rate_int),flush=True)
    bad=controls(data,prefs);dump(out/'negative_controls.json',bad)
    forbidden=[m for m in ['prefactor_model','order_model','ad_kkt','decimal_interval','second_order_ad','sympy','numpy','scipy','mpmath'] if m in sys.modules]
    assert not forbidden,forbidden
    report={
        'reference_tests': tests,
        'negative_controls_count': bad['count'],
        'sides': summary,
        'forbidden_modules_imported': forbidden,
    }
    dump(out/'summary.json',report)
    print('All Stage-20 independent prefactor-layer audit checks passed.',flush=True)
if __name__=='__main__':main()
