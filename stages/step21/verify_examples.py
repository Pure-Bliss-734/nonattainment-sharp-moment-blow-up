#!/usr/bin/env python3
"""Exact instance/identity checks accompanying a separately written theorem.

All order decisions use Fraction. The examples test hypotheses and catch
invalid generalizations; they are not a formal proof of the general theorem.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, json, random
from exact_algebra import *


def enc(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, dict): return {k:enc(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [enc(v) for v in x]
    return x


def check_allocation(weights, crosses, alphas):
    if len(weights) != len(crosses) or len(weights) != len(alphas):
        raise ValueError('length mismatch')
    for w, q, a in zip(weights, crosses, alphas):
        if a <= 0 or len(w) != len(q): raise ValueError('invalid completion')
        if psd(w)['rank'] != len(w): raise ValueError('weights must be PD in this test')
    C = sum((dot(q,mv(w,q))/a for w,q,a in zip(weights,crosses,alphas)), F(0))
    if C <= 0: raise ValueError('positive coefficient required')
    H = [scale(outer(q), 1/(a*C)) for q,a in zip(crosses,alphas)]
    if sum((pairing(w,h) for w,h in zip(weights,H)),F(0)) != 1:
        raise AssertionError('budget normalization')
    certs = [psd(border(h,q,a*C)) for h,q,a in zip(H,crosses,alphas)]
    # Exact identity for equality construction: border = vv^T/(a*C).
    for h,q,a in zip(H,crosses,alphas):
        v = list(q)+[a*C]
        if border(h,q,a*C) != scale(outer(v),1/(a*C)):
            raise AssertionError('rank-one identity')
    return {'C':C,'H0':H,'border_certificates':certs,'budget':F(1)}


def scalar_case(name, h, q, alpha, t, delta):
    cert = [psd(border([[a]],[b],c*t)) for a,b,c in zip(h,q,alpha)]
    thresholds = [b*b/(c*a) for a,b,c in zip(h,q,alpha)]
    expected = max([F(0)]+thresholds)
    if t != expected: raise AssertionError('not the exact optimum in this example')
    witness = None
    if t:
        i = thresholds.index(t); v = [-q[i]/h[i],F(1)]
        form = dot(v,mv(border([[h[i]]],[q[i]],alpha[i]*(t-1)),v))
        if form != -alpha[i]: raise AssertionError('failed exact lower witness')
        witness = {'block':i,'vector':v,'quadratic_at_t_minus_1':form}
    return {'name':name,'delta':delta,'t':t,'scaled_t':delta*t,
            'PSD_certificates':cert,'negative_witness':witness}


def examples():
    ds = [F(1,2),F(1,10),F(1,100),F(1,10**6)]
    samples = []; descriptors = {}
    for d in ds:
        samples.append(scalar_case('free_allocation', [d/5,4*d/5], [F(1),F(2)], [F(1)]*2, 5/d,d))
        samples.append(scalar_case('locked_unbalanced',[d/4,3*d/4],[F(1)]*2,[F(1)]*2,4/d,d))
        samples.append(scalar_case('locked_balanced',[d/2,d/2],[F(1)]*2,[F(1)]*2,2/d,d))
        samples.append(scalar_case('vanishing_crossing',[d],[d],[F(1)],d,d))
        samples.append(scalar_case('quadratic_budget',[d*d],[F(1)],[F(1)],1/(d*d),d))
        # Global branch z=0 is feasible; a chosen endpoint at z=1 has q*=1.
        samples.append(scalar_case('nonconcentrating_branch',[d],[F(0)],[F(1)],F(0),d))
        samples.append(scalar_case('duplicate_same_constraint',[d,d],[F(1)]*2,[F(1)]*2,1/d,d))
        # Curved Schur chart S1=z,S2=x-z-x^2, exact optimal allocation.
        budget=d-d*d; val=5/budget
        samples.append(scalar_case('curved_exact_chart',[budget/5,4*budget/5],[F(1),F(2)],[F(1)]*2,val,d))
        if val-5/d != 5/(1-d) or not 5 <= val-5/d <= 10:
            raise AssertionError('bounded constructive excess')
    A=matrix([[F(1,4)],[F(3,4)]])
    B=matrix([[F(1,2)],[F(1,2)]])
    if in_range(A,[F(1,2),F(1,2)]) or not in_range(B,[F(1,2),F(1,2)]):
        raise AssertionError('target accessibility')
    if rank(A)!=1 or rank(B)!=1: raise AssertionError('rank check')
    descriptors={
        'free_allocation': {'C':5,'true_scaled_limit':5,'formula_valid':True},
        'locked_unbalanced':{'C_from_valid_lower_certificate':2,'true_scaled_limit':4,
             'failed_upper_hypothesis':'H0 is not in range(J); rank J=1<2'},
        'locked_balanced':{'C':2,'true_scaled_limit':2,
             'observation':'rank J=1<2 but both needed targets are attainable; full row rank is not necessary'},
        'vanishing_crossing':{'C':0,'exact_tau':'delta','failed_growth_hypothesis':'C>0'},
        'quadratic_budget':{'exact_tau':'delta^-2','failed_hypothesis':'T/delta -> 1 and first-order normalization'},
        'nonconcentrating_branch':{'local_q_at_selected_endpoint':1,'global_tau':0,
             'failed_hypothesis':'uniform endpoint crossing concentration; z may remain zero'},
        'duplicate_same_constraint':{'incorrect_double_counted_C':2,'true_scaled_limit':1,
             'failed_hypothesis':'same PSD contribution cannot be counted twice'},
        'curved_exact_chart':{'C':5,'exact_tau':'5/(delta*(1-delta))',
             'bounded_excess':'tau-5/delta=5/(1-delta) in [5,10] for 0<delta<=1/2'} }
    # Coefficient identities clearing all denominators: not an asymptotic fit.
    # 5*x - (5 - 5*x)*x == 5*x^2, equivalent to the stated exact remainder.
    polycheck=padd([F(0),F(5)],pscale(pmul([F(5),F(-5)],[F(0),F(1)]),-1))
    if polycheck != [F(0),F(0),F(5)]: raise AssertionError('curved remainder identity')
    return {'examples':descriptors,'rational_sample_count':len(samples),'samples':samples,
            'symbolic_curved_remainder_coefficients':polycheck}


def high_rank_and_invariance():
    W1=matrix([[2,1],[1,3]]);W2=matrix([[1,0,0],[0,2,0],[0,0,4]])
    q1=[F(1),F(2)];q2=[F(1),F(-1),F(0)];alphas=[F(1),F(2)]
    out=check_allocation([W1,W2],[q1,q2],alphas)
    if out['C']!=F(39,2): raise AssertionError('higher-dimensional kernel prefactor')
    # Nonorthogonal change of kernel coordinates: q'=R^T q, W'=R^-1 W R^-T.
    R=matrix([[1,2],[-1,1]]);Ri=inverse(R)
    qp=mv(transpose(R),q1);Wp=mul(mul(Ri,W1),transpose(Ri))
    if dot(qp,mv(Wp,qp))!=dot(q1,mv(W1,q1)):raise AssertionError('coordinate invariance')
    regularized=[]
    for r in [F(1,10),F(1,1000),F(1,10**8)]:
        omega=trace(W1)+trace(W2); H=[add(scale(h,r),scale(eye(len(h)),r*r/omega)) for h in out['H0']]
        delta=r+r*r; t=out['C']/r
        for a,q,al in zip(H,[q1,q2],alphas): psd(border(a,q,al*t))
        if sum((pairing(w,h) for w,h in zip([W1,W2],H)),F(0))!=delta:raise AssertionError('regularized budget')
        if t-out['C']/delta != out['C']/(1+r):raise AssertionError('bounded excess identity')
        regularized.append({'r':r,'delta':delta,'t':t,'scaled_t':delta*t,
                            'excess_over_C_div_delta':t-out['C']/delta})
    out.update({'kernel_dimensions':[2,3],'alphas':alphas,'weights':[W1,W2],'crosses':[q1,q2],
                'kernel_coordinate_change':R,'changed_weight':Wp,'changed_cross':qp,
                'regularized_family_samples':regularized})
    return out


def random_identities():
    rng=random.Random(210907);n_test=48; inverse_coefficients=0
    for ix in range(n_test):
        n=1+ix%4
        M=matrix([[rng.randint(-3,3) for _ in range(n)] for _ in range(n)])
        P=add(mul(transpose(M),M),eye(n)); psd(P);Pi=inverse(P)
        q=[F(rng.randint(-3,3),rng.randint(1,3)) for _ in range(n)]
        if not any(q):q[0]=F(1)
        e=[F(rng.randint(-3,3),rng.randint(1,3)) for _ in range(n)]
        a=F(rng.randint(1,7),rng.randint(1,5))
        v=mv(Pi,q);Q=dot(q,v);b=dot(e,v);d=dot(e,mv(Pi,e))
        if Q<=0:raise AssertionError('positive Q')
        # S(r)[(Q+a*r)P^-1-vv^T] = r^2(Q+a*r)I.
        # Exact coefficient arrays, not evaluation at selected r.
        s_coeff=[zeros(n),scale(outer(q),1/a),P]
        n_coeff=[add(scale(Pi,Q),scale(outer(v),-1)),scale(Pi,a)]
        prod=[zeros(n) for _ in range(4)]
        for i,A in enumerate(s_coeff):
            for j,B in enumerate(n_coeff):prod[i+j]=add(prod[i+j],mul(A,B))
        want=[zeros(n),zeros(n),scale(eye(n),Q),scale(eye(n),a)]
        if prod!=want:raise AssertionError('cleared matrix inverse identity')
        inverse_coefficients+=4*n*n
        # Scalar expansion for f=q+r e:
        # f^T S^-1 f = a/r + (2ab+Qd-b^2-a^2+a d r)/(Q+a r).
        lhs=padd(pmul([Q,2*b,d],[Q,a]),pscale(pmul([Q,b],[Q,b]),-1))
        rhs=padd(pmul([F(0),a],[Q,a]),pmul([F(0),F(0),F(1)], [2*a*b+Q*d-b*b-a*a,a*d]))
        if lhs!=rhs:raise AssertionError('bounded remainder numerator identity')
        for r in [F(1,10),F(1,100),F(1,10**7)]:
            S=add(scale(outer(q),r/a),scale(P,r*r));f=[qi+r*ei for qi,ei in zip(q,e)]
            exact=dot(f,mv(inverse(S),f))
            rhsval=a/r+(2*a*b+Q*d-b*b-a*a+a*d*r)/(Q+a*r)
            if exact!=rhsval:raise AssertionError('inverse quadratic reference')
        # q=0 has no pole; evaluate independently instead of dividing by Q=0.
        for r in [F(1,7),F(1,10**5)]:
            f=[r*x for x in e];S=scale(P,r*r)
            if dot(f,mv(inverse(S),f))!=d:raise AssertionError('zero-crossing branch')
    return {'random_seed':210907,'matrix_cases':n_test,
            'cleared_inverse_coefficient_entries_checked':inverse_coefficients,
            'scalar_polynomial_identities':n_test,'inverse_quadratic_references':3*n_test,
            'zero_crossing_references':2*n_test}


def controls():
    rejected=[]
    def expect(name,fn):
        try: fn()
        except (AssertionError,ValueError,ZeroDivisionError):rejected.append(name)
        else:raise AssertionError('bad input passed: '+name)
    expect('negative_PSD_diagonal',lambda:psd([[1,0],[0,-1]]))
    expect('zero_diagonal_nonzero_crossing',lambda:psd([[0,1],[1,2]]))
    expect('singular_inverse',lambda:inverse([[1,2],[2,4]]))
    expect('nonpositive_completion_slope',lambda:check_allocation([[[F(1)]]],[[F(1)]],[F(0)]))
    expect('indefinite_weight',lambda:check_allocation([matrix([[1,2],[2,1]])],[[F(1),F(2)]],[F(1)]))
    expect('zero_C_in_positive_prefactor_branch',lambda:check_allocation([[[F(1)]]],[[F(0)]],[F(1)]))
    def bad_allocation():
        W=matrix([[2,1],[1,3]]);q=[F(1),F(2)];H=scale(outer(q),F(1,17))
        if pairing(W,H)!=1:raise AssertionError('false budget')
    expect('wrong_prefactor_breaks_budget',bad_allocation)
    def bad_rotation():
        W=matrix([[2,1],[1,3]]);q=[F(1),F(2)];R=matrix([[1,2],[0,1]])
        if dot(mv(transpose(R),q),mv(W,mv(transpose(R),q)))!=dot(q,mv(W,q)):
            raise AssertionError('forgot to transform W')
    expect('kernel_rescaling_without_weight_rescaling',bad_rotation)
    def bad_target():
        if not in_range(matrix([[F(1,4)],[F(3,4)]]),[F(1,2),F(1,2)]):
            raise AssertionError('unattainable target')
    expect('rank_deficient_unattainable_target',bad_target)
    def bad_count():
        delta=F(1,10);t=1/delta;psd([[delta,1],[1,t]])
        if delta*t<2: raise AssertionError('duplicated lower budget')
    expect('double_counting_single_PSD_block',bad_count)
    return {'rejected_count':len(rejected),'rejected':rejected}


def run(output_dir):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    reports={'examples':examples(),'higher_rank':high_rank_and_invariance(),
             'identity_checks':random_identities(),'negative_controls':controls()}
    for name,report in reports.items():
        (out/(name+'.json')).write_text(json.dumps(enc(report),indent=2,ensure_ascii=False)+'\n')
    summary={'all_step21_instance_checks_passed':True,
             'proof_boundary':'General statements are proved in the paper. These finite exact checks are regression and counterexample certificates, not formal theorem-prover certification.',
             'scenario_count':len(reports['examples']['examples'])+1,
             'scalar_example_points':reports['examples']['rational_sample_count'],
             'higher_rank_prefactor':reports['higher_rank']['C'],
             'identity_checks':reports['identity_checks'],
             'negative_controls':reports['negative_controls']['rejected_count'],
             'third_party_dependencies':[]}
    (out/'new_checks_summary.json').write_text(json.dumps(enc(summary),indent=2,ensure_ascii=False)+'\n')
    print('All Step-21 exact instance and polynomial-identity checks passed.')
    print('Nine scenarios; eight scalar scenarios at four rational parameters plus a higher-rank example.')
    print('The universal theorem is a written proof, not a claim proved by these samples.')
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',default='results/verified')
    run(parser.parse_args().output_dir)
