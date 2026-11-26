#!/usr/bin/env python3
"""Exact Schur completion, a nonextendible rational projection, and blowup.
No floating eigenvalues, SDP solver, SymPy, NumPy or external CAS are used.
"""
from __future__ import annotations
import argparse,json,subprocess,sys
from pathlib import Path
from fractions import Fraction as F
from completion_model import *
from verify_step12 import verify_layout,coordinates,rank_kernel,verify_state_equations
from verify_step11 import reconstruct_base
from verify_step10 import decimal

def dump(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
def form_json(form):return [[str(c) for c in p] for p in form]
def direction_pol(form,d):
    return linear([(v,(p,)) for v,p in zip(d,form[1:])],1)[0]

def verify_structure():
    b,P,C=model();state=reconstruct_base(UP/'inputs')
    layout,details=verify_layout(b,state)
    result={}
    for name,c in C.items():
        raw,ops,U,B=details[name]
        v=coordinates(c['W'],raw['basis']);basis=B+[v]
        rank,_=rank_kernel(transpose(U+basis),len(U+basis))
        check(rank==len(raw['basis'])==len(U)+len(basis),'Augmented basis does not span full quotient')
        check(direction_pol(c['diagonal'],RDIR)==(F(1),),'Completion diagonal does not grow as t')
        check(all(not any(direction_pol(f,RDIR)) for f in c['cross']),'Cross depends on t')
        for row in c['projected']['mat']:
            for f in row:check(not any(direction_pol(f,RDIR)),'Projected block depends on t')
        # Rebuild every augmented matrix entry from the raw Gram via basis coordinates.
        allops=ops+[c['W']]
        for i,vi in enumerate(basis):
            for j,vj in enumerate(basis):
                f=linear([(vi[k]*vj[l],en) for k,row in enumerate(raw['mat']) for l,en in enumerate(row) if vi[k]*vj[l]],22)
                exact=b['old_form'](wa.inner(allops[i],allops[j]))
                check(f==exact,'Raw congruence failed')
        result[name]={'W':operator_json(c['W']),'projected_basis':[operator_json(op) for op in ops],
            'raw_dimension':len(raw['basis']),'common_kernel_dimension':len(U),'projected_dimension':len(ops),
            'full_augmented_dimension':len(basis),'cross_forms':[form_json(f) for f in c['cross']],
            'diagonal_form':form_json(c['diagonal']),'raw_basis_change':[[str(v) for v in row] for row in basis]}
    for raw in b['blocks']:
        if raw['name'] in C:continue
        if raw['name']=='M4_None':
            check(raw['mat'][0][0][21]==(F(1),) and not any(any(p) for p in raw['mat'][0][0][:21]),'M4 is not [t]')
        else:
            for row in raw['mat']:
                for f in row:check(not any(direction_pol(f,RDIR)),'Unexpected varying raw block')
    print('Exact raw-to-completed block congruences passed.',flush=True)
    return result

def verify_interior(side):
    b,P,C=model();data=json.loads((UP/f'inputs/primal_{side}.json').read_text());x=F(data['x']);z=list(map(F,data['quotient_z']))
    cuts=thresholds(x,z);t=cuts['minimum'];raw=full_psd(x,z,t)
    active=next(name for name in C if cuts[name]['value']==t)
    y=cuts[active]['solution'];v=[-q for q in y]+[F(1)]
    G=augmented(active,x,z,t-1)
    check(dot(v,mv(G,v))==-1,'Below-threshold negative vector identity failed')
    check(t<F(data['completion']),'No improvement over previous completion')
    return {'side':side,'x':str(x),'old_completion':data['completion'],'minimum_completion':str(t),
            'minimum_display':decimal(t,35),'ceiling_completion':-(-t.numerator//t.denominator),
            'thresholds':{k:str(cuts[k]['value']) for k in C},'active_block':active,
            'exact_minimum_raw_blocks':raw,'total_rank_at_minimum':sum(a['rank'] for a in raw),
            'negative_vector_below_minimum':list(map(str,v)),
            'scope':'Minimum over the one completion parameter for this fixed projected moment point, not a minimum over all moments at this x.'}

def verify_singular_family():
    b,P,C=model();data=json.loads((UP/'inputs/primal_lower.json').read_text());x=F(data['x']);z0=list(map(F,data['quotient_z']))
    name='M2_None';c=C[name];A=evaluate_matrix(P[name],x,z0);n=len(A);u=[F(0)]*n;u[1]=1;u[2]=-1
    kappa=solve(A,u);a=dot(u,kappa);check(a>0,'Rank-one normalization not positive')
    theta=1/(3*a);direction=[F(-1)]+[F(0)]*20
    zstar=[v+theta*d for v,d in zip(z0,direction)]
    D=[[peval(en[1],x)*(-1) for en in row] for row in P[name]['mat']]
    check(D==[[-3*u[i]*u[j] for j in range(n)] for i in range(n)],'Target direction not -3 u u^T')
    B=evaluate_matrix(P[name],x,zstar)
    check(mv(B,kappa)==[F(0)]*n,'Exact null relation failed')
    pivs=[]
    for bl in P.values():
        pp=exact_psd_rank(evaluate_matrix(bl,x,zstar))
        check(len(pp)==len(bl['mat'])-int(bl['name']==name),'Unexpected boundary rank')
        pivs.append({'name':bl['name'],'rank':len(pp),'positive_pivots':list(map(str,pp))})
    v0=[evaluate_form(f,x,z0) for f in c['cross']];vstar=[evaluate_form(f,x,zstar) for f in c['cross']]
    check(v0==vstar,'Cross is not constant along this family')
    c0=evaluate_form(c['diagonal'],x,z0);check(c0==evaluate_form(c['diagonal'],x,zstar),'Constant diagonal changed')
    q=dot(kappa,vstar);check(q!=0,'No exact completion obstruction')
    try:solve(B,vstar)
    except ValueError:pass
    else:raise AssertionError('Singular point unexpectedly in the matrix range')
    # Independent direct-word zero norm and nonzero cross check.
    E=operator_linear(kappa,P[name]['basis']);ee=b['old_form'](wa.inner(E,E));ew=b['old_form'](wa.inner(E,c['W']))
    check(evaluate_form(ee,x,zstar)==0 and evaluate_form(ew,x,zstar)==q,'Original-word obstruction disagrees')
    check(not any(direction_pol(ee,RDIR)) and not any(direction_pol(ew,RDIR)),'Obstruction depends on completion t')
    # Raw nonlinear equality verification, including the whole free completion ray.
    gv=[evaluate_form(f,x,zstar) for f in b['gcoords']]
    state={w:md.affine({w:F(1)},b['g'])[0]+dot(md.affine({w:F(1)},b['g'])[1:],gv) for w in b['g']['forms']}
    delta={w:F(int(len(w)==8 and all(i<2 for i in w))) for w in state}
    eqcount=verify_state_equations(state,delta)
    coeff=q*q/a
    # Check the inverse formula as a POLYNOMIAL IDENTITY after clearing epsilon.
    ident=[[F(int(i==j)) for j in range(n)] for i in range(n)]
    invA=transpose([solve(A,col) for col in ident])
    U=[[u[i]*u[j]/a for j in range(n)] for i in range(n)]
    K=[[kappa[i]*kappa[j]/a for j in range(n)] for i in range(n)]
    # B(e) = (A-U)+e*U; e*inverse(B(e)) = K+e*(invA-K).
    A0=[[A[i][j]-U[i][j] for j in range(n)] for i in range(n)]
    V=[[invA[i][j]-K[i][j] for j in range(n)] for i in range(n)]
    H0=mm(A0,K);H1a=mm(A0,V);H1b=mm(U,K);H2=mm(U,V)
    check(all(v==0 for row in H0+H2 for v in row),'Cleared inverse constant/quadratic identity failed')
    check([[H1a[i][j]+H1b[i][j] for j in range(n)] for i in range(n)]==ident,'Cleared inverse linear coefficient failed')
    # The whole two-parameter (epsilon, completion) family satisfies original equations.
    gv0=[evaluate_form(f,x,z0) for f in b['gcoords']]
    de={w:dot(md.affine({w:F(1)},b['g'])[1:],[r-s for r,s in zip(gv0,gv)]) for w in state}
    for _,_,eq in wa.equations(9):
        total={}
        for mon,c_eq in eq.items():
            pol={(0,0):c_eq}
            for w in mon:
                new={}
                for (e,t),val in pol.items():
                    for e_delta,t_delta,fac in ((0,0,state[w]),(1,0,de[w]),(0,1,delta[w])):
                        key=(e+e_delta,t+t_delta);new[key]=new.get(key,F(0))+val*fac
                pol=new
            for key,val in pol.items():total[key]=total.get(key,F(0))+val
        check(not any(total.values()),'Two-parameter original equality failed')
    protected={}
    for key in ('x','p','a','d','j'):
        check(peval(b['obs'][key][1],x)==0,'Protected observable changes along the singular family')
        protected[key]=str(evaluate_form(b['obs'][key],x,zstar))
    check(evaluate_form(b['obs']['NQ'],x,zstar)>0,'Singular witness lost Q positivity')
    tbase=thresholds(x,z0)[name]['value'];regular=tbase-coeff
    # For both endpoint representatives, the OTHER augmented block is strict at t=478.
    # Affine interpolation and monotonicity in t give the same sufficient bound
    # throughout 0<=epsilon<=1. Thus the target controls the exact full minimum
    # whenever regular + coeff/epsilon exceeds 478.
    for zz in (z0,zstar):
        check(len(exact_psd_rank(augmented('M0_1',x,zz,F(478))))==7,'Other-block uniform completion bound failed')
    eps_cut=F(1,10**11)
    check(regular+coeff/eps_cut>478 and coeff>0,'Pole dominance bound failed')
    # Sherman-Morrison coefficients verified above, plus independent sample solves.
    samples=[]
    for estr in ('1','1/10000','1/100000000','1/1000000000000','1/10000000000000000','1/100000000000000000000'):
        eps=F(estr);z=[eps*v+(1-eps)*w for v,w in zip(z0,zstar)]
        ts=thresholds(x,z)
        check(ts[name]['value']==regular+coeff/eps,'Exact 1/epsilon completion formula failed')
        # Norm identity in original words, independently of the inverse formula.
        check(evaluate_form(ee,x,z)==eps*a and evaluate_form(ew,x,z)==q,'Family zero-norm scaling failed')
        t=ts['minimum'];raw=full_psd(x,z,t)
        samples.append({'epsilon':estr,'minimum_completion':str(t),'minimum_display':decimal(t,32),
             'target_threshold':str(ts[name]['value']),'target_display':decimal(ts[name]['value'],32),
             'active_block':max(C,key=lambda k:ts[k]['value']),
             'all_12_raw_blocks_psd':True,'total_positive_pivots':sum(r['rank'] for r in raw)})
    # Use a unit-last-component null vector for a manageable displayed cross value.
    kn=[v/kappa[-1] for v in kappa];qn=q/kappa[-1]
    generated={'x':str(x),'theta':str(theta),'zstar':list(map(str,zstar)),
        'null_vector':list(map(str,kappa)),'cross':str(q),'coefficient':str(coeff)}
    fixed=json.loads((ROOT/'inputs/singular_family.json').read_text())
    check(generated==fixed,'Frozen construction does not reproduce')
    # The omitted block can be made epsilon-close to PSD at the singular fibre:
    # its diagonal is t+c0 and b fixed. Choose eta=||b||^2/(t+c0)>0.
    # Adding eta I to B makes the Schur complement nonnegative; eta -> 0.
    check(dot(vstar,vstar)>0,'Unexpected zero cross vector')
    return {'x':str(x),'x_display':decimal(x,12),'theta':str(theta),'theta_display':decimal(theta,45),
        'target_block':name,'zstar':list(map(str,zstar)),'direction':list(map(str,direction)),
        'A':[[str(v) for v in row] for row in A],'u':list(map(str,u)),
        'kappa':list(map(str,kappa)),'a':str(a),'b':list(map(str,vstar)),
        'q_nonzero':str(q),'normalized_null_vector':list(map(str,kn)),
        'normalized_cross':str(qn),'normalized_cross_display':decimal(qn,40),
        'boundary_projected_blocks':pivs,'projected_total_positive_pivots':sum(v['rank'] for v in pivs),
        'original_factorized_equalities_for_fibre':eqcount,
        'two_parameter_original_equalities_checked':eqcount,
        'inverse_identity_coefficients_checked_exactly':True,
        'protected_observables':protected,
        'Q_at_singular_point':str(evaluate_form(b['obs']['NQ'],x,zstar)),
        'uniform_other_block_completion_cap':'478',
        'exact_pole_controls_full_minimum_for_epsilon_le':str(eps_cut),
        'null_operator_E':operator_json(E),'extension_operator_W':operator_json(c['W']),
        'completion_pole_coefficient':str(coeff),'pole_coefficient_display':decimal(coeff,40),
        'completion_regular_part':str(regular),'regular_part_display':decimal(regular,40),
        'family_formula':'target_threshold(epsilon) = regular_part + pole_coefficient/epsilon, 0 < epsilon <= 1',
        'samples':samples,
        'claims':{'all_projected_blocks_psd_at_star':True,'no_finite_completion_at_star':True,
                  'full_feasible_approximating_family_exists':True,'fixed_projection_fibre_weakly_infeasible':True,
                  'exact_projection_nonclosed_at_this_fixed_x':True,'is_radius_extremum':False},
        'witness_is_truncated_pseudomoments_not_physical_state':True}

def controls():
    rejected=0
    for A,b in [([[F(0)]],[F(1)]),([[F(1),0],[0,0]],[0,F(1)])]:
        try:solve(A,b)
        except ValueError:rejected+=1
    try:exact_psd_rank([[F(0),F(1)],[F(1),F(10000000000000000000000)]])
    except AssertionError:rejected+=1
    else:raise AssertionError('Zero diagonal/nonzero cross incorrectly accepted')
    check(rejected==3,'Negative controls failed')
    check(solve([[F(1),0],[0,0]],[F(2),0])==[F(2),0],'Consistent singular solve failed')
    exact_psd_rank([[F(1),0,F(2)],[0,0,0],[F(2),0,F(4)]])
    return {'negative_controls_rejected':3,'singular_range_positive_control':True}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    geometry=verify_structure();dump(args.output_dir/'completion_structure.json',geometry)
    interiors={side:verify_interior(side) for side in ('lower','upper')};dump(args.output_dir/'interior_minima.json',interiors)
    for side,item in interiors.items():print(side,'minimum eighth moment =',item['minimum_display'],flush=True)
    family=verify_singular_family();dump(args.output_dir/'nonclosed_projection.json',family)
    tests=controls();dump(args.output_dir/'negative_controls.json',tests)
    summary={
        'all_projected_blocks_psd': True,
        'certified_nonextendible_projection_x': family['x'],
        'exact_nonzero_null_cross': family['normalized_cross_display'],
        'exact_projection_nonclosed': True,
        'explicit_positive_1_over_epsilon_coefficient': family['pole_coefficient_display'],
        'finite_full_family_samples_verified': len(family['samples']),
        'fixed_projected_interior_completion_minima': {side:d['minimum_display'] for side,d in interiors.items()},
        'two_exact_completion_range_conditions': True,
    }
    dump(args.output_dir/'summary.json',summary)
    print('CERTIFIED: finite completion iff both cross vectors lie in the corresponding ranges.',flush=True)
    print('CERTIFIED: nonclosed moment projection at fixed x = 0.3855596.',flush=True)
    print('CERTIFIED: positive exact 1/epsilon divergence of the eighth moment.',flush=True)
    print('Radius-extremum attainment remains unproved; no numerical KKT candidate is used as proof.',flush=True)
    print('All Step-13 exact checks passed.',flush=True)
if __name__=='__main__':main()
