"""Exact interval consequences of the certified KKT roots.

The roots are algebraic implicitly: isolated zeros of the displayed rational
KKT system. No univariate eliminant or minimal-polynomial degree is claimed.
"""
from fractions import Fraction as F
import json
import load_frozen
from dyadic_intervals import *
from interval_kkt import KKTIntervals
from symmetric_model import KEEP,small_blocks
from step12_model import build
from step12_layout import projected
from completion_model import model

def eval_form(form,x,z):
    xx=[I(1),x,x*x];co=[I(1)]+z
    return sum((co[k]*sum((I(q)*xx[p] for p,q in enumerate(poly)),I(0)) for k,poly in enumerate(form)),I(0))
def reconstruct_z(y):
    z=[I(0)]*21
    for i,j in enumerate(KEEP):z[j]=y[i+1]
    z[2]=-eval_form(build()['obs']['a'],y[0],z)/2
    return z

def certify_consequences(side,data,prob,box,infos):
    y=box[:20];x=y[0];b=build();z=reconstruct_z(y);obs={k:eval_form(f,x,z) for k,f in b['obs'].items()}
    gamma=next(d['multiplier'][0][0] for d in infos if d['item']['name']=='single_Q_plus')
    L=F(7,2)*obs['p']+3*obs['a'];Lmax=F(109,20)
    # In all non-Q blocks the variable coefficients are x independent;
    # in Q the x-derivative of these coefficients equals L_z.
    const=[I(0)]*3
    for d in infos:
        P=prob.polys[d['item']['index']]
        Y=d['dual']
        for power in range(3):
            const[power]=const[power]+ip(Y,[[I(v) for v in row] for row in P[power][0]])
    # The constant part of L in the section z20=0 is (39/2)*X.
    L0slope=F(39,2)
    F2=const[2]-gamma*L0slope
    sig=1 if side=='lower' else -1
    Fprime=I(sig)-gamma*L  # exact at the root, by the x stationarity equation
    broad=I(F(3,8),F(2,5))
    derivative=Fprime+2*F2*(broad-x)
    if side=='lower':
        assert derivative.lo>0
        globalderivative=derivative
    else:
        globalderivative=derivative+gamma*Lmax
        assert globalderivative.hi<0
    lipschitz=derivative.absupper()+gamma.upper()*Lmax
    assert lipschitz<2
    glob={'gamma':record(gamma),'L_at_root':record(L),'F_quadratic_coefficient':record(F2),'F_derivative_on_radius_interval':record(derivative),'exclusion_derivative':record(globalderivative),'sum_bound_lipschitz_upper':str(lipschitz),'sum_bound_lipschitz_upper_decimal':dec(lipschitz,25),'sum_bound_lipschitz_less_than_two':True}
    # An interior complete point was certified at Step 12. Retaining its other
    # coordinates and reflecting z2 gives a strictly positive symmetric point.
    interior=json.loads((load_frozen.UP12/f'inputs/primal_{side}.json').read_text())
    xin=F(interior['x']);zin=list(map(F,interior['quotient_z']))
    yi=[I(xin)]+[I(zin[j]) for j in KEEP]
    if side=='lower':assert I(xin).lo>x.hi
    else:assert I(xin).hi<x.lo
    dx=yi[0]-x;dz=[yi[j]-y[j] for j in range(1,20)];pivots=[]
    for ib,bl in enumerate(prob.blocks):
        P=prob.polys[ib];n=len(bl['mat']);Mi=prob.matrices(ib,yi)[0]
        for k in range(1,20):assert all(not q for row in P[2][k] for q in row)
        H=zmat(n)
        for i in range(n):
            for j in range(n):H[i][j]=dx*dx*P[2][0][i][j]+dx*sum((dz[k-1]*P[1][k][i][j] for k in range(1,20)),I(0))
        pi=pd_pivots(Mi);pc=pd_pivots(sub(Mi,H))
        pivots.append({'block':bl['name'],'interior_pivots':[record(p) for p in pi],'corrected_pivots':[record(p) for p in pc]})
    reachable={'interior_x':str(xin),'all_18_blocks_strictly_positive_after_curvature_correction':True,'positive_pivots_each':sum(len(v['interior_pivots']) for v in pivots),'blocks':pivots}
    # Lift every small active kernel to its unsymmetrized parent block.
    # Positive definite kernel multipliers force the FULL parent matrix to
    # annihilate these columns at every optimal point.
    raw={a['name']:a for a in projected()};coeffrows=[];labels=[];K_by_name={}
    for d in infos:
        sb=prob.blocks[d['item']['index']];ob=raw[sb['parent']];T=[[I(c) for c in row] for row in sb['T']];K=mm(T,d['K']);K_by_name[sb['name']]=K
        nr=len(ob['mat']);nu=len(K[0]);A=[]
        for k in range(20):
            A.append([[sum((I(q)*x**power for power,q in enumerate(ob['mat'][i][j][k+1])),I(0)) for j in range(nr)] for i in range(nr)])
        AK=[mm(Ak,K) for Ak in A]
        for i in range(nr):
            for j in range(nu):coeffrows.append([a[i][j] for a in AK]);labels.append([sb['name'],i,j])
    selection=json.loads((load_frozen.ROOT/f'inputs/face_rows_{side}.json').read_text());assert labels==selection['row_labels']
    rows=selection['indices'];assert len(set(rows))==20
    pp=nonsingular_pivots([coeffrows[i] for i in rows]);assert len(pp)==20
    unique={'annihilation_rows':len(coeffrows),'rank':20,'selected_rows':rows,'labels':[labels[i] for i in rows],'nonsingular_pivots':[record(v) for v in pp],'unique_full_projected_optimal_point':True}
    # The M0_1 block is one of the two finite-completion blocks.
    _,_,comp=model();K=K_by_name['M0_1_plus'];assert len(K[0])==1
    nu=[row[0] for row in K];cross=[eval_form(f,x,z) for f in comp['M0_1']['cross']]
    q=dot(nu,cross);assert q.lower()>F(1,250)
    lam=next(d['multiplier'][0][0] for d in infos if d['item']['name']=='M0_1_plus');assert lam.lower()>F(1,250)
    # 2*|x-xi| >= sum >= lambda*rho, and near the unique optimizer
    # s > 1/500, so t+c0 > 1/(125000000*|x-xi|).
    nonattain={'block':'M0_1','kernel_vector':[record(a) for a in nu],'cross':record(q),'cross_exceeds':'1/250','kernel_dual_weight':record(lam),'kernel_dual_weight_exceeds':'1/250','finite_completion_impossible':True,'eventual_radius_gap_divergence_constant':'1/125000000','divergence_is_for_every_full_feasible_sequence_approaching_this_radius_endpoint':True,'neighborhood_radius_not_quantified':True}
    return {'global_exclusion':glob,'inner_approximation':reachable,'optimal_face':unique,'nonattainment':nonattain,'observables':{k:record(v) for k,v in obs.items()}}
