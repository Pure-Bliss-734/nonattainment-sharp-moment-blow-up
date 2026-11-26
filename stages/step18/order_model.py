"""Construct and certify a 1/distance upper family, without numerical solvers.

All intervals enclose the same isolated algebraic endpoint. Equalities are
checked as rational polynomial identities; interval bounds are only used for
strict inequalities. The Step-16 model and interval library are shared premises.
"""
from fractions import Fraction as F
from decimal_interval import Box,ZERO,ONE,ldlt,solve_spd,dot,mul,display
from independent_model import (ROOT,KEEP,from_form,evaluate,derivative,
                               pc,pv,plus,times,scaling,substitute)
import json


def square(a):
    """An enclosure for a*a that preserves nonnegativity."""
    lo,hi=a.lo(),a.hi()
    return Box(0 if lo<=0<=hi else min(lo*lo,hi*hi),max(lo*lo,hi*hi))


def curvature(poly,direction):
    if any(len(m)>2 for m in poly):raise AssertionError('polynomial degree exceeds two')
    return sum((c*direction[m[0]]*direction[m[1]] for m,c in poly.items() if len(m)==2),ZERO)


def path_coefficients(poly,star,direction):
    return (evaluate(poly,star),
            sum((evaluate(derivative(poly,i),star)*direction[i] for i in range(len(star))),ZERO),
            curvature(poly,direction))


def bernstein(a):return [a[0],a[0]+a[1]/2,sum(a,ZERO)]

def eval_path(a,s):return a[0]+s*(a[1]+s*a[2])


def fraction_raw_coordinates(model,y):
    return [evaluate(model.to_symmetric[k],list(y)+[F(0)]) for k in range(21)]


def verify_formal_identities(model,interior):
    """Check secant/Bernstein identities with the endpoint still a formal variable."""
    S=pv(20);O=plus(pc(1),scaling(S,-1))
    substitutions={i:plus(times(O,pv(i)),scaling(S,c)) for i,c in enumerate(interior)}
    d=[plus(pc(c),scaling(pv(i),-1)) for i,c in enumerate(interior)]
    secant=0
    for b in model.small:
        for row in b['mat']:
            for p in row:
                if any(len(m)>2 for m in p):raise AssertionError('nonquadratic pencil')
                h={}
                for mon,c in p.items():
                    if len(mon)==2:h=plus(h,scaling(times(d[mon[0]],d[mon[1]]),c))
                rhs=plus(plus(times(O,p),scaling(S,evaluate(p,interior))),scaling(times(times(S,O),h),-1))
                if substitute(p,substitutions)!=rhs:raise AssertionError('secant identity failed')
                secant+=1
    raw_to_y={k:{m:c for m,c in p.items() if 20 not in m} for k,p in model.to_symmetric.items()}
    bc=0
    for cp in model.data['completion'].values():
        for form in cp['cross']+[cp['diagonal']]:
            p=substitute(from_form(form),raw_to_y)
            path=substitute(p,substitutions);powers=[{}, {}, {}]
            for mon,c in path.items():
                j=mon.count(20)
                if j>2:raise AssertionError('nonquadratic crossing')
                mm=tuple(k for k in mon if k!=20)
                powers[j]=plus(powers[j],{mm:c})
            b0=powers[0];b1=plus(b0,scaling(powers[1],F(1,2)))
            b2=plus(plus(b0,powers[1]),powers[2])
            rhs=plus(plus(times(times(O,O),b0),scaling(times(times(S,O),b1),2)),times(times(S,S),b2))
            if rhs!=path:raise AssertionError('Bernstein identity failed')
            bc+=1
    return {'exact_secant_entry_identities':secant,'exact_Bernstein_entry_identities':bc,
            'formal_endpoint_variables':20,'interval_equality_not_used_as_identity_proof':True}


def get_family(model,data,box):
    side=data['side'];star=box[:20]
    record=json.loads((ROOT/f'inputs/interior_{side}.json').read_text())
    inside_f=[F(record['x'])]+[F(record['quotient_z'][k]) for k in KEEP]
    inside=[Box(c) for c in inside_f];direction=[v-u for u,v in zip(star,inside)]
    gap=direction[0] if side=='lower' else -direction[0]
    if gap.lo()<=0:raise AssertionError('interior is not on the required side')
    raw_star=model.raw_coordinates(star)
    raw_in=fraction_raw_coordinates(model,inside_f)
    raw_dir=[Box(v)-u for u,v in zip(raw_star,raw_in)]
    return {'star':star,'inside':inside,'inside_f':inside_f,'direction':direction,
            'gap':gap,'raw_star':raw_star,'raw_inside':raw_in,'raw_direction':raw_dir}


def certify_family(model,data,box,infos,parameters):
    side=data['side'];family=get_family(model,data,box)
    theta=F(parameters['theta']);K=F(parameters[side]['K'])
    radius=F(parameters['radius']);rounded=F(parameters[side]['rounded_C'])
    if not (0<theta<1 and K>0 and 0<radius<family['gap'].lo()):
        raise AssertionError('invalid upper-family parameters')
    exact=verify_formal_identities(model,family['inside_f'])
    pd=[]
    # This proves M(s) >= theta*s*M_in for ALL s in (0,1], not only samples.
    for i,b in enumerate(model.small):
        a=model.matrix(i,family['inside']);h=[[curvature(p,family['direction']) for p in row] for row in b['mat']]
        _,pa=ldlt(a)
        margin=[[(1-theta)*q-v for q,v in zip(row,hr)] for row,hr in zip(a,h)]
        _,pm=ldlt(margin)
        pd.append({'name':b['name'],'interior_positive_pivots':pa,'common_lower_matrix_positive_pivots':pm})
    raw_by_name={b['name']:b for b in model.raw}
    completion={};poles={}
    for name,cp in model.data['completion'].items():
        raw=raw_by_name[name]
        projected=[[from_form(f) for f in row] for row in cp['projected']['mat']]
        if projected!=raw['mat']:raise AssertionError('completion pencil does not match raw pencil')
        a=[[Box(evaluate(p,family['raw_inside'])) for p in row] for row in raw['mat']]
        bp=[path_coefficients(from_form(f),family['raw_star'],family['raw_direction']) for f in cp['cross']]
        c=path_coefficients(from_form(cp['diagonal']),family['raw_star'],family['raw_direction'])
        controls=[[bernstein(coef)[j] for coef in bp] for j in range(3)]
        cc=bernstein(c)
        if min(q.lo() for q in cc)<=0:raise AssertionError('completion diagonal is not uniformly positive')
        qrecords=[];qmax=F(0)
        for v in controls:
            sol,piv=solve_spd(a,[[q] for q in v])
            quad=dot(v,[row[0] for row in sol]);qmax=max(qmax,quad.hi())
            qrecords.append({'control_vector':v,'inverse_quadratic_enclosure':quad,'LDL_pivots':piv})
        cap=F(parameters[side]['metric_caps'][name])
        if qmax>=cap:raise AssertionError('displayed rational metric cap is not proved')
        if theta*K<=cap:raise AssertionError('K too small for completion bound')
        completion[name]={'dimension':len(a),'inverse_quadratic_controls':qrecords,
                          'quadratic_upper':qmax,'rounded_metric_cap':cap,'theta_K_minus_metric_cap':theta*K-cap,
                          'theta_K_minus_quadratic_upper':theta*K-qmax,
                          'diagonal_Bernstein_controls':cc,'cross_power_coefficients':bp,
                          'diagonal_power_coefficients':c}
        matches=[inf for inf in infos if model.small[inf['item']['index']]['parent']==name]
        if len(matches)!=1 or matches[0]['item']['nu']!=1:raise AssertionError('simple nullspace required for pole formula')
        inf=matches[0];sb=model.small[inf['item']['index']]
        nu=[r[0] for r in mul([[Box(q) for q in row] for row in sb['T']],inf['kernel'])]
        slope_matrix=[[path_coefficients(p,family['raw_star'],family['raw_direction'])[1] for p in row] for row in raw['mat']]
        slope=dot(nu,[dot(row,nu) for row in slope_matrix])
        if slope.lo()<=0:raise AssertionError('kernel eigenvalue does not open to first order')
        q0=dot(nu,[coef[0] for coef in bp]);residue=square(q0)/slope
        # At the exact root, the leading PSD Schur block has one-dimensional
        # kernel. The scalar kernel Schur complement has derivative slope.
        # Consequently s*B(s)^(-1) -> nu*nu^T/slope, including q0=0 cases.
        poles[name]={'kernel':nu,'kernel_dimension':1,'kernel_slope':slope,
                     'cross_at_endpoint':q0,'s_times_minimal_completion_limit':residue,
                     'distance_times_minimal_completion_limit':residue*family['gap']}
    coeff=K*family['gap']
    if coeff.hi()>=rounded:raise AssertionError('rounded upper coefficient is too small')
    if poles['M0_1']['distance_times_minimal_completion_limit'].lo()<=max(F(0),poles['M2_None']['distance_times_minimal_completion_limit'].hi()):
        raise AssertionError('dominant completion block not separated')
    return {'side':side,'theta':theta,'K':K,'radius':radius,'interior_parameter':family['inside_f'][0],
            'root_to_interior_distance':family['gap'],'y_inside':family['inside_f'],
            'root_y_enclosures':family['star'],'exact_identities':exact,
            'all_small_common_positive_matrices':pd,'completion_certificates':completion,
            'explicit_family':'y(s)=(1-s)*y_star+s*y_inside, 0<s<=1; t=K/s',
            'exact_upper_coefficient_enclosure':coeff,'rounded_upper_coefficient':rounded,
            'eighth_moment_upper':'t = K*abs(x_inside-xi)/abs(x-xi) < rounded_C/abs(x-xi)',
            'minimal_completion_along_this_path_poles':poles,'asymptotic_dominant_block':'M0_1',
            'path_optimal_prefactor_is_not_claimed_globally_optimal':True}


def certify_samples(model,data,box,params):
    """Algebraic sample coordinates are intervals, not asserted rational points."""
    fam=get_family(model,data,box);K=F(params[data['side']]['K']);out=[]
    for text in params['sample_distances']:
        distance=F(text)
        if not 0<distance<fam['gap'].lo():raise AssertionError('sample outside the interpolation interval')
        s=Box(distance)/fam['gap'];y=[a+s*d for a,d in zip(fam['star'],fam['direction'])]
        raw=model.raw_coordinates(y);t=K/s;blocks=[]
        for bl in model.raw:
            name=bl['name'];mat=[[evaluate(p,raw) for p in row] for row in bl['mat']]
            if name in model.data['completion']:
                cp=model.data['completion'][name];b=[evaluate(from_form(f),raw) for f in cp['cross']]
                c=evaluate(from_form(cp['diagonal']),raw)
                mat=[row+[b[i]] for i,row in enumerate(mat)]+[b+[t+c]]
            _,pp=ldlt(mat);blocks.append({'name':name,'positive_pivots':pp})
        _,pp=ldlt([[t]]);blocks.append({'name':'M4_None','positive_pivots':pp})
        if len(blocks)!=12 or sum(len(bl['positive_pivots']) for bl in blocks)!=54:
            raise AssertionError('full quotient dimension mismatch')
        out.append({'distance':distance,'s_enclosure':s,'x_enclosure':y[0],'t_enclosure':t,
                    'projected_coordinates':raw,'complete_quotient_blocks':blocks,
                    'algebraic_points_enclosed_by_rational_intervals':True})
    return out
