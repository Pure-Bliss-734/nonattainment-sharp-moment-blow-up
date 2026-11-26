from __future__ import annotations

import json

from decimal import Decimal, localcontext

from fractions import Fraction as Q

import model as md

from model import wa, check

def decimal(q,digits=45):
    q=Q(q)
    with localcontext() as context:
        context.prec=digits
        return str(Decimal(q.numerator)/Decimal(q.denominator))

def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')

def polynomial_value(coeffs,x):
    result=Q(0)
    for c in reversed(coeffs):result=result*x+c
    return result

def exact_psd_rank(matrix):
    A=[[Q(v) for v in row] for row in matrix];n=len(A);pivots=[]
    check(all(A[i][j]==A[j][i] for i in range(n) for j in range(n)),'Nonsymmetric matrix')
    while A:
        check(all(A[i][i]>=0 for i in range(len(A))),'Negative PSD diagonal')
        pivot=next((i for i in range(len(A)) if A[i][i]>0),None)
        if pivot is None:
            check(not any(v for row in A for v in row),'Zero diagonal with a nonzero off-diagonal')
            break
        p=A[pivot][pivot];pivots.append(p)
        rest=[i for i in range(len(A)) if i!=pivot]
        A=[[A[i][j]-A[i][pivot]*A[pivot][j]/p for j in rest] for i in rest]
    return pivots

def op_poly(terms,ground):
    op=md.operator_from_json(terms)
    check(op and max(map(wa.level,op))<=4,'Operator outside the declared cutoff')
    check(len({wa.charge(w) for w in op})==1,'Operator mixes charge sectors')
    return wa.inner(op,op,ground)

def add_word_polynomial(out,poly,weight):
    for word,c in poly.items():wa.add_term(out,word,c*weight)

def foundation(inputs):
    # These checks use the old level-8 quotient, not the new level-9 lift.
    s8=wa.quotient(8);n=s8['n'];forms=s8['forms']
    def reduce(poly):return wa.evaluate_trace_polynomial(poly,forms)
    def combo(*terms):return wa.linear_combine(*terms)
    x={n+1:Q(1)};one={n:Q(1)}
    p=reduce({(2,3):Q(-1)});a=reduce({(0,0,1,1):Q(1)})
    z={(0,):Q(1)};pi={(2,):Q(1)}
    identities=[
        (wa.inner(z,z),x),(wa.inner(z,pi),combo((Q(1,2),one))),
        (wa.inner(pi,pi),p),(wa.inner(z,z,True),combo((Q(1,2),one))),
        (wa.inner(z,pi,True),p),(wa.inner(pi,pi,True),combo((1,x),(Q(1,2),one))),
        (wa.inner({(0,1):Q(1),(1,0):Q(-1)},{(0,1):Q(1),(1,0):Q(-1)}),combo((1,p),(-1,x))),
        (wa.inner({(0,0):Q(1)},{(0,0):Q(1)}),a)]
    for polynomial,expected in identities:check(reduce(polynomial)==expected,'Low-order Gram identity failed')
    # M2 and N2 imply x>0, p>=x, 4x²-2x-1<=0, and 8x³+4x²-1>=0.
    lo,hi=Q(3,8),Q(2,5)
    check(8*lo**3+4*lo**2-1<0,'Lower rational radius comparison failed')
    check(4-2-1>0,'Elementary upper x<1 comparison failed')
    rays={}
    for side in ('lower','upper'):
        data=json.loads((inputs/f'radius_{side}_ray.json').read_text())
        total={};c_weight=None
        for atom in data['atoms']:
            weight=Q(atom['weight']);norm=Q(atom['norm2'])
            check(weight>0 and norm>0,'Nonpositive prior ray atom')
            check(norm==sum(v*v for v in atom['vector']),'Prior atom normalization mismatch')
            op=md.operator_from_json(atom['operator_terms'])
            add_word_polynomial(total,wa.inner(op,op,atom['kind']=='ground'),weight/norm)
            if atom['kind']=='ordinary' and op==md.C_OPERATOR:
                check(c_weight is None,'Repeated distinguished atom')
                c_weight=2*weight/norm
        coeffs=[Q(v) for v in data['coefficients_ascending']]
        check(reduce(wa.clean(total))=={n+i:v for i,v in enumerate(coeffs) if v},'Prior ray word identity failed')
        check(c_weight is not None and c_weight>0,'Missing positive C-norm atom')
        rays[side]={'coeffs':coeffs,'j_weight':c_weight}
    upper=rays['upper']['coeffs']
    check(upper[2]>0 and upper[1]+2*upper[2]<0,'Upper ray is not decreasing on [0,1]')
    check(polynomial_value(upper,hi)<0,'Cannot deduce x<2/5')
    return rays,{'level8_nonempty_words':n,'low_order_Gram_identities':len(identities),
        'prior_radius_ray_identities':2,'strict_radius_bracket':['3/8','2/5'],
        'proof_note':'Low Gram bounds give x>3/8 and x<1. The upper ray decreases on [0,1] and is negative at 2/5. Both rays were rechecked from raw words.'}

def verify_moment_certificate(data,model,obs):
    length=len(obs['one']);total=(Q(0),)*length
    scalars=md.scalar_vectors(obs);physical=[];scalar_weights={k:Q(0) for k in scalars}
    maxcoef=0
    for atom in data['atoms']:
        weight,norm=Q(atom['weight']),Q(atom['norm2'])
        check(weight>0 and norm>0,'Certificate has a nonpositive weight or norm')
        check(norm==sum(int(v)**2 for v in atom['vector']),'Atom norm does not match integer vector')
        maxcoef=max(maxcoef,*map(abs,atom['vector']))
        if atom['kind']=='scalar':
            check(atom['block'] in scalars and norm==1,'Unsupported scalar inequality')
            value=scalars[atom['block']];scalar_weights[atom['block']]+=weight
        else:
            check(atom['kind'] in ('ordinary','ground'),'Unknown operator kind')
            poly=op_poly(atom['operator'],atom['kind']=='ground')
            check(max(map(wa.level,poly),default=0)<=9,'Gram entry above level 9')
            value=md.affine(poly,model)
            physical.append(atom['kind'])
        total=md.add((1,total),(weight/norm,value))
    K=Q(data['K']);expected=md.add((K,obs['one']),(-1,obs['d']))
    check(total==expected,'Exact certificate identity K-d=sum(nonnegative terms) failed')
    check(K<Q(8,25),'Certificate does not prove the requested simple upper bound')
    return {'K':str(K),'K_decimal':decimal(K,65),'simple_K':'8/25',
            'stationarity_rows_zero':length-1,'constant_identity_verified':True,
            'positive_weights':len(data['atoms']),'ordinary_atoms':physical.count('ordinary'),
            'ground_atoms':physical.count('ground'),'scalar_atoms':len(data['atoms'])-len(physical),
            'maximum_integer_atom_coefficient':maxcoef,
            'scalar_weights':{k:str(v) for k,v in scalar_weights.items() if v},
            'solver_or_floating_point_used_in_verification':False}
