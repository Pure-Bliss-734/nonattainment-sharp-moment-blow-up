from __future__ import annotations

import json

from fractions import Fraction as F

import model as md

from model import wa, check

from verify_step10 import decimal

QOP={(0,1,0):F(3),(0,1,2):F(-1),(2,1,0):F(-1)}

def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))

def verify_certificate(data, model, obs):
    total=(F(0),)*len(obs['one']);scalars=md.scalar_vectors(obs)
    counts={'ordinary':0,'ground_old':0,'ground_new':0,'scalar':0}
    coefficients=[]
    for term in data['terms']:
        c=F(term['coefficient']);check(c>0,'Nonpositive coefficient')
        kind=term['kind']
        if kind=='scalar':
            check(term['block'] in scalars,'Unknown scalar inequality')
            val=scalars[term['block']];counts['scalar']+=1
        else:
            check(kind in ('ordinary','ground'),'Unknown positivity type')
            op=md.operator_from_json(term['operator'])
            check(bool(op),'Empty operator')
            level=max(map(wa.level,op))
            check(level<=4,'Operator exceeds cutoff')
            check(len({wa.charge(w) for w in op})==1,'Mixed charge operator')
            if kind=='ground' and level==4:
                check(op==QOP,'Unexpected additional ground-state constraint')
                counts['ground_new']+=1
            elif kind=='ground':
                check(level<=3,'Old ground operator exceeds level 3')
                counts['ground_old']+=1
            else:counts['ordinary']+=1
            poly=wa.inner(op,op,kind=='ground')
            check(max(map(wa.level,poly),default=0)<=9,'Gram polynomial exceeds level 9')
            val=md.affine(poly,model)
        total=md.add((1,total),(c,val));coefficients.append(str(c))
    K=F(data['K']);simple=F(data['simple_bound'])
    check(total==md.add((K,obs['one']),(-1,obs['d'])),'Certificate identity fails')
    check(counts['ground_new']==1,'Exactly one new scalar form required')
    check(K<simple,'Advertised rational bound does not follow')
    return {'K':str(K),'K_decimal':decimal(K,55),'simple_bound':str(simple),
            'term_count':len(data['terms']),'counts':counts,
            'exact_affine_stationarity_rows':len(obs['one'])-1,
            'constant_identity':True,'all_coefficients_strictly_positive':True,
            'new_operator':[{ 'word':list(w),'coefficient':str(c)} for w,c in QOP.items()],
            'max_integer_basis_coefficient':max(abs(v) for a in data['selected_integer_basis_vectors'] for v in a['vector']),
            'claim':'A sufficient bound in the old Gram baseline plus level-9 equalities and the single scalar N(Q)>=0; not an optimal bound.'}

def reconstruct_base(inputs):
    data=read_json(inputs/'primal_point.json');x=F(data['x']);fixed=md.build_fixed(9,x)
    values={tuple(t['word']):F(t['value']) for t in data['free_values']}
    check(set(values)=={fixed['words'][j] for j in fixed['free']},'Primal basis mismatch')
    v={j:values[fixed['words'][j]] for j in fixed['free']};v[fixed['n']]=F(1)
    state={w:sum((c*v[j] for j,c in form.items()),F(0)) for w,form in fixed['forms'].items()}
    return state
