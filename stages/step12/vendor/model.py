"""Exact level-9 word quotient, using only the Python standard library.

The physical planar products are q=x*x, u=x*p, v=x*a. We keep them as
independent affine variables ONLY while checking the certificate identity.
All inequalities used for them are subsequently justified after substitution.
"""
from __future__ import annotations
from fractions import Fraction as Q
from functools import lru_cache
import word_algebra as wa

@lru_cache(None)
def weighted_words(level: int, charge: int | None = None):
    out, stack = [], [()]
    while stack:
        word = stack.pop()
        lev = wa.level(word)
        if charge is None or wa.charge(word) == charge:
            out.append(word)
        for letter in range(4):
            if lev + wa.WEIGHT[letter] <= level:
                stack.append(word + (letter,))
    return tuple(sorted(out, key=lambda w: (wa.level(w), len(w), w)))

# Same words as the retained reference implementation, without enumerating 4**L.
wa.words = weighted_words

def check(condition, message):
    if not condition:
        raise AssertionError(message)

def finish(words, rows, n, **extra):
    pivots = wa.sparse_rref(rows, n)
    free = [i for i in range(n) if i not in pivots]
    forms = {(): {n: Q(1)}}
    def variable_form(i):
        return ({j: -v for j, v in pivots[i].items() if j != i}
                if i in pivots else {i: Q(1)})
    forms.update({w: variable_form(i) for i, w in enumerate(words)})
    return dict(words=words, rows=rows, n=n, rank=len(pivots), free=free,
                forms=forms, variable_forms=[variable_form(i) for i in range(n)], **extra)

def prove_virial_coefficient():
    low = wa.quotient(5)
    relation = {(2,3): Q(-1), (0,1): Q(-1),
                (0,1,0,1): Q(-2), (0,0,1,1): Q(2)}
    residual = wa.evaluate_trace_polynomial(relation, low['forms'])
    check(not residual, 'Independent low-order virial identity failed')
    # p=x+2(b-a), hence xb=xa+(xp-x^2)/2, NOT /4.
    return True

def build_global9():
    prove_virial_coefficient()
    words = [w for w in weighted_words(9, 0) if w]
    index = {w:i for i,w in enumerate(words)}
    nraw = len(words)
    q,u,v = nraw,nraw+1,nraw+2
    n = nraw+3
    low = wa.low_moments()
    rows=[]
    for _,_,equation in wa.equations(9):
        row={}
        for mon,c in equation.items():
            if not mon:
                wa.add_term(row,n,c)
            elif len(mon)==1:
                wa.add_term(row,index[mon[0]],c)
            elif len(mon)==2:
                w1,w2=mon
                if wa.level(w1)>wa.level(w2):w1,w2=w2,w1
                check(w1 in low,'Unexpected nonlinear product at level 9')
                known=low[w1]
                if 0 in known:
                    check(len(known)==1,'Unexpected low-order constant')
                    wa.add_term(row,index[w2],c*known[0])
                else:
                    check(known=={1:Q(1)},'Expected a factor x')
                    if w2 in low:
                        for power,value in low[w2].items():
                            wa.add_term(row,index[(0,1)] if power==0 else q,c*value)
                    else:
                        check(wa.level(w2)==4,'Unexpected high product factor')
                        if len(w2)==2:
                            check(set(w2)=={2,3},'Expected momentum pair')
                            wa.add_term(row,u,-c)
                        else:
                            check(len(w2)==4 and all(i<2 for i in w2),'Expected coordinates')
                            cyclic=min(w2[i:]+w2[:i] for i in range(4))
                            check(cyclic in ((0,0,1,1),(0,1,0,1)),'Unexpected fourth moment')
                            wa.add_term(row,v,c)
                            if cyclic==(0,1,0,1):
                                wa.add_term(row,u,c/2)
                                wa.add_term(row,q,-c/2)
            else:
                raise AssertionError('Products of more than two traces are not supported')
        rows.append(wa.clean(row))
    return finish(words,rows,n,lifts=(q,u,v),level=9)

def build_fixed(level=9,x=Q('0.38900458')):
    if not 3<=level<=10:raise ValueError('This fixed-x generator supports levels 3 through 10')
    words=[w for w in weighted_words(level,0) if w]
    index={w:i for i,w in enumerate(words)};n=len(words)
    low={w:sum(c*x**power for power,c in f.items()) for w,f in wa.low_moments().items()}
    rows=[]
    for _,_,equation in wa.equations(level):
        row={}
        for mon,c in equation.items():
            if not mon:wa.add_term(row,n,c)
            elif len(mon)==1:wa.add_term(row,index[mon[0]],c)
            elif len(mon)==2:
                u,v=mon
                if u in low:wa.add_term(row,index[v],c*low[u])
                elif v in low:wa.add_term(row,index[u],c*low[v])
                else:raise AssertionError('Nonlinear term not affine at fixed x')
            else:raise AssertionError('Unexpected product')
        rows.append(wa.clean(row))
    rows.append({index[(0,1)]:Q(1),n:-x})
    return finish(words,rows,n,x=x,level=level)

def reduced(poly, model):
    return wa.evaluate_trace_polynomial(poly,model['forms'])

def affine(poly,model):
    form=reduced(poly,model);free_index={j:i+1 for i,j in enumerate(model['free'])}
    result=[Q(0)]*(1+len(free_index))
    for j,c in form.items():result[0 if j==model['n'] else free_index[j]]+=c
    return tuple(result)

def variable_affine(index,model):
    form=model['variable_forms'][index];fk={j:i+1 for i,j in enumerate(model['free'])}
    result=[Q(0)]*(1+len(fk))
    for j,c in form.items():result[0 if j==model['n'] else fk[j]]+=c
    return tuple(result)

def add(*terms):
    """Linear combination of (rational coefficient, affine vector) pairs."""
    if not terms:raise ValueError('At least one vector required')
    out=[Q(0)]*len(terms[0][1])
    for c,vector in terms:
        for i,v in enumerate(vector):out[i]+=c*v
    return tuple(out)

def operator_from_json(terms):
    out={}
    for item in terms:wa.add_term(out,tuple(item['word']),Q(item['coefficient']))
    return wa.clean(out)

C_OPERATOR={(0,0,1):Q(1),(1,0,0):Q(-1)}
D_OPERATOR={(2,1,0):Q(1),(3,0,0):Q(1)}

def observables(model):
    one=affine({():Q(1)},model)
    out={'one':one,'x':affine({(0,1):Q(1)},model),
         'p':affine({(2,3):Q(-1)},model),
         'a':affine({(0,0,1,1):Q(1)},model),
         'b':affine({(0,1,0,1):Q(1)},model),
         'd':affine(wa.inner(D_OPERATOR,D_OPERATOR),model),
         'j':affine({w:c/2 for w,c in wa.inner(C_OPERATOR,C_OPERATOR).items()},model)}
    if 'lifts' in model:
        for name,i in zip(('q','u','v'),model['lifts']):out[name]=variable_affine(i,model)
    return out

def scalar_vectors(obs):
    lo,hi=Q(3,8),Q(2,5)
    return {
      'xu':add((hi,obs['one']),(-1,obs['x'])),
      'square_chord':add((lo+hi,obs['x']),(-lo*hi,obs['one']),(-1,obs['q'])),
      'xp_ul':add((hi,obs['p']),(-1,obs['u'])),
      'xa_ul':add((hi,obs['a']),(-1,obs['v']))}

def verify_all_rows(model):
    n=model['n'];vs=model['variable_forms']
    for row in model['rows']:
        total={}
        for j,c in row.items():
            f={n:Q(1)} if j==n else vs[j]
            for k,v in f.items():wa.add_term(total,k,c*v)
        check(not wa.clean(total),'A generated equality was not satisfied identically')
    return len(model['rows'])
