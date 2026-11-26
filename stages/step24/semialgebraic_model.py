"""Structural polynomial audit; NOT a quantifier-elimination implementation.

The accompanying written proof invokes real-closed-field quantifier elimination
and univariate semialgebraic monotonicity. This module checks the concrete finite
rational input to which those theorems are applied and exports its full matrices.
"""
from fractions import Fraction as F
from itertools import combinations

def canonical(p):
    d={}
    for m,c in p:
        m=tuple(m);c=F(c)
        if tuple(sorted(m))!=m:raise ValueError('unordered monomial')
        if any(not isinstance(i,int) or not 0<=i<21 for i in m):raise ValueError('unexpected variable')
        d[m]=d.get(m,F(0))+c
    return {m:c for m,c in d.items() if c}

def export_model(data):
    blocks=[];raw_count=0;count=0;terms=0;maxdegree=0
    for b in data['raw']:
        A=[[canonical(v) for v in row] for row in b['mat']];n=len(A)
        if any(len(row)!=n for row in A):raise ValueError('nonsquare')
        raw_count+=n*n
        for i in range(n):
            for j in range(n):
                if A[i][j]!=A[j][i]:raise ValueError('not exactly symmetric')
        if b['name'] in data['completion']:
            cp=data['completion'][b['name']];v=[canonical(x) for x in cp['cross']]
            if len(v)!=n:raise ValueError('wrong completion dimension')
            c=canonical(cp['diagonal']);c[(21,)]=c.get((21,),F(0))+1
            A=[row+[q] for row,q in zip(A,v)]+[v+[c]]
        blocks.append({'name':b['name'],'matrix':A})
    blocks.append({'name':'M4_None','matrix':[[{(21,):F(1)}]]})
    for b in blocks:
        A=b['matrix'];n=len(A)
        for i in range(n):
            for j in range(n):
                if A[i][j]!=A[j][i]:raise ValueError('bad completed symmetry')
                p=A[i][j];count+=1;terms+=len(p)
                for m in p:
                    if len(m)>2 or sum(i!=0 for i in m)>1:raise ValueError('not affine in witnesses / quadratic in joint coordinates')
                    maxdegree=max(maxdegree,len(m))
    if raw_count!=343 or len(blocks)!=12:raise ValueError('wrong finite model layout')
    # All principal minors, NOT just leading minors, can encode PSD.
    counts=[2**len(b['matrix'])-1 for b in blocks]
    encoded=[{'name':b['name'],'dimension':len(b['matrix']),
              'matrix':[[[[list(m),str(c)] for m,c in sorted(p.items())] for p in row] for row in b['matrix']]} for b in blocks]
    return {'variables':['X']+['z'+str(i) for i in range(20)]+['t'],
            'coefficient_field':'Q; substitute X=xi+sigma*delta with real algebraic xi',
            'blocks':encoded,'block_dimensions':[b['dimension'] for b in encoded],
            'raw_entries':raw_count,'complete_entries':count,'nonzero_terms':terms,
            'max_entry_total_degree':maxdegree,'all_principal_minors_count':sum(counts),
            'principal_minor_counts':counts,
            'equivalent_PSD_predicate':'For each block G, all nonempty principal minors are >=0. Equivalently exists real L: G=L*L^T.',
            'minors_expanded':False,'quantifier_elimination_executed':False,
            'meaning':'These are full completion matrices, not just projected PSD blocks. The invisible ninth-level coordinates have the inherited affine lift.'}

def formulas():
    return {
      'feasible':'Phi_s(delta,z,t): delta>0 AND t>=0 AND ALL full G_b(xi_s+sigma_s*delta,z,t) are PSD',
      'graph_infimum':'v=tau_s(delta) IFF [FORALL(z,t): Phi_s(delta,z,t) => v<=t] AND [FORALL eps>0 EXISTS(z,t): Phi_s(delta,z,t) AND t<v+eps]',
      'graph_residue':'(delta,v) belongs to graph_tau AND delta*(v-r)=C_s',
      'limit_singleton':'d=D_s IFF FORALL eps>0 EXISTS h>0 FORALL(delta,v): (0<delta<min(h,delta0) AND graph_tau(delta,v)) => -eps*delta < delta*(v-d)-C_s < eps*delta',
      'proof_steps':['Tarski-Seidenberg gives semialgebraicity of these finite quantified predicates.',
                     'Previous bounds make tau finite on every sufficiently small positive delta and residue bounded.',
                     'Univariate semialgebraic monotonicity gives a finite limit.',
                     'The singleton limit is definable over the real algebraic numbers, hence is algebraic.',
                     'Growth dichotomy gives a positive rational leading exponent for the nonzero remaining error, or the error vanishes identically near zero.'],
      'scope':'Written mathematical implication using standard theorems. This JSON does not claim a machine execution of general quantifier elimination.'}

def examples():
    count=0;records=[]
    for j in (2,3,5,7,11,101):
        d=F(1,j);b=F(2,3)
        t=1/d+b
        if d*(t-b)!=1:raise AssertionError('shifted scalar model')
        count+=1
        t=1/(d*(1-d))
        if t-1/d-1!=d/(1-d):raise AssertionError('curvature remainder')
        count+=1
        # Delta=q^2 makes the fractional-power model rational for reference checks.
        q=F(1,j);d=q*q;t=1/d+b+q
        if d*(t-b-q)!=1 or q*q!=d:raise AssertionError('fractional example')
        count+=1
        # Nonattained infimum: y>0, y*s=1, lower threshold C/delta+b+y.
        d=F(1,j);y=F(1,j*j);s=1/y;t=1/d+b+y
        if y*s!=1 or t<=1/d+b or d*(t-b-y)!=1:raise AssertionError('unattained example')
        count+=1
    records=[{'name':'constant shift','tau':'1/delta+2/3','D':'2/3'},
             {'name':'curvature','tau':'1/(delta*(1-delta))','D':'1','next_power':'1'},
             {'name':'fractional power','constraints':'q>=0, q^2=delta, [[delta,1],[1,t-2/3-q]]>=0','tau':'1/delta+2/3+sqrt(delta)','D':'2/3','next_power':'1/2'},
             {'name':'unattained infimum','constraints':'y>=0,y*s=1,[[delta,1],[1,t-2/3-y]]>=0','tau':'1/delta+2/3','attained':False}]
    # Zero leading principal minors alone do not certify PSD.
    bad=[[F(0),F(0)],[F(0),F(-1)]]
    assert bad[0][0]==0 and bad[0][0]*bad[1][1]-bad[0][1]**2==0 and bad[1][1]<0
    return {'exact_rational_reference_points':count,'cases':records,
            'not_a_proof_by_samples':'Continuous examples are established by the scalar determinant identities shown; generic semialgebraic theorem is a written proof.',
            'leading_minors_only_negative_control_rejected':True}
