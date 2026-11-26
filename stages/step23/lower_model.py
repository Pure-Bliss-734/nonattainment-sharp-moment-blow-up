"""Stage23: a completed-square certificate for the universal bounded deficit.

Algebraic coefficients are exact functions of the already-certified KKT root.
I objects enclose those coefficients; they are NOT definitions by midpoints.
Zero-containing enclosures are regression diagnostics, not identity proofs.
The identities rely on the explicit inverse definitions and the inherited root
and stationarity equations, as detailed in the proof in the cited paper.
"""
from fractions import Fraction as F
from stage23_sources import prepare
OLD,bm,DATA,INTEGRITY=prepare()
ap=bm.ap
from audit_interval import I,Z,O,zeros,eye,tr,dot,mm,mv,add,sub,scale,solve,spd


def poly(raw):return {tuple(m):I(c) for m,c in raw}
def plus(*polys):
    out={}
    for p in polys:
        for mon,c in p.items():out[mon]=out.get(mon,Z)+c
    return out

def sc(p,s):return {m:c*s for m,c in p.items()}
def mul(p,q):
    out={}
    for m,c in p.items():
        for n,d in q.items():
            k=tuple(sorted(m+n));out[k]=out.get(k,Z)+c*d
    return out

def dx(p,k):return dict(ap.deriv(list(p.items()),k))
def val(p,y):return ap.peval(list(p.items()),y)
def trace_poly(A,P):
    if len(A)!=len(P):raise ValueError('trace dimensions')
    return plus(*(sc(P[i][j],A[j][i]) for i in range(len(P)) for j in range(len(P))))

def sym(A):return scale(add(A,tr(A)),F(1,2))
def zero_diagnostic(v,label):
    if not v.contains(0):raise AssertionError('identity inconsistent: '+label)
def affine_shape(p):
    # All polynomials are affine in z and at most quadratic in (X,z).
    for m in p:
        if len(m)>2 or sum(k!=0 for k in m)>1:
            raise AssertionError('unsupported polynomial structure')

def at_x(p,xi):
    out={}
    for m,c in p.items():
        rest=tuple(k for k in m if k!=0);v=c
        for _ in range(m.count(0)):v=v*xi
        out[rest]=out.get(rest,Z)+v
    return out

def quotient_at_root(p,xi,sigma):
    """Return p/[sigma*(X-xi)] given the separately proved p(xi,z)=0.

    For p=P0(z)+X P1(z)+X² P2, the result is
    sigma*[P1(z)+(X+xi)P2]. No division by a small distance is used.
    """
    affine_shape(p);parts=[{}, {},{}]
    for m,c in p.items():
        k=m.count(0);rest=tuple(a for a in m if a!=0)
        parts[k][rest]=parts[k].get(rest,Z)+c
    return sc(plus(parts[1],sc(parts[2],xi),
                   {tuple(sorted((0,)+m)):c for m,c in parts[2].items()}),sigma)


def reconstruct(side,metric_shift,data=None,half=F(1,2)):
    """Construct a correction by a positive weighted right inverse.

    H has one row for every e_i^T B_b K_b e_j in the original 20-coordinate
    space, not merely the 19-coordinate symmetry slice. P is block diagonal
    W_b tensor (B_b(y*)+eta I)^-1. eta is an algorithmic metric, not a physical
    regularizer or additional constraint.
    """
    if data is None:data=DATA
    eta=F(metric_shift)
    if eta<=0:raise ValueError('positive metric required')
    pref=ap.reconstruct(data,side)
    inherited=ap.certify_lower(data,pref) # Recheck parameter mixed coefficients.
    C=pref['C'];sigma=pref['sigma'];gamma=inherited['gamma'];F2=inherited['F2']
    star=pref['rawstar'];xi=star[0]
    raw={name:[[sc(dict(pol),O) for pol in row] for row in P] for name,P in pref['raw'].items()}
    L=plus(sc(poly(data['observables']['p']),F(7,2)),sc(poly(data['observables']['a']),3))
    Lstar=val(L,star)
    A0={};rows=[]
    for num,inf in enumerate(pref['infos']):
        K=inf['rawK'];B=raw[inf['parent']];n=len(K)
        if inf['parent'] in data['completion']:
            if inf['nu']!=1:raise AssertionError('this implementation uses simple completion kernels')
            cp=data['completion'][inf['parent']]
            q=pref['completion'][inf['parent']]['cross'];w=inf['W'][0][0]
            s=plus(*(sc(poly(b),v[0]) for v,b in zip(K,cp['cross'])))
            A0=plus(A0,sc(s,w*q))
        for j in range(inf['nu']):
            for i in range(n):
                h=plus(*(sc(B[i][a],K[a][j]) for a in range(n)))
                affine_shape(h)
                zero_diagnostic(val(h,star),'root kernel h')
                row=[val(dx(h,k),star) for k in range(1,21)]
                rows.append({'info':num,'column':j,'row':i,'h':h,'H':row})
    if len(rows)!=(76 if side=='lower' else 82):raise AssertionError('full face equation count')
    H=[r['H'] for r in rows]
    PH=[];metric_pivots=[]
    for num,inf in enumerate(pref['infos']):
        n=len(inf['rawK']);W=inf['W'];rank=inf['nu']
        B=ap.matval(pref['raw'][inf['parent']],star)
        Areg=add(B,scale(eye(n),eta));pivots=spd(sym(Areg))
        metric_pivots.append({'name':inf['name'],'pivots':pivots})
        Ainv=solve(Areg,eye(n))
        Hb=[[v['H'] for v in rows if v['info']==num and v['column']==j] for j in range(rank)]
        AH=[mm(Ainv,h) for h in Hb]
        for j in range(rank):
            PH.extend([[sum((W[j][k]*AH[k][i][col] for k in range(rank)),Z)
                        for col in range(20)] for i in range(n)])
    Gram=sym(mm(tr(H),PH));gram_pivots=spd(Gram)
    target=[half*C*sigma*gamma*val(dx(L,k),star)-val(dx(A0,k),star) for k in range(1,21)]
    solution=solve(Gram,[[v] for v in target]);alpha=[r[0] for r in mm(PH,solution)]
    # Exact relation H^T alpha=target is by definition alpha=P H (H^T P H)^-1 target.
    for a,b in zip(mv(tr(H),alpha),target):zero_diagnostic(a-b,'weighted right inverse')
    Es=[zeros(len(i['rawK']),i['nu']) for i in pref['infos']]
    A=A0
    for a,row in zip(alpha,rows):
        Es[row['info']][row['row']][row['column']]+=a
        A=plus(A,sc(row['h'],a))
    kpoly={};partners=[]
    for inf,E in zip(pref['infos'],Es):
        V=tr(solve(inf['W'],tr(E))) # V W=E, W symmetric positive definite.
        VVW=mm(mm(V,inf['W']),tr(V))
        kpoly=plus(kpoly,trace_poly(VVW,raw[inf['parent']]))
        bottom=None
        if inf['parent'] in data['completion']:
            cp=data['completion'][inf['parent']]
            q=pref['completion'][inf['parent']]['cross'];w=inf['W'][0][0]
            vb=plus(*(sc(poly(b),v[0]) for v,b in zip(V,cp['cross'])))
            kpoly=plus(kpoly,sc(vb,2*w*q),sc(poly(cp['diagonal']),w*q.sq()))
            bottom=[q]
        partners.append({'name':inf['name'],'parent':inf['parent'],'K':inf['rawK'],
                         'W':inf['W'],'V':V,'completion_bottom':bottom})
    Dpoly=plus({():O-sigma*gamma*Lstar-F2*sigma*xi,(0,):F2*sigma},sc(L,sigma*gamma))
    P=plus(A,sc(plus({():O},Dpoly),-C/2))
    # P(xi,z)=0 follows from the coefficient equation and its value at z*.
    # These checks are diagnostics, not replacing the root/inverse definitions.
    at=at_x(P,xi)
    for m,c in at.items():zero_diagnostic(c,'anchored polynomial coefficient '+str(m))
    rpoly=quotient_at_root(P,xi,sigma)
    deficit=plus(sc(kpoly,O/C),sc(rpoly,-2))
    for pol in [A,kpoly,Dpoly,P,rpoly,deficit]:affine_shape(pol)
    delta_poly={():-sigma*xi,(0,):I(sigma)}
    Tpoly={}
    for inf in pref['infos']:
        Y=mm(mm(inf['rawK'],inf['W']),tr(inf['rawK']))
        Tpoly=plus(Tpoly,trace_poly(Y,raw[inf['parent']]))
    supporting=plus(Tpoly,sc(mul(delta_poly,Dpoly),-1))
    for m,c in supporting.items():zero_diagnostic(c,'support identity '+str(m))
    quotient=plus(P,sc(mul(delta_poly,rpoly),-1))
    for m,c in quotient.items():zero_diagnostic(c,'anchored quotient '+str(m))
    return {'side':side,'eta':eta,'pref':pref,'inherited_lower':inherited,'C':C,'xi':xi,
            'sigma':sigma,'gamma':gamma,'F2':F2,'L':L,'Lstar':Lstar,'raw':raw,
            'rows':rows,'H':H,'PH':PH,'Gram':Gram,'gram_pivots':gram_pivots,
            'metric_pivots':metric_pivots,'target':target,'alpha':alpha,'partners':partners,
            'A0':A0,'A':A,'D':Dpoly,'k':kpoly,'r':rpoly,'deficit':deficit,'delta':delta_poly,
            'anchored_coefficient_diagnostics':at,'support_identity_diagnostics':supporting,
            'quotient_identity_diagnostics':quotient}


def bound_deficit(model,delta0,claim,data=None):
    if data is None:data=DATA
    delta0=F(delta0);claim=F(claim)
    side=model['side'];par=data['neighborhood_parameters'][side]
    if not 0<delta0<=F(data['neighborhood_parameters']['delta']):raise ValueError('outside inherited neighborhood')
    rad0=F(par['coordinate_A'])*bm.sqrt_interval(I(delta0)).hi()+F(par['coordinate_B'])*delta0
    rad=min(F(par['R']),rad0)
    star=model['pref']['rawstar'];xi=star[0];sigma=model['sigma']
    # One-sided X interval, including endpoint. Bound every coefficient of z.
    X=I(xi.lo(),xi.hi()+delta0) if sigma==1 else I(xi.lo()-delta0,xi.hi())
    center=[X]+star[1:]
    S=model['deficit'];at_center=val(S,center)
    grad=[val(dx(S,k),center) for k in range(1,21)]
    Lipschitz=sum(g.norm() for g in grad)
    upper=at_center.hi()+rad*Lipschitz
    if not upper<claim:raise AssertionError('claimed bounded deficit is too small')
    return {'side':side,'distance_limit':delta0,'metric_shift':model['eta'],
            'lower_remainder_bound':claim,'computed_deficit_upper':upper,
            'deficit_at_endpoint':val(S,star),'X_interval':X,'coordinate_radius':rad,
            'coordinate_bound_before_min':rad0,'center_deficit':at_center,
            'gradient_intervals':grad,'gradient_l1_upper':Lipschitz,
            'strict_margin':claim-upper,
            'conclusion':'For EVERY finite complete feasible moment sequence: t > C/delta - lower_remainder_bound.',
            'identity':'sum_b tr W_b (U_b-C/delta*K_b)^T G_b (U_b-C/delta*K_b) = C*(t-C/delta+deficit(X,z)) >= 0',
            'not_an_optimized_subleading_coefficient':True}


def serialize_poly(p):return [[list(m),v] for m,v in sorted(p.items())]
def coefficient_report(model):
    return {'side':model['side'],'C':model['C'],'xi':model['xi'],'sigma':model['sigma'],
            'gamma':model['gamma'],'F2':model['F2'],'metric_shift':model['eta'],
            'row_layout':[{k:v for k,v in r.items() if k in ['info','column','row']} for r in model['rows']],
            'H':model['H'],'weighted_Gram':model['Gram'],'Gram_positive_pivots':model['gram_pivots'],
            'metric_positive_pivots':model['metric_pivots'],'target':model['target'],'alpha':model['alpha'],
            'partners':model['partners'],
            'polynomials':{name:serialize_poly(model[name]) for name in ['A0','A','D','k','r','deficit']},
            'identity_semantics':'Exact algebraic functions of the inherited certified root. Inverse definitions and KKT equations prove the identities; interval-zero diagnostics are not used as standalone equality proofs.',
            'full_original_coordinate_count':20,'full_annihilation_row_count':len(model['rows'])}
