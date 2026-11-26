"""Quantitative second-order lifting on the inherited, explicit matrix pencil.

Proof computations: integers/Fraction plus the pinned outward interval engine.
An interval containing zero is a diagnostic, never an identity certificate.
Exact cancellations used below follow from the certified root equations,
Schur differentiation, and the definition of the right-inverse directions.
"""
from fractions import Fraction as F
import json,sys
from source_tools import prepare
_,_old,_integrity=prepare()
sys.path.insert(0,str(_old))
import audit_prefactor as ap
from audit_interval import I,Z,O,zeros,eye,tr,dot,mm,mv,add,sub,scale,congr,part,spd,solve,max_row_norm

def encode(x):
    if isinstance(x,I):return x.record()
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x

def load_data():return json.loads((_old/'inputs/shared_pencils.json').read_text())
def matrix_sym(a):return scale(add(a,tr(a)),F(1,2))
def identity_diagnostic(A,B,label):
    if len(A)!=len(B) or any(len(r)!=len(s) for r,s in zip(A,B)):raise ValueError(label+' shape')
    if not all((x-y).contains(0) for r,s in zip(A,B) for x,y in zip(r,s)):
        raise AssertionError(label+' enclosure inconsistent with exact identity')

def build_quadratic(data,side):
    p=ap.reconstruct(data,side);star=p['star'];d0=p['d0'];g2=[];curves=[]
    for inf,layout in zip(p['infos'],data['sides'][side]['layout']):
        P=p['small'][inf['index']]['mat'];rank=inf['rank'];nu=inf['nu'];lead=layout['perm'][:rank]
        M=ap.matval(P,star)
        Md=[[sum((ap.peval(ap.deriv(pol,j),star)*d0[j] for j in range(20)),Z) for pol in row] for row in P]
        Mdd=[[sum((2*I(c)*d0[m[0]]*d0[m[1]] for m,c in pol if len(m)==2),Z) for pol in row] for row in P]
        A=part(M,lead,lead)
        if rank:
            E=part(mm(Md,inf['K']),lead,list(range(nu)))
            cur=sub(congr(Mdd,inf['K']),scale(mm(tr(E),solve(A,E)),2))
        else:cur=congr(Mdd,inf['K'])
        for i,j,_ in inf['locations']:g2.append(cur[i][j])
        curves.append({'name':inf['name'],'g_second_d0':cur})
    target=[h-v/2 for h,v in zip(p['hI'],g2)];d2=mv(p['R'],target)
    if not all((v-w).contains(0) for v,w in zip(mv(p['J'],d2),target)):
        raise AssertionError('right inverse diagnostic')
    coords=[[a,b,c] for a,b,c in zip(star,d0,d2)]
    rawcoords=[[a,b,c] for a,b,c in zip(p['rawstar'],ap.raw_dir(data,d0),ap.raw_dir(data,d2))]
    a=p['sigma']*d2[0]
    if a.lo()<=0:raise AssertionError('this specific construction requires a>0')
    active={inf['index']:inf for inf in p['infos']};projected=[]
    for i,b in enumerate(p['small']):
        ms=ap.matrix_series(b['mat'],coords)
        if i in active:
            inf=active[i];cs=[congr(m,inf['T']) for m in ms]
            projected.append({'name':b['name'],'index':i,'active':True,'cs':cs,'rank':inf['rank'],
                              'H0':p['H0'][i],'HI':p['HI'][i]})
        else:projected.append({'name':b['name'],'index':i,'active':False,'cs':ms})
    completed={}
    for name,cp in data['completion'].items():
        inf=next(i for i in p['infos'] if i['parent']==name);n=len(p['raw'][name])
        # Other parity block first, then the same positive-range columns as
        # the active Schur chart. This keeps its second coefficient exactly HI.
        cols=[]
        for j,b in enumerate(p['small']):
            if b['parent']==name and j!=inf['index']:cols+=tr(b['T'])
        own=mm(p['small'][inf['index']]['T'],inf['T'])
        cols+=tr(own)[:inf['rank']];cols+=[[r[0] for r in inf['rawK']]]
        T=tr(cols)
        if len(cols)!=n:raise AssertionError('incomplete completion chart')
        chart_pivots=spd(mm(tr(T),T))
        cs=[congr(m,T) for m in ap.matrix_series(p['raw'][name],rawcoords)]
        bc0=[ap.peval_series(ap.poly(pol),rawcoords) for pol in cp['cross']]
        bc=[mv(tr(T),[v[k] for v in bc0]) for k in range(7)]
        cc=ap.peval_series(ap.poly(cp['diagonal']),rawcoords)
        completed[name]={'name':name,'index':inf['index'],'cs':cs,'rank':n-1,'bc':bc,'cc':cc,
                         'T':T,'chart_pivots':chart_pivots,'q':p['completion'][name]['cross'],
                         'H0':p['H0'][inf['index']],'HI':p['HI'][inf['index']]}
    return {'pref':p,'side':side,'d2':d2,'a':a,'coords':coords,'rawcoords':rawcoords,
            'curvatures':curves,'projected':projected,'completed':completed}

def schur_remainder(item,U):
    """S=r H0+r² HI+r³ R(r), using cancellation-free divided differences.

    Exact identities: E=r(E1+r Ehat), A=A0+r Ahat,
    (A^-1 Etilde-A0^-1 E1)/r=A^-1(Ehat-Ahat A0^-1 E1).
    No division by an interval containing r=0 occurs.
    """
    cs=item['cs'];r=item['rank'];nu=len(cs[0])-r;ids=list(range(r));tail=list(range(r,r+nu))
    val=ap.value_series(cs,U)
    if r:
        A=part(val,ids,ids);pa=spd(A);A0=part(cs[0],ids,ids)
        E1=part(cs[1],ids,tail)
        Etilde=part(ap.value_series(cs,U,1),ids,tail)
        Ehat=part(ap.value_series(cs,U,2),ids,tail)
        Ahat=part(ap.value_series(cs,U,1),ids,ids)
        G0=solve(A0,E1);G=solve(A,Etilde)
        diff=solve(A,sub(Ehat,mm(Ahat,G0)))
        R=sub(part(ap.value_series(cs,U,3),tail,tail),add(mm(tr(Ehat),G),mm(tr(E1),diff)))
        second=sub(part(cs[2],tail,tail),mm(tr(E1),G0))
    else:
        A=[];pa=[];Etilde=[]
        R=part(ap.value_series(cs,U,3),tail,tail);second=part(cs[2],tail,tail)
    identity_diagnostic(part(cs[1],tail,tail),item['H0'],item['name']+' first Schur')
    identity_diagnostic(second,item['HI'],item['name']+' second Schur')
    # The off-block and kernel-constant zeros follow from endpoint feasibility
    # and the certified kernel definition; diagnostics cannot establish them.
    identity_diagnostic(part(cs[0],tail,tail),zeros(nu),item['name']+' kernel')
    if r:identity_diagnostic(part(cs[0],ids,tail),zeros(r,nu),item['name']+' off kernel')
    return {'A':A,'leading_pivots':pa,'Etilde':Etilde,'R':matrix_sym(R)}

def estimate_strip(curve,r0):
    r0=F(r0)
    if r0<=0:raise ValueError('positive strip required')
    U=I(0,r0);p=curve['pref'];OI=O/p['Omega'];projected=[]
    for it in curve['projected']:
        if not it['active']:
            piv=spd(ap.value_series(it['cs'],U));projected.append({'name':it['name'],'type':'inactive','pivots':piv});continue
        s=schur_remainder(it,U);norm=max_row_norm(s['R']);margin=OI.lo()-r0*norm
        if margin<=0:raise AssertionError(it['name']+' second-order remainder too large')
        projected.append({'name':it['name'],'type':'active','leading_pivots':s['leading_pivots'],
                          'R':s['R'],'remainder_row_norm_upper':norm,'HI_minus_r_normR_lower':margin})
    comp={}
    for name,it in curve['completed'].items():
        s=schur_remainder(it,U);n=it['rank'];E=[row[0] for row in s['Etilde']]
        b=[ap.horner([v[k] for v in it['bc']],U) for k in range(n+1)]
        f=b[:-1];ai_f=[row[0] for row in solve(s['A'],[[v] for v in f])]
        k=ap.horner([v[-1] for v in it['bc'][1:]],U)-dot(E,ai_f)
        z=ap.horner(it['cc'],U)-dot(f,ai_f)
        identity_diagnostic([[it['bc'][0][-1]]],[[it['q']]],name+' crossing normalization')
        comp[name]={'R':s['R'][0][0],'k':k,'Z':z,'q':it['q'],
                    'h0':it['H0'][0][0],'hI':it['HI'][0][0],
                    'leading_pivots':s['leading_pivots'],'chart_pivots':it['chart_pivots']}
    return {'side':curve['side'],'r_limit':r0,'projected':projected,'completion':comp}

def certify_constant(curve,est,m,M_target=None):
    """Certify t=C/r+m, then translate to t<=C/delta+M."""
    m=F(m);r0=est['r_limit'];U=I(0,r0);C=curve['pref']['C'];a=curve['a']
    if m<0:raise ValueError('nonnegative completion constant required')
    records={}
    for name,v in est['completion'].items():
        identity_diagnostic([[C*v['h0']]],[[v['q'].sq()]],name+' C*h0=q^2')
        Mt=m+v['Z'];bottom=C+U*Mt
        if bottom.lo()<=0:raise AssertionError(name+' scaled bottom not positive')
        # After scaling diag(r^-1/2,r^1/2), eliminating the bottom scalar:
        # S/r=hI+rR+[h0*(m+Z)-2qk-r*k²]/[C+r*(m+Z)].
        V=v['hI']+U*v['R']+(v['h0']*Mt-2*v['q']*v['k']-U*v['k'].sq())/bottom
        if V.lo()<=0:raise AssertionError(name+' scaled Schur coefficient not positive')
        records[name]={'scaled_bottom':bottom,'scaled_Schur_over_r':V,**v}
    # a>0 implies delta=r+a r² is increasing, delta(r0)>=r0.
    # Every 0<delta<=r0 has a unique 0<r<=r0. Also
    # C/r+m = C/delta + m + C*a/(1+a*r) <= C/delta + m + C_up*a_up.
    extra=m+C.hi()*a.hi();M=F(-((-extra.numerator)//extra.denominator))
    if M<extra:raise AssertionError('upward integer rounding failed')
    if M_target is not None:
        M=F(M_target)
        if not extra<M:raise AssertionError('requested M does not bound conversion cost')
    return {'side':curve['side'],'r_limit':r0,'distance_limit':r0,'completion_m':m,
            'C':C,'a':a,'M':M,'unrounded_M_upper':extra,
            'parameter_derivative_lower':F(1),'completion':records,'projected':est['projected'],
            'continuous_all_positive_r':True,
            'family':'y=y*+r d0+r^2 d2; delta=r+a r^2; t=C/r+m',
            'conclusion':'tau(delta)<=C/delta+M for EVERY 0<delta<=distance_limit',
            'no_two_sided_bounded_remainder_claim':True}


def sqrt_interval(x):
    """Enclose the nonnegative square root using integer arithmetic only."""
    from math import isqrt
    from audit_interval import SCALE
    x=I(x)
    if x.l<0:raise ValueError('negative square-root argument')
    lo=isqrt(x.l*SCALE);hi=isqrt(x.h*SCALE)
    if hi*hi<x.h*SCALE:hi+=1
    return I.raw(lo,hi)

def sample_quadratic(data,curve,cert,delta):
    """Pointwise regression: interval enclosures, not rational midpoints."""
    d=F(delta)
    if not 0<d<=cert['distance_limit']:raise ValueError('distance outside certified interval')
    a=curve['a'];C=curve['pref']['C']
    r=2*I(d)/(O+sqrt_interval(O+4*a*d))
    if r.lo()<=0 or r.hi()>cert['r_limit']:raise AssertionError('inverse parameter inclusion')
    y=[s+r*v+r.sq()*w for s,v,w in zip(curve['pref']['star'],curve['pref']['d0'],curve['d2'])]
    raw=ap.enc_raw(data,y);t=C/r+cert['completion_m'];records=[]
    # Stable identity for the excess; no subtraction of two huge poles.
    excess=cert['completion_m']+C*a/(O+a*r)
    if excess.hi()>=cert['M']:raise AssertionError('cost exceeds certified constant')
    for name,P in curve['pref']['raw'].items():
        A=ap.matval(P,raw)
        if name in data['completion']:
            cp=data['completion'][name];b=[ap.peval(ap.poly(pol),raw) for pol in cp['cross']]
            c=ap.peval(ap.poly(cp['diagonal']),raw)
            A=[row+[b[i]] for i,row in enumerate(A)]+[b+[t+c]]
        records.append({'name':name,'pivots':spd(A)})
    records.append({'name':'M4_None','pivots':spd([[t]])})
    if len(records)!=12 or sum(len(v['pivots']) for v in records)!=54:raise AssertionError('incomplete sample layout')
    return {'side':curve['side'],'delta':d,'r':r,'x':y[0],'t':t,'y':y,'raw_coordinates':raw,
            'excess_over_C_div_delta':excess,'blocks':records,
            'definition':'Algebraic point from the certified root and r=2*delta/(1+sqrt(1+4*a*delta)); intervals are enclosures, not midpoint certificates.'}
