"""Independent Stage-20 reconstruction of Stage-19 lower and upper arguments.

Shared: exported rational pencils; certified root equations and root boxes.
New: binary interval arithmetic; fixed-kernel Schur derivatives (no AD);
coordinate-series convolution; scalar Schur-complement positivity with an
absolute perturbation bound (no Stage-19 determinant expansion).
"""
from fractions import Fraction as F
from audit_interval import I,Z,O,zeros,eye,tr,dot,mm,mv,add,sub,scale,congr,part,spd,solve,max_row_norm,overlap


def poly(p):return [(tuple(m),F(c)) for m,c in p]
def peval(p,y):
    out=Z
    for mon,c in p:
        v=I(c)
        for k in mon:v=v*y[k]
        out=out+v
    return out

def deriv(p,k):
    out=[]
    for mon,c in p:
        if k in mon:
            m=list(mon);count=m.count(k);m.remove(k);out.append((tuple(m),c*count))
    return out

def pcoeff(p,mon):return sum((c for m,c in p if m==tuple(mon)),F(0))
def matval(p,y):return [[peval(q,y) for q in r] for r in p]
def pmatrix(p):return [[poly(q) for q in r] for r in p]
def peval_series(p,coords):
    """Direct substitution by convolution, not derivatives/Hessian formulas."""
    out=[Z for _ in range(7)]
    for mon,c in p:
        v=[I(c)]
        for k in mon:
            w=[Z for _ in range(len(v)+len(coords[k])-1)]
            for i,a in enumerate(v):
                for j,b in enumerate(coords[k]):w[i+j]=w[i+j]+a*b
            v=w
        if len(v)>7:raise AssertionError('degree larger than six')
        for i,a in enumerate(v):out[i]=out[i]+a
    return out

def matrix_series(p,coords):
    n=len(p);v=[[peval_series(q,coords) for q in r] for r in p]
    return [[[v[i][j][k] for j in range(n)] for i in range(n)] for k in range(7)]

def horner(a,u):
    v=Z
    for c in reversed(a):v=c+u*v
    return v

def value_series(ms,u,start=0):
    return [[horner([m[i][j] for m in ms[start:]],u) for j in range(len(ms[0]))] for i in range(len(ms[0]))]

def enc_raw(data,y):return [peval(poly(p),y) for p in data['to_raw']]
def raw_dir(data,d):
    return [sum((c*d[m[0]] for m,c in poly(p) if len(m)==1),Z) for p in data['to_raw']]

def reconstruct(data,side):
    sd=data['sides'][side];sigma=1 if side=='lower' else -1
    star=[I(F(c)-F(sd['radius']),F(c)+F(sd['radius'])) for c in sd['center']]
    y=star[:20];lamb=star[20:];raw=enc_raw(data,y)
    small=[dict(b,mat=pmatrix(b['mat']),T=[[I(q) for q in r] for r in b['T']]) for b in data['small']]
    rawmat={b['name']:pmatrix(b['mat']) for b in data['raw']}
    J=[];infos=[];offset=0;Omega=Z
    for it in sd['layout']:
        b=small[it['index']]; assert b['name']==it['name']
        n=len(b['mat']);r=it['rank'];nu=it['nu'];perm=it['perm']
        assert sorted(perm)==list(range(n)) and r+nu==n
        M=matval(b['mat'],y)
        A=part(M,perm[:r],perm[:r]);B=part(M,perm[:r],perm[r:]);piv=spd(A)
        S=solve(A,B) if r else []
        K=zeros(n,nu)
        for i,k in enumerate(perm[:r]):K[k]=[-v for v in S[i]]
        for i,k in enumerate(perm[r:]):K[k][i]=O
        W=zeros(nu);locations=[]
        for i in range(nu):
            for j in range(i,nu):
                W[i][j]=W[j][i]=lamb[offset]/(1 if i==j else 2)
                locations.append((i,j,offset));offset+=1
        wp=spd(W);Omega=Omega+sum((W[i][i] for i in range(nu)),Z)
        # At an exact root, d(Schur)=K^T (dM) K. Derive J this way,
        # without differentiating elimination or reusing the Stage-19 AD.
        tensors=[congr([[peval(deriv(p,k),y) for p in row] for row in b['mat']],K) for k in range(20)]
        for i,j,_ in locations:J.append([ten[i][j] for ten in tensors])
        T=zeros(n)
        for j,k in enumerate(perm[:r]):T[k][j]=O
        for i in range(n):
            for j in range(nu):T[i][r+j]=K[i][j]
        rawK=mm(b['T'],K)
        infos.append({'index':it['index'],'name':b['name'],'parent':b['parent'],'rank':r,'nu':nu,
                      'K':K,'rawK':rawK,'W':W,'T':T,'locations':locations,
                      'leading_pivots':piv,'weight_pivots':wp})
    assert offset==sd['nc']
    # Diagnostics only: exact stationarity is a shared root equation.
    res=[sum((lamb[k]*J[k][j] for k in range(len(J))),Z)-I(sigma if j==0 else 0) for j in range(20)]
    assert all(v.contains(0) for v in res),'KKT normalization inconsistent with independent derivatives'
    cps={};C=Z
    for name,cp in data['completion'].items():
        choices=[i for i in infos if i['parent']==name]
        assert len(choices)==1 and choices[0]['nu']==1,'completion requires simple kernel'
        inf=choices[0];v=[r[0] for r in inf['rawK']]
        cross=dot(v,[peval(poly(p),raw) for p in cp['cross']]);w=inf['W'][0][0]
        assert w.lo()>0
        term=w*cross.sq();C=C+term
        cps[name]={'kernel':v,'cross':cross,'weight':w,'contribution':term,'index':inf['index']}
    assert C.lo()>0 and Omega.lo()>0
    Gram=mm(J,tr(J));gp=spd(Gram)
    inv=solve(Gram,eye(len(J)));R=mm(tr(J),inv)
    id=mm(J,R);assert all(id[i][j].contains(int(i==j)) for i in range(len(J)) for j in range(len(J)))
    h0=[Z]*len(J);hi=[Z]*len(J);H0={};HI={}
    for inf in infos:
        v=zeros(inf['nu']);vi=scale(eye(inf['nu']),O/Omega)
        if inf['parent'] in cps:v[0][0]=cps[inf['parent']]['cross'].sq()/C
        for i,j,k in inf['locations']:h0[k]=v[i][j];hi[k]=vi[i][j]
        H0[inf['index']]=v;HI[inf['index']]=vi
    assert dot(lamb,h0).contains(1) and dot(lamb,hi).contains(1)
    d0=mv(R,h0);dI=mv(R,hi)
    assert d0[0].contains(sigma) and dI[0].contains(sigma)
    # The exact identity follows from lambda^T J=sigma e_x^T,
    # J R=I and lambda^T h=1; inclusion is NOT used to prove equality.
    d0[0]=I(sigma);dI[0]=I(sigma)
    v=[b-a for a,b in zip(d0,dI)];v[0]=Z
    return {'side':side,'sigma':sigma,'star':y,'rawstar':raw,'small':small,'raw':rawmat,
            'infos':infos,'C':C,'Omega':Omega,'completion':cps,'J':J,'R':R,
            'gram_pivots':gp,'h0':h0,'hI':hi,'H0':H0,'HI':HI,'d0':d0,'dI':dI,'v':v,
            'normalization_diagnostics':res}


def certify_upper(data,pref,u0):
    u0=F(u0)
    if not 0<u0<1:raise AssertionError('invalid positive strip')
    U=I(0,u0);C=pref['C'];OI=O/pref['Omega'];active={v['index']:v for v in pref['infos']}
    coords=[[a,Z,b,c] for a,b,c in zip(pref['star'],pref['d0'],pref['v'])]
    projected=[]
    for i,b in enumerate(pref['small']):
        ms=matrix_series(b['mat'],coords)
        if i not in active:
            piv=spd(value_series(ms,U));projected.append({'name':b['name'],'type':'inactive','pivots':piv});continue
        inf=active[i];cs=[congr(m,inf['T']) for m in ms];r=inf['rank'];n=inf['nu'];ids=list(range(r));tail=list(range(r,r+n))
        A=part(value_series(cs,U),ids,ids);pa=spd(A)
        E=part(value_series(cs,U,2),ids,tail)
        product=mm(tr(E),solve(A,E)) if r else zeros(n)
        rem=sub(part(value_series(cs,U,4),tail,tail),product)
        # A scalar norm bound is a different certificate from the old
        # matrix LDL of HI+uR: HI=I/Omega and ||R||2 <= ||R||infty.
        norm=max_row_norm(rem);margin=OI.lo()-u0*norm
        assert margin>0,b['name']+': uniform scalar remainder bound fails'
        projected.append({'name':b['name'],'type':'active','leading_pivots':pa,
                          'remainder_row_norm_upper':norm,'HI_minus_u_normR_lower':margin})
    # Both augmentations: Schur complement in the positive bottom scalar.
    rawcoords=[[a,Z,b,c] for a,b,c in zip(pref['rawstar'],raw_dir(data,pref['d0']),raw_dir(data,pref['v']))]
    completed={}
    for name,cp in data['completion'].items():
        v=pref['completion'][name]['kernel'];n=len(v)
        choices=[j for j in range(n) if not v[j].contains(0)]
        pivot=max(choices,key=lambda j:min(abs(v[j].lo()),abs(v[j].hi())))
        keep=[j for j in range(n) if j!=pivot];T=zeros(n)
        for j,k in enumerate(keep):T[k][j]=O
        for i in range(n):T[i][-1]=v[i]
        cs=[congr(m,T) for m in matrix_series(pref['raw'][name],rawcoords)]
        ids=list(range(n-1));A=part(value_series(cs,U),ids,ids);piv=spd(A)
        E=[row[-1] for row in value_series(cs,U,2)[:-1]]
        bc=[peval_series(poly(p),rawcoords) for p in cp['cross']]
        b=[horner(p,U) for p in bc];f=[b[i] for i in keep]
        so=solve(A,[[e,ff] for e,ff in zip(E,f)])
        ai_e=[r[0] for r in so];ai_f=[r[1] for r in so]
        rem=horner([m[-1][-1] for m in cs[4:]],U)-dot(E,ai_e)
        cross_rem=dot(v,[horner(p[2:],U) for p in bc])-dot(E,ai_f)
        c=horner(peval_series(poly(cp['diagonal']),rawcoords),U)
        diag_rem=c-dot(f,ai_f)
        q=pref['completion'][name]['cross'];q2=q.sq();h0=q2/C
        bottom=C*(O+2*U)+U.sq()*diag_rem
        assert bottom.lo()>0
        # Independent scalar identity:
        # (H - cross^2/bottom)/u = hI+h0 + u*B,
        # B=R+[h0*(1-2u)*Z-4*q^2-2*q*k-u^2*k^2]/bottom.
        correction=rem.norm()+(
            h0.norm()*(1+2*u0)*diag_rem.norm()+4*q2.hi()+
            2*q.norm()*cross_rem.norm()+u0*u0*cross_rem.norm()**2)/bottom.lo()
        margin=OI.lo()+h0.lo()-u0*correction
        assert margin>0,name+': scalar Schur lower bound fails'
        completed[name]={'pivot':pivot,'leading_pivots':piv,'R':rem,'k':cross_rem,'Z':diag_rem,
                         'bottom_lower':bottom.lo(),'correction_absolute_upper':correction,
                         'Schur_over_u_lower':margin,'cross_squared_divided_by_C':h0}
    return {'side':pref['side'],'u_limit':u0,'distance_limit':u0*u0,
            'projected':projected,'completion':completed,
            'criterion':'independent scalar Schur-complement lower bound; no determinant-over-u expansion',
            'covers_every_positive_u_in_strip':True}


def certify_lower(data,pref):
    side=pref['side'];par=data['neighborhood_parameters'];p=par[side]
    ar=F(p['coordinate_A']);br=F(p['coordinate_B']);rad=F(p['R']);delta=F(par['delta'])
    raw=pref['rawstar'];wide=[I(F(3,8),F(2,5))]+raw[1:]
    e1=e2=F(0);details={};cmax=F(0)
    for name,cp in data['completion'].items():
        item=pref['completion'][name];nu=item['kernel'];q=item['cross'];w=item['weight']
        polys=[poly(f) for f in cp['cross']]
        for pol in polys:
            assert all(len(m)<=2 and (len(m)<2 or m==(0,0)) for m,c in pol),'cross must be affine in z, no xz'
        gs=[sum((v*peval(deriv(pol,k),raw) for v,pol in zip(nu,polys)),Z) for k in range(1,21)]
        gz=sum(v.norm() for v in gs)
        gx=sum((v*peval(deriv(pol,0),wide) for v,pol in zip(nu,polys)),Z).norm()
        aa=gz*ar;bb=gz*br+gx
        # (q+e)^2 >= q^2 - 2|q||e|. No need to subtract e^2.
        e1+=2*w.hi()*q.norm()*aa;e2+=2*w.hi()*q.norm()*bb
        cp0=poly(cp['diagonal']);assert all(len(m)<=1 for m,c in cp0)
        c0=peval(cp0,raw)
        cb=c0.hi()+sum(abs(c) for m,c in cp0 if len(m)==1 and m[0]>0)*rad+abs(pcoeff(cp0,(0,)))*delta
        cmax=max(cmax,cb)
        details[name]={'cross_A':aa,'cross_B':bb,'diagonal_upper':cb,'weight':w,'q':q}
    assert cmax<1,'unit common diagonal upper bound not verified'
    # Recompute F2 and gamma from independently lifted dual matrices.
    L={}
    for name,mult in [('p',F(7,2)),('a',F(3))]:
        for m,c in poly(data['observables'][name]):L[m]=L.get(m,F(0))+c*mult
    L={m:c for m,c in L.items() if c};assert all(len(m)==1 for m in L)
    lx=abs(L.get((0,),F(0)));lz=sum(abs(c) for m,c in L.items() if m[0]>0)
    F2=Z;gamma=None;mixed_checks=0
    for inf in pref['infos']:
        D=mm(mm(inf['rawK'],inf['W']),tr(inf['rawK']));pols=pref['raw'][inf['parent']]
        for i,row in enumerate(pols):
            for j,pol in enumerate(row):F2=F2+D[i][j]*pcoeff(pol,(0,0))
        if inf['parent']=='single_Q':
            assert len(D)==1;gamma=D[0][0]
    if gamma is None:
        # Exported naming is explicit; do not guess silently.
        names=[i['parent'] for i in pref['infos']]
        raise AssertionError('Q scalar block name missing: '+str(names))
    assert gamma.lo()>0
    for name,pols in pref['raw'].items():
        for i,row in enumerate(pols):
            for j,pol in enumerate(row):
                for k in range(1,21):
                    expected=L.get((k,),F(0)) if name=='single_Q' else F(0)
                    assert pcoeff(pol,(0,k))==expected,'wrong mixed coefficient'
                    mixed_checks+=1
    F2=F2-gamma*L.get((0,),F(0))
    b1=gamma.hi()*lz*ar;b2=gamma.hi()*(lz*br+lx)+F2.norm()
    # Explicit O(sqrt(delta)) modulus on a smaller strip u<=u0 from upper proof.
    return {'radius':delta,'e_coefficients':[F(0),e1,e2],
            'denominator_coefficients':[F(1),b1,b2],'common_diagonal_upper':F(1),
            'computed_diagonal_upper':cmax,'gamma':gamma,'F2':F2,'blocks':details,
            'raw_mixed_coefficient_checks':mixed_checks,
            'inequality':'delta*t >= [C-e1*sqrt(delta)-e2*delta]/[1+b1*sqrt(delta)+b2*delta]-delta',
            'covers_arbitrary_feasible_sequences_without_smoothness_assumption':True}


def sample_full(data,pref,delta):
    from math import isqrt
    d=F(delta);u=F(isqrt(d.numerator),isqrt(d.denominator))
    assert u*u==d and d>0
    y=[s+d*a+d*u*b for s,a,b in zip(pref['star'],pref['d0'],pref['v'])]
    raw=enc_raw(data,y);t=pref['C']*(1+2*u)/d;records=[]
    for name,p in pref['raw'].items():
        A=matval(p,raw)
        if name in data['completion']:
            cp=data['completion'][name];b=[peval(poly(p),raw) for p in cp['cross']];c=peval(poly(cp['diagonal']),raw)
            A=[row+[b[i]] for i,row in enumerate(A)]+[b+[t+c]]
        records.append({'name':name,'pivots':spd(A)})
    records.append({'name':'M4_None','pivots':spd([[t]])})
    assert len(records)==12 and sum(len(b['pivots']) for b in records)==54
    return {'delta':d,'u':u,'x':y[0],'t':t,'scaled_cost':d*t,'blocks':records}
