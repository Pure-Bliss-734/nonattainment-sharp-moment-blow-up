"""Optimal asymptotic completion coefficient and a Puiseux feasible family.

Every Box encloses a quantity at the *same* certified algebraic endpoint.
Identities J R = I and d_x=sigma follow algebraically from the definitions and
KKT equations. Interval residual inclusion is only a regression check, never
an identity proof. No numerical optimizer is used.
"""
from fractions import Fraction as F
from decimal_interval import Box, ZERO, ONE, zeros, eye, mul, transpose, dot, solve_spd, ldlt
from independent_model import evaluate, derivative, from_form, pc, pv, plus, times, scaling, substitute
from ad_kkt import active_constraints
from order_model import square


def matrix_add(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def matrix_scale(a,s):return [[x*s for x in r] for r in a]
def evaluate_matrix(p,coords):return [[evaluate(x,coords) for x in row] for row in p]
def congruence(a,t):return mul(mul(transpose(t),a),t)
def matrix_columns(a,rows,cols):return [[a[i][j] for j in cols] for i in rows]
def poly_path(p,star,d,v):
    """Coefficients for p(star+u^2 d+u^3 v); p has total degree <= 2."""
    if any(len(mon)>2 for mon in p):raise AssertionError('nonquadratic polynomial')
    c=[ZERO for _ in range(7)];c[0]=evaluate(p,star)
    for i in range(len(star)):
        g=evaluate(derivative(p,i),star)
        c[2]=c[2]+g*d[i];c[3]=c[3]+g*v[i]
    for mon,a in p.items():
        if len(mon)==2:
            i,j=mon;c[4]=c[4]+a*d[i]*d[j]
            c[5]=c[5]+a*(d[i]*v[j]+v[i]*d[j]);c[6]=c[6]+a*v[i]*v[j]
    return c

def path_matrices(polys,star,d,v):
    c=[[poly_path(p,star,d,v) for p in row] for row in polys]
    return [[[p[k] for p in row] for row in c] for k in range(7)]
def eval_coeffs(c,u,start=0):
    result=c[-1]
    for a in reversed(c[start:-1]):result=a+u*result
    return result

def eval_matrices(ms,u,start=0):
    n=len(ms[0]);return [[eval_coeffs([m[i][j] for m in ms],u,start) for j in range(n)] for i in range(n)]
def matvec(a,v):return [dot(row,v) for row in a]

def formal_path_checks(model):
    """Check the degree-six substitution as a rational polynomial identity.

    The endpoint and both tangent directions are still formal variables.
    Root nullspace cancellations are justified separately by the KKT equations.
    """
    U=pv(60);U2=times(U,U);U3=times(U2,U)
    sub={i:plus(plus(pv(i),times(U2,pv(20+i))),times(U3,pv(40+i))) for i in range(20)}
    powers=[pc(1)]
    for _ in range(6):powers.append(times(powers[-1],U))
    def check(p):
        if any(len(m)>2 for m in p):raise AssertionError('unexpected degree')
        c=[{} for _ in range(7)];c[0]=p
        for i in range(20):
            dp=derivative(p,i)
            c[2]=plus(c[2],times(dp,pv(20+i)))
            c[3]=plus(c[3],times(dp,pv(40+i)))
        for m,a in p.items():
            if len(m)==2:
                i,j=m
                c[4]=plus(c[4],scaling(times(pv(20+i),pv(20+j)),a))
                c[5]=plus(c[5],scaling(plus(times(pv(20+i),pv(40+j)),times(pv(40+i),pv(20+j))),a))
                c[6]=plus(c[6],scaling(times(pv(40+i),pv(40+j)),a))
        rhs={}
        for k in range(7):rhs=plus(rhs,times(c[k],powers[k]))
        if rhs!=substitute(p,sub):raise AssertionError('formal Puiseux coefficient identity failed')
    nsmall=0
    for b in model.small:
        for row in b['mat']:
            for p in row:check(p);nsmall+=1
    raw_sub={k:{m:c for m,c in p.items() if 20 not in m} for k,p in model.to_symmetric.items()}
    nraw=ncross=0
    for b in model.raw:
        for row in b['mat']:
            for p in row:check(substitute(p,raw_sub));nraw+=1
    for cp in model.data['completion'].values():
        for f in cp['cross']+[cp['diagonal']]:check(substitute(from_form(f),raw_sub));ncross+=1
    # Exact scalar identity used to cancel the constant determinant term.
    # Symbols: u,q,C,hI,R,Z,k. h0=q^2/C is cleared by multiplication by C.
    import itertools
    def pvar(i):return pv(i)
    u,q,C,hI,R,Z,k=[pvar(i) for i in range(7)]
    q2=times(q,q);u2=times(u,u)
    # C*H=q^2+u(C*hI-q^2)+u^2*C*R
    CH=plus(plus(q2,times(u,plus(times(C,hI),scaling(q2,-1)))),times(u2,times(C,R)))
    bottom=plus(times(C,plus(pc(1),scaling(u,2))),times(u2,Z))
    cross=plus(q,times(u2,k))
    lhs=plus(times(CH,bottom),scaling(times(C,times(cross,cross)),-1))
    h1C=plus(times(C,hI),scaling(q2,-1))
    # C times [D(u)/u], expressed without inverse C.
    r0=times(C,plus(times(C,hI),q2))
    r1=plus(plus(scaling(times(C,h1C),2),times(times(C,C),R)),plus(times(q2,Z),scaling(times(C,times(q,k)),-2)))
    r2=plus(scaling(times(times(C,C),R),2),times(h1C,Z))
    r3=times(C,plus(times(R,Z),scaling(times(k,k),-1)))
    rhs=times(u,plus(plus(r0,times(u,r1)),plus(times(u2,r2),times(times(u2,u),r3))))
    if lhs!=rhs:raise AssertionError('scaled completion determinant identity failed')
    return {'small_matrix_entries':nsmall,'raw_matrix_entries':nraw,'cross_and_diagonal_entries':ncross,
            'formal_endpoint_and_direction_variables':60,'scalar_determinant_identity':True}


def build_prefactor(model,data,box,infos):
    y=box[:20];raw=model.raw_coordinates(y);m=data['nc'];sigma=1 if data['side']=='lower' else -1
    completions={};C=ZERO;offset=0;locations={};Omega=ZERO
    for inf in infos:
        it=inf['item'];nu=it['nu'];ind=it['index'];W=inf['multiplier']
        ids=[]
        for i in range(nu):
            for j in range(i,nu):ids.append((i,j,offset));offset+=1
        locations[ind]=ids
        Omega=Omega+sum((W[i][i] for i in range(nu)),ZERO)
    if offset!=m or Omega.lo()<=0:raise AssertionError('invalid multiplier trace')
    for name,cp in model.data['completion'].items():
        matching=[inf for inf in infos if model.small[inf['item']['index']]['parent']==name]
        if len(matching)!=1 or matching[0]['item']['nu']!=1:
            raise AssertionError('need exactly one simple active kernel per completion block')
        inf=matching[0];sb=model.small[inf['item']['index']]
        kernel=[r[0] for r in mul([[Box(q) for q in row] for row in sb['T']],inf['kernel'])]
        q=dot(kernel,[evaluate(from_form(f),raw) for f in cp['cross']])
        w=inf['multiplier'][0][0]
        if w.lo()<=0:raise AssertionError('nonpositive completion dual weight')
        term=w*square(q);C=C+term
        completions[name]={'index':inf['item']['index'],'kernel':kernel,'cross':q,'weight':w,'contribution':term}
    if C.lo()<=0:raise AssertionError('positive prefactor missing')
    # J is the derivative of every independent active Schur entry.
    g,_=active_constraints(model,data['layout'],y)
    J=[q.gradient for q in g];Gram=mul(J,transpose(J))
    inv,piv=solve_spd(Gram,eye(m));right=mul(transpose(J),inv)
    # SPD Gram proves full row rank, hence R=J^T (J J^T)^-1 is a genuine right inverse.
    residual=mul(J,right)
    if any(not residual[i][j].contains(int(i==j)) for i in range(m) for j in range(m)):
        raise AssertionError('right inverse enclosure failed regression check')
    H0={inf['item']['index']:zeros(inf['item']['nu']) for inf in infos}
    HI={inf['item']['index']:matrix_scale(eye(inf['item']['nu']),ONE/Omega) for inf in infos}
    for name,cp in completions.items():H0[cp['index']][0][0]=square(cp['cross'])/C
    targets0=[ZERO]*m;targetsI=[ZERO]*m
    for idx,loc in locations.items():
        for i,j,k in loc:targets0[k]=H0[idx][i][j];targetsI[k]=HI[idx][i][j]
    d0=matvec(right,targets0);dI=matvec(right,targetsI)
    for d in (d0,dI):
        if not d[0].contains(sigma):raise AssertionError('normalization enclosure disagrees with KKT')
    # Exact: sigma*e_x^T=lambda^T J, lambda^T targets=1, so d_x=sigma.
    # This substitution uses that proved identity, not interval zero-inclusion.
    d0[0]=Box(sigma);dI[0]=Box(sigma)
    v=[b-a for a,b in zip(d0,dI)];v[0]=ZERO
    # Check derivative normalizations numerically enclosed, diagnostic only.
    for target,d in ((targets0,d0),(targetsI,dI)):
        if any(not (a-b).contains(0) for a,b in zip(matvec(J,d),target)):
            raise AssertionError('direction derivative enclosure failed')
    return {'C':C,'Omega':Omega,'completion':completions,'J':J,'right_inverse':right,'gram_pivots':piv,
            'H0':H0,'HI':HI,'target0':targets0,'targetI':targetsI,'d0':d0,'dI':dI,'v':v,
            'tangent_rank':m,'sigma':sigma,'normalization_proof':'lambda^T J=sigma e_x^T, lambda^T target0=lambda^T targetI=1',
            'lower_second_cross_zero_NOT_assumed':True}


def active_transform(model,inf):
    n=len(inf['kernel']);perm=inf['item']['perm'];r=inf['item']['rank'];nu=inf['item']['nu']
    T=zeros(n)
    for j,k in enumerate(perm[:r]):T[k][j]=ONE
    for i in range(n):
        for j in range(nu):T[i][r+j]=inf['kernel'][i][j]
    return T,r,nu


def full_to_raw_directions(model,pref):
    """to_symmetric is affine; its derivative transforms both path directions."""
    out=[]
    for direction in (pref['d0'],pref['v']):
        coords=[]
        for idx in range(21):
            p={m:c for m,c in model.to_symmetric[idx].items() if 20 not in m}
            if any(len(mon)>1 for mon in p):raise AssertionError('raw coordinate map not affine')
            coords.append(sum((c*direction[mon[0]] for mon,c in p.items() if mon),ZERO))
        out.append(coords)
    return out


def continuous_family(model,data,box,infos,pref,u_limit):
    """Certify every 0<u<=u_limit without sampling or division by [0,u_limit].

    Path: y=y*+u^2 d0+u^3(dI-d0); x=xi+sigma*u^2.
    Choose t=C*(1+2u)/u^2. Root and tangent equalities cancel the
    constant/linear terms before interval evaluation.
    """
    h=F(u_limit)
    if not 0<h<1:raise AssertionError('invalid u strip')
    U=Box(0,h);y=box[:20];d0=pref['d0'];v=pref['v'];C=pref['C']
    active={inf['item']['index']:inf for inf in infos};actrep=[];inactrep=[]
    for idx,b in enumerate(model.small):
        ms=path_matrices(b['mat'],y,d0,v)
        if idx not in active:
            _,piv=ldlt(eval_matrices(ms,U));inactrep.append({'name':b['name'],'pivots':piv});continue
        inf=active[idx];T,r,nu=active_transform(model,inf);cs=[congruence(m,T) for m in ms]
        RR=list(range(r));NN=list(range(r,r+nu))
        if r:
            A=matrix_columns(eval_matrices(cs,U),RR,RR)
            # E/u^2; its constant zero comes from the exact root kernel.
            E=matrix_columns(eval_matrices(cs,U,2),RR,NN)
            sol,pa=solve_spd(A,E);prod=mul(transpose(E),sol)
        else:pa=[];prod=zeros(nu)
        Drem=matrix_columns(eval_matrices(cs,U,4),NN,NN)
        R=[[a-b for a,b in zip(ar,br)] for ar,br in zip(Drem,prod)]
        # Schur/u^2=(1-u)H0+u(HI+u R). H0 is PSD by construction.
        margin=matrix_add(pref['HI'][idx],matrix_scale(R,U))
        _,pm=ldlt(margin)
        actrep.append({'name':b['name'],'rank':r,'nullity':nu,'leading_pivots':pa,
                       'HI_plus_u_R_positive_pivots':pm,'remainder':R})
    # Two augmented raw completion matrices and their crossing vectors.
    rawstar=model.raw_coordinates(y);rd0,rv=full_to_raw_directions(model,pref)
    rawby={b['name']:b for b in model.raw};reports={}
    for name,cp in model.data['completion'].items():
        item=pref['completion'][name];nu=item['kernel'];n=len(nu)
        pivot=max(range(n),key=lambda i:min(abs(nu[i].lo()),abs(nu[i].hi())) if not nu[i].contains(0) else F(0))
        if nu[pivot].contains(0):raise AssertionError('no certified nonzero kernel coordinate')
        keep=[i for i in range(n) if i!=pivot];T=zeros(n)
        for j,k in enumerate(keep):T[k][j]=ONE
        for i in range(n):T[i][-1]=nu[i]
        ms=path_matrices(rawby[name]['mat'],rawstar,rd0,rv);cs=[congruence(m,T) for m in ms]
        inds=list(range(n-1));last=[n-1]
        A=matrix_columns(eval_matrices(cs,U),inds,inds)
        E=matrix_columns(eval_matrices(cs,U,2),inds,last)
        bco=[poly_path(from_form(f),rawstar,rd0,rv) for f in cp['cross']]
        b=[eval_coeffs(co,U) for co in bco]
        f=[[b[i]] for i in keep]
        allrhs=[E[i]+f[i] for i in range(n-1)]
        sol,piv=solve_spd(A,allrhs);Einv=[row[0] for row in sol];finv=[row[1] for row in sol]
        e=[row[0] for row in E];ff=[row[0] for row in f]
        Drem=eval_coeffs([c[-1][-1] for c in cs],U,4)
        R=Drem-dot(e,Einv)
        ghat=dot(nu,[eval_coeffs(co,U,2) for co in bco])
        k=ghat-dot(e,finv)
        cco=poly_path(from_form(cp['diagonal']),rawstar,rd0,rv)
        Z=eval_coeffs(cco,U)-dot(ff,finv)
        q=item['cross'];q2=square(q);h0=q2/C;hI=ONE/pref['Omega'];h1=hI-h0
        # After elimination and diag(1/u,u) scaling, the determinant/u is:
        Dscaled=C*hI+q2 + U*(2*h1*C+R*C+h0*Z-2*q*k) + square(U)*(2*R*C+h1*Z) + U*square(U)*(R*Z-square(k))
        bottom=C*(ONE+2*U)+square(U)*Z
        if Dscaled.lo()<=0:raise AssertionError(name+': positive scaled determinant not certified')
        if bottom.lo()<=0:raise AssertionError(name+': positive bottom diagonal not certified')
        reports[name]={'kernel_pivot':pivot,'leading_pivots':piv,'R_interval':R,'cross_remainder_k':k,
                        'diagonal_remainder_Z':Z,'scaled_determinant_over_u':Dscaled,'scaled_bottom_diagonal':bottom,
                        'identity':'det(scaled completed 2x2)/u=C/Omega+q^2+O(u); explicit coefficients verified',
                        'all_positive_u_in_strip_full_PSD':True}
    return {'side':data['side'],'u_limit':h,'distance_limit':h*h,'active_blocks':actrep,'inactive_blocks':inactrep,
            'completion':reports,'exact_path':'y(u)=y_star+u^2*d0+u^3*(dI-d0), t(u)=C*(1+2u)/u^2',
            'all_continuous_parameters_certified':True,'curve_not_rounded_to_rational_midpoints':True}


def quantitative_lower(model,box,pref,consequences,parameters,side):
    """A bound valid for every full feasible point, including sqrt(delta) paths.

    It uses the re-certified Stage-17 uniform coordinate error, not any tangent
    assumption on a minimizing sequence.
    """
    y=box[:20];raw=model.raw_coordinates(y);row=parameters[side]
    AR=F(row['coordinate_A']);BR=F(row['coordinate_B']);R=F(row['R']);delta0=F(parameters['delta'])
    err=[F(0)]*5;details={};cbound=F(0)
    for name,cp in model.data['completion'].items():
        item=pref['completion'][name];nu=item['kernel'];w=item['weight'];q=item['cross']
        # Cross is affine in z, with no X*z term. Build coefficients as intervals
        # because the fixed null vector has algebraic entries.
        forms=[from_form(f) for f in cp['cross']]
        for p in forms:
            for k in range(1,21):
                if derivative(derivative(p,k),0):raise AssertionError('mixed Xz crossing not supported')
                for j in range(1,21):
                    if derivative(derivative(p,k),j):raise AssertionError('nonaffine z crossing')
        gradients=[sum((nu[i]*evaluate(derivative(forms[i],k),raw) for i in range(len(nu))),ZERO) for k in range(1,21)]
        gz=sum(c.norm() for c in gradients)
        broad=[Box(F(3,8),F(2,5))]+raw[1:]
        gx=sum((nu[i]*evaluate(derivative(forms[i],0),broad) for i in range(len(nu))),ZERO).norm()
        A=gz*AR;B=gz*BR+gx
        e=[F(0),2*q.norm()*A,2*q.norm()*B+A*A,2*A*B,B*B]
        for i in range(5):err[i]+=w.hi()*e[i]
        cpoly=from_form(cp['diagonal'])
        if any(len(m)>1 for m in cpoly):raise AssertionError('expected affine completion diagonal')
        cstar=evaluate(cpoly,raw);cz=sum(abs(cpoly.get((k,),F(0))) for k in range(1,21));cx=abs(cpoly.get((0,),F(0)))
        cbound=max(cbound,cstar.hi()+cz*R+cx*delta0)
        details[name]={'cross_A_sqrt_delta':A,'cross_B_delta':B,'root_cross_absolute_bound':q.norm(),
                       'weight_upper':w.hi(),'diagonal_upper':cstar.hi()+cz*R+cx*delta0}
    # Reconstruct L = 7p/2+3a in all original coordinates.
    L=plus(scaling(from_form(model.data['observables']['p']),F(7,2)),scaling(from_form(model.data['observables']['a']),3))
    if any(len(m)!=1 for m in L):raise AssertionError('L must be linear')
    Lz=sum(abs(L.get((k,),F(0))) for k in range(1,21));Lx=abs(L.get((0,),F(0)))
    gamma=F(consequences['global']['gamma']['upper']);f2=consequences['global']['F2'];f2norm=max(abs(F(f2['lower'])),abs(F(f2['upper'])))
    b1=gamma*Lz*AR;b2=gamma*(Lz*BR+Lx)+f2norm
    cbound=max(F(1),F((cbound.numerator+cbound.denominator-1)//cbound.denominator))
    return {'radius':delta0,'coordinate_A':AR,'coordinate_B':BR,'coordinate_R':R,'E_coefficients_in_u':err,
            'T_over_delta_upper_coefficients':[F(1),b1,b2],'common_diagonal_upper':cbound,'blocks':details,
            'universal_bound':'delta*t >= (C-E(sqrt(delta)))/(1+b1*sqrt(delta)+b2*delta)-cbar*delta',
            'all_full_feasible_sequences_not_only_differentiable_paths':True,
            'limit_inferior_at_least_C':True}


def full_sample(model,box,pref,delta):
    """Independent pointwise evaluation of all twelve completed quotient blocks.

    Samples are regression checks. Continuous validity is proved separately.
    """
    from math import isqrt
    delta=F(delta);a=isqrt(delta.numerator);b=isqrt(delta.denominator)
    if a*a!=delta.numerator or b*b!=delta.denominator or delta<=0:
        raise AssertionError('sample requires a positive rational square')
    u=F(a,b);y=[v+delta*d+delta*u*e for v,d,e in zip(box[:20],pref['d0'],pref['v'])]
    raw=model.raw_coordinates(y);t=pref['C']*(1+2*u)/delta;out=[]
    for bl in model.raw:
        M=evaluate_matrix(bl['mat'],raw);name=bl['name']
        if name in model.data['completion']:
            cp=model.data['completion'][name];cross=[evaluate(from_form(f),raw) for f in cp['cross']]
            c=evaluate(from_form(cp['diagonal']),raw)
            M=[row+[cross[i]] for i,row in enumerate(M)]+[cross+[t+c]]
        _,piv=ldlt(M);out.append({'name':name,'positive_pivots':piv})
    _,piv=ldlt([[t]]);out.append({'name':'M4_None','positive_pivots':piv})
    if len(out)!=12 or sum(len(b['positive_pivots']) for b in out)!=54:
        raise AssertionError('original completed dimensions changed')
    return {'delta':delta,'u':u,'x':y[0],'t':t,'delta_times_t':delta*t,'blocks':out,
            'algebraic_sample_enclosed_not_rational_midpoint':True}
