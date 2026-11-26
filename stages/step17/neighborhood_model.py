"""Exact interval constants for a quantitative optimal-face error bound.

This module imports only the verified Step-16 interval and polynomial layers.
All decisions use Fraction or outward-rounded rational intervals.
"""
from fractions import Fraction as F
import json
from prepare_source import prepare
SRC,_ = prepare()
from decimal_interval import (Box,ZERO,ONE,eye,zeros,mul,transpose,dot,solve_spd,
                              display)
from independent_model import Model, derivative, evaluate, from_form
from ad_kkt import primal_dual_data

def inverse(a):
    """Outward interval Gauss-Jordan inverse with nonzero pivot checks."""
    n=len(a)
    if any(len(row)!=n for row in a):raise ValueError('not square')
    w=[list(row)+eye(n)[i] for i,row in enumerate(a)]
    pivots=[]
    for j in range(n):
        choices=[i for i in range(j,n) if w[i][j].excludes_zero()]
        if not choices:raise ValueError('inverse pivot could be zero')
        k=max(choices,key=lambda i:min(abs(w[i][j].lo()),abs(w[i][j].hi())))
        w[j],w[k]=w[k],w[j]
        p=w[j][j];pivots.append(p)
        w[j]=[q/p for q in w[j]]
        for i in range(n):
            if i==j:continue
            fac=w[i][j]
            if fac.zero():continue
            w[i]=[q-fac*r for q,r in zip(w[i],w[j])]
        # The known exact eliminated part is the identity. Resetting these
        # entries discards only interval dependency, not actual uncertainty.
        for i in range(n):w[i][j]=ONE if i==j else ZERO
    return [row[n:] for row in w],pivots

def upward(q,scale=1):
    """Round a nonnegative rational upward to a multiple of 1/scale."""
    q=F(q)
    if q<0:raise ValueError('negative upper bound')
    a=q*scale
    return F(-(-a.numerator//a.denominator),scale)

def bound_data(model,data,values,infos):
    raw=model.raw_coordinates(values[:20]); broad=Box(F(3,8),F(2,5))
    xr=[broad]+raw[1:]
    byname={b['name']:b for b in model.raw}
    rows=[];H=[];all_k={}
    for inf in infos:
        small=model.small[inf['item']['index']]
        k=mul([[Box(q) for q in row] for row in small['T']],inf['kernel'])
        all_k[small['name']]=k
        if small['name']=='single_Q_plus':
            if len(k)!=1 or len(k[0])!=1 or k[0][0].lo()!=1 or k[0][0].hi()!=1:
                raise AssertionError('scalar Q kernel normalization changed')
        winv,_=solve_spd(inf['multiplier'],eye(len(inf['multiplier'])))
        parent=byname[small['parent']]['mat'];n=len(parent);nu=len(k[0])
        for i in range(n):
            for j in range(nu):
                vec=[row[j] for row in k]
                hj=[dot([evaluate(derivative(p,v),raw) for p in parent[i]],vec)
                    for v in range(1,21)]
                H.append(hj)
                diag=parent[i][i]
                u=evaluate(diag,xr).norm()
                v=sum(evaluate(derivative(diag,a),xr).norm() for a in range(1,21))
                # h_l(X,z)-h_l(xi,z) <= |X-xi| (T0_l + T1_l ||z-z*||inf)
                t0=dot([evaluate(derivative(p,0),xr) for p in parent[i]],vec).norm()
                t1=sum(dot([evaluate(derivative(derivative(p,0),a),xr)
                             for p in parent[i]],vec).norm() for a in range(1,21))
                alpha=winv[j][j].hi()
                if alpha<=0:raise AssertionError('positive inverse diagonal missing')
                rows.append({'label':[small['name'],i,j], 'S0':2*alpha*u,'S1':2*alpha*v,
                             'T0':t0,'T1':t1,'inverse_weight':alpha,
                             'diagonal_root_radius_bound':u,'diagonal_coordinate_slope_bound':v})
    frozen=json.loads((SRC/f"inputs/face_rows_{data['side']}.json").read_text())
    if [r['label'] for r in rows]!=frozen['row_labels']:raise AssertionError('row order')
    idx=frozen['indices'];h=[H[i] for i in idx]
    inv,piv=inverse(h)
    for prod in (mul(h,inv),mul(inv,h)):
        if any(not prod[i][j].contains(int(i==j)) for i in range(20) for j in range(20)):
            raise AssertionError('inverse identity enclosure check failed')
    selected=[rows[i] for i in idx]
    mu=max(sum(q.norm() for q in r) for r in inv)
    s0=max(r['S0'] for r in selected);s1=max(r['S1'] for r in selected)
    t0=max(r['T0'] for r in selected);t1=max(r['T1'] for r in selected)
    nu=[r[0] for r in all_k['M0_1_plus']]
    cp=model.data['completion']['M0_1']
    polys=[from_form(f) for f in cp['cross']]
    atroot=dot(nu,[evaluate(p,raw) for p in polys])
    lz=sum(dot(nu,[evaluate(derivative(p,a),xr) for p in polys]).norm() for a in range(1,21))
    lx=dot(nu,[evaluate(derivative(p,0),xr) for p in polys]).norm()
    mixed=sum(dot(nu,[evaluate(derivative(derivative(p,0),a),xr)
                     for p in polys]).norm() for a in range(1,21))
    # Cross vectors have no X*z, and are affine in z. Checking this makes the
    # bound |s(X,z)-s(xi,z*)| <= Lz*r+Lx*distance exact globally in z.
    if mixed:raise AssertionError('unexpected cross parameter-coordinate product')
    for p in polys:
        if any(derivative(derivative(p,i),j) for i in range(1,21) for j in range(1,21)):
            raise AssertionError('cross not affine')
    weight=next(inf['multiplier'][0][0] for inf in infos if inf['item']['name']=='M0_1_plus')
    return {'raw_root':raw,'nu':nu,'cross_at_root':atroot,'cross_Lz':lz,'cross_Lx':lx,
            'weight':weight,'selected_indices':idx,'selected_rows':selected,'H':h,
            'H_inverse':inv,'inverse_pivots':piv,'mu':mu,'S0':s0,'S1':s1,'T0':t0,'T1':t1}

def main():
    model=Model()
    for side in ['lower','upper']:
        data=json.loads((SRC/f'inputs/kkt_{side}.json').read_text())
        r=F(data['radius']);values=[Box(F(q)-r,F(q)+r) for q in data['center']]
        infos,_=primal_dual_data(model,data,values)
        b=bound_data(model,data,values,infos)
        print(side,flush=True)
        for key in ['mu','S0','S1','T0','T1','cross_Lz','cross_Lx']:
            print(key,display(b[key],15), 'ceil=',upward(b[key]),flush=True)
        print('cross',b['cross_at_root'],'weight',b['weight'],'c0',b['raw_root'][16]/2,flush=True)
if __name__=='__main__':main()

def sqrt_up(q,scale=10**6):
    """Rational upper bound on sqrt(q), with an exact squared check."""
    from math import isqrt
    q=F(q)
    if q<0:raise ValueError('negative radicand')
    n=q.numerator*scale*scale;d=q.denominator
    k=isqrt(n//d)
    if k*k*d<n:k+=1
    out=F(k,scale)
    if out*out<q:raise AssertionError('sqrt enclosure failure')
    return out

def rowwise_bounds(model,b,R):
    inv=b['H_inverse'];rows=b['selected_rows'];raw=b['raw_root'];nu=b['nu']
    # Q itself is scalar and satisfies 0<=Q<=support/gamma, avoiding a
    # needless square-root loss for its annihilation row.
    ss=[];tt=[]
    for row in rows:
        if row['label'][0]=='single_Q_plus':
            ss.append(F(0));tt.append(2*row['inverse_weight']+row['T0']+row['T1']*R)
        else:
            ss.append(sqrt_up(row['S0']+row['S1']*R))
            tt.append(row['T0']+row['T1']*R)
    coordA=max(sum(q.norm()*s for q,s in zip(row,ss)) for row in inv)
    coordB=max(sum(q.norm()*t for q,t in zip(row,tt)) for row in inv)
    cp=model.data['completion']['M0_1'];polys=[from_form(f) for f in cp['cross']]
    grad=[dot(nu,[evaluate(derivative(p,k),raw) for p in polys]) for k in range(1,21)]
    left=mul([grad],inv)[0]
    crossA=sum(q.norm()*s for q,s in zip(left,ss))
    crossB=b['cross_Lx']+sum(q.norm()*t for q,t in zip(left,tt))
    return {'R':R,'row_sqrt_bounds':ss,'row_linear_bounds':tt,'coord_A':coordA,'coord_B':coordB,
            'cross_left_inverse':left,'cross_A':crossA,'cross_B':crossB}

def conditions(b,q,delta):
    # Prove coordinate norm < R by contradiction for norm>=R; then use the
    # affine cross identity to obtain the positive cross margin.
    R=q['R'];a=q['coord_A'];bb=q['coord_B'];A=q['cross_A'];B=q['cross_B']
    first=R-bb*delta
    ok1=first>0 and a*a*delta<first*first
    gap=b['cross_at_root'].lo()-F(1,500)-B*delta
    ok2=gap>0 and A*A*delta<gap*gap
    ok3=(b['raw_root'][16]/2).hi()+R/2<F(1,10)
    return ok1,ok2,ok3
