"""KKT evaluation by generic automatic differentiation of scalar elimination."""
from decimal_interval import Box,ZERO,ONE,zeros,eye,mul,transpose,solve_spd,ldlt
from second_order_ad import Jet,N

def active_constraints(model,layout,y,with_jets=True):
    coords=[Jet.variable(v,i) for i,v in enumerate(y)] if with_jets else y
    equations=[];info=[]
    for it in layout:
        idx=it['index'];perm=it['perm'];r=it['rank'];nu=it['nu'];b=model.small[idx]
        if b['name']!=it['name']:raise AssertionError('active layout name mismatch')
        mat=model.matrix(idx,coords);n=len(mat)
        if sorted(perm)!=list(range(n)) or r+nu!=n:raise AssertionError('bad Schur partition')
        work=[[mat[i][j] for j in perm] for i in perm]
        # Differentiate Gaussian elimination itself; no prederived derivative.
        for k in range(r):
            p=work[k][k]
            if (p.value if with_jets else p).low<=0:raise ValueError('active leading block not positive')
            inv=p.reciprocal()
            for i in range(k+1,n):
                factor=work[i][k]*inv
                for j in range(i,n):
                    val=work[i][j]-factor*work[k][j]
                    work[i][j]=val;work[j][i]=val
        for i in range(nu):
            for j in range(i,nu):equations.append(work[r+i][r+j])
        info.append((it,mat))
    return equations,info

def evaluate_kkt(model,data,values):
    nc=data['nc'];dim=20+nc
    if len(values)!=dim:raise ValueError('KKT input dimension mismatch')
    g,_=active_constraints(model,data['layout'],values[:20]);lam=values[20:]
    if len(g)!=nc:raise AssertionError('wrong constraint count')
    sigma=1 if data['side']=='lower' else -1
    f=[q.value for q in g]
    f.extend(Box(sigma if j==0 else 0)-sum((lam[k]*q.gradient[j] for k,q in enumerate(g)),ZERO) for j in range(20))
    jac=[q.gradient+[ZERO]*nc for q in g]
    jac.extend([-sum((lam[k]*q.h(i,j) for k,q in enumerate(g)),ZERO) for j in range(20)]+[-q.gradient[i] for q in g] for i in range(20))
    return f,jac

def primal_dual_data(model,data,values):
    """Rebuild kernel columns by positive LDL triangular solves, not old inverse."""
    y=values[:20];lam=values[20:];off=0;info=[];prank=drank=0;seen=set()
    for it in data['layout']:
        idx=it['index'];seen.add(idx);p=it['perm'];r=it['rank'];nu=it['nu']
        m=model.matrix(idx,y);n=len(m)
        if r:
            a=[[m[i][j] for j in p[:r]] for i in p[:r]]
            b=[[m[i][j] for j in p[r:]] for i in p[:r]]
            w,piv=solve_spd(a,b);e=[[-v for v in row] for row in w]+eye(nu)
        else:piv=[];e=eye(nu)
        k=zeros(n,nu)
        for i,j in enumerate(p):k[j]=e[i]
        mult=zeros(nu)
        for i in range(nu):
            for j in range(i,nu):
                mult[i][j]=mult[j][i]=lam[off]/(1 if i==j else 2);off+=1
        _,dual_piv=ldlt(mult)
        prank+=len(piv);drank+=len(dual_piv)
        info.append({'item':it,'kernel':k,'multiplier':mult,'dual':mul(mul(k,mult),transpose(k)),
                     'leading_pivots':piv,'dual_pivots':dual_piv})
    rest=[]
    for i,b in enumerate(model.small):
        if i in seen:continue
        _,piv=ldlt(model.matrix(i,y));prank+=len(piv)
        rest.append({'name':b['name'],'pivots':[p.record() for p in piv]})
    if off!=data['nc'] or prank+drank!=51:raise AssertionError('rank partition wrong')
    return info,{'primal_rank':prank,'dual_rank':drank,'active_blocks':[{'name':d['item']['name'],'positive_primal_pivots':[p.record() for p in d['leading_pivots']],'positive_dual_pivots':[p.record() for p in d['dual_pivots']]} for d in info],'inactive_blocks':rest}
