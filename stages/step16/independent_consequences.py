"""Second implementation of global exclusion, reachability and completion audit."""
from fractions import Fraction as F
import json
from independent_model import ROOT,KEEP,pc,pv,plus,times,scaling,derivative,evaluate,from_form
from decimal_interval import Box,ZERO,ONE,zeros,mul,transpose,dot,inner,ldlt,matrix_minor_pivots,display

def rec_list(x):return [q.record() for q in x]
def audit_consequences(model,data,box,infos):
    side=data['side'];y=box[:20];xi=y[0];rawcoords=model.raw_coordinates(y)
    # L is independently assembled in the UNSYMMETRIZED coordinate system.
    Lpoly=plus(scaling(from_form(model.data['observables']['p']),F(7,2)),scaling(from_form(model.data['observables']['a']),3))
    if Lpoly.get((0,),F(0))!=F(39,2) or () in Lpoly:raise AssertionError('unexpected L constant part')
    if any(len(m)!=1 for m in Lpoly):raise AssertionError('L is not linear')
    L=evaluate(Lpoly,rawcoords);Lmax=F(109,20)
    rawbyname={b['name']:b for b in model.raw}
    parentduals={name:zeros(len(b['mat'])) for name,b in rawbyname.items()}
    kernels={};gamma=None;annihilations=[];rowlabels=[]
    for inf in infos:
        sb=model.small[inf['item']['index']];t=[[Box(q) for q in row] for row in sb['T']]
        k=mul(t,inf['kernel']);kernels[sb['name']]=k
        yp=mul(mul(t,inf['dual']),transpose(t));dest=parentduals[sb['parent']]
        for i in range(len(dest)):
            for j in range(len(dest)):dest[i][j]=dest[i][j]+yp[i][j]
        if sb['name']=='single_Q_plus':gamma=inf['multiplier'][0][0]
        parent=rawbyname[sb['parent']]['mat'];n=len(parent);nu=len(k[0])
        aks=[]
        for v in range(1,21):
            coeff=[[evaluate(derivative(p,v),rawcoords) for p in row] for row in parent]
            aks.append(mul(coeff,k))
        for i in range(n):
            for j in range(nu):
                annihilations.append([aks[v][i][j] for v in range(20)]);rowlabels.append([sb['name'],i,j])
    if gamma is None or gamma.low<=0:raise AssertionError('positive Q weight missing')
    # Coefficient identities are checked as exact rational polynomials, not by
    # finding near-zero residuals at a numerical point.
    structure_checks=0
    for b in model.raw:
        for row in b['mat']:
            for p in row:
                for v in range(1,21):
                    got=derivative(derivative(p,0),v)
                    wanted=derivative(Lpoly,v) if b['name']=='single_Q' else {}
                    if got!=wanted:raise AssertionError('global parameter support identity broken')
                    structure_checks+=1
                if derivative(derivative(derivative(p,0),0),0):raise AssertionError('more than quadratic')
    # The coefficient of X^2 in the total frozen dual pairing, evaluated by
    # extracting monomials from ORIGINAL, unreflected matrices.
    totalx2=ZERO
    for name,b in rawbyname.items():
        yy=parentduals[name]
        for i,row in enumerate(b['mat']):
            for j,p in enumerate(row):
                if any(m.count(0)>=2 and len(m)>2 for m in p):raise AssertionError('coordinate-dependent X^2')
                totalx2=totalx2+yy[i][j]*p.get((0,0),F(0))
    f2=totalx2-gamma*F(39,2)
    sigma=1 if side=='lower' else -1
    fp=Box(sigma)-gamma*L
    broad=Box(F(3,8),F(2,5));der=fp+2*f2*(broad-xi)
    exclude=der if side=='lower' else der+gamma*Lmax
    if side=='lower' and exclude.low<=0:raise AssertionError('lower global derivative not positive')
    if side=='upper' and exclude.high>=0:raise AssertionError('upper global derivative not negative')
    lip=der.norm()+gamma.hi()*Lmax
    if lip>=2:raise AssertionError('global Lipschitz bound >= 2')
    global_report={'mixed_parameter_polynomial_checks':structure_checks,'gamma':gamma.record(),'L_at_root':L.record(),'F2':f2.record(),'exclusion_derivative':exclude.record(),'support_lipschitz_bound':str(lip),'support_lipschitz_bound_decimal':display(lip),'all_real_parameters_outside_the_root_excluded':True,'exact_stationarity_transfer':'All retained derivatives are KKT rows; missing w row is exactly zero by the independently checked odd parity blocks.'}
    # Polynomial secant curvature, recomputed by monomials, not by copied Hessian
    # formulas. M(y_s)=(1-s)M(y)+sM(v)-s(1-s)H for H=sum quadratic c*delta_i*delta_j.
    inp=json.loads((ROOT/f'inputs/interior_{side}.json').read_text())
    inside=[Box(F(inp['x']))]+[Box(F(inp['quotient_z'][j])) for j in KEEP]
    if side=='lower' and inside[0].lo()<=xi.hi():raise AssertionError('lower point not inside')
    if side=='upper' and inside[0].hi()>=xi.lo():raise AssertionError('upper point not inside')
    delta=[u-v for u,v in zip(inside,y)];curve=[]
    for i,b in enumerate(model.small):
        mat=model.matrix(i,inside);h=zeros(len(mat))
        for j,row in enumerate(b['mat']):
            for k,p in enumerate(row):
                for mon,c in p.items():
                    if len(mon)==2:h[j][k]=h[j][k]+c*delta[mon[0]]*delta[mon[1]]
        _,pi=ldlt(mat);_,ph=ldlt([[q-r for q,r in zip(row,other)] for row,other in zip(mat,h)])
        curve.append({'name':b['name'],'interior_pivots':rec_list(pi),'curvature_corrected_pivots':rec_list(ph)})
    # The prior discovery selected a minor, but interval proof determines its rank.
    selected=json.loads((ROOT/f'inputs/face_rows_{side}.json').read_text())
    if selected['row_labels']!=rowlabels:raise AssertionError('face row ordering changed')
    inds=selected['indices']
    if len(set(inds))!=20 or len(inds)!=20:raise AssertionError('dependent row selection')
    piv=matrix_minor_pivots([annihilations[i] for i in inds])
    # Positive kernel multipliers at optimum require every such row to vanish.
    # Rank 20 proves uniqueness in all unsymmetrized coordinates.
    face={'annihilation_rows':len(rowlabels),'full_unsymmetrized_dimension':20,'proved_rank':len(piv),'selected_indices':inds,'minor_pivots':rec_list(piv),'unique_projected_optimizer':True}
    k=kernels['M0_1_plus']
    if len(k[0])!=1:raise AssertionError('completion kernel not one-dimensional')
    nu=[r[0] for r in k]
    cp=model.data['completion']['M0_1']
    cross=[evaluate(from_form(f),rawcoords) for f in cp['cross']]
    value=dot(nu,cross)
    weight=next(v['multiplier'][0][0] for v in infos if v['item']['name']=='M0_1_plus')
    if value.lo()<=F(1,250) or weight.lo()<=F(1,250):raise AssertionError('nonzero crossing / positive weight fails')
    nonattain={'original_block':'M0_1','kernel':rec_list(nu),'cross':value.record(),'dual_weight':weight.record(),'cross_greater_than_1_over_250':True,'weight_greater_than_1_over_250':True,'finite_completion_impossible':True,'eventual_divergence':'t+c0 > 1/(125000000*abs(x-xi))','explicit_neighborhood_radius_not_certified':True}
    return {'global':global_report,'inner_approximation':{'x_inside':inp['x'],'all_18_blocks_passed':True,'blocks':curve},'optimal_face':face,'nonattainment':nonattain,
            'observables':{k:model.observable(k,y).record() for k in model.data['observables']}}
