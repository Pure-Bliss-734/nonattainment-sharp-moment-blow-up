#!/usr/bin/env python3
"""Coefficientwise checks of asymptotic identities used in the paper.

Only finite polynomial identities are proved by this script. No endpoint
interval, asymptotic quantifier, or analytic theorem is established here.
"""
from __future__ import annotations
from fractions import Fraction
from pathlib import Path
import argparse, json

class Poly:
    """Sparse commutative polynomial ring over the rationals."""
    def __init__(self, terms=None):
        self.terms={k:Fraction(v) for k,v in (terms or {}).items() if v}
    @staticmethod
    def scalar(a): return Poly({():Fraction(a)})
    @staticmethod
    def var(name): return Poly({(name,):Fraction(1)})
    def __add__(self,other):
        if not isinstance(other,Poly): other=Poly.scalar(other)
        q=dict(self.terms)
        for k,v in other.terms.items():q[k]=q.get(k,Fraction(0))+v
        return Poly(q)
    __radd__=__add__
    def __neg__(self):return Poly({k:-v for k,v in self.terms.items()})
    def __sub__(self,other):return self+-other
    def __rsub__(self,other):return -self+other
    def __mul__(self,other):
        if not isinstance(other,Poly): other=Poly.scalar(other)
        q={}
        for k,v in self.terms.items():
            for l,w in other.terms.items():
                m=tuple(sorted(k+l));q[m]=q.get(m,Fraction(0))+v*w
        return Poly(q)
    __rmul__=__mul__
    def __pow__(self,p):
        if not isinstance(p,int) or p<0:raise ValueError('nonnegative integer powers only')
        q=Poly.scalar(1)
        for _ in range(p):q=q*self
        return q
    def __eq__(self,other):return isinstance(other,Poly) and self.terms==other.terms

def main(output:Path):
    v=lambda names:tuple(Poly.var(n) for n in names.split())
    checks=[];controls=[]
    def same(name,a,b):
        if a!=b:raise RuntimeError('failed identity: '+name)
        checks.append({'name':name,'expanded_terms':len(a.terms),'coefficientwise_equal':True})
    def different(name,a,b):
        if a==b:raise RuntimeError('negative control failed: '+name)
        controls.append({'name':name,'nonzero_residual_terms':len((a-b).terms)})
    C,d,th,r,k,t=v('C delta theta r k t')
    lhs=d*C*t+d*k-C*C*(1+th)-2*C*d*r+C*C*th
    rhs=C*d*t-C*C+d*k-2*C*d*r
    same('completed square, after multiplication by delta',lhs,rhs)
    bad=d*C*t+d*k-Fraction(3,2)*C*C*(1+th)-2*C*d*r+C*C*th
    different('incorrect three-quarter partner coefficient',bad,rhs)

    a,r,Q0,Q1,Q2=v('a r Q0 Q1 Q2')
    left=(Q0+2*r*Q1+r*r*Q2)*(Q0+a*r)-(Q0+r*Q1)**2
    right=a*r*(Q0+a*r)+r*r*(2*a*Q1+Q0*Q2-Q1*Q1-a*a+a*r*Q2)
    same('rank-one inverse, after clearing r^2(Q0+ar)',left,right)
    different('missing negative a-square correction',left,right+r*r*a*a)

    C,r,q,h,R,mu,zeta,k=v('C r q h R mu zeta k')
    bottom=C+r*(mu+zeta)
    left=(q*q+C*r*h+C*r*r*R)*bottom-C*(q+r*k)**2
    right=r*(C*bottom*(h+r*R)+q*q*(mu+zeta)-2*C*q*k-C*r*k*k)
    same('quadratic completion Schur numerator',left,right)
    different('missing quadratic cross correction',left,right+C*r*r*k*k)

    C,u,q,h,R,zeta,k=v('C u q h R zeta k')
    left=(q*q+u*(C*h-q*q)+C*u*u*R)*(C*(1+2*u)+u*u*zeta)-C*(q+u*u*k)**2
    right=u*(C*(C*h+q*q)+u*(2*C*(C*h-q*q)+C*C*R+q*q*zeta-2*C*q*k)
             +u*u*(2*C*C*R+(C*h-q*q)*zeta)+u*u*u*C*(R*zeta-k*k))
    same('sharp-family scaled determinant',left,right)
    different('missing constant positive h contribution',left,right-u*C*C*h)

    C,r,kappa,mu=v('C r kappa mu')
    left=C*(1+kappa*r)+mu*r*(1+kappa*r)-C
    right=mu*r*(1+kappa*r)+C*kappa*r
    same('radius-distance cost conversion',left,right)
    different('omitted radius curvature',left,mu*r*(1+kappa*r))
    report={'status':'passed','scope':'Finite exact polynomial identities and altered-formula controls only',
            'identity_count':len(checks),'rejected_perturbations':len(controls),'identities':checks,'negative_controls':controls}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path('formal_identity_report.json'))
    main(p.parse_args().output)
