"""Generic second-order forward automatic differentiation over intervals.

A jet stores value, gradient and upper triangular Hessian. Schur complements
are differentiated by ordinary scalar elimination using these jets, rather
than any prederived Schur derivative or inverse derivative formula.
"""
from decimal_interval import Box, ZERO, ONE
N=20
PAIRS=tuple((i,j) for i in range(N) for j in range(i,N))
INDEX={p:k for k,p in enumerate(PAIRS)}
class Jet:
    __slots__=('value','gradient','hessian')
    def __init__(self,value,gradient=None,hessian=None):
        self.value=Box.co(value)
        self.gradient=[ZERO]*N if gradient is None else gradient
        self.hessian=[ZERO]*len(PAIRS) if hessian is None else hessian
    @staticmethod
    def variable(value,index):
        g=[ZERO]*N;g[index]=ONE;return Jet(value,g)
    @staticmethod
    def co(v):return v if isinstance(v,Jet) else Jet(v)
    def __add__(self,b):
        b=Jet.co(b)
        return Jet(self.value+b.value,[a+c for a,c in zip(self.gradient,b.gradient)],[a+c for a,c in zip(self.hessian,b.hessian)])
    __radd__=__add__
    def __neg__(self):return Jet(-self.value,[-q for q in self.gradient],[-q for q in self.hessian])
    def __sub__(self,b):return self+-Jet.co(b)
    def __rsub__(self,b):return Jet.co(b)+-self
    def __mul__(self,b):
        if not isinstance(b,Jet):
            b=Box.co(b)
            return Jet(self.value*b,[q*b for q in self.gradient],[q*b for q in self.hessian])
        g=[self.gradient[i]*b.value+self.value*b.gradient[i] for i in range(N)]
        h=[self.hessian[k]*b.value+self.value*b.hessian[k]+self.gradient[i]*b.gradient[j]+self.gradient[j]*b.gradient[i] for k,(i,j) in enumerate(PAIRS)]
        return Jet(self.value*b.value,g,h)
    __rmul__=__mul__
    def reciprocal(self):
        rv=self.value.reciprocal();c1=-(rv*rv);c2=2*rv*rv*rv
        return Jet(rv,[c1*q for q in self.gradient],[c1*self.hessian[k]+c2*self.gradient[i]*self.gradient[j] for k,(i,j) in enumerate(PAIRS)])
    def __truediv__(self,b):return self*Jet.co(b).reciprocal()
    def h(self,i,j):return self.hessian[INDEX[min(i,j),max(i,j)]]
