"""Independent outward-rounded fixed-decimal intervals.

Arithmetic decisions are integer/Fraction only, not Decimal or binary float.
This does not import or adapt the previous dyadic interval implementation.
Decimal is used solely to print rational endpoints.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
DIGITS=180
DEN=10**DIGITS
class Box:
    __slots__=('low','high')
    def __init__(self,low=0,high=None):
        a=F(low); b=F(low if high is None else high)
        if a>b:raise ValueError('reversed interval')
        self.low=a.numerator*DEN//a.denominator
        self.high=-((-b.numerator*DEN)//b.denominator)
    @classmethod
    def integers(cls,l,h):
        if l>h:raise ValueError('reversed integer interval')
        o=object.__new__(cls);o.low=l;o.high=h;return o
    def zero(self):return self.low==0 and self.high==0
    @staticmethod
    def co(v):return v if isinstance(v,Box) else Box(v)
    def __add__(self,b):
        b=Box.co(b);return Box.integers(self.low+b.low,self.high+b.high)
    __radd__=__add__
    def __neg__(self):return Box.integers(-self.high,-self.low)
    def __sub__(self,b):return self+-Box.co(b)
    def __rsub__(self,b):return Box.co(b)+-self
    def __mul__(self,b):
        b=Box.co(b)
        if self.zero() or b.zero():return ZERO
        endpoints=[a*c for a in (self.low,self.high) for c in (b.low,b.high)]
        return Box.integers(min(endpoints)//DEN,-((-max(endpoints))//DEN))
    __rmul__=__mul__
    def reciprocal(self):
        if self.low<=0<=self.high:raise ZeroDivisionError('interval divisor contains zero')
        # 1/x is decreasing on each connected component of its domain.
        a=F(DEN,self.high); b=F(DEN,self.low)
        return Box(a,b)
    def __truediv__(self,b):return self*Box.co(b).reciprocal()
    def __rtruediv__(self,b):return Box.co(b)*self.reciprocal()
    def __pow__(self,n):
        if not isinstance(n,int):raise TypeError('integer powers required')
        if n<0:return self.reciprocal()**(-n)
        result=ONE
        for _ in range(n):result=result*self
        return result
    def lo(self):return F(self.low,DEN)
    def hi(self):return F(self.high,DEN)
    def norm(self):return F(max(abs(self.low),abs(self.high)),DEN)
    def contains(self,q):return self.lo()<=F(q)<=self.hi()
    def excludes_zero(self):return self.low>0 or self.high<0
    def record(self):return {'lower':str(self.lo()),'upper':str(self.hi()),'lower_decimal':display(self.lo()),'upper_decimal':display(self.hi())}
    def __repr__(self):return f'Box({display(self.lo(),8)},{display(self.hi(),8)})'
ZERO=Box(0);ONE=Box(1)
def display(q,n=45):
    q=F(q)
    with localcontext() as c:
        c.prec=n;return str(Decimal(q.numerator)/Decimal(q.denominator))
def zeros(n,m=None):return [[ZERO for _ in range(n if m is None else m)] for _ in range(n)]
def eye(n):return [[ONE if i==j else ZERO for j in range(n)] for i in range(n)]
def transpose(a):return [list(r) for r in zip(*a)]
def dot(a,b):
    if len(a)!=len(b):raise ValueError('dot dimension mismatch')
    s=ZERO
    for u,v in zip(a,b):
        u=Box.co(u);v=Box.co(v)
        if not u.zero() and not v.zero():s=s+u*v
    return s
def mul(a,b):
    if not a or not b:raise ValueError('empty matrix product')
    if len(a[0])!=len(b):raise ValueError('matrix product dimension mismatch')
    return [[dot(r,c) for c in transpose(b)] for r in a]
def add(a,b):return [[u+v for u,v in zip(r,s)] for r,s in zip(a,b)]
def sub(a,b):return [[u-v for u,v in zip(r,s)] for r,s in zip(a,b)]
def inner(a,b):return sum((dot(r,s) for r,s in zip(a,b)),ZERO)
def ldlt(a):
    """Nonpivoted LDL^T with interval pivots. Positive pivots prove PD."""
    n=len(a);L=eye(n);D=[]
    if any(len(r)!=n for r in a):raise ValueError('not square')
    for j in range(n):
        p=a[j][j]-sum((L[j][k]*L[j][k]*D[k] for k in range(j)),ZERO)
        if p.low<=0:raise ValueError(f'nonpositive/uncertain LDL pivot {j}: {p}')
        D.append(p)
        for i in range(j+1,n):
            L[i][j]=(a[i][j]-sum((L[i][k]*D[k]*L[j][k] for k in range(j)),ZERO))/p
    return L,D
def solve_spd(a,b):
    """Solve A X=B through triangular LDL solves, not a matrix inverse."""
    L,D=ldlt(a);n=len(a);m=len(b[0]);u=zeros(n,m);w=zeros(n,m);out=zeros(n,m)
    for i in range(n):
        for j in range(m):u[i][j]=b[i][j]-sum((L[i][k]*u[k][j] for k in range(i)),ZERO)
    for i in range(n):
        for j in range(m):w[i][j]=u[i][j]/D[i]
    for i in reversed(range(n)):
        for j in range(m):out[i][j]=w[i][j]-sum((L[k][i]*out[k][j] for k in range(i+1,n)),ZERO)
    return out,D

def matrix_minor_pivots(a):
    """Row-pivoted interval elimination with first usable, not largest, pivot."""
    a=[r[:] for r in a];n=len(a);piv=[]
    if any(len(r)!=n for r in a):raise ValueError('minor not square')
    for j in range(n):
        k=next((i for i in range(j,n) if a[i][j].excludes_zero()),None)
        if k is None:raise ValueError('minor pivot includes zero')
        a[j],a[k]=a[k],a[j];p=a[j][j];piv.append(p)
        for i in range(j+1,n):
            fac=a[i][j]/p
            for k in range(j+1,n):a[i][k]=a[i][k]-fac*a[j][k]
            a[i][j]=ZERO
    return piv
