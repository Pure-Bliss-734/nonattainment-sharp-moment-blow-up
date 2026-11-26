"""Outward-rounded dyadic intervals. All arithmetic is integer/Fraction only.

Endpoints are integers divided by 2**BITS. Rounding always expands an interval.
No hardware floating-point operation is used in a proof decision.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
BITS=720
SCALE=1<<BITS

def ceildiv(a,b):
    if b<=0:raise ValueError('positive denominator required')
    return -((-a)//b)
class I:
    __slots__=('lo','hi')
    def __init__(self,lo=0,hi=None,*,raw=False):
        if raw:
            self.lo=int(lo);self.hi=int(hi)
        else:
            a=F(lo);b=F(lo if hi is None else hi)
            self.lo=(a.numerator*SCALE)//a.denominator
            self.hi=ceildiv(b.numerator*SCALE,b.denominator)
        if self.lo>self.hi:raise ValueError('reversed interval')
    @staticmethod
    def co(v):return v if isinstance(v,I) else I(v)
    def __add__(self,o):
        o=I.co(o);return I(self.lo+o.lo,self.hi+o.hi,raw=True)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo,raw=True)
    def __sub__(self,o):return self+-I.co(o)
    def __rsub__(self,o):return I.co(o)+-self
    def __mul__(self,o):
        o=I.co(o)
        if self.lo==self.hi==0 or o.lo==o.hi==0:return I(0)
        p=(self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi)
        return I(min(p)//SCALE,ceildiv(max(p),SCALE),raw=True)
    __rmul__=__mul__
    def reciprocal(self):
        if self.lo<=0<=self.hi:raise ZeroDivisionError('interval contains zero')
        if self.hi<0:return -(-self).reciprocal()
        return I((SCALE*SCALE)//self.hi,ceildiv(SCALE*SCALE,self.lo),raw=True)
    def __truediv__(self,o):return self*I.co(o).reciprocal()
    def __rtruediv__(self,o):return I.co(o)*self.reciprocal()
    def __pow__(self,n):
        if not isinstance(n,int):raise ValueError('integer power only')
        if n<0:return self.reciprocal()**(-n)
        out=I(1);a=self
        while n:
            if n&1:out=out*a
            n//=2
            if n:a=a*a
        return out
    def iszero(self):return self.lo==self.hi==0
    def excludes_zero(self):return self.lo>0 or self.hi<0
    def absupper(self):return F(max(abs(self.lo),abs(self.hi)),SCALE)
    def lower(self):return F(self.lo,SCALE)
    def upper(self):return F(self.hi,SCALE)
    def contains(self,v):
        q=F(v);return self.lo*q.denominator<=q.numerator*SCALE<=self.hi*q.denominator
    def __repr__(self):return f'I({dec(self.lower(),10)}, {dec(self.upper(),10)})'

def dec(v,digits=35):
    q=F(v)
    with localcontext() as c:
        c.prec=digits
        return str(Decimal(q.numerator)/Decimal(q.denominator))

def record(a):
    a=I.co(a)
    return {'lo':str(a.lower()),'hi':str(a.upper()),'lo_decimal':dec(a.lower(),65),'hi_decimal':dec(a.upper(),65)}

def zmat(n,m=None):return [[I(0) for _ in range(n if m is None else m)] for __ in range(n)]
def eye(n):return [[I(int(i==j)) for j in range(n)] for i in range(n)]
def tr(A):return [list(row) for row in zip(*A)]
def add(A,B):return [[a+b for a,b in zip(r,s)] for r,s in zip(A,B)]
def sub(A,B):return [[a-b for a,b in zip(r,s)] for r,s in zip(A,B)]
def scale(A,c):return [[a*c for a in row] for row in A]
def dot(a,b):
    s=I(0)
    for x,y in zip(a,b):
        x=I.co(x);y=I.co(y)
        if not x.iszero() and not y.iszero():s=s+x*y
    return s

def mm(A,B):
    BT=tr(B)
    return [[dot(r,c) for c in BT] for r in A]
def mv(A,v):return [dot(row,v) for row in A]
def trace(A):return sum((A[i][i] for i in range(len(A))),I(0))
def ip(A,B):return sum((dot(a,b) for a,b in zip(A,B)),I(0))
def inverse(A):
    n=len(A);R=[list(row)+eye(n)[i] for i,row in enumerate(A)]
    for k in range(n):
        # Exact range test. Pivot selection itself has no proof significance.
        choices=[i for i in range(k,n) if R[i][k].excludes_zero()]
        if not choices:raise ValueError('interval pivot does not exclude zero')
        j=max(choices,key=lambda i:min(abs(R[i][k].lo),abs(R[i][k].hi)))
        R[k],R[j]=R[j],R[k];p=R[k][k];R[k]=[v/p for v in R[k]]
        for i in range(n):
            if i==k:continue
            q=R[i][k]
            R[i]=[v-q*w for v,w in zip(R[i],R[k])]
    return [row[n:] for row in R]

def pd_pivots(A):
    """Enclosures of positive LDL pivots, without taking square roots."""
    R=[list(row) for row in A];out=[]
    for k in range(len(R)):
        p=R[k][k]
        if p.lo<=0:raise ValueError(f'Cannot certify positive pivot {k}: {p}')
        out.append(p)
        for i in range(k+1,len(R)):
            for j in range(i,len(R)):
                v=R[i][j]-R[i][k]*R[k][j]/p
                R[i][j]=R[j][i]=v
    return out

def nonsingular_pivots(A):
    R=[list(row) for row in A];out=[]
    for k in range(len(R)):
        cand=[i for i in range(k,len(R)) if R[i][k].excludes_zero()]
        if not cand:raise ValueError('Cannot certify nonsingular minor')
        j=max(cand,key=lambda i:min(abs(R[i][k].lo),abs(R[i][k].hi)))
        R[k],R[j]=R[j],R[k];p=R[k][k];out.append(p)
        for i in range(k+1,len(R)):
            q=R[i][k]/p
            for j in range(k+1,len(R)):R[i][j]=R[i][j]-q*R[k][j]
    return out

def self_test():
    # Contains exact reference values, including negative operands and reciprocals.
    vals=[F(-7,3),F(-1,10**12),F(0),F(1,7),F(19,4)]
    count=0
    for a in vals:
        for b in vals:
            assert (I(a)+I(b)).contains(a+b);assert (I(a)*I(b)).contains(a*b);count+=2
            if b:assert (I(a)/I(b)).contains(a/b);count+=1
    assert (I(-2,3)*I(-4,5)).contains(-12)
    assert (I(-2,3)*I(-4,5)).contains(15)
    try:I(-1,1).reciprocal()
    except ZeroDivisionError:pass
    else:raise AssertionError('zero-containing denominator accepted')
    p=pd_pivots([[I(2),I(1)],[I(1),I(1)]])
    assert p[0].contains(2) and p[1].contains(F(1,2))
    try:pd_pivots([[I(0),I(1)],[I(1),I(3)]])
    except ValueError:pass
    else:raise AssertionError('indefinite matrix accepted')
    return count+5
