"""Stage 20: outward-rounded binary integer intervals, independent of Stage 19.

Each endpoint is an integer over 2**768. No hardware floats are used for a
proof decision. Matrix solves below use row-oriented Gauss-Jordan elimination;
positivity uses symmetric Schur elimination, not an eigensolver.
"""
from fractions import Fraction

BITS=768
SCALE=1 << BITS

def ceil_div(a,b):
    if b <= 0: raise ValueError('positive denominator required')
    return -((-a)//b)

class I:
    __slots__=('l','h')
    def __init__(self,a=0,b=None):
        if isinstance(a,I):
            if b is not None: raise ValueError('ambiguous interval')
            self.l,self.h=a.l,a.h; return
        a=Fraction(a); b=a if b is None else Fraction(b)
        if a>b: raise ValueError('reversed interval')
        self.l=(a.numerator*SCALE)//a.denominator
        self.h=ceil_div(b.numerator*SCALE,b.denominator)
    @classmethod
    def raw(cls,l,h):
        if l>h: raise ValueError('invalid interval')
        x=object.__new__(cls);x.l=l;x.h=h;return x
    def lo(self): return Fraction(self.l,SCALE)
    def hi(self): return Fraction(self.h,SCALE)
    def norm(self): return Fraction(max(abs(self.l),abs(self.h)),SCALE)
    def contains(self,x):
        x=Fraction(x);return self.lo()<=x<=self.hi()
    def __neg__(self):return I.raw(-self.h,-self.l)
    def __add__(self,o):
        o=I(o);return I.raw(self.l+o.l,self.h+o.h)
    __radd__=__add__
    def __sub__(self,o):return self+-I(o)
    def __rsub__(self,o):return I(o)+-self
    def __mul__(self,o):
        o=I(o);v=[self.l*o.l,self.l*o.h,self.h*o.l,self.h*o.h]
        return I.raw(min(v)//SCALE,ceil_div(max(v),SCALE))
    __rmul__=__mul__
    def inv(self):
        if self.l<=0<=self.h:raise ZeroDivisionError('interval contains zero')
        if self.h<0:return -(-self).inv()
        return I.raw((SCALE*SCALE)//self.h,ceil_div(SCALE*SCALE,self.l))
    def __truediv__(self,o):return self*I(o).inv()
    def __rtruediv__(self,o):return I(o)*self.inv()
    def sq(self):
        mx=max(self.l*self.l,self.h*self.h)
        mn=0 if self.l<=0<=self.h else min(self.l*self.l,self.h*self.h)
        return I.raw(mn//SCALE,ceil_div(mx,SCALE))
    def record(self):return {'lower':str(self.lo()),'upper':str(self.hi())}
    def __repr__(self):return 'I('+str(self.lo())+', '+str(self.hi())+')'

Z=I(0);O=I(1)
def zeros(n,m=None):return [[Z for _ in range(n if m is None else m)] for _ in range(n)]
def eye(n):return [[I(int(i==j)) for j in range(n)] for i in range(n)]
def tr(a):return [list(v) for v in zip(*a)]
def dot(a,b):
    if len(a)!=len(b):raise ValueError('dot dimensions')
    return sum((x*y for x,y in zip(a,b)),Z)
def mm(a,b):
    if not a or not b:return []
    if len(a[0])!=len(b):raise ValueError('matrix dimensions')
    return [[dot(row,col) for col in tr(b)] for row in a]
def mv(a,b):return [dot(r,b) for r in a]
def add(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def sub(a,b):return [[x-y for x,y in zip(r,s)] for r,s in zip(a,b)]
def scale(a,s):return [[x*s for x in r] for r in a]
def congr(a,t):return mm(mm(tr(t),a),t)
def part(a,rows,cols):return [[a[i][j] for j in cols] for i in rows]
def overlap(a,b):return max(a.lo(),b.lo())<=min(a.hi(),b.hi())

def spd(a):
    """Certify symmetric PD by scalar Schur elimination. Return all pivots."""
    n=len(a);w=[[I(v) for v in r] for r in a];p=[]
    if any(len(r)!=n for r in w):raise ValueError('square required')
    for k in range(n):
        d=w[k][k]
        if d.l<=0:raise AssertionError('nonpositive pivot at '+str(k))
        p.append(d)
        for i in range(k+1,n):
            for j in range(i,n):
                v=w[i][j]-w[i][k]*w[k][j]/d
                w[i][j]=v;w[j][i]=v
    return p

def solve(a,b):
    """Independent interval Gauss-Jordan solver. Exact input solution is enclosed."""
    n=len(a);m=len(b[0]);w=[[I(v) for v in r]+[I(v) for v in s] for r,s in zip(a,b)]
    for k in range(n):
        choices=[i for i in range(k,n) if not w[i][k].contains(0)]
        if not choices:raise AssertionError('no certified pivot')
        i=max(choices,key=lambda j:min(abs(w[j][k].lo()),abs(w[j][k].hi())))
        w[i],w[k]=w[k],w[i]
        inv=w[k][k].inv()
        w[k]=[v*inv for v in w[k]]
        # The exact normalized pivot is one; it is an identity of row operations.
        w[k][k]=O
        for i in range(n):
            if i==k:continue
            f=w[i][k]
            w[i]=[v-f*u for v,u in zip(w[i],w[k])];w[i][k]=Z
    return [r[n:n+m] for r in w]

def max_row_norm(a):return max(sum(x.norm() for x in row) for row in a)
