"""Small exact rational matrix routines used only by the new examples.

The general theorem is proved in the paper; finite tests do not prove a
universally quantified theorem. No floating point or third-party import.
"""
from fractions import Fraction as F


def matrix(rows):
    out = [[F(x) for x in row] for row in rows]
    if not out or not out[0] or any(len(row) != len(out[0]) for row in out):
        raise ValueError('nonempty rectangular matrix required')
    return out


def zeros(n, m=None):
    return [[F(0) for _ in range(n if m is None else m)] for _ in range(n)]


def eye(n):
    a = zeros(n)
    for i in range(n): a[i][i] = F(1)
    return a


def transpose(a): return [list(row) for row in zip(*a)]


def add(a, b):
    if len(a) != len(b) or any(len(r) != len(s) for r, s in zip(a, b)):
        raise ValueError('shape mismatch')
    return [[x+y for x, y in zip(r, s)] for r, s in zip(a, b)]


def scale(a, c): return [[F(c)*v for v in row] for row in a]


def mul(a, b):
    if len(a[0]) != len(b): raise ValueError('shape mismatch')
    return [[sum((x*y for x, y in zip(r, s)), F(0)) for s in transpose(b)] for r in a]


def mv(a, v): return [sum((x*y for x, y in zip(r, v)), F(0)) for r in a]


def dot(u, v):
    if len(u) != len(v): raise ValueError('shape mismatch')
    return sum((F(a)*F(b) for a, b in zip(u, v)), F(0))


def outer(u, v=None):
    if v is None: v = u
    return [[F(a)*F(b) for b in v] for a in u]


def trace(a): return sum((a[i][i] for i in range(len(a))), F(0))


def pairing(a, b):
    # The trace pairing for symmetric matrices counts both off-diagonal entries.
    return trace(mul(a, b))


def rref(a):
    a = matrix(a); m, n = len(a), len(a[0]); piv = []; i = 0
    for j in range(n):
        k = next((k for k in range(i, m) if a[k][j]), None)
        if k is None: continue
        a[i], a[k] = a[k], a[i]
        q = a[i][j]; a[i] = [v/q for v in a[i]]
        for k in range(m):
            if k != i and a[k][j]:
                q = a[k][j]; a[k] = [x-q*y for x, y in zip(a[k], a[i])]
        piv.append(j); i += 1
        if i == m: break
    return a, piv


def rank(a): return len(rref(a)[1])


def inverse(a):
    n = len(a)
    if len(a[0]) != n: raise ValueError('square matrix required')
    e = eye(n); rr, pp = rref([r+s for r, s in zip(a, e)])
    if pp[:n] != list(range(n)) or [r[:n] for r in rr] != e:
        raise ValueError('singular matrix')
    inv = [r[n:] for r in rr]
    if mul(a, inv) != e or mul(inv, a) != e: raise AssertionError('inverse identity')
    return inv


def in_range(a, v):
    if len(a) != len(v): raise ValueError('shape mismatch')
    return rank(a) == rank([row+[F(x)] for row, x in zip(a, v)])


def psd(a):
    """Return a rational Schur-pivot certificate, or reject a non-PSD input."""
    a = matrix(a)
    if len(a) != len(a[0]) or a != transpose(a):
        raise ValueError('symmetric square matrix required')
    pivots = []; labels = list(range(len(a))); permutation = []
    while a:
        n = len(a)
        if any(a[i][i] < 0 for i in range(n)):
            raise ValueError('negative Schur diagonal')
        k = next((i for i in range(n) if a[i][i] > 0), None)
        if k is None:
            if any(v for row in a for v in row):
                raise ValueError('zero diagonal with nonzero crossing')
            return {'rank':len(pivots),'pivots':pivots,'pivot_indices':permutation,
                    'zero_dimension':n}
        permutation.append(labels.pop(k)); ids = [i for i in range(n) if i != k]
        p = a[k][k]; pivots.append(p)
        a = [[a[i][j]-a[i][k]*a[k][j]/p for j in ids] for i in ids]
    return {'rank':len(pivots),'pivots':pivots,'pivot_indices':permutation,
            'zero_dimension':0}


def border(a, q, c):
    if len(a) != len(q): raise ValueError('shape mismatch')
    return [row+[F(v)] for row, v in zip(a, q)]+[[F(v) for v in q]+[F(c)]]


# Univariate polynomials stored low coefficient first; used for exact identities.
def padd(a, b):
    c = [F(0)]*max(len(a), len(b))
    for i, x in enumerate(a): c[i] += x
    for i, x in enumerate(b): c[i] += x
    while len(c)>1 and c[-1] == 0: c.pop()
    return c


def pscale(a, k): return [F(k)*F(x) for x in a]


def pmul(a, b):
    c = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): c[i+j] += x*y
    while len(c)>1 and c[-1] == 0: c.pop()
    return c


def peval(a, r):
    out = F(0)
    for v in reversed(a): out = out*r+v
    return out
