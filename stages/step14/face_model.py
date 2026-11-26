"""Reconstruct the exposed face from original words; proof path has no third-party imports."""
from __future__ import annotations
import json
from fractions import Fraction as F
from functools import lru_cache
from load_upstream import UP13, UP12
import completion_model as cm
from verify_step12 import rank_kernel


def require(value: bool, message: str) -> None:
    if not value:
        raise ValueError(message)


def at_x(form, x):
    return cm.evaluate_form(form, x)[:21]  # constant + twenty coordinates; z20=0


def combine_forms(weights, forms):
    return [sum((c * f[k] for c, f in zip(weights, forms)), F(0)) for k in range(21)]


def affine_value(f, z):
    return f[0] + cm.dot(f[1:], z[:20])


@lru_cache(None)
def reconstruct():
    b, blocks, completion = cm.model()
    base_input = json.loads((UP12 / 'inputs/primal_lower.json').read_text())
    x = F(base_input['x'])
    require(x == F(963899, 2500000), 'Wrong fixed radius')
    zbase = list(map(F, base_input['quotient_z']))
    require(len(zbase) == 21 and zbase[-1] == 0, 'Bad quotient representative')
    A = cm.evaluate_matrix(blocks['M2_None'], x, zbase)
    u = [F(0), F(1), F(-1), F(0), F(0), F(0), F(0)]
    k = cm.solve(A, u)
    theta = 1 / (3 * cm.dot(u, k))
    zstar = list(zbase)
    zstar[0] -= theta
    require(k[-1] != 0, 'Kernel chart is singular')
    nu = [c / k[-1] for c in k]
    frozen = json.loads((UP13 / 'inputs/singular_family.json').read_text())
    require(x == F(frozen['x']) and zstar == list(map(F, frozen['zstar'])),
            'Independent reconstruction differs from the frozen singular point')
    expected = list(map(F, frozen['null_vector']))
    require(nu == [c / expected[-1] for c in expected], 'Kernel reconstruction differs')

    # E_i(z)=(B2 nu)_i and q(z)=nu^T b2. Every entry is rebuilt from word forms.
    rows = [combine_forms(nu, [at_x(f, x) for f in row])
            for row in blocks['M2_None']['mat']]
    q = combine_forms(nu, [at_x(f, x) for f in completion['M2_None']['cross']])
    rho = combine_forms(nu, rows)
    c2 = at_x(completion['M2_None']['diagonal'], x)
    H = [f[1:] for f in rows]
    rank, face_kernel = rank_kernel(H, 20)
    require(rank == 7 and len(face_kernel) == 13, 'Unexpected face affine dimension')
    require(all(affine_value(f, zstar) == 0 for f in rows), 'Point is outside the face')
    require(affine_value(q, zstar) < 0, 'Wrong signed crossing convention')

    # A necessary condition for full completion on this face is q=0.
    Hc = H + [q[1:]]
    rhs = [-f[0] for f in rows] + [-q[0]]
    rankc, kernelc = rank_kernel(Hc, 20)
    require(rankc == 8 and len(kernelc) == 12, 'Wrong compatible affine dimension')
    error = [cm.dot(row, zstar[:20]) - c for row, c in zip(Hc, rhs)]
    multipliers = cm.solve(cm.mm(Hc, cm.transpose(Hc)), error)
    correction = cm.mv(cm.transpose(Hc), multipliers)
    center = [c - a for c, a in zip(zstar[:20], correction)]
    require(cm.mv(Hc, center) == rhs, 'Compatible affine center is wrong')
    require(all(not any(cm.mv(Hc, v)) for v in kernelc), 'Wrong tangent basis')

    return dict(x=x, b=b, blocks=blocks, completion=completion, zbase=zbase,
                zstar=zstar, nu=nu, theta=theta, equations=rows, cross=q, rho=rho, c2=c2,
                face_kernel=face_kernel, compatibility_rows=Hc,
                compatibility_rhs=rhs, center=center, compatible_kernel=kernelc)


def atom_form(data, atom):
    blocks = data['blocks']
    name = atom['block']
    require(name in blocks, 'Unknown positive block')
    v = list(map(F, atom['vector']))
    require(v and any(v), 'Zero atom')
    n = len(blocks[name]['mat'])
    # The M2 compatible face can be represented by its first six rows/columns.
    require(len(v) == n or (name == 'M2_None' and len(v) == n - 1), 'Wrong atom dimension')
    norm = cm.dot(v, v)
    require(norm == F(atom['normalizer']) and norm > 0, 'Wrong normalizer')
    weight = F(atom['weight'])
    require(weight > 0, 'Nonpositive atom weight')
    total = [F(0)] * 21
    for i in range(len(v)):
        for j in range(len(v)):
            f = at_x(blocks[name]['mat'][i][j], data['x'])
            c = weight * v[i] * v[j] / norm
            total = [a + c * t for a, t in zip(total, f)]
    return total
