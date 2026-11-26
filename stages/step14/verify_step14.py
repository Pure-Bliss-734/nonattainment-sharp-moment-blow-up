#!/usr/bin/env python3
"""Exact whole-face obstruction certificate. Uses Python standard library only.
No numerical KKT output is an input to this proof.
"""
from __future__ import annotations
import argparse, copy, json, sys, subprocess
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
from load_upstream import UP13
import face_model as fm
from face_model import cm, require
from verify_step12 import verify_state_equations
ROOT = Path(__file__).resolve().parent
sys.set_int_max_str_digits(0)


def decimal(q, prec=55):
    with localcontext() as ctx:
        ctx.prec = prec
        return str(Decimal(q.numerator) / Decimal(q.denominator))


def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')


def verify_identity(data, cert):
    require(F(cert['x']) == data['x'], 'Certificate radius mismatch')
    require(len(cert['atoms']) == 13, 'Wrong certificate support size')
    forms = [fm.atom_form(data, atom) for atom in cert['atoms']]
    total = [sum((f[k] for f in forms), F(0)) for k in range(21)]
    require(sum((F(a['weight']) for a in cert['atoms']), F(0)) == 1, 'Weights are not normalized')
    for v in data['compatible_kernel']:
        require(cm.dot(total[1:], v) == 0, 'Farkas stationarity fails')
    value = fm.affine_value(total, data['center'])
    require(value == F(cert['value']) and value < -F(3, 100000), 'Not a strict exact separator')

    # Recover a global affine identity, without assuming E=0 or q=0.
    generators = data['equations'] + [data['cross']]
    lam = cm.solve(cm.transpose([r[1:] for r in generators]), total[1:])
    constant = total[0] - cm.dot(lam, [r[0] for r in generators])
    expected = [constant] + [F(0)] * 20
    for c, row in zip(lam, generators):
        expected = [a + c * b for a, b in zip(expected, row)]
    require(total == expected and constant == value, 'Global face identity mismatch')
    beta = -lam[-1]
    require(beta > 0, 'Wrong signed crossing multiplier')
    delta = -constant
    gap = delta / beta
    require(gap > F(1, 5000000), 'The stated simple cross-gap is not certified')
    require(gap * gap / 4 > F(1, 10**14), 'Wrong neighborhood product bound')
    return dict(total=total, value=value, alpha=lam[:-1], delta=delta, beta=beta, gap=gap)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'results/verified')
    args = parser.parse_args()
    data = fm.reconstruct()
    ranks = []
    strict_count = 0
    for name, block in data['blocks'].items():
        S = cm.evaluate_matrix(block, data['x'], data['zstar'])
        pivots = cm.exact_psd_rank(S)
        require(len(pivots) == len(S) - (name == 'M2_None'), 'Unexpected relative-interior rank')
        strict_count += len(pivots)
        ranks.append(dict(name=name, size=len(S), rank=len(pivots), positive_pivots=list(map(str, pivots))))
        base = cm.evaluate_matrix(block, data['x'], data['zbase'])
        require(len(cm.exact_psd_rank(base)) == len(base), 'Base is not strictly projected feasible')
    B2 = cm.evaluate_matrix(data['blocks']['M2_None'], data['x'], data['zstar'])
    require(not any(cm.mv(B2, data['nu'])), 'The claimed kernel does not vanish')
    require(len(cm.exact_psd_rank([r[:6] for r in B2[:6]])) == 6, 'Complementary principal block is not positive')

    # Return to original matrix words for the exposing polynomial and crossing.
    basis = data['blocks']['M2_None']['basis']
    operator = cm.operator_linear(data['nu'], basis)
    W = data['completion']['M2_None']['W']
    rho_from_words = fm.at_x(data['b']['old_form'](cm.wa.inner(operator, operator)), data['x'])
    cross_from_words = fm.at_x(data['b']['old_form'](cm.wa.inner(operator, W)), data['x'])
    require(rho_from_words == data['rho'], 'Exposing operator word identity fails')
    require(cross_from_words == data['cross'], 'Crossing word identity fails')
    gvals = [cm.evaluate_form(v, data['x'], data['zstar']) for v in data['b']['gcoords']]
    state = {}
    for word in data['b']['g']['forms']:
        f = cm.md.affine({word: F(1)}, data['b']['g'])
        state[word] = f[0] + cm.dot(f[1:], gvals)
    count = verify_state_equations(state)
    require(count == 5417, 'Unexpected count of original factorized equations')

    cert = json.loads((ROOT / 'inputs/whole_face_atomic_certificate.json').read_text())
    proof = verify_identity(data, cert)
    require(data['c2'] == [F(0)] * 20 + [F(1)], 'The omitted diagonal is not z19')
    original_atoms = []
    for atom in cert['atoms']:
        block = data['blocks'][atom['block']]
        op = cm.operator_linear(list(map(F, atom['vector'])), block['basis'][:len(atom['vector'])])
        global_form = cm.md.affine(cm.wa.inner(op, op, block['ground']), data['b']['g'])
        lifted_form = cm.linear([(global_form[0], data['b']['obs']['one'])] + list(zip(global_form[1:], data['b']['gcoords'])), 22)
        raw = fm.at_x(lifted_form, data['x'])
        coefficient = F(atom['weight']) / F(atom['normalizer'])
        require([coefficient * c for c in raw] == fm.atom_form(data, atom), 'Atom word reconstruction mismatch')
        original_atoms.append(dict(atom, ground=block['ground'], operator=cm.operator_json(op)))
    negative = {}
    variants = {}
    bad = copy.deepcopy(cert); bad['atoms'][0]['weight'] = '-' + bad['atoms'][0]['weight']; variants['negative_weight'] = bad
    bad = copy.deepcopy(cert); bad['value'] = str(-F(bad['value'])); variants['wrong_separator_sign'] = bad
    bad = copy.deepcopy(cert); bad['atoms'][0]['normalizer'] = str(F(bad['atoms'][0]['normalizer']) + 1); variants['wrong_normalizer'] = bad
    bad = copy.deepcopy(cert); bad['atoms'][0]['weight'] = str(F(bad['atoms'][0]['weight']) + F(1, 10**10)); variants['perturbed_weight'] = bad
    for name, bad in variants.items():
        try:
            verify_identity(data, bad)
        except (ValueError, AssertionError):
            negative[name] = True
        else:
            raise ValueError('Negative control incorrectly passed: ' + name)

    geometry = dict(x=str(data['x']), face_affine_dimension=13, face_equation_rank=7,
                    compatible_affine_dimension=12, compatible_equation_rank=8,
                    normalized_kernel=list(map(str, data['nu'])), theta=str(data['theta']),
                    face_equations=list(map(lambda r:list(map(str,r)), data['equations'])),
                    face_tangent_basis=list(map(lambda r:list(map(str,r)), data['face_kernel'])),
                    strict_relative_interior_point=list(map(str,data['zstar'])),
                    crossing_at_point=str(fm.affine_value(data['cross'],data['zstar'])),
                    original_operator=cm.operator_json(operator), added_coordinate_word=cm.operator_json(W),
                    ranks=ranks, original_equations_checked=count,
                    proof='B2 nu=0 has rank 7; at zstar B2 has rank 6 and every other block is positive definite. Hence the exposed face has dimension 13.')
    identity = dict(constant=str(proof['value']), delta=str(proof['delta']), beta=str(proof['beta']),
                    alpha=list(map(str,proof['alpha'])), total_affine_form=list(map(str,proof['total'])),
                    cross_form=list(map(str,data['cross'])), rho_form=list(map(str,data['rho'])),
                    cross_gap=str(proof['gap']), cross_gap_decimal=decimal(proof['gap']),
                    constant_decimal=decimal(proof['value']), beta_decimal=decimal(proof['beta']),
                    term_count=len(cert['atoms']), exact_identity='sum positive atoms = -delta + alpha.T*(B2 nu) - beta*(nu.T b2)',
                    whole_face_conclusion='nu.T b2 <= -delta/beta < -1/5000000 on the entire exposed face')
    neighborhood = dict(condition='abs(alpha.T*(B2 nu)) <= delta/2',
                        exact_product_bound=str(proof['gap']**2/4), simplified_product_bound='1/100000000000000',
                        inequality='rho*(t+c2) > 1/10^14, where rho=nu.T*B2*nu and c2=z19',
                        conclusion='Every full feasible sequence with projected limit in this face has t -> +infinity.',
                        does_not_assert_a_rate_in_radius_distance=True)
    summary = {'all_points_of_exposed_face_lack_finite_completion': True, 'cross_gap_decimal': decimal(proof['gap']), 'every_point_of_face_is_in_closure_of_completable_projection': True, 'exactly_checked_stationarity_coordinates': 12, 'face_affine_equation_rank': 7, 'face_dimension': 13, 'fixed_x': str(data['x']), 'necessary_compatibility_affine_dimension': 12, 'original_factorized_equations': count, 'positive_atoms': 13, 'positive_pivots_at_relative_interior': strict_count, 'whole_face_uniform_neighborhood_divergence': True}
    for name, obj in [('face_geometry.json',geometry), ('whole_face_identity.json',identity),
                      ('neighborhood_bound.json',neighborhood), ('negative_controls.json',negative), ('original_word_atoms.json',original_atoms), ('summary.json',summary)]:
        save(args.output_dir / name, obj)
    print('EXACT: a 13-dimensional exposed face has no finite completion anywhere.')
    print('EXACT: 13 positive atoms; 12 stationarity coordinates and 21 affine coefficients checked.')
    print('SEPARATION:', decimal(proof['value']))
    print('UNIFORM CROSS GAP:', decimal(proof['gap']), '> 1/5000000')
    print('EXACT: rho*(t+c2) > 1/10^14 in the stated face neighborhood.')
    print('Radius-extremum attainment remains unproved. KKT files are discovery only.')
    print('All Step-14 exact checks passed.')

if __name__ == '__main__':
    main()
