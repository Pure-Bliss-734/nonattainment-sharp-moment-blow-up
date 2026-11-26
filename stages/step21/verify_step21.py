#!/usr/bin/env python3
"""Check exact Schur examples and the finite hypotheses of the completion theorem."""
from pathlib import Path
import argparse,json,subprocess,sys
from source_tools import ROOT,prepare,SOURCE_SHA
from verify_examples import run as run_examples,enc
from quadratic_lift_tests import run as run_quadratic


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output-dir',default='results/verified')
    p.add_argument('--examples-only',action='store_true',help='do not read or audit the prior matrix model')
    args=p.parse_args()
    out=Path(args.output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    # The top-level runner checks the preceding stage separately.
    ex=run_examples(out);quad=run_quadratic()
    (out/'quadratic_lift_checks.json').write_text(json.dumps(enc(quad),indent=2,ensure_ascii=False)+'\n')
    bridge=None
    if not args.examples_only:
        from condition_bridge import audit_bridge
        bridge=audit_bridge(out)
    summary={
        'counterexample_and_positive_scenarios': ex['scenario_count'],
        'scalar_exact_sample_points': ex['scalar_example_points'],
        'matrix_identity_checks': ex['identity_checks'],
        'second_Schur_derivative_reference_cases': quad['second_Schur_derivative_reference_cases'],
        'second_Schur_derivative_matrix_entries': quad['matrix_entries_compared'],
        'negative_controls_rejected': ex['negative_controls'],
        'source_model_bridge_checked': bridge is not None,
        'source_model_ranks': None if bridge is None else {s:bridge['sides'][s]['certified_J_rank'] for s in ['lower','upper']},
        'curved_quadratic_path_continuous_interval_proved': '0<r<=1/100, delta=r+2*r^2, t=5/delta+10',
    }
    (out/'summary.json').write_text(json.dumps(enc(summary),indent=2,ensure_ascii=False)+'\n')
    print('All Stage-21 executable checks passed.')
    print('See the paper for the general proof.')
    print('No explicit M or full feasibility neighborhood is certified for the new physical-model quadratic path.')


if __name__=='__main__':main()
