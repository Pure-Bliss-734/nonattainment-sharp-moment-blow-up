#!/usr/bin/env python3
"""Rebuild frozen Step-15 inputs and export coefficients, not verification answers.

This is the only bridge importing the earlier implementation. The independent
verifier reads its JSON and never imports old symmetry, interval or KKT code.
"""
from pathlib import Path
from fractions import Fraction
import sys, json, hashlib, shutil, argparse
ROOT=Path(__file__).resolve().parent
SHA='d7191cf97e1860f35364164b230f3bf51cc58f95aae70df20de28d2e07e0fbc6'
def plain(v):
    if isinstance(v,Fraction): return str(v)
    if isinstance(v,dict): return {str(k):plain(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [plain(x) for x in v]
    return v
def main():
    old=ROOT.parent/'step15'
    sys.path.insert(0,str(old))
    import load_frozen
    from step12_model import build
    from step12_layout import projected
    from completion_model import model
    from symmetric_model import encode_blocks
    b=build();_,_,comp=model()
    data={'source_sha256':SHA,'raw_blocks':projected(),'observables':b['obs'],'completion':comp,
          'legacy_small_blocks_for_comparison':encode_blocks()}
    path=ROOT/'inputs';path.mkdir(exist_ok=True)
    (path/'pencils.json').write_text(json.dumps(plain(data),indent=2)+'\n')
    for name in ['symmetry.json','kkt_lower.json','kkt_upper.json','face_rows_lower.json','face_rows_upper.json']:
        shutil.copyfile(old/'inputs'/name,path/name)
    for side in ['lower','upper']:
        shutil.copyfile(load_frozen.UP12/f'inputs/primal_{side}.json',path/f'interior_{side}.json')
    print('Exported original rational pencils and frozen candidate data; no old certificate booleans used.')
if __name__=='__main__':main()
