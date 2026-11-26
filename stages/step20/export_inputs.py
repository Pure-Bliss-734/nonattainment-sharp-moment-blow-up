"""Export the rational pencil for the independent prefactor checks."""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import sys,json
ROOT=Path(__file__).resolve().parent
SHA19='8701173566eb0730e2a2a740af2f2dad3770a6e9a33a899279dbe0e74bc0c819'
NAME19='matrix_bootstrap_step19_optimal_prefactor'
def prepare():
    p19 = ROOT.parent/'step19'
    sys.path.insert(0,str(p19))
    import frozen_source
    p18,p17,p16,_=frozen_source.prepare()
    return p19,p18,p17,p16

def ep(p):return [[list(k),str(v)] for k,v in sorted(p.items())]
def export():
    p19,p18,p17,p16=prepare()
    from independent_model import Model,from_form,ROOT as SRC
    model=Model()
    inp={'source19_sha256':SHA19,'scope':'rational matrices shared with audited word quotient',
         'small':[{'name':b['name'],'parent':b['parent'],'T':[[str(v) for v in r] for r in b['T']],
                   'mat':[[ep(p) for p in r] for r in b['mat']]} for b in model.small],
         'raw':[{'name':b['name'],'mat':[[ep(p) for p in r] for r in b['mat']]} for b in model.raw],
         'to_raw':[ep({m:c for m,c in model.to_symmetric[k].items() if 20 not in m}) for k in range(21)],
         'completion':{k:{'cross':[ep(from_form(f)) for f in v['cross']], 'diagonal':ep(from_form(v['diagonal']))} for k,v in model.data['completion'].items()},
         'observables':{k:ep(from_form(v)) for k,v in model.data['observables'].items()},
         'sides':{s:json.loads((SRC/f'inputs/kkt_{s}.json').read_text()) for s in ['lower','upper']},
         'neighborhood_parameters':json.loads((p17/'inputs/neighborhood_parameters.json').read_text())}
    dst=ROOT/'inputs/shared_pencils.json';text=json.dumps(inp,indent=2)+'\n'
    if dst.exists():assert dst.read_text()==text,'export not deterministic'
    dst.write_text(text)
    print('Exported shared rational pencils, 18 small / 11 raw blocks.',flush=True)
    print('input_sha256='+sha256(text.encode()).hexdigest(),flush=True)
    return p19
if __name__=='__main__':export()
