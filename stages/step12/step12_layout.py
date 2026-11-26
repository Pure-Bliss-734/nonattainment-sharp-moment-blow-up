"""Rebuild the exact projected Gram blocks, using the standard library only."""
import json
from pathlib import Path
from step12_model import *
ROOT=Path(__file__).resolve().parent

@lru_cache(None)
def projected():
 b=build();data=json.loads((ROOT/'inputs/projected_layout.json').read_text());out=[]
 for layout,raw in zip(data['blocks'],b['blocks']):
  check(layout['name']==raw['name'],'Layout ordering mismatch')
  ops=[md.operator_from_json(o) for o in layout['basis']]
  if not ops:continue
  mat=[[b['obs']['NQ']]] if raw['name']=='single_Q' else [[b['old_form'](wa.inner(u,v,raw['ground'])) for v in ops] for u in ops]
  out.append(dict(raw,basis=ops,mat=mat))
 return out

