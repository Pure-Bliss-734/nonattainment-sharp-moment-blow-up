#!/usr/bin/env python3
"""Check source-level references, labels, and integration invariants."""
from __future__ import annotations
import argparse,json,re
from collections import Counter
from pathlib import Path

def main(source:Path,output:Path)->None:
    text=source.read_text(encoding='utf-8')
    text=re.sub(r'(?<!\\)%[^\n]*','',text)
    labels=re.findall(r'\\label\{([^}]+)\}',text)
    refs=re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',text)
    bib=re.findall(r'\\bibitem\{([^}]+)\}',text)
    cites=[]
    for m in re.finditer(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}',text):cites.extend(k.strip() for k in m[1].split(','))
    duplicates=lambda values:sorted(k for k,n in Counter(values).items() if n>1)
    errors={
      'duplicate_labels':duplicates(labels),
      'missing_reference_labels':sorted(set(refs)-set(labels)),
      'duplicate_bibliography_keys':duplicates(bib),
      'missing_citation_keys':sorted(set(cites)-set(bib)),
      'legacy_fallback_commands':re.findall(r'\\(?:PriorEq|IntroEq|PriorRef)\b',text),
      'figure_environments':re.findall(r'\\begin\{figure\*?\}',text),
      'temporary_signed_urls':re.findall(r'\?token=|github_pat_|ghp_[A-Za-z0-9]{20,}',text),
    }
    errors={k:v for k,v in errors.items() if v}
    report={'status':'passed' if not errors else 'failed','labels':len(labels),'reference_uses':len(refs),
            'bibliography_entries':len(bib),'citation_uses':len(cites),
            'main_theorems':text.count(r'\begin{maintheorem}'),
            'general_theorems':text.count(r'\begin{generaltheorem}'),
            'propositions':text.count(r'\begin{setupprop}'),
            'corollaries':text.count(r'\begin{paperendpointcor}'),
            'proof_environments':text.count(r'\begin{proof}'),'errors':errors,
            'scope':'Syntactic consistency; does not establish mathematical correctness.'}
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if errors:raise RuntimeError('manuscript consistency check failed')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path,help='Path to the manuscript LaTeX source')
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.source,a.output)
