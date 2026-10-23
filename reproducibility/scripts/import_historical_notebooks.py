#!/usr/bin/env python3
"""Stage historical paper notebooks from the author's source archive, without modifying it."""
import argparse
import hashlib
import json
import re
from pathlib import Path

HERE=Path(__file__).resolve().parents[2]
DEFAULT=Path.home()/'cerebellum/cerebellum/ref_sims/paperModelingSZ1_reproducibility'
NOTEBOOKS={
    'autor':'META_PAPER_MODELING_AUTOR_PREFIXED.ipynb',
    'basal':'META_PAPER_MODELING_BASAL_PREFIXED.ipynb',
    'evoked':'META_PAPER_MODELING_EVOKED_v5_PREFIXED.ipynb',
}


def digest(b):
    return hashlib.sha256(b).hexdigest()


def source_text(cell):
    val=cell.get('source','')
    return ''.join(val) if isinstance(val,list) else str(val)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,default=DEFAULT)
    ap.add_argument('--check-only',action='store_true')
    ap.add_argument('--overwrite',action='store_true')
    a=ap.parse_args()
    source=a.source.expanduser().resolve()
    if not source.is_dir():
        ap.error('Notebook archive not found: '+str(source))
    dest=HERE/'notebooks/historical'
    report={'source_dir':str(source),'notebooks':{},'note':'Original code cell sources preserved; execution outputs and execution counts cleared for sharing.'}
    prepared={}
    for mode,name in NOTEBOOKS.items():
        item=source/name
        if not item.is_file():
            ap.error('Historical notebook missing: '+str(item))
        raw=item.read_bytes()
        notebook=json.loads(raw.decode('utf-8'))
        if not isinstance(notebook,dict) or notebook.get('nbformat')!=4 or not isinstance(notebook.get('cells'),list):
            ap.error('Invalid Jupyter Notebook v4: '+str(item))
        sources=[]
        num_outputs=0
        paths=[]
        pkl_mentions=0
        for cell in notebook['cells']:
            text=source_text(cell)
            if cell.get('cell_type')=='code':
                sources.append(text)
                pkl_mentions+=len(re.findall(r'\.pkl\b',text,re.I))
                for pattern in [r'/home/[^\s\"\']+',r'/media/[^\s\"\']+',r'[^\s\"\']+\.pkl']:
                    matches=re.findall(pattern,text)
                    for m in matches[:5]:
                        cleaned=re.sub(r'/home/[^/]+','/home/<USER>',m)
                        cleaned=re.sub(r'/media/[^/]+','/media/<USER>',cleaned)
                        paths.append(cleaned[:180])
                num_outputs+=len(cell.get('outputs',[]))
                cell['outputs']=[]
                cell['execution_count']=None
        csha=digest('\n'.join(sources).encode('utf-8'))
        prepared[mode]=(name,notebook,csha)
        report['notebooks'][mode]={
            'filename':name,'source_sha256':digest(raw),'code_cells':len(sources),
            'code_sha256':csha,'outputs_removed':num_outputs,
            'pkl_references':pkl_mentions,
            'path_examples':sorted(set(paths))[:15],
        }
        print(f'CHECK {mode:7s} {name}: {len(sources)} code cells; {num_outputs} outputs to clear; {pkl_mentions} .pkl references')
    if a.check_only:
        print('CHECK-ONLY OK; no notebook files written')
        return
    dest.mkdir(parents=True,exist_ok=True)
    for mode,(name,notebook,csha) in prepared.items():
        target=dest/name
        txt=json.dumps(notebook,ensure_ascii=False,indent=1)+'\n'
        if target.exists() and not a.overwrite:
            current=json.loads(target.read_text(encoding='utf-8'))
            old='\n'.join(source_text(cell) for cell in current.get('cells',[]) if cell.get('cell_type')=='code')
            if digest(old.encode('utf-8'))!=csha:
                ap.error(f'Target has different code; use --overwrite deliberately: {target}')
            print('SKIP existing, same code:',target.name)
        else:
            target.write_text(txt,encoding='utf-8')
            print('COPIED:',target)
    report['source_dir']='Historical local paperModelingSZ1_reproducibility archive (absolute path withheld)'
    (dest/'NOTEBOOK_PROVENANCE.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('Notebook source code retained, execution results cleared; legacy notebooks require historical intermediate files')

if __name__=='__main__':
    main()
