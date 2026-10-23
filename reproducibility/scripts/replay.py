#!/usr/bin/env python3
"""Portable replay: default smoke test (three k levels) or --full (all 16)."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime,timezone
from pathlib import Path
import yaml

HERE=Path(__file__).resolve().parents[2]
REG=HERE/'configurations/SZ_repro_2025/protocols.yaml'


def hashfile(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for s in iter(lambda:f.read(1024*1024),b''):
            h.update(s)
    return h.hexdigest()


def write_json(p,d):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n',encoding='utf-8')


def command(argv,cwd,log):
    print('EXEC:', ' '.join(map(str,argv)),flush=True)
    with log.open('w',encoding='utf-8') as f:
        p=subprocess.Popen(argv,cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
        for line in p.stdout:
            print(line,end='',flush=True)
            f.write(line)
        code=p.wait()
    if code:
        raise RuntimeError('Comando fallito; vedere '+str(log))


def copy_inputs(case,step,entry,out):
    source=HERE/'configurations/SZ_repro_2025/experiments'/case/step
    morph=HERE/'morphologies/SZ_repro_2025'/step
    if not source.is_dir() or not morph.is_dir():
        raise RuntimeError('Sorgenti mancanti: '+str(source))
    out.mkdir(parents=True,exist_ok=False)
    shutil.copytree(source/'configurations',out/'configurations')
    shutil.copytree(morph,out/'morphologies')
    shutil.copytree(HERE/'cerebellum',out/'cerebellum')
    k=step.split('_deg_')[-1]
    local=out/'morphologies'/f'morphologies_{k}.yaml'
    for file in (out/'configurations').rglob('*.yaml'):
        text=file.read_text(encoding='utf-8')
        data=yaml.safe_load(text)
        def refs(o):
            if isinstance(o,dict):
                for key,val in o.items():
                    if key=='ref' and isinstance(val,str):
                        yield val
                    yield from refs(val)
            elif isinstance(o,list):
                for val in o:
                    yield from refs(val)
        for r in set(refs(data)):
            path,sep,suffix=r.partition('#')
            old=(source/file.relative_to(out).parent/path).resolve() if path else None
            expected=(HERE/'morphologies/SZ_repro_2025'/step/f'morphologies_{k}.yaml').resolve()
            if old==expected:
                new=Path(os.path.relpath(local,file.parent)).as_posix()+(sep+suffix if sep else '')
                text=text.replace(r,new)
        file.write_text(text,encoding='utf-8')
    morph_cfg=yaml.safe_load(local.read_text(encoding='utf-8'))
    for item in morph_cfg.get('morphologies',[]):
        fp=(local.parent/item['file']).resolve()
        if not fp.is_file() or not local.parent.resolve() in fp.parents:
            raise RuntimeError(f'SWC runtime non valido: {item["file"]}')
        item['file']=str(fp)
    local.write_text(yaml.safe_dump(morph_cfg,sort_keys=False,allow_unicode=True),encoding='utf-8')
    main=out/Path(entry['main'])
    if not main.is_file():
        raise RuntimeError('Main config non trovato: '+str(main))
    for f in (out/'configurations').rglob('*.yaml'):
        d=yaml.safe_load(f.read_text(encoding='utf-8'))
        def check(o):
            if isinstance(o,dict):
                for k,v in o.items():
                    if k=='ref' and isinstance(v,str) and v.split('#')[0].endswith('.yaml'):
                        if not (f.parent/v.split('#')[0]).is_file():
                            raise RuntimeError(f'Ref mancante {v} in {f}')
                    check(v)
            elif isinstance(o,list):
                for v in o: check(v)
        check(d)
    for f in (out/'morphologies').rglob('*.yaml'):
        d=yaml.safe_load(f.read_text(encoding='utf-8'))
        for item in d.get('morphologies',[]):
            target=Path(item['file'])
            target=target if target.is_absolute() else f.parent/target
            if not target.is_file():
                raise RuntimeError(f'SWC mancante in {f}: '+item['file'])
    return main


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',default='NR_basal4hz')
    selection=p.add_mutually_exclusive_group()
    selection.add_argument('--full',action='store_true',help='Run all 16 archived k-levels for each selected case')
    selection.add_argument('--levels',type=str,help='Explicit comma-separated k levels, e.g. 1.000,0.760,0.600')
    selection.add_argument('--limit',type=int,default=None,help='Legacy compatibility: first N levels, 0 means all 16')
    p.add_argument('--mpi',type=int,default=8)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--prepare-only',action='store_true',help='Copia e verifica gli input, senza eseguire NEST')
    a=p.parse_args()
    if (a.limit is not None and a.limit<0) or a.mpi<1:
        p.error('limit>=0 e mpi>=1 richiesti')
    default_smoke=('1.000','0.760','0.600')
    if a.levels:
        user_levels=tuple(v.strip() for v in a.levels.split(',') if v.strip())
        if not user_levels or len(set(user_levels))!=len(user_levels):
            p.error('The --levels list must contain distinct k values')
        for value in user_levels:
            try:
                normalized=f'{float(value):.3f}'
            except ValueError:
                p.error(f'Invalid --levels value: {value}')
            if value!=normalized:
                p.error(f'Use 3 decimal places for --levels: {normalized}')
    else:
        user_levels=default_smoke
    reg=yaml.safe_load(REG.read_text(encoding='utf-8'))['experiments']
    targets=list(reg) if a.case=='all' else [a.case]
    if any(c not in reg for c in targets):
        p.error('Caso sconosciuto. Disponibili: '+', '.join(reg))
    target=a.output.expanduser().resolve()
    if target==HERE or HERE in target.parents or target in HERE.parents:
        p.error('Output deve essere esterno al checkout GitHub')
    if not a.prepare_only and (not os.getenv('VIRTUAL_ENV') or shutil.which('mpirun') is None or shutil.which('bsb') is None):
        p.error('Attiva il venv e installa MPI, BSB, NEST 3.7 e cerebmodule')
    for c in targets:
        entries=reg[c]['levels']
        if a.full or a.limit==0:
            selected=entries
        elif a.limit is not None:
            selected=entries[:a.limit]
        else:
            by_k={v['k']:v for v in entries}
            missing=[k for k in user_levels if k not in by_k]
            if missing:
                p.error(f'{c}: requested k not in archive: {missing}')
            selected=[by_k[k] for k in user_levels]
        print(f'CASE {c}: {len(selected)}/{len(entries)} levels: {", ".join(v["k"] for v in selected)}',flush=True)
        for item in selected:
            step=item['step']
            out=target/c/step
            status=out/'run_status.json'
            if status.is_file():
                if json.loads(status.read_text(encoding='utf-8')).get('status')=='completed':
                    print('SKIP completato',c,step)
                    continue
            if out.exists():
                raise RuntimeError(f'Output già esistente e non completato: {out}. Non sovrascrivo.')
            main=copy_inputs(c,step,item,out)
            print('INPUT OK',c,step,main,flush=True)
            if a.prepare_only:
                continue
            cfg=yaml.safe_load(main.read_text(encoding='utf-8'))
            storage=(cfg.get('storage') or {}).get('root')
            if not storage or Path(str(storage)).is_absolute() or '..' in Path(str(storage)).parts:
                raise RuntimeError(f'storage.root non sicuro in {main}')
            hdf5=out/str(storage)
            if hdf5.exists():
                raise RuntimeError(f'HDF5 preesistente: {hdf5}')
            write_json(status,{'status':'running','case':c,'step':step,'main':str(main.relative_to(out))})
            try:
                command(['mpirun','-np',str(a.mpi),'bsb','compile',str(main)],out,out/'compile.log')
                if not hdf5.is_file():
                    raise RuntimeError(f'HDF5 assente: {hdf5}')
                command(['mpirun','-np',str(a.mpi),'bsb','reconfigure',str(hdf5),str(main)],out,out/'reconfigure.log')
                sim='kase_stim_vitro' if reg[c]['protocol']=='kase' else 'basal_activity'
                command(['mpirun','-np',str(a.mpi),'bsb','simulate',str(hdf5),sim],out,out/'simulate.log')
                nio=sorted(p for p in out.glob('*.nio') if p.stat().st_size)
                if not nio:
                    raise RuntimeError('Nessun NIO non vuoto')
                write_json(status,{'status':'completed','case':c,'step':step,'mpi':a.mpi,'seed':item['nest_seed'],
                                   'hdf5':{'path':hdf5.name,'sha256':hashfile(hdf5)},
                                   'nio':[{'path':f.name,'sha256':hashfile(f)} for f in nio],
                                   'timestamp':datetime.now(timezone.utc).isoformat()})
                print('OK',c,step,len(nio),'file NIO',flush=True)
            except Exception as exc:
                write_json(status,{'status':'failed','error':str(exc),'case':c,'step':step})
                raise
    print('OPERAZIONE TERMINATA')

if __name__=='__main__':
    main()
