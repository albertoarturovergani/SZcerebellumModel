#!/usr/bin/env python3
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import yaml

BASE=Path('/media/alberto/MyPassportAlbertoColoreROSSO/PCDESKTOPUNIPV')
DEFAULT_SOURCE=Path('/home/alberto/cerebellum/cerebellum/ref_sims/paperModelingSZ1_reproducibility')
DEFAULT_DEST=BASE/'SZ_repro'
DEFAULT_REPO=BASE/'SZcerebellumModel'
CASES={
    'HR_basal':'20251120193505_basal_activity_15_0_glutamate_all',
    'NR_basal':'20251120172014_basal_activity_2_0_glutamate_all',
    'HR_kase':'20251113214858_kase_stim_vitro_15_0_glutamate_all',
    'NR_kase':'20251113161643_kase_stim_vitro_2_0_glutamate_all'
}
CELLS=('GranuleCell','GolgiCell','PurkinjeCell','StellateCell','BasketCell')
SUFFIXES={'.yaml','.yml','.json','.swc'}
OLD_PREFIXES=('/home/alberto/cerebellum/cerebellum/','/home/alberto/cerebellum/')


def digest(file):
    h=hashlib.sha256()
    with file.open('rb') as f:
        for part in iter(lambda:f.read(1024*1024),b''):
            h.update(part)
    return h.hexdigest()


def copy_checked(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists():
        if digest(src)!=digest(dst):
            raise RuntimeError(f'File esistente diverso: {dst}. Non viene sovrascritto.')
        return
    shutil.copy2(src,dst)


def run_subprocess(cmd,cwd,log):
    print(f'\nRUN: {" ".join(map(str,cmd))}',flush=True)
    log.parent.mkdir(parents=True,exist_ok=True)
    with log.open('w',encoding='utf-8') as f:
        p=subprocess.Popen(list(map(str,cmd)),cwd=str(cwd),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
        for line in p.stdout:
            f.write(line)
            f.flush()
            print(line,end='',flush=True)
        code=p.wait()
    if code:
        raise RuntimeError(f'Comando fallito (exit {code}). Log: {log}')


def json_write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf-8')


def step_dirs(run_dir):
    return sorted((p for p in run_dir.iterdir() if p.is_dir() and re.match(r'^\d+_deg_\d+\.\d+$',p.name)),key=lambda x:x.name)


def destination_for_yaml(file):
    name=file.name
    if re.fullmatch(r'morphologies_\d+\.\d+\.ya?ml',name):
        return Path('morphologies')/name
    if name.startswith(('dcn_','dcn_io_','dcn_vitro_nest_','dcn_io_vitro_nest_')):
        return Path('configurations/mouse/dcn-io')/name
    if name.startswith(('basal_vitro_','kase_stim_vitro_','pawan_stim_vitro_','vbt_stim_protocol_vitro_')):
        return Path('configurations/mouse/nest')/name
    if name.startswith('mouse_cerebellar_cortex_') and 'structure_summary' not in name:
        return Path('configurations/mouse')/name
    return None


def rewrite_local_paths(file,root):
    original=file.read_text(encoding='utf-8')
    rewritten=original
    for prefix in OLD_PREFIXES:
        rewritten=rewritten.replace(prefix,str(root)+'/')
    if rewritten!=original:
        file.write_text(rewritten,encoding='utf-8')


def _all_refs(data):
    if isinstance(data,dict):
        for key,val in data.items():
            if key in ('ref','file') and isinstance(val,str):
                yield str(key),val
            yield from _all_refs(val)
    elif isinstance(data,list):
        for item in data:
            yield from _all_refs(item)


def check_references(root):
    missing=[]
    for file in (*root.rglob('*.yaml'),*root.rglob('*.yml')):
        if file.relative_to(root).parts[0]=='cerebellum':
            continue
        data=yaml.safe_load(file.read_text(encoding='utf-8'))
        for field,ref in _all_refs(data):
            name=ref.split('#',1)[0]
            if not name or '://' in name and not name.startswith('file://'):
                continue
            if name.startswith('file://'):
                name=name[7:]
            if field=='file' and not name.lower().endswith(('.swc','.yaml','.yml','.json')):
                continue
            if field=='ref' and not name.lower().endswith(('.yaml','.yml','.json','.swc')):
                continue
            path=Path(name)
            candidate=path if path.is_absolute() else (file.parent/path)
            if not candidate.exists():
                missing.append((str(file.relative_to(root)),f'{field}={ref}'))
    return missing


def preview_seed(root,main,kind,k):
    cfg=yaml.safe_load(main.read_text(encoding='utf-8'))
    simulations=cfg.get('simulations') or {}
    seeds={s:v.get('seed') for s,v in simulations.items() if isinstance(v,dict) and 'seed' in v}
    sources={s:str(main.relative_to(root)) for s in seeds}
    if kind=='basal' and 'basal_activity' not in seeds:
        basal=root/'configurations/mouse/nest'/f'basal_vitro_{k}.yaml'
        if basal.is_file():
            data=yaml.safe_load(basal.read_text(encoding='utf-8'))
            sim=(data.get('simulations') or {}).get('basal_activity') or {}
            if 'seed' in sim:
                seeds['basal_activity']=sim['seed']
                sources['basal_activity']=str(basal.relative_to(root))
    return {'seed_values':seeds,'seed_sources':sources,'storage_root':(cfg.get('storage') or {}).get('root')}


def choose_main(root,k,kind):
    if kind=='kase':
        return root/'configurations/mouse/nest'/f'kase_stim_vitro_{k}.yaml'
    dcn=root/'configurations/mouse/dcn-io'/f'dcn_io_vitro_nest_{k}.yaml'
    basal=root/'configurations/mouse/nest'/f'basal_vitro_{k}.yaml'
    return dcn if dcn.exists() else basal


def prepare_step(source,run_name,step,destination,repo):
    k=step.name.split('_deg_',1)[1]
    kind='kase' if 'kase_stim_vitro' in run_name else 'basal'
    archive=destination/'originals'/run_name/step.name
    root=destination/'runs'/run_name/step.name
    source_hashes={}
    paths=[]
    for file in step.rglob('*'):
        if not file.is_file() or file.suffix.lower() not in SUFFIXES:
            continue
        rel=file.relative_to(step)
        copy_checked(file,archive/rel)
        source_hashes[str(rel)]=digest(file)
        if file.suffix.lower()=='.swc':
            parts=rel.parts
            if 'morphologies' not in parts:
                continue
            pos=parts.index('morphologies')
            sub=Path(*parts[pos+1:])
            if not sub.parts:
                continue
            target=Path('morphologies')/sub
        elif file.suffix.lower() in ('.yaml','.yml'):
            target=destination_for_yaml(file)
        else:
            target=None
        if target is not None:
            active=root/target
            active.parent.mkdir(parents=True,exist_ok=True)
            if file.suffix.lower() in ('.yaml','.yml'):
                expected=file.read_text(encoding='utf-8')
                for prefix in OLD_PREFIXES:
                    expected=expected.replace(prefix,str(root)+'/')
                if active.exists():
                    if active.read_text(encoding='utf-8')!=expected:
                        raise RuntimeError(f'Configurazione già presente ma diversa: {active}')
                else:
                    active.write_text(expected,encoding='utf-8')
            else:
                copy_checked(file,active)
            paths.append(target)
    if not root.exists():
        raise RuntimeError(f'Nessuna configurazione trovata in: {step}')
    pkg=repo/'cerebellum'
    if not pkg.is_dir():
        raise RuntimeError(f'Package cerebellum mancante: {pkg}')
    local_pkg=root/'cerebellum'
    if not local_pkg.exists():
        shutil.copytree(pkg,local_pkg,ignore=shutil.ignore_patterns('__pycache__','*.pyc','.git','*.hdf5','*.nio'))
    required=[root/'configurations/mouse'/f'mouse_cerebellar_cortex_{k}.yaml',root/'morphologies'/f'morphologies_{k}.yaml']
    required.extend(root/'morphologies'/c/f'{c}_{k}_atr.swc' for c in CELLS)
    main=choose_main(root,k,kind)
    required.append(main)
    missing=[str(p.relative_to(root)) for p in required if not p.is_file() or p.stat().st_size==0]
    unresolved=check_references(root)
    meta={'source_run':run_name,'source_step':step.name,'k':k,'kind':kind,'main_configuration':str(main.relative_to(root)),'original_file_hashes':source_hashes,'missing_required':missing,'missing_imports':unresolved,'prepared_at':datetime.now().isoformat()}
    if main.is_file():
        meta.update(preview_seed(root,main,kind,k))
    json_write(root/'repro_manifest.json',meta)
    return meta


def prepare(args):
    if not args.source.is_dir():
        raise RuntimeError(f'Cartella archivio non trovata: {args.source}')
    if not args.repo.is_dir():
        raise RuntimeError(f'Repository sul MyPassport non trovato: {args.repo}')
    args.dest.mkdir(parents=True,exist_ok=True)
    report=[]
    for run in sorted((d for d in args.source.iterdir() if d.is_dir() and d.name[:4].isdigit()),key=lambda p:p.name):
        steps=step_dirs(run)
        if not steps:
            continue
        info=run/'info.yaml'
        if info.is_file():
            copy_checked(info,args.dest/'originals'/run.name/'info.yaml')
        for step in steps:
            meta=prepare_step(args.source,run.name,step,args.dest,args.repo)
            report.append({'run':run.name,'step':step.name,'missing_required':meta['missing_required'],'missing_imports':meta['missing_imports'],'seed':meta.get('seed_values',{})})
            mark='OK' if not meta['missing_required'] and not meta['missing_imports'] else 'CHECK'
            print(f'{mark} {run.name}/{step.name} seed={meta.get("seed_values",{})}',flush=True)
            if mark!='OK':
                print('  missing:',meta['missing_required'],meta['missing_imports'],flush=True)
    json_write(args.dest/'prepare_report.json',report)
    bad=[x for x in report if x['missing_required'] or x['missing_imports']]
    print(f'\nPreparati: {len(report)}, con dipendenze mancanti: {len(bad)}')
    print(f'Configurazioni copiate in: {args.dest}')
    print(f'Verifica: {args.dest}/prepare_report.json')


def run_step(root,mpi,resume):
    meta=json.loads((root/'repro_manifest.json').read_text(encoding='utf-8'))
    if meta['missing_required'] or meta['missing_imports']:
        raise RuntimeError(f'Configurazione incompleta: {root}. Controllare repro_manifest.json')
    status_file=root/'run_status.json'
    if status_file.exists():
        prev=json.loads(status_file.read_text(encoding='utf-8'))
        if prev.get('status')=='completed':
            print(f'SKIP già completata: {root}')
            return
    if not os.getenv('VIRTUAL_ENV'):
        raise RuntimeError('Attivare prima BSB_venv_repro e nest_vars.sh')
    yaml_file=root/meta['main_configuration']
    storage=meta.get('storage_root')
    if not storage:
        raise RuntimeError(f'Manca storage.root in {yaml_file}')
    storage=Path(str(storage))
    if storage.is_absolute():
        raise RuntimeError(f'storage.root non relativo: {storage}. Correggere la copia in SZ_repro, non gli originali.')
    hdf5=root/storage
    if '..' in storage.parts:
        raise RuntimeError(f'storage.root deve restare nella cartella di lavoro: {storage}')
    before=set(root.glob('*.nio'))
    json_write(status_file,{'status':'running','started_at':datetime.now().isoformat(),'config':str(yaml_file.relative_to(root)),'mpi':mpi})
    try:
        if hdf5.exists():
            if not resume:
                raise RuntimeError(f'HDF5 già esistente: {hdf5}. Usa --resume soltanto se vuoi riutilizzare questa ricostruzione.')
            print(f'Riutilizzo connettoma esistente: {hdf5}')
        else:
            run_subprocess(['mpirun','-np',str(mpi),'bsb','compile',str(yaml_file)],root,root/'compile.log')
        if not hdf5.is_file() or hdf5.stat().st_size==0:
            raise RuntimeError(f'Compilazione terminata ma HDF5 mancante/vuoto: {hdf5}')
        run_subprocess(['mpirun','-np',str(mpi),'bsb','reconfigure',str(hdf5),str(yaml_file)],root,root/'reconfigure.log')
        sim='basal_activity' if meta['kind']=='basal' else 'kase_stim_vitro'
        run_subprocess(['mpirun','-np',str(mpi),'bsb','simulate',str(hdf5),sim],root,root/'simulate.log')
        files=sorted(p for p in root.glob('*.nio') if p not in before and p.stat().st_size>0)
        if not files:
            raise RuntimeError('Simulazione finita senza nuovi file .nio non vuoti')
        outputs=[{'file':str(p.relative_to(root)),'bytes':p.stat().st_size,'sha256':digest(p)} for p in files]
        json_write(status_file,{'status':'completed','finished_at':datetime.now().isoformat(),'mpi':mpi,'simulation':sim,'seed_values':meta.get('seed_values',{}),'hdf5':{'file':str(hdf5.relative_to(root)),'bytes':hdf5.stat().st_size,'sha256':digest(hdf5)},'nio':outputs})
        print(f'OK {root}: {len(files)} file NIO',flush=True)
    except Exception as e:
        json_write(status_file,{'status':'failed','finished_at':datetime.now().isoformat(),'error':str(e),'mpi':mpi})
        raise


def run(args):
    if args.case=='all':
        cases=list(CASES)
    else:
        cases=[args.case]
    for case in cases:
        name=CASES[case]
        folder=args.dest/'runs'/name
        if not folder.is_dir():
            raise RuntimeError(f'Configurazioni per {case} non preparate. Eseguire prima il comando prepare.')
        steps=sorted((s for s in folder.iterdir() if s.is_dir() and '_deg_' in s.name),key=lambda p:p.name)
        if args.limit:
            steps=steps[:args.limit]
        print(f'\n=== {case} | {len(steps)} livelli ===',flush=True)
        for i,step in enumerate(steps,1):
            print(f'\n[{i}/{len(steps)}] {step.name}',flush=True)
            run_step(step,args.mpi,args.resume)
    print('\nTutte le esecuzioni richieste sono terminate.')


def main():
    parser=argparse.ArgumentParser(description='Copia e riproduce configurazioni SZ originali senza usare HDF5 archiviati')
    parser.add_argument('--source',type=Path,default=DEFAULT_SOURCE)
    parser.add_argument('--dest',type=Path,default=DEFAULT_DEST)
    parser.add_argument('--repo',type=Path,default=DEFAULT_REPO)
    subs=parser.add_subparsers(dest='action',required=True)
    subs.add_parser('prepare')
    x=subs.add_parser('run')
    x.add_argument('--case',choices=['all',*CASES],default='NR_basal')
    x.add_argument('--limit',type=int,default=0)
    x.add_argument('--mpi',type=int,default=8)
    x.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    try:
        if args.action=='prepare':
            prepare(args)
        else:
            run(args)
    except (RuntimeError,OSError,yaml.YAMLError) as exc:
        print('ERROR:',exc,file=sys.stderr)
        sys.exit(1)


if __name__=='__main__':
    main()
