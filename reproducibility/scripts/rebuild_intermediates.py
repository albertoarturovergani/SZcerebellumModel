#!/usr/bin/env python3
"""Stage the original 2025 analysis layout and optionally recreate its intermediate PKLs."""
import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
import yaml

ROOT=Path(__file__).resolve().parents[2]
REG=ROOT/'configurations/SZ_repro_2025/protocols.yaml'
NOTEBOOKS=('META_PAPER_MODELING_AUTOR_PREFIXED.ipynb','META_PAPER_MODELING_BASAL_PREFIXED.ipynb','META_PAPER_MODELING_EVOKED_v5_PREFIXED.ipynb')


def json_write(path,data):
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf-8')


def link_readonly(src,dst):
    src=src.resolve()
    if not src.is_file():
        raise FileNotFoundError(src)
    if dst.is_symlink():
        if dst.resolve()!=src:
            raise RuntimeError(f'Collegamento esistente diverso: {dst}')
        return
    if dst.exists():
        if os.path.samefile(src,dst):
            return
        raise RuntimeError(f'File esistente diverso: {dst}')
    try:
        dst.symlink_to(src)
    except OSError:
        try:
            os.link(src,dst)
        except OSError as exc:
            raise RuntimeError(f'Impossibile collegare senza duplicare {src}. Usa un output Linux nativo.') from exc


def copy_checked(src,dst):
    if not src.is_file():
        raise FileNotFoundError(src)
    if dst.is_file():
        if src.read_bytes()!=dst.read_bytes():
            raise RuntimeError(f'Input YAML diverso già presente: {dst}')
        return
    dst.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(src,dst)


def completed_step(src):
    status=src/'run_status.json'
    if not status.is_file():
        return None
    entry=json.loads(status.read_text(encoding='utf-8'))
    if entry.get('status')!='completed':
        return None
    hname=(entry.get('hdf5') or {}).get('path')
    if not isinstance(hname,str) or Path(hname).name!=hname:
        raise RuntimeError(f'Status HDF5 inatteso: {status}')
    hdf=src/hname
    nio=sorted(src.glob('*.nio'))
    if not hdf.is_file() or len(nio)!=8 or any(x.stat().st_size==0 for x in nio):
        raise RuntimeError(f'Risultati incompleti: {src}')
    return hdf,nio


def stage_case(case,reg,input_root,workspace):
    meta=reg['experiments'][case]
    legacy=workspace/meta['archive_2025']
    legacy.mkdir(parents=True,exist_ok=True)
    info={
        'experiment':{'name':meta['archive_2025'],'n_steps':0,'stim_name':meta['protocol']},
        'modeling_parameters':{'reserve':{'k0':0.5,'sharpness':15 if meta['resilience']=='HR' else 2}}
    }
    available=[]
    for level in meta['levels']:
        step=level['step']
        src=input_root/case/step
        result=completed_step(src)
        if result is None:
            continue
        hdf,nio=result
        dest=legacy/step
        struct=dest/'1.structure'
        stim=dest/'2.stimulus'
        function=dest/'3.function'
        for folder in (struct,stim,function):
            folder.mkdir(parents=True,exist_ok=True)
        link_readonly(hdf,struct/hdf.name)
        for f in nio:
            link_readonly(f,function/f.name)
        cfg=src/'configurations/mouse'
        for f in cfg.rglob('*.yaml'):
            copy_checked(f,stim/f.name)
        morph=src/'morphologies'/f'morphologies_{level["k"]}.yaml'
        if morph.is_file():
            copy_checked(morph,struct/morph.name)
        available.append(step)
        print('STAGED',case,step,'NIO=',len(nio),flush=True)
    if available:
        info['experiment']['n_steps']=len(available)
        path=legacy/'info.yaml'
        if path.exists():
            old=yaml.safe_load(path.read_text(encoding='utf-8'))
            if old.get('modeling_parameters')!=info['modeling_parameters']:
                raise RuntimeError('Reserve parameters are inconsistent in '+str(path))
        path.write_text(yaml.safe_dump(info,sort_keys=False),encoding='utf-8')
    return legacy,available


def execute(cmd,cwd,log):
    log.parent.mkdir(parents=True,exist_ok=True)
    print('EXEC',' '.join(str(x) for x in cmd),'LOG',log,flush=True)
    with log.open('w',encoding='utf-8') as stream:
        code=subprocess.run([str(x) for x in cmd],cwd=cwd,stdout=stream,stderr=subprocess.STDOUT).returncode
    if code:
        tail=log.read_text(encoding='utf-8',errors='replace').splitlines()[-25:]
        raise RuntimeError(f'Comando fallito (exit {code}): {log}\n'+'\n'.join(tail))


def rebuild_case(case,reg,input_root,legacy,steps,sigma_baseline,sigma_kase):
    scripts=ROOT/'scripts_py'
    for step in steps:
        struct=legacy/step/'1.structure'
        summaries=list(struct.glob('*_structure_summary.yaml'))
        if not summaries:
            files=list(struct.glob('*.hdf5'))
            if len(files)!=1:
                raise RuntimeError(f'HDF5 ambiguo in {struct}: {len(files)} file')
            source=input_root/case/step
            original_hdf5=source/files[0].name
            if not original_hdf5.is_file():
                raise FileNotFoundError(f'HDF5 originale del replay assente: {original_hdf5}')
            execute([sys.executable,scripts/'do_structure_yaml.py','--path',original_hdf5,'--output-dir',struct],source,legacy/f'structure_{step}.log')
            if not list(struct.glob('*_structure_summary.yaml')):
                raise RuntimeError('Il riepilogo di struttura non è stato generato in '+str(struct))
    execute([sys.executable,scripts/'do_network_analysis.py','--path',legacy],legacy,legacy/'network_analysis.log')
    spkl=legacy/'B_collective_results/1.structures/df_summary.pkl'
    if not spkl.is_file():
        raise RuntimeError('PKL strutturale assente: '+str(spkl))
    if reg['experiments'][case]['protocol']=='kase':
        for sigma in sigma_kase:
            execute([sys.executable,scripts/'do_completeAnalysis_vbtkase_v2.py','--path',legacy,'--sigma',str(sigma),'--stimName','kase_stim_vitro'],legacy,legacy/f'kase_sigma_{sigma}.log')
            out=legacy/f'B_collective_results/2.functions/sigma_{sigma}/{sigma}_df.pkl'
            if not out.is_file():
                raise RuntimeError('PKL Kase assente: '+str(out))
    else:
        execute([sys.executable,scripts/'do_completeAnalysis_v2.py','--path',legacy,'--sigma',str(sigma_baseline)],legacy,legacy/'basal_analysis.log')
        out=legacy/'B_collective_results/2.functions/df_summary.pkl'
        if not out.is_file():
            raise RuntimeError('PKL funzionale assente: '+str(out))
    json_write(legacy/'postprocess_status.json',{'case':case,'steps':steps,'sigma_baseline':sigma_baseline,'sigma_kase':sigma_kase,'status':'completed'})
    print('INTERMEDIATES GENERATED',case,len(steps),'levels',flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,required=True,help='Root of completed replay.py results')
    ap.add_argument('--output',type=Path,help='Historical analysis workspace (default INPUT/legacy_analysis)')
    ap.add_argument('--case',default='all')
    ap.add_argument('--execute',action='store_true',help='Run the original structural and functional post-processing scripts')
    ap.add_argument('--sigma-baseline',type=int,default=100)
    ap.add_argument('--sigma-kase',default='10',help='Comma-separated sigma values, e.g. 10,25,50,75,100')
    ap.add_argument('--require-smoke',action='store_true',help='Require k=1.000,0.760,0.600 for each selected case')
    a=ap.parse_args()
    if a.sigma_baseline<=0:
        ap.error('sigma-baseline must be positive')
    try:
        sigmas=[int(v.strip()) for v in a.sigma_kase.split(',')]
        assert sigmas and all(v>0 for v in sigmas) and len(sigmas)==len(set(sigmas))
    except (AssertionError,ValueError):
        ap.error('sigma-kase must be distinct positive integers, comma-separated')
    root=a.input.expanduser().resolve()
    output=(a.output or root/'legacy_analysis').expanduser().resolve()
    if not root.is_dir():
        ap.error('Input directory missing: '+str(root))
    if (ROOT==output or ROOT in output.parents or output in ROOT.parents) or (root==output or output in root.parents):
        ap.error('Use output outside the checkout and distinct from the input root')
    if output!=root/'legacy_analysis' and (output==root or root in output.parents):
        ap.error('Custom --output must be outside --input')
    reg=yaml.safe_load(REG.read_text(encoding='utf-8'))
    cases=list(reg['experiments']) if a.case=='all' else [a.case]
    if any(c not in reg['experiments'] for c in cases):
        ap.error('Unknown case; available: '+', '.join(reg['experiments']))
    output.mkdir(parents=True,exist_ok=True)
    (output/'tempPlots').mkdir(exist_ok=True)
    notebook_source=ROOT/'notebooks/historical'
    for fname in NOTEBOOKS:
        source=notebook_source/fname
        if source.is_file():
            copy_checked(source,output/fname)
    coverage={}
    for case in cases:
        dest,steps=stage_case(case,reg,root,output)
        coverage[case]=steps
        if a.require_smoke:
            missing=sorted({'000_deg_1.000','009_deg_0.760','015_deg_0.600'}-set(steps))
            if missing:
                raise RuntimeError(f'{case}: missing smoke levels {missing}')
        if a.execute and steps:
            rebuild_case(case,reg,root,dest,steps,a.sigma_baseline,sigmas)
    json_write(output/'layout_report.json',{'source_replay':str(root),'available_levels':coverage,'mode':'historical legacy adapter','executed':a.execute,'notes':'Intermediate PKLs use original postprocessing and may differ from direct raw NIO diagnostic estimates; historical notebooks still need execution in correct order.'})
    print('WORKSPACE',output,'cases',len(coverage),'levels',sum(len(v) for v in coverage.values()))
    if not a.execute:
        print('STAGING ONLY. To generate .pkl intermediates, rerun with --execute.')
    else:
        print('PKL rebuilding attempted for each available experiment; inspect postprocess logs and statuses.')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
