#!/usr/bin/env python3
"""Analyze raw output from the portable SZ replay, without historical intermediates."""
import argparse
import csv
import gc
import json
import math
from collections import defaultdict
from pathlib import Path

CASES={
    'HR_auto':'auto','NR_auto':'auto',
    'HR_basal4hz':'basal4hz','NR_basal4hz':'basal4hz',
    'HR_kase':'kase','NR_kase':'kase',
}
RECORD_TO_HDF5={
    'basket':'basket_cell','golgi':'golgi_cell','granule':'granule_cell',
    'purkinje':'purkinje_cell','stellate':'stellate_cell',
    'mossy_fibers':'mossy_fibers','glomerulus':'glomerulus',
    'dcn_i':'dcn_i','dcn_p':'dcn_p','io':'io',
}


def analyze_step(folder,case):
    import h5py
    from neo.io import NixIO
    from bsb import from_storage
    status_file=folder/'run_status.json'
    if not status_file.is_file():
        raise ValueError(f'Manca run_status.json: {folder}')
    status=json.loads(status_file.read_text(encoding='utf-8'))
    if status.get('status')!='completed':
        raise ValueError(f'Risultato non completato: {folder}')
    hdf=status.get('hdf5') or {}
    hname=hdf.get('path') or hdf.get('file')
    if not hname or Path(hname).is_absolute() or '..' in Path(hname).parts:
        raise ValueError(f'HDF5 non valido nel run_status: {folder}')
    hdf_path=folder/hname
    if not hdf_path.is_file():
        raise FileNotFoundError(hdf_path)
    with h5py.File(hdf_path,'r') as f:
        populations=tuple(sorted(f['placement'].keys()))
    scaffold=from_storage(str(hdf_path))
    sizes={p:int(len(scaffold.get_placement_set(p))) for p in populations}
    del scaffold
    nio=sorted(folder.glob('*.nio'))
    if len(nio)!=8:
        raise ValueError(f'Attesi 8 NIO MPI in {folder}, trovati {len(nio)}')
    spikes=defaultdict(int)
    counts=defaultdict(int)
    durations=defaultdict(set)
    for file in nio:
        with NixIO(str(file),mode='ro') as reader:
            blocks=reader.read_all_blocks()
            if len(blocks)!=1 or len(blocks[0].segments)!=1:
                raise ValueError(f'Layout NIO inatteso: {file}')
            seen=set()
            for st in blocks[0].segments[0].spiketrains:
                dev=str(st.annotations.get('device',''))
                if not dev.endswith('_record'):
                    continue
                if dev in seen:
                    raise ValueError(f'Device ripetuto in {file}: {dev}')
                seen.add(dev)
                population=dev.removesuffix('_record')
                if population not in RECORD_TO_HDF5:
                    raise ValueError(f'Device sconosciuto: {dev}')
                seconds=float((st.t_stop-st.t_start).rescale('s').magnitude)
                if not math.isfinite(seconds) or seconds<=0:
                    raise ValueError(f'Durata NIO non valida: {file}')
                spikes[population]+=len(st)
                counts[population]+=1
                durations[population].add(round(seconds,9))
        del blocks
        gc.collect()
    k=float(folder.name.split('_deg_')[-1])
    rows=[]
    for pop,num in sorted(spikes.items()):
        hpop=RECORD_TO_HDF5[pop]
        if hpop not in sizes or sizes[hpop]<1:
            raise ValueError(f'Popolazione HDF5 assente o vuota: {hpop}')
        if counts[pop]!=8 or len(durations[pop])!=1:
            raise ValueError(f'Numero file/durate incoerente per {pop}: {folder}')
        seconds=next(iter(durations[pop]))
        rows.append({
            'case':case,'protocol':CASES[case],'step':folder.name,
            'k':k,'sz_pct':round((1-k)*100,4),
            'population':pop,'cell_count':sizes[hpop],
            'n_spikes':num,'duration_s':seconds,
            'firing_rate_hz':round(num/(sizes[hpop]*seconds),9),
        })
    if not rows:
        raise ValueError(f'Nessun registratore di popolazione in {folder}')
    return rows


def plot_figures(rows,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for mode,prefix in [('auto','autor'),('basal4hz','basal'),('kase','evoked')]:
        fig,axes=plt.subplots(1,2,figsize=(11,4.2))
        for ax,pop,title in zip(axes,['purkinje','dcn_p'],['Purkinje','DCN-P']):
            for state in ('NR','HR'):
                entries=sorted((r for r in rows if r['protocol']==mode and r['population']==pop and r['case'].startswith(state+'_')),key=lambda r:r['sz_pct'])
                if entries:
                    ax.plot([r['sz_pct'] for r in entries],[r['firing_rate_hz'] for r in entries],marker='o',label=state)
            ax.set_xlabel('SZ factor (%)')
            ax.set_ylabel('Mean raw firing rate (Hz)')
            ax.set_title(title)
            ax.legend()
        fig.suptitle(f'{mode}: diagnostic raw firing rates (not manuscript estimates)')
        fig.tight_layout()
        fig.savefig(out/f'{prefix}_smoke_firing.png',dpi=170)
        plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,required=True,help='Root passed as --output to replay.py')
    ap.add_argument('--output',type=Path,default=None)
    ap.add_argument('--strict',action='store_true',help='Require all six cases x three smoke levels')
    a=ap.parse_args()
    root=a.input.expanduser().resolve()
    out=(a.output or root/'analysis_smoke').expanduser().resolve()
    if not root.is_dir():
        ap.error(f'Input directory missing: {root}')
    if root==out or root in out.parents and out==root:
        ap.error('Choose a dedicated output directory')
    rows=[]
    processed=[]
    missing=[]
    for case in CASES:
        case_root=root/case
        if not case_root.is_dir():
            missing.append(case)
            continue
        for step in sorted(case_root.glob('[0-9][0-9][0-9]_deg_*')):
            if not step.is_dir():
                continue
            try:
                rr=analyze_step(step,case)
                rows.extend(rr)
                processed.append((case,step.name))
                print('OK',case,step.name,len(rr),'populations',flush=True)
            except Exception as e:
                raise RuntimeError(f'Error in {case}/{step.name}: {e}') from e
    if not processed:
        ap.error(f'No completed replay outputs found below {root}')
    if a.strict:
        required={'000_deg_1.000','009_deg_0.760','015_deg_0.600'}
        actual=set(processed)
        notfound=sorted((case,step) for case in CASES for step in required if (case,step) not in actual)
        if notfound:
            ap.error('Smoke test incomplete: '+repr(notfound))
    out.mkdir(parents=True,exist_ok=True)
    csv_path=out/'smoke_firing_rates.csv'
    with csv_path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    plot_figures(rows,out)
    summary={
        'scope':'Raw population firing rates from newly generated NIO and HDF5 only',
        'limitations':'Diagnostic plots only; not full manuscript smoothed FR, transfer entropy, or Kase stimulus-locked analyses',
        'mpi':'8 files are rank partitions, not independent biological/statistical replicates',
        'steps_processed':[{'case':c,'step':s} for c,s in processed],
        'missing_cases':missing,
        'rows':len(rows),
    }
    (out/'smoke_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(f'{len(processed)} simulations analyzed; {len(rows)} population measurements')
    print('CSV:',csv_path)
    print('Figures:',out)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
