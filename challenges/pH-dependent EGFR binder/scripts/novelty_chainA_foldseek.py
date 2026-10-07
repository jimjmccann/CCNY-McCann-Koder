#!/usr/bin/env python3
"""Chain-A novelty: match sequences to folded structures, extract chain A, foldseek, score Level-4 halves."""
import sys,os,csv,glob,collections,subprocess,json
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
'TYR':'Y','VAL':'V'}
# SCRATCH: volume holding the folded-prediction tree and the foldseek DB.
# REPO: root of the design checkout.  Neither tree is deposited (size).
SCRATCH=os.environ['DESIGN_DATA_ROOT']; REPO='/PATH/TO/CHECKOUT'
csvin,outdir,tag=sys.argv[1],sys.argv[2],sys.argv[3]
def cseq(p,ch='A'):
    seen,seq=set(),[]
    for ln in open(p):
        if not ln.startswith('ATOM') or ln[21]!=ch: continue
        k=ln[22:27]
        if k in seen: continue
        seen.add(k); seq.append(AA3.get(ln[17:20].strip(),'X'))
    return ''.join(seq)
want={r['name']:r['sequence'].strip() for r in csv.DictReader(open(csvin))}
idx=os.path.join(outdir,'seqidx.json')
if os.path.exists(idx):
    byseq={k:v for k,v in json.load(open(idx)).items()}
else:
    byseq=collections.defaultdict(list)
    for d in glob.glob(os.path.join(SCRATCH,'*','**','predictions','*'),recursive=True):
        if not os.path.isdir(d): continue
        m=sorted(glob.glob(os.path.join(d,'*_model_0.pdb')))
        if not m: continue
        try: byseq[cseq(m[0])].append(m[0])
        except Exception: pass
    json.dump(byseq,open(idx,'w'))
q=os.path.join(outdir,'q_'+tag); os.makedirs(q,exist_ok=True)
for f in glob.glob(os.path.join(q,'*.pdb')): os.remove(f)
man=[]
for name,seq in want.items():
    hits=byseq.get(seq,[])
    if not hits: print("  MISSING STRUCTURE:",name); continue
    t=name.rsplit('-',1)[-1]
    pref=[h for h in hits if h.split('/')[-2].endswith('_'+t)]
    src=(pref or hits)[0]; dst=os.path.join(q,name+'.pdb')
    with open(dst,'w') as fh:
        for ln in open(src):
            if ln.startswith(('ATOM','TER')) and (len(ln)<22 or ln[21]=='A'): fh.write(ln)
        fh.write('END\n')
    assert cseq(dst)==seq, name
    man.append((name,len(seq),src.replace(SCRATCH,'$DESIGN_DATA_ROOT'),len(hits)))
print(f"  {len(man)} chain-A PDBs written")
res=os.path.join(outdir,f'hits_{tag}.tsv')
subprocess.run(['foldseek','easy-search',q,'/PATH/TO/foldseek-pdb-db',res,os.path.join(outdir,'tmp_'+tag),
 '--format-output','query,target,fident,alnlen,qlen,tlen,qcov,tcov,alntmscore,qtmscore,ttmscore,evalue,bits',
 '-e','10','--max-seqs','2000','--threads','8'],check=True,stdout=open(os.path.join(outdir,f'fs_{tag}.log'),'w'),stderr=subprocess.STDOUT)
cols="query target fident alnlen qlen tlen qcov tcov alntmscore qtmscore ttmscore evalue bits".split()
by=collections.defaultdict(list)
for ln in open(res):
    f=ln.rstrip("\n").split("\t")
    if len(f)!=len(cols): continue
    d=dict(zip(cols,f))
    for k in cols[2:]: d[k]=float(d[k])
    d['cov']=d['alnlen']/d['qlen']; by[d['query']].append(d)
print(f"\n{'design':24} {'len':>4} | {'STRUCT':>9} {'bestTM@c>.7':>11} {'HIGH(TM>=.8)':>12} | {'SEQ<=30%':>8} {'maxSeqID@c>.5':>13} | nearest")
print("-"*118)
out={}
for name,L,src,nh in sorted(man,key=lambda x:(x[1],x[0])):
    hs=by.get(name,[])
    ok=[h for h in hs if h['cov']>0.70]
    best=max(ok,key=lambda h:h['alntmscore']) if ok else None
    mod=bool(best and best['alntmscore']>=0.50); high=bool(best and best['alntmscore']>=0.80)
    sig=[h for h in hs if h['cov']>0.50]
    sid=max((h['fident'] for h in sig),default=0.0)*100
    nn=best['target'].split('-')[0] if best else '-'
    out[name]={'len':L,'moderate_struct':mod,'high_struct':high,
      'best_tm_cov70':(best['alntmscore'] if best else None),'max_seqid_cov50_pct':sid,
      'nearest':(best['target'] if best else None),'n_hits':len(hs),'source':src}
    print(f"{name:24} {L:>4} | {('DEFEATED' if mod else 'clears'):>9} {(f'{best[chr(39)+chr(39)] if False else best['alntmscore']:.3f}' if best else '  n/a'):>11} {('YES' if high else 'no'):>12} | {('FAIL' if sid>30 else 'pass'):>8} {sid:>12.1f}% | {nn}")
json.dump(out,open(os.path.join(outdir,f'verdict_{tag}.json'),'w'),indent=1)
