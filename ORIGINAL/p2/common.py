import csv,hashlib,json,os,resource,subprocess,sys,time
from pathlib import Path
W=Path(__file__).resolve().parents[1];C=json.loads((W/'config.json').read_text());P1=Path(C['p1_root']);OLD=Path(C['prior_ldh']);A=Path(C['project_root'])/'analysis/prokka_annotations'
resource.setrlimit(resource.RLIMIT_AS,(C['memory_limit_MiB']*1024**2,)*2)
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
for d in ['tables','sequences','logs','manifests','references','qc']:(W/d).mkdir(exist_ok=True)
SOURCES={};RESOURCE=[];START=time.monotonic()
prior_cache={}
with (W/'preflight/p1_source_verification.tsv').open() as f:
 for x in csv.DictReader(f,delimiter='\t'):prior_cache[x['path']]=x

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def track(p,role):
 p=Path(p).resolve();k=str(p)
 if k not in SOURCES:
  cached=prior_cache.get(k);st=p.stat()
  h=cached['observed_sha256'] if cached and int(cached['mtime_ns'])==st.st_mtime_ns and int(cached['bytes'])==st.st_size else sha(p)
  SOURCES[k]={'path':k,'role':role,'bytes':st.st_size,'sha256_before':h,'mtime_ns_before':st.st_mtime_ns}
 return p
def read(p,sep='\t',role='source evidence'):
 with track(p,role).open(newline='') as f:return list(csv.DictReader(f,delimiter=sep))
def table(name,rows,fields=None):
 rows=list(rows);fields=fields or list(dict.fromkeys(k for r in rows for k in r)) or ['no_records']
 p=Path(name);p=p if p.is_absolute() else W/p
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter=',' if p.suffix=='.csv' else '\t');w.writeheader();w.writerows(rows)
def dump(name,x):(W/name).write_text(json.dumps(x,indent=2)+'\n')
def fasta(p,role='original sequences'):
 p=track(p,role);k=None;s=[]
 with p.open() as f:
  for line in f:
   if line.startswith('>'):
    if k is not None:yield k,''.join(s)
    k=line[1:].split()[0];s=[]
   else:s.append(line.strip())
  if k is not None:yield k,''.join(s)
def emit(name,d):
 with (W/name).open('w') as f:
  for k,s in d.items():f.write('>'+k+'\n'+s+'\n')
def stage_end(stage):
 table('manifests/'+stage+'_inputs.tsv',SOURCES.values())
 RESOURCE.append({'command':stage+'_python','argv':json.dumps(sys.argv),'wall_seconds':time.monotonic()-START,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'exit_code':0,'threads_requested':1,'memory_limit_MiB':C['memory_limit_MiB']})
 table('logs/'+stage+'_resource_usage.tsv',RESOURCE)
def run(argv,label,timeout=None):
 if time.monotonic()-START>C['analysis_budget_seconds']:raise RuntimeError('Declared stage runtime budget reached')
 argv=list(map(str,argv));t=time.monotonic();lim=timeout or C['command_timeout_seconds'];timelog=W/'logs'/(label+'.time.txt');env=dict(os.environ)
 with (W/'logs'/(label+'.stdout')).open('w') as out,(W/'logs'/(label+'.stderr')).open('w') as err:
  try:r=subprocess.run(['/usr/bin/time','-v','-o',str(timelog),*argv],stdout=out,stderr=err,timeout=lim,env=env);code=r.returncode
  except subprocess.TimeoutExpired:code='TIMEOUT'
 rec={'command':label,'argv':json.dumps(argv),'wall_seconds':time.monotonic()-t,'peak_RSS_KiB':'','exit_code':code,'threads_requested':1,'memory_limit_MiB':C['memory_limit_MiB']}
 if timelog.exists():
  for line in timelog.read_text().splitlines():
   if 'Maximum resident set size' in line:rec['peak_RSS_KiB']=line.split(':')[-1].strip()
 RESOURCE.append(rec);table('logs/commands_resource_usage.tsv',RESOURCE)
 if code!=0:raise RuntimeError(f'{label} failed: {code}; see logs')
 return rec
# NCBI genetic code 11 amino acids, with bacterial alternative starts handled separately.
BASES='TCAG';AAS='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODONS={a+b+c:aa for (a,b,c),aa in zip(((a,b,c) for a in BASES for b in BASES for c in BASES),AAS)}
STARTS={'TTG','CTG','ATT','ATC','ATA','ATG','GTG'}
def trans(s,start=True):
 aa=''.join(CODONS.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
 if start and s[:3] in STARTS and aa:aa='M'+aa[1:]
 return aa
def rc(s):return s.translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))[::-1]
def feature_extract(contig,x):
 s=contig[x['start']-1:x['end']]
 return rc(s) if x['strand']=='-' else s
