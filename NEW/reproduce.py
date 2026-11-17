#!/usr/bin/env python3
"""NEW self-contained KJS release workflow. No network access or original-host data required."""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,csv,json,hashlib,subprocess,shutil,sys,time,os,re
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data';COMMANDS=[]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def rows(p,sep='\t'):
 with Path(p).open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter=sep))
def table(p,rr,fields=None):
 rr=list(rr);fields=fields or list(rr[0]);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rr)
def fasta(p):
 data={};k=None
 with Path(p).open() as f:
  for line in f:
   if line.startswith('>'):k=line[1:].split()[0];assert k not in data;data[k]=''
   elif k:data[k]+=line.strip()
 return data
def emit(p,d):
 with p.open('w') as f:
  for k,s in d.items():f.write('>'+k+'\n'+s+'\n')
def relative(p):return os.path.relpath(p,ROOT)
def run(cmd,out,label,env=None):
 argv=[relative(x) if isinstance(x,Path) else str(x) for x in cmd];t=time.monotonic()
 with (out/(label+'.stdout')).open('w') as a,(out/(label+'.stderr')).open('w') as b:
  p=subprocess.run(argv,cwd=ROOT,stdout=a,stderr=b,env=env,timeout=600)
 COMMANDS.append(dict(label=label,argv=argv,exit_code=p.returncode,seconds=time.monotonic()-t))
 if p.returncode:raise RuntimeError(f'{label} failed; inspect {relative(out/(label+".stderr"))}')
def parse_hsp(p):
 rr=[]
 with Path(p).open() as f:
  for i,line in enumerate(f,1):
   v=line.rstrip().split('\t');assert len(v)>=12
   rr.append(dict(q=v[0],s=v[1],identity=float(v[2]),length=int(v[3]),qs=int(v[4]),qe=int(v[5]),ss=int(v[6]),se=int(v[7]),qlen=int(v[8]),slen=int(v[9]),evalue=float(v[10]),bitscore=float(v[11]),line=i))
 return rr

def filter_stage(out):
 run(['python3',ROOT/'RECONSTRUCTED/rebuild_historical_filter.py','--raw',D/'archived/global_blastp.tsv','--expected',D/'archived/filtered_hits.tsv','--out',out/'filtered_hits.tsv','--report',out/'filter.json'],out,'filter')
 r=json.loads((out/'filter.json').read_text());assert r['byte_identical'] and r['reconstructed_rows']==7373
 return dict(records=7373,comparison='byte-identical observable reconstruction, not recovered original generator')

def p1_stage(out):
 ann={x['genome_id']+'|'+x['locus_tag']:x for x in rows(D/'p1/annotations.tsv')};seq=fasta(D/'p1/screen_loci.faa');hist={x['protein_id'] for x in rows(D/'archived/filtered_hits.tsv') if x['reference']=='D_LDH'};assert set(ann)==set(seq)==hist
 best={}
 for i,line in enumerate((D/'archived/global_blastp.tsv').open(),1):
  v=line.rstrip().split('\t')
  if not v[0].startswith('D_LDH|') or float(v[2])<30 or int(v[3])<80:continue
  rank=(float(v[11]),-float(v[10]),float(v[2]),int(v[3]),-i)
  if v[1] not in best or rank>best[v[1]][0]:best[v[1]]=rank,v,i
 result=[];generic=[];ec=set()
 for k in sorted(hist):
  a=ann[k];rank,v,line=best[k];target='1.1.1.28' in re.findall(r'\d+\.\d+\.\d+\.[\d-]+',a['EC_number'])
  if target:ec.add(k)
  if not target and a['product']=='Putative 2-hydroxyacid dehydrogenase':generic.append(k)
  result.append(dict(genome_id=k.split('|')[0],stable_protein_id=k,product=a['product'],EC_number=a['EC_number'],protein_aa=len(seq[k]),protein_sha256=hashlib.sha256(seq[k].encode()).hexdigest(),concordance='both_sources' if target else 'BLAST_only',identity_pct=float(v[2]),query_span_coverage=(abs(int(v[7])-int(v[6]))+1)/int(v[12]),subject_span_coverage=(abs(int(v[9])-int(v[8]))+1)/int(v[13]),alignment_length=int(v[3]),query_length=int(v[12]),subject_length=int(v[13]),bitscore=float(v[11]),evalue=float(v[10]),raw_source_line=line))
 table(out/'p1_locus_evidence.tsv',result);assert (out/'p1_locus_evidence.tsv').read_bytes()==(D/'expected/p1_locus_evidence.tsv').read_bytes()
 census=rows(D/'p1/genome_census.tsv');unique=defaultdict(list)
 for k in generic:unique[hashlib.sha256(seq[k].encode()).hexdigest()].append(k)
 queries={'V5U%03d'%i:seq[ks[0]] for i,(h,ks) in enumerate(sorted(unique.items()),1)};emit(out/'generic_unique.faa',queries)
 table(out/'unique_sequences.tsv',[dict(query='V5U%03d'%i,protein_sha256=h,loci=len(ks),locus_ids=';'.join(sorted(ks))) for i,(h,ks) in enumerate(sorted(unique.items()),1)])
 discord=hist-ec;generic_ec={k for k in discord if '1.1.1.-' in re.findall(r'\d+\.\d+\.\d+\.[\d-]+',ann[k]['EC_number'])}
 r=dict(genomes=len(census),frozen_CDS_census_sum=sum(int(x['CDS_records']) for x in census),historical_screen_loci=len(hist),target_EC_loci=len(ec),discordant_loci=len(discord),generic_EC_discordant=len(generic_ec),generic_loci=len(generic),generic_genomes=len({k.split('|')[0] for k in generic}),generic_unique_sequences=len(unique),largest_sequence_multiplicities=sorted([len(v) for v in unique.values()],reverse=True)[:3],archived_bilateral85=sum(min(x['query_span_coverage'],x['subject_span_coverage'])>=.85 for x in result if x['stable_protein_id'] in discord),oracle='byte-identical projection of archived original P1 scientific table; full 14-table native P1 needs external source collection')
 assert [r[k] for k in ['genomes','frozen_CDS_census_sum','historical_screen_loci','target_EC_loci','discordant_loci','generic_EC_discordant','generic_loci','generic_genomes','generic_unique_sequences','archived_bilateral85']]==[507,1558611,1483,969,514,472,468,467,75,7]
 return r

def panel_stage(out,blastp):
 policy=json.loads((ROOT/'RECONSTRUCTED/panel_policy.json').read_text());qq=fasta(out/'generic_unique.faa');old={x['protein_sha256']:x['query'] for x in rows(D/'archived/P2_query_members.tsv')};alias={q:old[hashlib.sha256(s.encode()).hexdigest()] for q,s in qq.items()};oldbest={}
 for line in (D/'archived/P2_panel.raw.tsv').open():
  v=line.rstrip().split('\t');k=v[0],v[1];oldbest[k]=max(oldbest.get(k,float('-inf')),float(v[11]))
 summaries=[]
 for go,ge in policy['gap_settings']:
  raw=out/f'panel_{go}_{ge}.tsv';run([blastp,'-task','blastp','-query',out/'generic_unique.faa','-subject',D/'references/panel.faa','-matrix','BLOSUM62','-comp_based_stats','2','-evalue','1e-5','-seg','no','-num_threads','1','-max_target_seqs','10000','-gapopen',go,'-gapextend',ge,'-outfmt','6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore','-out',raw],out,f'panel_{go}_{ge}')
  best={}
  for x in parse_hsp(raw):best[x['q'],x['s']]=max(best.get((x['q'],x['s']),float('-inf')),x['bitscore'])
  prefs=sum(max(best[q,r] for r in policy['alternative_external'] if (q,r) in best)>max(best[q,r] for r in policy['D_LDH_associated_external'] if (q,r) in best) for q in qq)
  mismatches=[(q,r) for (q,r),v in best.items() if go==11 and oldbest.get((alias[q],r))!=v];assert prefs==75 and not mismatches
  summaries.append(dict(gap=f'{go}/{ge}',alternative_higher_unique=prefs,represented_loci=468,baseline_archived_score_mismatches=len(mismatches) if go==11 else None))
 return summaries

def seed_stage(out,blastp):
 native=out/'native_seed';
 for d in ['scripts','tables','sequences','preflight','references','logs','manifests','qc']:(native/d).mkdir(parents=True,exist_ok=True)
 for n in ['common.py','seed_RBH.py']:shutil.copyfile(ROOT/'ORIGINAL/p2'/n,native/'scripts'/n)
 source=(native/'scripts/seed_RBH.py').read_bytes();original_hash=hashlib.sha256(source).hexdigest();resolved=shutil.which(blastp) or blastp
 if resolved!='/usr/bin/blastp':
  t=source.decode().replace("'/usr/bin/blastp'",repr(resolved));(native/'scripts/seed_RBH.py').write_text(t)
 project=native/'project';shutil.copytree(D/'project/analysis',project/'analysis');mp=rows(D/'project/metadata/processed/genome_metadata.tsv');(project/'metadata/processed').mkdir(parents=True)
 with (project/'metadata/processed/genome_metadata.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(mp[0]));w.writeheader();w.writerows(mp)
 p1=native/'p1';(p1/'results_v003/qc').mkdir(parents=True);shutil.copyfile(D/'p1/genome_census.tsv',p1/'results_v003/qc/genome_census.tsv')
 (native/'preflight/p1_source_verification.tsv').write_text('path\tobserved_sha256\tmtime_ns\tbytes\n')
 config=dict(p1_root=str(p1.resolve()),project_root=str(project.resolve()),prior_ldh=str(native.resolve()),memory_limit_MiB=2867,analysis_budget_seconds=1800,command_timeout_seconds=180)
 (native/'config.json').write_text(json.dumps(config));shutil.copyfile(D/'seeds/policy.json',native/'seed_RBH_policy.json');shutil.copyfile(D/'seeds/priority_7_candidates.tsv',native/'tables/priority_7_candidates.tsv');shutil.copyfile(D/'seeds/candidate_loci.faa',native/'sequences/candidate_loci.faa')
 env=dict(os.environ);env['PYTHONDONTWRITEBYTECODE']='1';run(['python3',native/'scripts/seed_RBH.py'],out,'native_seed',env)
 actual=native/'tables/seed_complete_proteome_RBH.tsv';assert actual.read_bytes()==(D/'expected/seed_RBH.tsv').read_bytes();shutil.copyfile(actual,out/'seed_RBH.tsv')
 summary=json.loads((native/'tables/seed_RBH_summary.json').read_text());return dict(seeds=summary['seeds_tested'],reciprocal=summary['reciprocal_co_best_seeds'],bilateral85_both_directions=summary['reciprocal_and_bilateral85_seeds'],original_seed_source_sha256=original_hash,executed_seed_source_sha256=sha(native/'scripts/seed_RBH.py'),byte_identical_scientific_table=True,relocation='none' if source==(native/'scripts/seed_RBH.py').read_bytes() else 'RECONSTRUCTED runtime executable-path relocation only')

def domain_stage(out,hmmsearch):
 p=out/'domains.domtbl';run([hmmsearch,'--cpu','1','--cut_ga','--noali','--domtblout',p,D/'hmmer/Pfam_targeted_panel.hmm',D/'hmmer/queries.faa'],out,'hmmsearch')
 tokens=lambda p:[line.split() for line in p.read_text().splitlines() if line.strip() and not line.startswith('#')]
 a=tokens(p);b=tokens(D/'expected/targeted_Pfam.original.domtbl');assert a==b and len(a)==374
 (out/'domain_data.tsv').write_text('\n'.join('\t'.join(x) for x in a)+'\n');return dict(domains=len(a),data_tokens_identical=True,scope='targeted frozen model panel, not full Pfam discovery')

def consistency_stage(out):
 nb=rows(D/'context/neighbors.tsv');by=defaultdict(dict)
 for x in nb:
  if int(x['transcription_CDS_offset']) in [-1,1]:by[x['stable_protein_id']][int(x['transcription_CDS_offset'])]=x
 pat=[]
 for k,v in by.items():
  u=v.get(-1);d=v.get(1)
  if u and d and re.sub(r'_\d+$','',u['neighbor_gene'])=='clpP' and re.sub(r'_\d+$','',d['neighbor_gene'])=='dppA':pat.append((u,d))
 assert len(by)==468 and len(pat)==459 and all(u['same_strand']=='True' and d['same_strand']=='False' for u,d in pat)
 H={x['protein_id'] for x in rows(D/'archived/filtered_hits.tsv') if x['reference']=='D_LDH'};ann={x['genome_id']+'|'+x['locus_tag']:x for x in rows(D/'p1/annotations.tsv')};Dids={k for k in H if ann[k]['EC_number']!='1.1.1.28'};counts=[]
 for label in ['11_1','9_2']:
  best={}
  for x in parse_hsp(D/f'archived/v005_full_{label}.tsv'):
   if x['identity']<30 or x['length']<80:continue
   rank=(-x['bitscore'],x['evalue'],-x['identity'],x['line'])
   if x['s'] not in best or rank<best[x['s']][0]:best[x['s']]=rank,x
  b85=sum(min((abs(x['qe']-x['qs'])+1)/x['qlen'],(abs(x['se']-x['ss'])+1)/x['slen'])>=.85 for k,(_,x) in best.items() if k in Dids);counts.append(dict(gap=label,archived_v005_passing=len(best),archived_v005_bilateral85=b85))
 assert [x['archived_v005_passing'] for x in counts]==[1486,1488] and [x['archived_v005_bilateral85'] for x in counts]==[43,17]
 return dict(neighborhood=459,neighborhood_denominator=468,dppA_opposite_strand=459,context_verification_level='recalculation from frozen source-derived coordinate table; original GFF parsing not rerun',recorded_full_searches=counts,full_search_verification_level='archived output recomputation; no 507-proteome search rerun from this compact package')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['all','filter','p1','panel','seeds','domains','consistency']);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--blastp',default=shutil.which('blastp') or 'blastp');ap.add_argument('--hmmsearch',default=shutil.which('hmmsearch') or 'hmmsearch');a=ap.parse_args();out=a.out.resolve()
 if not out.is_relative_to(ROOT/'runs') or out.exists():raise SystemExit('Use a fresh directory under repository/runs; existing outputs are never overwritten')
 out.mkdir(parents=True);before={relative(p):sha(p) for p in sorted(D.rglob('*')) if p.is_file()};table(out/'input_sha256.tsv',[dict(path=k,sha256=v) for k,v in before.items()])
 stages=['filter','p1','panel','seeds','domains','consistency'] if a.mode=='all' else ['p1','panel'] if a.mode=='panel' else [a.mode]
 result={};versions={}
 try:
  for tool,arg in [(a.blastp,'-version'),(a.hmmsearch,'-h')]:
   if (tool==a.blastp and any(s in stages for s in ['panel','seeds'])) or (tool==a.hmmsearch and 'domains' in stages):
    v=subprocess.check_output([tool,arg],text=True);versions[Path(tool).name]=dict(version_lines=v.splitlines()[:4],binary_sha256=sha(shutil.which(tool) or tool))
    if tool==a.blastp and '2.16.0+' not in v:raise RuntimeError('Exact archived comparisons require tested BLASTP 2.16.0+')
    if tool==a.hmmsearch and 'HMMER 3.4' not in v:raise RuntimeError('Exact domain comparisons require tested HMMER 3.4')
  for stage in stages:
   result[stage]=filter_stage(out) if stage=='filter' else p1_stage(out) if stage=='p1' else panel_stage(out,a.blastp) if stage=='panel' else seed_stage(out,a.blastp) if stage=='seeds' else domain_stage(out,a.hmmsearch) if stage=='domains' else consistency_stage(out)
  assert all(sha(ROOT/k)==v for k,v in before.items()),'Immutable included input changed'
  result.update(status='PASS',python_version=sys.version,tools=versions,input_files=len(before),input_hashes_unchanged=True)
 finally:
  (out/'commands.json').write_text(json.dumps(COMMANDS,indent=2));(out/'results.json').write_text(json.dumps(result,indent=2))
 scientific=[p for p in out.iterdir() if p.is_file() and p.suffix in ['.tsv','.faa','.domtbl']]
 table(out/'output_sha256.tsv',[dict(path=p.name,sha256=sha(p)) for p in sorted(scientific)])
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
