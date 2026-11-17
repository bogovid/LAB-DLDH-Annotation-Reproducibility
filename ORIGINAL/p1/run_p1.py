#!/usr/bin/env python3
"""P1: stdlib only, read-only inputs, new output directory required."""
import argparse, csv, hashlib, json, platform, re, sys
from pathlib import Path
from collections import defaultdict, Counter
from datetime import datetime, timezone

ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
C=json.loads(args.config.read_text()); O=args.out.resolve()
O.mkdir(parents=True,exist_ok=False)
for d in ['tables','qc','manifests']: (O/d).mkdir()
R=Path(C['project_root']); B=Path(C['global_screen']); P=Path(C['prior_ldh']); A=R/'analysis/prokka_annotations'
inputs={}; qc=[]
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
 return h.hexdigest()
def track(p,role):
 p=Path(p).resolve()
 if str(p) not in inputs: inputs[str(p)]={'path':str(p),'role':role,'bytes':p.stat().st_size,'sha256_before':digest(p),'mtime_ns_before':p.stat().st_mtime_ns}
 return p
def read(p,sep='\t',role='evidence'):
 p=track(p,role)
 with p.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter=sep))
def table(name,rows,fields=None):
 rows=list(rows)
 if not fields:fields=list(dict.fromkeys(k for row in rows for k in row)) if rows else ['no_records']
 sep=',' if name.endswith('.csv') else '\t'
 with (O/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter=sep);w.writeheader();w.writerows(rows)
def check(name,condition,detail):qc.append({'check':name,'status':'PASS' if condition else 'FAIL','detail':str(detail)})
def ecset(s):return set(re.findall(r'\d+\.\d+\.\d+\.[\d-]+',s))
def fasta(p):
 track(p,'original protein sequences or reference');name=None;seq=[]
 with p.open() as f:
  for line in f:
   if line.startswith('>'):
    if name is not None:yield name,''.join(seq)
    name=line[1:].split()[0];seq=[]
   else:seq.append(line.strip())
  if name is not None:yield name,''.join(seq)
def sid(g,l):return g+'|'+l
track(args.config,'analysis configuration');track(Path(__file__),'analysis script')
meta=read(R/'metadata/processed/genome_metadata.csv',',','genome census');md={x['Genome_ID']:x for x in meta};check('507 unique genome metadata records',len(md)==len(meta)==507,len(meta))
cat=read(R/'metadata/processed/gene_catalog_analysis_ready.csv',',','catalog used by source builder')
validated=read(R/'metadata/processed/gene_catalog_validated.csv',',','validated catalog comparison')
track(R/'scripts/python/build_fermentation_matrix.py','historical matching method')
cat_audit=[]
for ec in sorted(set(e for x in cat for e in ecset(x['EC_number']))):
 rr=[x for x in cat if ec in ecset(x['EC_number'])];vv=[x for x in validated if ec in ecset(x['EC_number'])]
 cat_audit.append({'EC':ec,'analysis_ready_rows':len(rr),'gene_ids':';'.join(x['Gene_ID'] for x in rr),'modules':';'.join(x['Module_ID'] for x in rr),'validated_rows':len(vv)})
table('tables/catalog_ec_reuse.tsv',cat_audit)
marker=read(R/'analysis/functional_modules/fermentation/fermentation_marker_hits.tsv',role='historical catalog assignments')
assign=defaultdict(list)
for x in marker:assign[sid(x['Genome_ID'],x['locus_tag'])].append(x)
catalog_loci=set(assign)
table('tables/catalog_locus_assignments.tsv',({'stable_protein_id':k,**r} for k in sorted(assign) for r in assign[k]))
filtered=read(B/'global_507_filtered_hits.tsv',role='historical filtered BLASTP')
F={(x['reference'],x['protein_id']):x for x in filtered}; D={x['protein_id'] for x in filtered if x['reference']=='D_LDH'}
raw=[]; raw_by=defaultdict(list); rawpath=track(B/'global_blastp.tsv','raw global BLASTP HSPs')
fields=['qseqid','sseqid','identity','alignment_length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','query_length','subject_length']
with rawpath.open() as f:
 for line_no,line in enumerate(f,1):
  v=line.rstrip().split('\t');assert len(v)==14,(line_no,len(v));r=dict(zip(fields,v));r['source_line']=line_no;r['reference']=r['qseqid'].split('|')[0]
  for k in ['identity','evalue','bitscore']:r[k]=float(r[k])
  for k in ['alignment_length','mismatch','gapopen','qstart','qend','sstart','send','query_length','subject_length']:r[k]=int(r[k])
  r['query_span_coverage']=(abs(r['qend']-r['qstart'])+1)/r['query_length'];r['subject_span_coverage']=(abs(r['send']-r['sstart'])+1)/r['subject_length']
  r['alignment_over_query']=r['alignment_length']/r['query_length'];r['alignment_over_subject']=r['alignment_length']/r['subject_length']
  r['historical_filter_pass']=r['identity']>=C['historical_filter']['identity_min'] and r['alignment_length']>=C['historical_filter']['alignment_length_min']
  raw.append(r);raw_by[(r['reference'],r['sseqid'])].append(r)
regen={k for k,rr in raw_by.items() if any(x['historical_filter_pass'] for x in rr)}
check('historical filtered unique keys reproduced',regen==set(F),f'expected={len(F)} reproduced={len(regen)} missing={len(set(F)-regen)} extra={len(regen-set(F))}')
check('filtered unique rows',len(F)==len(filtered),len(filtered))
filter_mismatches=[]
for k,x in F.items():
 matches=[r for r in raw_by.get(k,[]) if r['historical_filter_pass'] and all(float(x[a])==float(r[b]) for a,b in [('identity','identity'),('alignment_length','alignment_length'),('bitscore','bitscore'),('query_length','query_length'),('subject_length','subject_length'),('evalue','evalue')])]
 if not matches:filter_mismatches.append(k)
check('every filtered record matches a passing raw HSP',not filter_mismatches,filter_mismatches[:10])
table('tables/global_blast_hsp_evidence.tsv',raw)
# Prior domains are joined only by independently checked original-protein hashes.
domains=defaultdict(list);prior_hash={}; priorclass={x['id']:x for x in read(P/'phase2c1/tables/candidate_classification.tsv',role='prior functional comparison, not orthology')}
domain_long=[]
for phase,memberfile in [('phase2c0','query_members.tsv'),('phase2c1','candidate_members.tsv')]:
 members=read(P/phase/'tables'/memberfile,role='prior sequence-to-query map'); dd=defaultdict(list)
 for x in read(P/phase/'tables/pfam_domains.tsv',role='prior Pfam profile evidence'):dd[x['query']].append(x)
 for m in members:
  k=sid(m['genome'],m['locus']);prior_hash[(phase,k)]=m['protein_sha256']
  for d in dd[m['query']]:
   rec={'stable_protein_id':k,'source_phase':phase,'query':m['query'],'expected_protein_sha256':m['protein_sha256'],**d};domain_long.append(rec);domains[k].append(rec)
chim=read(B/'chimera_dual_candidates_corrected.tsv',role='historical chimera list'); chimids={x['protein'] for x in chim}
arch=read(P/'phase2c0/tables/chimera_architecture.tsv',role='prior segment comparison')
for x in arch:
 k=next((k for k in chimids if k.endswith('|'+x['locus'])),None)
 if k:
  prior_hash[('phase2c0',k)]=x['sha256']
  for name,acc,sp in [('Aminotran_1_2','PF00155',x['Pfam_Aminotran_1_2']),('2-Hacid_dh_C','PF02826',x['Pfam_2_Hacid_dh_C'])]:
   if sp:
    domains[k].append({'source_phase':'phase2c0_chimera_architecture','expected_protein_sha256':x['sha256'],'profile':name,'accession':acc,'aa_start':sp.split('-')[0],'aa_end':sp.split('-')[1],'independent_evalue':'not_in_summary'})
reverse=read(P/'phase2c1/tables/reciprocal_orientation_checks.tsv',role='panel orientation check, not genome-wide RBH')
reverse_by=defaultdict(list)
for x in reverse:reverse_by[x['id']].append(x)
# Inventory known prior LDH groups; these labels are operational sequence inventories.
prior_groups=defaultdict(set)
for g in ['D330','L309','L320']:
 for x in read(P/'tables'/f'{g}_candidate_inventory.tsv',role='prior LDH inventory'):
  # schema varies: retain raw inventory separately and recover IDs below
  vals='|'.join(str(v) for v in x.values()); tags=re.findall(r'LABG\d{4}_\d+',vals)
  for l in tags: prior_groups[sid('LABG_'+l[4:8],l)].add(g)
# Build a locus-level table for catalog/BLAST/prior candidates, and a census of ALL original CDS.
selected=set(assign)|{r['sseqid'] for r in raw}|set(priorclass)|set(domains)|set(prior_groups)|chimids
annotations={};sequences={};allcounts=Counter();genome_qc=[];ec_ids=set();lec_ids=set();dup_tags=[]
for ix,(g,m) in enumerate(sorted(md.items()),1):
 t=A/g/(g+'.tsv');f=A/g/(g+'.faa');seen=set();cds=0;stats=Counter();target=set()
 if t.exists():
  for x in read(t,role='original Prokka annotation'):
   if x.get('ftype')!='CDS':continue
   k=sid(g,x['locus_tag']);cds+=1
   if k in seen:dup_tags.append(k)
   seen.add(k);ecs=ecset(x.get('EC_number',''));ed=C['modules']['D_LDH'] in ecs
   if ed:ec_ids.add(k);selected.add(k)
   if C['modules']['L_LDH'] in ecs:lec_ids.add(k);selected.add(k)
   es='target_EC_recorded' if ed else ('EC_field_empty' if not ecs else 'other_EC_recorded')
   bs='filtered_hit' if k in D else ('raw_only_hit' if ('D_LDH',k) in raw_by else 'no_reported_hit')
   allcounts[(es,bs)]+=1;stats[(es,bs)]+=1
   if k in selected:annotations[k]=x;target.add(x['locus_tag'])
 if f.exists():
  fn=0;fs=set()
  for name,s in fasta(f):
   fn+=1
   if name in target:sequences[sid(g,name)]=s;fs.add(name)
 else:fn=0;fs=set()
 genome_qc.append({'genome':g,'species':m['Species_current'],'accession':m['Accession'],'duplicate_accession_flag':m['Duplicate_accession'],'annotation_tsv':str(t),'annotation_available':t.exists(),'faa_available':f.exists(),'CDS_records':cds,'FAA_records':fn,'selected_loci':len(target),'selected_FAA_missing':len(target-fs)})
 if ix%50==0:print(f'Prokka audit {ix}/507; selected={len(annotations)}',flush=True)
check('all 507 annotation TSVs and FAAs available',all(x['annotation_available'] and x['faa_available'] for x in genome_qc),len(genome_qc))
check('no duplicate CDS locus tags within genome',not dup_tags,dup_tags[:10])
check('all selected IDs map to original annotations',not(selected-set(annotations)),sorted(selected-set(annotations))[:20])
check('all selected original FAA sequences found',set(annotations)<=set(sequences),len(set(annotations)-set(sequences)))
check('selected FAA lengths match Prokka CDS length',all(int(annotations[k]['length_bp']) in (3*len(s),3*(len(s)+1)) for k,s in sequences.items()),'CDS length permits terminal stop codon')
hashes={k:hashlib.sha256(s.encode()).hexdigest() for k,s in sequences.items()}
check('all prior sequence mappings match original FAA',all(hashes.get(k)==h for (_,k),h in prior_hash.items()),sum(hashes.get(k)!=h for (_,k),h in prior_hash.items()))
# Verify raw subject/reference lengths, without treating score/similarity as function.
check('raw subject lengths match original FAA',all(r['subject_length']==len(sequences.get(r['sseqid'],'')) for r in raw),'all 10,292 HSPs')
refs=dict(fasta(B/'reference_proteins.faa'))
check('raw query lengths match original reference FAA',all(r['query_length']==len(refs.get(r['qseqid'],'')) for r in raw),list(refs))
check('EC catalog D-LDH IDs equal direct Prokka EC IDs',{k for k,rr in assign.items() if any(C['modules']['D_LDH'] in ecset(x['EC_number']) for x in rr)}==ec_ids,len(ec_ids))
# Exact copy identity never collapses distinct loci.
hash_ids=defaultdict(list)
for k in selected:
 if k in hashes:hash_ids[hashes[k]].append(k)
evidence=[];discord=[]
for k in sorted(selected):
 g,l=k.split('|');a=annotations.get(k,{});seq=sequences.get(k);rr=raw_by.get(('D_LDH',k),[]);passing=[x for x in rr if x['historical_filter_pass']]
 best=max(passing or rr,key=lambda x:(x['bitscore'],-x['evalue'],x['identity'],x['alignment_length'],-x['source_line'])) if rr else None
 ecs=ecset(a.get('EC_number',''));ec_state='target_EC_recorded' if k in ec_ids else ('annotation_unavailable' if not a else ('EC_field_empty' if not ecs else 'other_EC_recorded'))
 blast_state='filtered_hit' if k in D else ('raw_only_hit' if rr else 'no_reported_hit')
 label='both_sources' if k in ec_ids and k in D else 'EC_only' if k in ec_ids else 'BLAST_only' if k in D else 'neither_target_source'
 valid_domains=[x for x in domains.get(k,[]) if hashes.get(k)==x['expected_protein_sha256']]
 significant=[x for x in valid_domains if x['independent_evalue']=='not_in_summary' or float(x['independent_evalue'])<=C['domain_independent_evalue_max']]
 ds=';'.join(sorted(set(x['accession'].split('.')[0]+'@'+str(x['aa_start'])+'-'+str(x['aa_end']) for x in significant)))
 # Preserve positional architecture AND a coordinate-free profile set for stratification.
 dset=';'.join(sorted(set(x['accession'].split('.')[0] for x in significant)))
 dstate='profiles_recorded' if significant else ('screened_no_profile_below_P1_cutoff' if valid_domains else 'unavailable')
 pc=priorclass.get(k,{});gap=pc.get('canonical_minus_alternative_relative_score','');priorok=hashes.get(k)==pc.get('protein_sha256')
 if not priorok:pc={};gap=''
 flags=[]
 if k in chimids:flags.append('P2_chimera_segments_and_genomic_context')
 if label in ['EC_only','BLAST_only']:flags.append('P2_discordant_annotation_homology')
 if best and min(best['query_span_coverage'],best['subject_span_coverage'])<0.85:flags.append('P2_partial_or_architecture')
 if dstate=='unavailable' and (k in ec_ids or k in D):flags.append('P2_domain_scan_unavailable')
 if gap and float(gap)<0.10:flags.append('P2_competing_substrate_controls')
 if k in ec_ids or k in D:flags.extend(['P2_genome_wide_reciprocal_homology','P2_orthogroups_and_gene_tree'])
 row={'genome_id':g,'stable_protein_id':k,'locus_tag':l,'species':md[g]['Species_current'],'accession':md[g]['Accession'],'gene':a.get('gene',''),'product':a.get('product',''),'EC_number':a.get('EC_number',''),'COG':a.get('COG',''),'annotation_state':ec_state,'BLAST_state':blast_state,'concordance':label,'protein_aa':len(seq) if seq is not None else '', 'protein_sha256':hashes.get(k,''),'identical_selected_locus_copies':len(hash_ids.get(hashes.get(k),[])),'catalog_assignment_rows':len(assign.get(k,[])),'catalog_gene_ids':';'.join(sorted(set(x['Gene_ID'] for x in assign.get(k,[])))),'catalog_modules':';'.join(sorted(set(x['Module_ID'] for x in assign.get(k,[])))),'prior_inventory_groups':';'.join(sorted(prior_groups[k])),'D_LDH_raw_HSPs':len(rr),'D_LDH_passing_HSPs':len(passing),'D_LDH_query':best['qseqid'] if best else '', 'identity_pct':best['identity'] if best else '', 'query_span_coverage':best['query_span_coverage'] if best else '', 'subject_span_coverage':best['subject_span_coverage'] if best else '', 'alignment_length':best['alignment_length'] if best else '', 'query_length':best['query_length'] if best else '', 'subject_length':best['subject_length'] if best else '', 'bitscore':best['bitscore'] if best else '', 'evalue':best['evalue'] if best else '', 'raw_source_line':best['source_line'] if best else '', 'domain_state':dstate,'domain_architecture':ds,'domain_profile_set':dset,'domain_sources':';'.join(sorted(set(x['source_phase'] for x in valid_domains))),'orthology_state':'unavailable_no_validated_locus_orthogroup_join','panel_orientation_state':'available_not_genome_wide_RBH' if reverse_by[k] else 'unavailable','prior_similarity_bin':pc.get('primary_similarity_bin',''),'prior_review_reasons':pc.get('review_reasons',''),'canonical_control':pc.get('best_canonical_DLDH_control',''),'alternative_control':pc.get('best_alternative_substrate_control',''),'canonical_minus_alternative_relative_score':gap,'chimera_flag':k in chimids,'validation_flags':';'.join(flags)}
 evidence.append(row)
 if label in ['EC_only','BLAST_only']:discord.append(row)
table('tables/protein_locus_evidence.tsv',evidence);table('tables/dldh_evidence.tsv',[x for x in evidence if x['stable_protein_id'] in ec_ids|D]);table('tables/dldh_discordant.tsv',discord,list(evidence[0]))
table('tables/domain_hits_prior.tsv',domain_long);table('tables/panel_orientation_prior.tsv',reverse)
table('qc/genome_census.tsv',genome_qc)
table('tables/chimera_verification.tsv',({'stable_protein_id':k,'protein_aa':len(sequences[k]),'sha256':hashes[k],'EC':annotations[k]['EC_number'],'gene':annotations[k]['gene'],'D_LDH_filtered':k in D,'prior_hash_matches':all(hashes[k]==h for (_,kk),h in prior_hash.items() if kk==k)} for k in sorted(chimids)))
check('five identical 543-aa chimera loci',len(chimids)==5 and len({hashes[k] for k in chimids})==1 and all(len(sequences[k])==543 for k in chimids),sorted(chimids))
table('tables/concordance_all_CDS.tsv',({'EC_evidence_state':a,'BLAST_evidence_state':b,'loci':n} for (a,b),n in sorted(allcounts.items())))
conc=[]
for scope,universe in [('DLDH_evidence_union',ec_ids|D),('selected_catalog_blast_prior_union',selected)]:
 for e in [True,False]:
  for b in [True,False]:conc.append({'scope':scope,'target_EC_recorded':e,'filtered_DLDH_hit':b,'loci':sum((k in ec_ids)==e and (k in D)==b for k in universe)})
# Full-proteome neither cell is explicit, without mistaking it for confirmed absence.
for e in [True,False]:
 for b in [True,False]:conc.append({'scope':'all_original_CDS','target_EC_recorded':e,'filtered_DLDH_hit':b,'loci':sum(n for (es,bs),n in allcounts.items() if (es=='target_EC_recorded')==e and (bs=='filtered_hit')==b)})
table('tables/concordance.tsv',conc);table('tables/concordance.csv',conc)
# Multiple orthogonal summaries; do not sum across dimensions.
strata=[]
for dim in ['product','gene','genome_id','domain_state','domain_profile_set','orthology_state','prior_similarity_bin']:
 cc=Counter((x['concordance'],x[dim] or 'unavailable_or_blank') for x in discord)
 strata.extend({'dimension':dim,'concordance':a,'stratum':b,'loci':n} for (a,b),n in sorted(cc.items()))
for dim in ['identity_pct','query_span_coverage','subject_span_coverage']:
 cc=Counter()
 for x in discord:
  v=x[dim]
  if v=='':bucket='unavailable_no_reported_alignment'
  elif dim=='identity_pct':bucket='<40' if v<40 else '40-<70' if v<70 else '70-<90' if v<90 else '>=90'
  else:bucket='<0.50' if v<0.5 else '0.50-<0.85' if v<0.85 else '>=0.85'
  cc[x['concordance'],bucket]+=1
 strata.extend({'dimension':dim,'concordance':a,'stratum':b,'loci':n} for (a,b),n in sorted(cc.items()))
table('tables/discordance_strata.tsv',strata)
table('tables/priority_candidates.tsv',[x for x in evidence if x['stable_protein_id'] in ec_ids|D or x['chimera_flag']])
# Extensible target configuration: L-LDH awaits a matching genome-wide reference-panel search.
module_rows=[]
for label,ec in C['modules'].items():
 ids={k for k,a in annotations.items() if ec in ecset(a.get('EC_number',''))}
 module_rows.append({'module':label,'EC':ec,'annotated_distinct_loci':len(ids),'catalog_assignment_rows':sum(ec in ecset(x['EC_number']) for x in marker),'sequence_search_state':'historical_global_search_audited' if label=='D_LDH' else 'not_in_original_global_three_reference_search','interpretation':'no activity assignment'})
table('tables/module_readiness.tsv',module_rows)
# Audit manifests for every consumed prior output; do not imply a full-tree integrity audit.
manifest_checks=[]
for phase in ['phase2c0','phase2c1']:
 mm=read(P/phase/'manifests/output_manifest_sha256.tsv',role='prior output manifest')
 for x in mm:
  p=(P/x['path']).resolve()
  if str(p) in inputs:manifest_checks.append({'path':str(p),'manifest':str(P/phase/'manifests/output_manifest_sha256.tsv'),'expected':x['sha256'],'observed':inputs[str(p)]['sha256_before'],'matches':x['sha256']==inputs[str(p)]['sha256_before']})
check('consumed prior outputs match published checksums',all(x['matches'] for x in manifest_checks),len(manifest_checks))
table('qc/prior_manifest_checks.tsv',manifest_checks)
summary={'original_CDS_records':sum(allcounts.values()),'genomes':len(md),'catalog_rows':len(marker),'catalog_distinct_loci':len(catalog_loci),'catalog_repeated_locus_rows':len(marker)-len(catalog_loci),'DLDH_EC_loci':len(ec_ids),'DLDH_catalog_EC_rows':sum(C['modules']['D_LDH'] in ecset(x['EC_number']) for x in marker),'filtered_DLDH_loci':len(D),'both_sources':len(ec_ids&D),'BLAST_only':len(D-ec_ids),'EC_only':len(ec_ids-D),'DLDH_union':len(ec_ids|D),'chimera_loci':len(chimids),'chimera_sequence_hashes':len({hashes[k] for k in chimids}),'selected_loci':len(selected),'discordant_loci':len(discord),'discordant_domain_available':sum(x['domain_state']!='unavailable' for x in discord),'prior_manifest_records_checked':len(manifest_checks),'historical_filter_reproduced':regen==set(F),'filter_origin':'empirically exactly reconstructed from raw/filtered records; original filtering script not located','software':{'python':sys.version,'platform':platform.platform()},'config':C}
for name,actual,expected in [('historical EC count',len(ec_ids),969),('historical filtered count',len(D),1483),('historical BLAST-only count',len(D-ec_ids),514)]:check(name,actual==expected,f'{actual} expected {expected}')
# Rehash every consumed source at completion; no source writes anywhere in this script.
for ix,x in enumerate(inputs.values(),1):
 p=Path(x['path']);x['sha256_after']=digest(p);x['unchanged']=x['sha256_after']==x['sha256_before'];x['mtime_ns_after']=p.stat().st_mtime_ns
 if ix%200==0:print(f'Provenance rehash {ix}/{len(inputs)}',flush=True)
check('all consumed inputs unchanged',all(x['unchanged'] for x in inputs.values()),len(inputs))
table('manifests/source_manifest_sha256.tsv',inputs.values());table('qc/checks.tsv',qc)
summary['QC_failures']=[x for x in qc if x['status']=='FAIL']
(O/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(O/'manifests/run_metadata.json').write_text(json.dumps({'completed_utc':datetime.now(timezone.utc).isoformat(),'argv':sys.argv,'config_sha256':digest(args.config),'script_sha256':digest(Path(__file__))},indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
if summary['QC_failures']:sys.exit(2)
