#!/usr/bin/env python3
from common import *
from collections import defaultdict
track(Path(__file__),'targeted seed reciprocal full-proteome search methods');track(W/'scripts/common.py','shared methods');policy=json.loads(track(W/'seed_RBH_policy.json','predeclared complete-proteome seed RBH policy').read_text())
seeds=read(W/'tables/priority_7_candidates.tsv');assert len(seeds)==7
bgen=policy['comparator_genome'];bpath=track(A/bgen/(bgen+'.faa'),'complete original comparator proteome')
seqs=dict(fasta(W/'sequences/candidate_loci.faa','prepared candidate protein sequences'));census={x['genome']:x for x in read(P1/'results_v003/qc/genome_census.tsv')};metadata={x['Genome_ID']:x for x in read(Path(C['project_root'])/'metadata/processed/genome_metadata.csv',',')}
def parse(p):
 rr=[]
 with track(p,'new complete-proteome seed raw search').open() as f:
  for line in f:
   q,s,ident,alen,qs,qe,ss,se,ql,sl,ev,bs=line.rstrip().split('\t');rr.append({'query':q,'subject':s,'identity':float(ident),'alignment_length':int(alen),'qstart':int(qs),'qend':int(qe),'sstart':int(ss),'send':int(se),'qlen':int(ql),'slen':int(sl),'evalue':float(ev),'bitscore':float(bs),'query_coverage':(abs(int(qe)-int(qs))+1)/int(ql),'subject_coverage':(abs(int(se)-int(ss))+1)/int(sl)})
 return rr
def bestset(rr,q=None):
 h=[x for x in rr if q is None or x['query']==q]
 if not h:return []
 score=max(x['bitscore'] for x in h);d={}
 for x in h:
  if x['bitscore']==score:
   if x['subject'] not in d or x['alignment_length']>d[x['subject']]['alignment_length']:d[x['subject']]=x
 return [d[k] for k in sorted(d)]
fmt='6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore';base=['/usr/bin/blastp','-evalue',str(policy['evalue']),'-seg','no','-num_threads','1','-max_target_seqs',str(policy['max_target_seqs']),'-outfmt',fmt]
summary=[];forwardall=[];reverseall=[];seqcounts={bgen:sum(1 for k,s in fasta(bpath))};bannotations={x['locus_tag']:x for x in read(A/bgen/(bgen+'.tsv'))}
for i,row in enumerate(seeds,1):
 k=row['stable_protein_id'];ag,al=k.split('|');ap=track(A/ag/(ag+'.faa'),'complete original source proteome');seqcounts[ag]=sum(1 for k,s in fasta(ap));assert max(seqcounts[ag],seqcounts[bgen])<policy['max_target_seqs'];assert seqcounts[ag]==int(census[ag]['FAA_records'])
 qfile=f'sequences/RBH_seed_{i:02}.faa';emit(qfile,{k:seqs[k]});fo=W/'tables'/f'RBH_seed_{i:02}_forward.raw.tsv'
 run([*base,'-query',W/qfile,'-subject',bpath,'-out',fo],f'RBH_seed_{i:02}_forward')
 fr=parse(fo);forwardall.extend({'seed':k,**x} for x in fr);top=bestset(fr);btags={x['subject'] for x in top}
 rseq={bgen+'|'+tag:s for tag,s in fasta(bpath) if tag in btags};rfile=f'sequences/RBH_seed_{i:02}_reverse_queries.faa';emit(rfile,rseq);reverse=[]
 if rseq:
  ro=W/'tables'/f'RBH_seed_{i:02}_reverse.raw.tsv';run([*base,'-query',W/rfile,'-subject',ap,'-out',ro],f'RBH_seed_{i:02}_reverse');reverse=parse(ro);reverseall.extend({'seed':k,**x} for x in reverse)
 for h in top:
  bkey=bgen+'|'+h['subject'];back=bestset(reverse,bkey);selfh=next((x for x in back if x['subject']==al),None);bilateral= min(h['query_coverage'],h['subject_coverage'])>=.85 and selfh is not None and min(selfh['query_coverage'],selfh['subject_coverage'])>=.85
  summary.append({'seed':k,'seed_genome':ag,'seed_assembly':metadata[ag]['Accession'],'comparator_genome':bgen,'comparator_assembly':metadata[bgen]['Accession'],'comparator_best_locus':bkey,'comparator_original_gene':bannotations[h['subject']]['gene'],'comparator_original_product':bannotations[h['subject']]['product'],'forward_best_identity':h['identity'],'forward_query_coverage':h['query_coverage'],'forward_subject_coverage':h['subject_coverage'],'forward_best_bitscore':h['bitscore'],'forward_co_best_loci':len(top),'reverse_co_best_loci':';'.join(ag+'|'+x['subject'] for x in back),'seed_in_reverse_co_best_set':selfh is not None,'both_orientations_bilateral85':bilateral,'scope':'targeted seed complete-proteome reciprocal-best-HSP sets; no orthology/function assertion'})
 if not top:summary.append({'seed':k,'seed_genome':ag,'seed_assembly':metadata[ag]['Accession'],'comparator_genome':bgen,'comparator_assembly':metadata[bgen]['Accession'],'comparator_best_locus':'','seed_in_reverse_co_best_set':False,'both_orientations_bilateral85':False,'scope':'no reported hit at declared settings; not protein absence'})
table('tables/seed_complete_proteome_RBH.tsv',summary);table('tables/seed_RBH_forward_HSPs.tsv',forwardall);table('tables/seed_RBH_reverse_HSPs.tsv',reverseall)
dump('tables/seed_RBH_summary.json',{'seeds_tested':len(seeds),'reciprocal_co_best_seeds':len({x['seed'] for x in summary if x['seed_in_reverse_co_best_set']}),'reciprocal_and_bilateral85_seeds':len({x['seed'] for x in summary if x['both_orientations_bilateral85']}),'complete_proteome_record_counts':seqcounts,'scope':policy['scope'],'complete_source_and_target_proteomes_used':True,'all_cohort_pairwise_RBH':False})
stage_end('seed_RBH');print((W/'tables/seed_RBH_summary.json').read_text())
