#!/usr/bin/env python3
"""RECONSTRUCTED: the archived filter's observable transformation, not its missing original script.
Evidence: P1 run_p1.py / config.json and archived global HSP/passing table.
Only identity>=30 and alignment length>=80 are identifiable as necessary here.
"""
import argparse,csv,json,hashlib
from pathlib import Path

def rebuild(raw_path,out_path):
 fields=['genome','reference','protein_id','bitscore','identity','alignment_length','query_length','subject_length','evalue']
 rows=[];seen=set()
 with Path(raw_path).open() as f:
  for line_no,line in enumerate(f,1):
   v=line.rstrip('\n').split('\t')
   if len(v)!=14:raise ValueError(('unexpected raw schema',line_no,len(v)))
   if float(v[2])<30 or int(v[3])<80:continue
   key=v[0].split('|')[0],v[1]
   if key in seen:raise ValueError('Duplicate passing key: original selection policy not recoverable; refusing to invent one')
   seen.add(key)
   rows.append(dict(genome=v[1].split('|')[0],reference=key[0],protein_id=v[1],bitscore=float(v[11]),identity=float(v[2]),alignment_length=int(v[3]),query_length=int(v[12]),subject_length=int(v[13]),evalue=float(v[10])))
 rows.sort(key=lambda x:(x['genome'],x['reference'],x['protein_id']))
 with Path(out_path).open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
 return rows

def compare(rows,expected):
 with Path(expected).open() as f:saved=list(csv.DictReader(f,delimiter='\t'))
 key=lambda x:(x['reference'],x['protein_id'])
 a={key(x):x for x in rows};b={key(x):x for x in saved};differences=[]
 for k in sorted(set(a)&set(b)):
  for field in ['bitscore','identity','alignment_length','query_length','subject_length','evalue']:
   if float(a[k][field])!=float(b[k][field]):differences.append(dict(reference=k[0],protein=k[1],field=field,new=a[k][field],archived=b[k][field]))
 return dict(reconstructed_rows=len(rows),archived_rows=len(saved),missing_keys=sorted(set(b)-set(a)),extra_keys=sorted(set(a)-set(b)),numeric_discrepancies=differences,semantic_equal=not differences and set(a)==set(b),interpretation='Observable archived transformation only; unknown upstream parameters and redundant filters remain unidentified')

def main():
 p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--expected',type=Path);p.add_argument('--report',type=Path);a=p.parse_args()
 root=Path(__file__).resolve().parents[1]
 if not a.out.resolve().is_relative_to(root/'runs'):raise SystemExit('Output must be a new file under this release/runs')
 a.out.parent.mkdir(parents=True,exist_ok=True);rows=rebuild(a.raw,a.out)
 if a.expected:
  r=compare(rows,a.expected);r['byte_identical']=a.out.read_bytes()==a.expected.read_bytes()
  if a.report:
   if not a.report.resolve().is_relative_to(root/'runs'):raise SystemExit('Report must be under release/runs')
   with a.report.open('x') as f:json.dump(r,f,indent=2)
  print(json.dumps(r,indent=2));return 0 if r['semantic_equal'] else 1
 return 0
if __name__=='__main__':raise SystemExit(main())
