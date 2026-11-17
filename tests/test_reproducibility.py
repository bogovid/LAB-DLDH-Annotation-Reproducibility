"""Self-contained public tests. Fresh searches use only the included frozen input set."""
from pathlib import Path
import unittest,subprocess,json,csv,hashlib,re,sys,time,os
ROOT=Path(__file__).resolve().parents[1]
class ReproducibilityTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.out=ROOT/'runs'/('test_'+str(time.time_ns()));env=dict(os.environ);env['PYTHONDONTWRITEBYTECODE']='1'
  p=subprocess.run(['python3',str(ROOT/'NEW/reproduce.py'),'all','--out',str(cls.out)],cwd=ROOT,capture_output=True,text=True,env=env,timeout=600)
  if p.returncode:raise RuntimeError(p.stdout+'\n'+p.stderr)
  cls.result=json.loads((cls.out/'results.json').read_text())
 def test_all_included_stages_pass(self):self.assertEqual(self.result['status'],'PASS');self.assertTrue(self.result['input_hashes_unchanged'])
 def test_historical_filter(self):self.assertEqual(self.result['filter']['records'],7373);self.assertTrue(json.loads((self.out/'filter.json').read_text())['byte_identical'])
 def test_P1_original_table_projection(self):self.assertEqual((self.out/'p1_locus_evidence.tsv').read_bytes(),(ROOT/'data/expected/p1_locus_evidence.tsv').read_bytes())
 def test_denominators_and_sequence_redundancy(self):
  p=self.result['p1'];self.assertEqual([p[k] for k in ['genomes','frozen_CDS_census_sum','historical_screen_loci','target_EC_loci','discordant_loci','generic_EC_discordant','generic_loci','generic_unique_sequences']],[507,1558611,1483,969,514,472,468,75]);self.assertEqual(p['largest_sequence_multiplicities'],[173,68,57])
 def test_fresh_panel_searches(self):self.assertEqual([x['alternative_higher_unique'] for x in self.result['panel']],[75,75]);self.assertEqual(self.result['panel'][0]['baseline_archived_score_mismatches'],0)
 def test_fresh_complete_proteome_reciprocal_searches(self):
  s=self.result['seeds'];self.assertEqual([s[k] for k in ['seeds','reciprocal','bilateral85_both_directions']],[7,6,5]);self.assertTrue(s['byte_identical_scientific_table'])
 def test_fresh_targeted_domains(self):self.assertEqual(self.result['domains']['domains'],374);self.assertTrue(self.result['domains']['data_tokens_identical'])
 def test_context_and_archived_new_counts(self):
  c=self.result['consistency'];self.assertEqual(c['neighborhood'],459);self.assertEqual(c['dppA_opposite_strand'],459);self.assertEqual([x['archived_v005_passing'] for x in c['recorded_full_searches']],[1486,1488]);self.assertEqual([x['archived_v005_bilateral85'] for x in c['recorded_full_searches']],[43,17])
 def test_original_code_integrity(self):
  with (ROOT/'provenance/file_origins.tsv').open() as f:records=list(csv.DictReader(f,delimiter='\t'))
  for x in records:
   if x['path'].startswith('ORIGINAL/'):
    b=(ROOT/x['path']).read_bytes();self.assertEqual(hashlib.sha256(b).hexdigest(),x['source_sha256']);self.assertEqual(hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest(),x['git_blob_oid'])
 def test_public_privacy_scope(self):
  forbidden=[re.compile(r'(?<![A-Za-z0-9:])/(?:home|mnt)/'),re.compile(r'[A-Za-z]:\\Users\\'),re.compile(r'ghp_[A-Za-z0-9]{30,}'),re.compile(r'github_pat_[A-Za-z0-9_]{30,}'),re.compile(r'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----')]
  for p in ROOT.rglob('*'):
   if not p.is_file() or any(x in p.relative_to(ROOT).parts for x in ['runs','.git','__pycache__']):continue
   self.assertNotEqual(p.suffix,'.pyc')
   text=p.read_text(errors='replace')
   for rule in forbidden:self.assertIsNone(rule.search(text),str(p.relative_to(ROOT)))
 def test_version_and_unapproved_licensing_metadata(self):
  d=json.loads((ROOT/'.zenodo.json').read_text());self.assertEqual(d['version'],'1.0.0');self.assertEqual(d['upload_type'],'software');self.assertNotIn('doi',d);
  if (ROOT/'LICENSE').exists():self.assertEqual(d.get('license'),'mit')
  else:self.assertNotIn('license',d);self.assertIn('license approval', (ROOT/'README.md').read_text().lower())
 def test_complete_native_P1_proof_not_mislabelled_self_contained(self):
  d=json.loads((ROOT/'validation/full_P1_external_verification.json').read_text());self.assertFalse(d['self_contained_public_package']);self.assertEqual(len(d['scientific_tables']),14);self.assertTrue(all(x['byte_identical'] for x in d['scientific_tables']))
if __name__=='__main__':unittest.main(verbosity=2)
