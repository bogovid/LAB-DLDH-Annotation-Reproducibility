# File-level licensing — active author grants, explicit third-party exclusions

The author has approved the licenses below. **There is no blanket license over every byte of this repository.** The authoritative scope is recorded in provenance/file_rights.tsv and the accession/sequence notices linked below.

| Material | Effective rights / license | Exclusions |
|---|---|---|
| Project-authored code and configuration in ORIGINAL/, RECONSTRUCTED/, NEW/ and tests/ | **MIT**, LICENSE and LICENSES/MIT.txt | Does not license datasets, underlying sequences, dependencies or third-party materials |
| Author-owned README, docs, metadata, provenance, validation reports and original result/derived-data contributions | **CC BY 4.0**, LICENSES/CC-BY-4.0.md | Only the author's contributions; underlying third-party material remains outside this grant |
| Eight UniProt external reference sequence copies and associated UniProt content | **UniProt Consortium CC BY 4.0**, retained attribution | Not a grant by the project author; no assumption of assay-construct identity or other-IP clearance |
| Targeted Pfam model file | **Upstream CC0-1.0** | Applies to model data, not the HMMER program or other files |
| NCBI-derived proteomes, selected source protein strings and local comparator sequences | **No new repository grant for underlying third-party sequence content**; NCBI molecular-data policy and any applicable submitter rights/notices remain | Neither MIT nor the author's CC BY grant transfers ownership or clears third-party rights |
| Mixed files (panel/query FASTA, raw alignments containing sequence strings, annotations/metadata) | Contribution/record-level split: author-generated results/selection/annotations CC BY 4.0; UniProt copies retain upstream CC BY; NCBI-derived content excluded from author grant | A single SPDX identifier must not be interpreted as licensing all components |

UniProt and Pfam attribution/terms are in docs/THIRD_PARTY_NOTICES.md. The accession review is in provenance/uniprot_accession_review.tsv and provenance/ncbi_accession_review.tsv; the latter includes all 507 source metadata accessions and distinguishes actually attempted checks from unreviewed records. HMMER sequence provenance is in provenance/HMMER_query_origins.tsv.

Accession-specific notices could not be fully inspected through the available record views. No specific prohibition on sequence redistribution was established, but complete third-party clearance is **not** claimed. P37666 is cited in patent publication EP3280794B1; that citation alone does not establish a restriction on this archive or the current patent scope/status. No patent freedom-to-operate conclusion is offered. See docs/ACCESSION_RIGHTS_REVIEW.md.

LICENSES/MIT.proposed.txt is retained only as a superseded preparation record. The active MIT grant is LICENSE. The original research data and analytical source files have not been rewritten to add license headers; scope is supplied by these notices and manifests.
