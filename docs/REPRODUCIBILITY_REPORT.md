# Reproducibility report — version 1.0.0 preparation

All included scientific stages were freshly executed with the curated data. No archived search output is labelled as a fresh search merely because its numbers were recomputed.

| Component | Actual validation | Result / limitation |
|---|---|---|
| Historical filtering | New reconstruction applied to unchanged original raw HSP table | 7,373 records and complete filtered file byte-identical; original generator remains unrecovered |
| Manuscript-specific P1 | New release code reads selected original annotations/proteins and original HSPs | All 1,483 projected rows byte-identical to the defined original P1 table projection |
| Complete native P1 | Original unchanged run_p1.py executed with the full external source collection during preparation | 14 scientific tables byte-identical; 1,037 original input pre/post checks pass; this full run is not self-contained in the compact package |
| P2 selected reciprocity | Original seed_RBH.py/common.py executed in a generated isolated tree using eight included complete proteomes | Seven seeds; 6 reciprocal co-best relationships; 5 bilateral85 in both directions; expected scientific table byte-identical |
| Targeted HMMER | Fresh HMMER 3.4 --cut_ga scan using the included original targeted models and 131 sequences | All 374 domain data records identical; path/time comments excluded from expected table |
| Reference panel | Fresh BLASTP 2.16.0+, 11/1 and 9/2 | 75/75 unique sequences representing 468 loci prefer the alternative group; baseline per-reference maximum scores match archive |
| Genomic context | Recalculation from the included original-derived neighbour coordinate table | 459/468, dppA opposite in all 459; this is not a new original-GFF parse |
| Recorded whole-collection searches | Recalculation from explicitly archived v005 11/1 and 9/2 outputs | 1,486/1,488 passing loci; 43/17 bilateral85 among historical 514; this is not a new 507-proteome search |
| CDS census | Summation of the included frozen per-genome census | 507 and 1,558,611; a raw annotation recount requires the omitted originals |

The two-stage distinction matters: exact downstream table reproduction cannot recover undocumented original search parameters. The original historical screen has 1,483 loci and seven selected-HSP bilateral85 cases; new declared searches must not silently replace those cohorts. Reference preference does not establish substrate activity, misannotation rate or cohort-wide orthology.

## Executed versions and parameters

Python 3.13.13, BLASTP 2.16.0+ (May 4 2025 build), HMMER 3.4. Binary hashes and version text are in validation/results.json. Commands are in validation/commands.json and native_seed_child_commands.json. Logical <repo> substitution removes host prefixes from native child commands without changing flags. Original native code keeps its default scoring settings; exact tool versions are enforced. The release runner explicitly records whether a nonstandard BLASTP path required a runtime-only relocation.

Data hashes are in validation/input_sha256.tsv, scientific output hashes in validation/output_sha256.tsv, and external full-P1 inputs in full_P1_external_input_hashes.tsv. Complete repository-file SHA-256 and archive sidecars accompany the release. Original source hashes and Git blob identities are independently checked.

## Tests and portability

`python -m unittest discover -s tests -v` executes the full included workflow in a fresh output directory. It checks numerical denominators, archived-file equality, genuine reciprocal/domain/panel reexecution, source-code identity, metadata version, licensing-pending state and private-path/credential exclusions. Dependencies must be installed at the documented tested versions; missing tools are not silently skipped.

The final version is also tested from a clean copy extracted from the release archive, with no private source configuration or original research data copied there. The final QA record identifies that test run and archive/commit hashes. A first documentation-wording check failed before README clarified the explicit pending-license status; scientific checks passed, the wording was corrected, and the final suite is rerun. No numerical oracle was changed to make a test pass.

## Not reproduced or not claimed

The earliest global BLASTP command, exact original reannotation/database environment, independent isolate provenance, full-Pfam discovery, all-genome orthology and biochemical candidate activities remain outside verified scope. No inaccessible preprint is declared verified. No full structural atlas, phylogeny, unpublished manuscript or irrelevant broader catalogue is included.
