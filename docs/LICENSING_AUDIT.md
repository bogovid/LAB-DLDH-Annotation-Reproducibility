# Licensing and publication audit

Reviewed 2026-11-17. This is a file-scope audit, not an author declaration or a claim to eliminate all possible third-party rights.

| Included group | Provenance / terms | Current disposition |
|---|---|---|
| ORIGINAL analytical Python | Unchanged project code, source hashes and Git identifiers retained | Proposed MIT, subject to author confirmation of rights and scope |
| RECONSTRUCTED, NEW and tests | Explicitly derived/new project code | Proposed MIT, not activated |
| README/docs/reports and author-derived tables | New or derived project material | No license assigned; author may approve a documentation/data license separately |
| Eight external UniProt reference sequences, plus their copies in alignments/query inputs | Official UniProt statement: copyrightable database parts CC BY 4.0 | Retain attribution, source URLs, exact accessions/hashes and transformations; do not replace by MIT |
| Targeted Pfam model file | Official InterPro/Pfam terms: CC0 | Retain existing model metadata and resource attribution; do not mislabel tool license |
| Eight local Prokka full proteomes and selected source records | Local derivatives of accession-identified NCBI molecular data | Necessary limited duplication; NCBI policy is not a transfer of all submitter rights; author redistribution/terms approval pending |
| Archived and expected computational outputs | Original project outputs or documented projections | Mark historical versus new-derived; no blanket third-party license assigned |
| Program dependencies | Not redistributed | Record versions/binary hashes; follow upstream terms on installation |

The official UniProt API documentation search result carried the CC BY 4.0 statement; direct help-page opens returned a JavaScript fallback. The InterPro team explicitly states downloadable Pfam data CC0. These findings support scoped third-party notices, not ownership of author-generated annotations. No inaccessible reference was declared a new verified biochemical construct.

The complete NCBI collection, private terminal histories, raw session transcripts, private configurations, host-specific logs, Git credentials, manuscripts, unrelated structural models and phylogenomic data were excluded. Provenance exports use logical source IDs; sanitization is explicit, while recovered source code remains unchanged.

## Required author action before publication

1. Approve the concrete MIT proposal and scope for project-authored code, or specify another license. Activate it as an actual license only after approval.
2. Approve terms and rights for author-generated data/documentation; retain third-party UniProt CC BY 4.0 and Pfam CC0 notices.
3. Add the approved record-level license/mixed-license information to .zenodo.json. Zenodo requires a license; the omission is intentional pending approval.
4. Authorize the actual GitHub/Zenodo publication separately. No approval is inferred from preparation of this package.

See [Zenodo licenses and rights](https://help.zenodo.org/docs/deposit/describe-records/licenses/). No public DOI, release date or publication identifier is invented.
