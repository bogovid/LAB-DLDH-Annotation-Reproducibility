# LAB-DLDH-Annotation-Reproducibility

Version **1.0.0**, prepared locally for `bogovid/LAB-DLDH-Annotation-Reproducibility`.

Reproducibility code and a compact frozen evidence set for the **KJS D330 annotation reproducibility study**. The package examines annotation-versus-similarity disagreement, exact-sequence redundancy, selected reciprocal comparisons and targeted profile-HMM evidence. It does **not** identify candidate enzyme activities.

**Publication status:** unpublished; author approval of the proposed code license and author-generated data/documentation terms is pending. No Zenodo DOI has been assigned. License approval is pending. See [licensing audit](docs/LICENSING_AUDIT.md).

## What runs from the included files

| Stage | Independently executed check | Expected result |
|---|---|---|
| `filter` | Reconstruct the observable archived filter | 7,373 rows; byte-identical output |
| `p1` | Recompute the manuscript-specific locus table from original annotation records, proteins and archived HSPs | 1,483 loci = 969 target EC + 514 other EC; table identical to a projection of the original P1 scientific table |
| `panel` | Fresh BLASTP searches of 75 unique sequences against the frozen 11-entry panel | Alternative-panel preference 75/75 sequences = 468 loci at 11/1 and 9/2 |
| `seeds` | Replay original P2 code against **eight complete proteomes included here** | Seven seed tests; 6 reciprocal co-best sets; 5 bilateral85 in both directions; archived table reproduced byte-for-byte |
| `domains` | Fresh targeted HMMER scan of 131 frozen P2 query sequences | All 374 domain data records identical |
| `consistency` | Recalculate neighbourhood and recorded full-search counts | 459/468 neighbour patterns; recorded v005 full-search counts 1,486/1,488 and coverage counts 43/17 |

The **earliest global BLASTP invocation remains unknown**. Reproducing the archived filter does not recover that command. The 507-genome full search is not rerun by the compact workflow: `consistency` reads explicitly labelled archived v005 outputs. Original GFF parsing and the full 1,558,611-CDS recount also require external originals; the public workflow validates the frozen source-derived context table and sums the included census.

The full native P1 implementation is preserved in `ORIGINAL/p1/run_p1.py`. Its complete 14-table execution was verified against the original external collection during preparation; see `validation/full_P1_external_verification.json`. **Only the manuscript-specific P1 projection is self-contained here.** Broader fermentation/catalogue datasets and unrelated structural/phylogenomic results are intentionally excluded.

## Installation

Linux/WSL, Python 3.10+ standard library, BLASTP **2.16.0+**, HMMER **3.4**, and GNU `time` (`/usr/bin/time`). Tested Python: 3.13.13. The exact-version checks protect comparison with archived scores.

An example environment, requiring an existing conda/mamba installation:

```bash
conda create -n lab-dldh -c conda-forge -c bioconda python=3.13 blast=2.16.0 hmmer=3.4
conda activate lab-dldh
```

This is an installation recipe, not a claimed solver lockfile. The release was tested with existing binaries; no environment was installed during preparation. Verify `blastp -version`, `hmmsearch -h` and `/usr/bin/time --version`. A filter/P1-only audit uses Python alone. Software versions, binary hashes and command-line parameters are recorded by each run.

## Run the complete included workflow

From repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 python NEW/reproduce.py all --out "$PWD/runs/run_001"
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -v
```

Every run must use a **new** directory under `runs/`; existing outputs and included input files are never overwritten. `runs/` is excluded from Git. No network access, original author-host paths or private histories are needed. `--blastp` and `--hmmsearch` can select executables explicitly.

Each run creates `results.json`, `commands.json`, `input_sha256.tsv`, `output_sha256.tsv` and scientific outputs. Original P2 source is copied into an isolated runtime tree. If BLASTP resides elsewhere, only the executable path is relocated in that runtime copy and both original/executed code hashes are recorded; preserved originals remain unchanged.

## Code classifications and provenance

- `ORIGINAL/`: three unchanged analytical source files, with original hashes and Git-source identifiers. These files are not the missing earliest global-search generator.
- `RECONSTRUCTED/`: the observable historical filter and the documented finite-panel policy. Any redundant original filtering conditions remain unidentified.
- `NEW/`: the release-specific portable execution and validation code.

`provenance/file_origins.tsv`, `provenance/original_annotation_inputs.tsv`, the original-code hash record and the data dictionary document sources and transformations. Private source prefixes were replaced with logical source identifiers **only in explicitly marked derived exports**. No private histories, credentials, manuscript drafts, unrelated models or large 507-genome database are included.

Only the eight full proteomes necessary for complete-proteome reciprocal searches are bundled. All other collection-wide proteins/GFFs are represented by selected records, accession metadata and input hashes. This avoids redistributing an unnecessary full NCBI collection.

## Interpretation and licensing

507 selected genomes are not independent population samples. The 468 generic candidates represent 75 exact sequences; 472 generic-EC discordant records form a broader group. A high BLAST score, reciprocal relationship, profile hit or neighbouring gene label is not an activity assay or proof of orthology/operon membership. Archived seven-case coverage and newly declared 43/17-case coverage results must remain distinct.

UniProt-derived reference material retains CC BY 4.0 attribution; the included Pfam models retain CC0 status. These terms are not replaced by the **proposed**, unapproved MIT grant for project-authored code. See `LICENSE.md`, `LICENSES/MIT.proposed.txt`, `docs/THIRD_PARTY_NOTICES.md` and `docs/LICENSING_AUDIT.md`. Documentation/author-derived data terms also need author approval.

Author: **Bogovid Živković**, Independent Researcher, Serbia. `CITATION.cff` and `.zenodo.json` contain supplied author/version metadata, without invented identifiers. `.zenodo.json` deliberately omits the still-unapproved top-level license; it must be completed before deposit. See `docs/PUBLICATION_STEPS.md` for approval-gated commands. Nothing has been pushed or uploaded.
