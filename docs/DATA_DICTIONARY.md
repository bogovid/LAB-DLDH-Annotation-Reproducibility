# Included evidence and verification levels

- data/archived/global_blastp.tsv: unchanged 14-column historical HSP output, without header. It is archived evidence, not a newly executed search.
- data/archived/filtered_hits.tsv: unchanged 7,373-row archived three-reference passing table. It includes 1,483 D330 loci. Other reference records are retained only to validate the same observed filter.
- data/p1/annotations.tsv and screen_loci.faa: selected original annotation fields/proteins for exactly the 1,483 historical D330 loci, with genome|locus identity retained. They are extracted from originals; no functional relabelling.
- data/p1/genome_census.tsv: freshly recounted source-derived per-genome census frozen during curation; the package can sum it but cannot independently recount omitted raw annotations.
- data/expected/p1_locus_evidence.tsv: defined field projection of the original P1 1,483-row scientific table; independent oracle for the release-specific P1 reconstruction. The broader 14-table P1 workflow uses external originals.
- data/project/: six nonprivate metadata fields for all 507 accessions, eight complete proteomes required by seed searches and the comparator's original TSV. The proteome subset is necessary for genuine full-proteome reciprocal tests, not an all-genome orthology claim.
- data/seeds/: seven predefined seed identifiers/sequences and original selected-seed policy. The seed code's original hard-coded comparator is retained.
- data/hmmer/: unchanged targeted model panel and 131-query frozen P2 set. Query sequences contain project and reference proteins; resource-specific attribution applies. No structural atlas is included.
- data/expected/targeted_Pfam.original.domtbl: original non-comment domain records; private command-path comments removed explicitly, with source hash retained. Domain comparison uses all data tokens, not timestamps/formatting.
- data/references/: unchanged 11-entry comparator panel; eight external accession-based references and three local operational comparators. Provenance omits only the private local-cache path column; construct limitations remain.
- data/context/neighbors.tsv: source-derived coordinates and transcription-relative offsets; original private GFF prefix replaced by logical source ID. Arithmetic/strand tests do not constitute a fresh GFF parse.
- data/archived/v005_full_11_1.tsv and v005_full_9_2.tsv: archived outputs of independently declared v005 searches, not original raw global search. Full 507-proteome searches are outside the compact execution; these outputs permit transparent count recalculation.
- provenance/original_annotation_inputs.tsv: logical source names and SHA-256 for the original 507 TSV/FAA pairs. No private host path is required.

No date of genome acquisition is asserted from an unverified statement in the earlier candidate archive. Literal inherited taxonomy is retained. Exact loci, exact sequences, genomes, HSPs and selected seeds are distinct units; repeated sequences do not establish independent sampling.
