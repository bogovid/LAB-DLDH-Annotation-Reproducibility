# Publication commands — not executed

All remote creation, pushes, visibility changes and deposit/release publication require explicit author approval. This local preparation is not that approval. No remote existence or authenticated account was verified: the locally configured account named bogovid failed authentication/network checks. No credential value is copied here.

## Licensing finalized; resolve remaining third-party notices before claiming complete clearance

The author-approved MIT and scoped CC BY 4.0 grants are active. Do not request the same author ownership/license approval again. UniProt CC BY and Pfam CC0 remain; underlying NCBI source material is excluded from author grants. Review provenance/file_rights.tsv and the accession-specific review before a public deposit. No complete third-party or patent clearance is asserted.

## Authenticate and inspect, without publishing

```bash
gh auth login --hostname github.com
gh auth status
gh api user --jq .login
gh repo view bogovid/LAB-DLDH-Annotation-Reproducibility --json name,owner,visibility,url
```

The authenticated login must be the intended bogovid account or an explicitly authorized collaborator. A failed repo lookup does not alone distinguish a nonexistent repository from missing permissions.

## Only after explicit authorization to create/push

The local origin URL is already configured. If authenticated inspection confirms the repository does not exist, create it **without push**:

```bash
gh repo create bogovid/LAB-DLDH-Annotation-Reproducibility --private --description 'Reproducibility evidence for the KJS D330 annotation study'
```

Creating even a private remote is an external action and was not performed. After separately approved push:

```bash
git push -u origin main
```

No visibility change is performed or implied by these instructions. If public visibility is desired, obtain explicit authorization for that action and use the GitHub settings/API accordingly.

## Only after release/publication approval

Following separate publication approval, review of remaining third-party notices and validation, create the version tag and release:

```bash
git tag -a v1.0.1 -m 'KJS D330 reproducibility software 1.0.1'
git push origin refs/tags/v1.0.1
gh release create v1.0.1 --title 'Version 1.0.1' --notes-file docs/RELEASE_NOTES_v1.0.1.md
```

Attach the approved archive and checksum sidecar if desired; do not attach local private audit/configuration files. None of these commands has been executed during preparation.

## Zenodo

The included .zenodo.json now has principal software license mit, version 1.0.1 and explicit file/contribution-level CC BY / CC0 / NCBI-source exceptions. It is not a blanket grant. No DOI, publication date, manuscript DOI or ORCID is invented.

Only after explicit publication approval and review of outstanding third-party notices should a GitHub/Zenodo version be deposited. If both CITATION.cff and .zenodo.json are present, Zenodo uses .zenodo.json for GitHub archiving: https://help.zenodo.org/docs/github/describe-software/ . No remote action has been executed.
