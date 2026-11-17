# Publication commands — not executed

All remote creation, pushes, visibility changes and deposit/release publication require explicit author approval. This local preparation is not that approval. No remote existence or authenticated account was verified: the locally configured account named bogovid failed authentication/network checks. No credential value is copied here.

## Complete the approved licensing decision first

Review LICENSES/MIT.proposed.txt and the file-scope licensing audit. After author approval, create LICENSE with the agreed MIT text (without the proposal header) and document the approved author-data/documentation terms. Retain UniProt CC BY 4.0 and Pfam CC0 notices. Add the approved software record license and mixed-license explanation to .zenodo.json. Re-run tests, privacy checks and SHA-256 manifests; commit the approved changes. Do not publish an incomplete license field merely to obtain a DOI.

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

Following final metadata/license approval and validation, create the version tag and release:

```bash
git tag -a v1.0.0 -m 'KJS D330 reproducibility software 1.0.0'
git push origin refs/tags/v1.0.0
gh release create v1.0.0 --title 'Version 1.0.0' --notes-file docs/RELEASE_NOTES_v1.0.0.md
```

Attach the approved archive and checksum sidecar if desired; do not attach local private audit/configuration files. None of these commands has been executed during preparation.

## Zenodo

The included .zenodo.json is a metadata draft with supplied author, affiliation, software type and version. Its license is intentionally pending. No DOI, publication date, manuscript DOI or ORCID has been invented. The related manuscript is identified by title only because no verified publication identifier was supplied.

After approval, either connect the intended GitHub repository to Zenodo and archive the approved versioned release, or create a manual software deposit using the approved archive/metadata. Review mixed file licenses in the deposit UI before publishing. Do not enable an automatic public archive while metadata/rights remain unapproved. If both CITATION.cff and .zenodo.json are present, Zenodo uses .zenodo.json for GitHub archiving: https://help.zenodo.org/docs/github/describe-software/ .
