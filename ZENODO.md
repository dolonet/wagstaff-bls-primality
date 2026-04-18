# Zenodo deposit — instructions for the author

This repository ships a `.zenodo.json` at its root. Once the GitHub
repository `dolonet/wagstaff-bls-primality` is linked to a Zenodo
account, publishing a GitHub **Release** will automatically trigger a
Zenodo deposit whose metadata is read from `.zenodo.json`; you click
**Publish** on the Zenodo side to mint a DOI.

## One-time setup (on Zenodo)

1. Log in to <https://zenodo.org> with the account you want to own the
   deposit.
2. Go to *Account → GitHub*: <https://zenodo.org/account/settings/github/>.
3. Flip the repository switch next to `dolonet/wagstaff-bls-primality`
   to **ON**. (If the repository is not listed, click *Sync now* at
   the top of the page.)

## Publishing a release

1. On GitHub, create a tag and release, e.g. `v1.5.0`, with release
   notes describing the deposit.
2. GitHub sends a webhook to Zenodo; Zenodo creates a draft deposit,
   fetches `.zenodo.json`, and attaches the tarball of the tagged
   commit.
3. Go to <https://zenodo.org/deposit> and open the draft. Review the
   metadata (title, authors, description, license). Click **Publish**.
4. Zenodo mints a DOI of the form `10.5281/zenodo.XXXXXXX` and sends
   you the deposit URL.

## Updating the manuscript DOI after publication

After the first publication, paste the resulting DOI into:

* `paper/wagstaff_bls_primality_v1.5.tex` --- replace "DOI forthcoming
  upon publication" in Section 5.2 with `\url{https://doi.org/10.5281/zenodo.XXXXXXX}`.
* `README.md` --- add a "Citation" badge / link.
* `CITATION.cff` --- bump `date-released` and add an `identifiers`
  block with the DOI, then tag a follow-up release so Zenodo mints a
  new version of the record with the DOI already wired in.

## Metadata overrides

The Zenodo deposit page lets you override any field set by
`.zenodo.json`. Override rather than re-tag if you only need a
metadata fix (e.g., a typo in the description).

## Which license is on the deposit

The top-level license declared in `.zenodo.json` is **MIT**, which
governs the code and the JSON certificates. The manuscript under
`paper/` carries CC-BY-4.0 via `LICENSE-PAPER`; Zenodo cannot hold two
licenses on one deposit, so the MIT choice reflects the primary
artifact (software). The description text calls out the CC-BY-4.0
paper split explicitly.

If you prefer two separate deposits (one software-only, one
manuscript-only under CC-BY-4.0), split the upload: create a second
repository `wagstaff-bls-primality-paper` containing just the `paper/`
directory, link that to Zenodo under CC-BY-4.0, and cross-reference
the two via `related_identifiers`.
