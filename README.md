# MS stromal transcriptomics

Reproducible code accompanying **Stromal transcriptional features and DSE network
perturbation in multiple sclerosis**, a manuscript revision prepared for PLOS ONE.
The manuscript is not represented here as an accepted or published article.

[Reproduction guide](docs/REPRODUCIBILITY.md) ·
[Code inventory](provenance/source_inventory.csv) · [License notices](THIRD_PARTY_NOTICES.md)

## Scope and analysis versions

`code/` contains the revised analysis used for the current manuscript.
`archive/legacy/` preserves the original Result1–Result6 analyses for provenance;
they are not alternate evidence for the revised conclusions. Document-production tools and earlier drafts are preserved in a separate local
companion folder; see `provenance/PRIVATE_TOOLS.md`. Third-party packages are
kept separately with their original licenses. Every discovered source file has a
mapping in `provenance/source_inventory.csv`; duplicate workstation/portable copies
map to one maintained revised script.

The revised study uses donor-level inference, nested classifier selection,
MS-only immune/drug reassessment, all 18 spatial sections, and three genuine DSE
virtual-network perturbations. Its reference results include OOF AUC 0.527,
new PBMC AUC 0.500, no BH-significant immune/drug group contrasts, and TIMP3 as
the shared nonself candidate across the three DSE perturbations. Expression and
model outputs do not establish disulfidptotic cell death or clinical treatment effects.

## Quick start

Python 3.11+ is sufficient for the offline source checks and data preparation:

```bash
python run.py check
python -m unittest discover -s tests -v
python run.py list
python run.py prepare --archive /path/to/S1_Data.zip
Rscript environment/install_dependencies.R
python run.py run validate
python run.py run figures
```

The companion `S1_Data.zip` is distributed separately. Its exact SHA-256 and size
are recorded in `config/data_archive.json`. Preparation verifies the archive and
each data-manifest entry, then extracts into `work/reproduction/`, which Git ignores.
Reference tables remain read-only in `results/reference/`; reruns write to the
workspace. No data repository URL or DOI has been invented for this local release.

For a different R installation or user library, put launcher options **before**
the command, e.g. `python run.py --rscript Rscript --r-library /path/to/library run validate`.
Detailed requirements and upstream-stage inputs are in the reproduction guide.

## Repository map

| Directory | Contents |
|---|---|
| `code/` | Revised analyses and optional Illustrator finishing script |
| `config/` | Data checksum, recorded seeds and model parameters |
| `resources/` | Frozen candidate genes, checkpoints and independent TF targets |
| `results/reference/` | All 40 CSV tables from the completed revision |
| `environment/` | Observed package versions and installation guidance |
| `tools/`, `tests/` | Safe data preparation and source/distribution checks |
| `archive/` | Historical analyses and one-time migrations |
| `third_party/` | Separately attributed third-party source and licenses |
| `provenance/` | Source coverage, data accessions and checksums |

The raw whole atlas, raw Visium files, model weights, R workspaces, manuscripts,
and production AI/TIFF files are excluded from the Git source tree. The data archive
contains processed/minimal inputs sufficient for the documented processed-data stages.

## License and citation

Project contributors' own code is MIT licensed. scTenifoldKnk/scTenifoldNet retain
GPL (>=2); the archived DrugReflector code retains its original noncommercial
terms. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

For now, cite the code repository and exact release/commit you use. The manuscript
title above identifies the companion study; a publication citation/DOI can be added
after actual publication or public archival.

Upload instructions: [GitHub 上传说明](docs/UPLOAD_GITHUB.zh-CN.md).
