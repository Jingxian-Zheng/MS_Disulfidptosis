# Reproduction guide

## 1. Prepare one isolated workspace

Use Python 3.11+ and R 4.3.2. The launcher uses the current interpreter for Python
stages and `Rscript` from PATH for R stages. Analysis Python dependencies are in
`environment/requirements-analysis.txt`; only raw bulk preparation needs them.
Install those with `python -m pip install -r environment/requirements-analysis.txt`.

Run `python run.py prepare --archive /path/to/S1_Data.zip` from the repository root.
The archive is supplied alongside this local release, not automatically uploaded
by this distribution. `config/data_archive.json` records its exact identity.
Use `python run.py --workspace work/second-run prepare --archive /path/to/S1_Data.zip`
for another independent workspace.

The launcher sets `REVISION_ROOT`, `LEGACY_STUDY_ROOT` and `MS_CODE_ROOT`. `--atlas`
sets `ATLAS_RDS`; `--bulk-inputs` sets `BULK_INPUT_ROOT`. Paths may contain spaces.
Use the launcher rather than sourcing scripts from an arbitrary working directory.

## 2. Processed-data workflow

| Stage | Required prepared inputs | Output |
|---|---|---|
| `validate` | Donor module scores, OOF predictions, drug matrix, spatial P values, KO tables | `reports/processed_validation.csv` |
| `nested-original` | Original-expression features, metadata | Primary 5-by-4 nested OOF predictions and metrics |
| `nested-rank` | Frozen rank features, all three cohort metadata | Separate rank-model OOF and locked external predictions |
| `virtual-ko` | Stromal counts/metadata and installed scTenifoldNet | Three virtual-perturbation results |
| `case-only` | Primary OOF, reconstructed gene expression, immune sets/caches, archived 198-drug matrix | S6, S11–S14 tables |
| `spatial-fdr` | All 144 planned spatial tests | Corrected S8 table, failed fits retained |
| `figures` | Included processed inputs and complete revision tables | Ten PDF/SVG/PNG statistical figures |

Each stage runs only on explicit selection. Full virtual KO, communication,
spatial GLS and nested selection can be computationally expensive. Reruns may
overwrite outputs within the chosen workspace; the repository reference files
remain unchanged. `validate` recomputes selected statistical quantities from
processed inputs and does not claim to rerun the complete machine-learning training.

`case-only` uses included cached immune estimates when present. Removing these
caches triggers GSVA/CIBERSORT recomputation and requires the documented software.
Drug predictions are archived oncoPredict outputs, not new patient IC50 measurements.

## 3. Upstream stages

`single-nucleus`, `communication` and `minimal-data` require `--atlas /path/to/atlas.rds`.
The exact original curated Seurat object must expose assay `originalexp` integer
UMI counts and metadata fields `patient_id`, `sample_id`, `celltype`, `subtype`,
`condition`, `lesion_type`, `age`, `sex`, and `batch_sn`. The original local RDS is
not bundled and no direct public RDS URL is asserted. Public source: GSE279180 and
Lerma-Martin et al., doi:10.1038/s41593-024-01796-z. Author/source curation is needed
if starting with a different public file format. Keep the donor mapping supplied
with the revision; do not count tissue sections as independent donors.

`download-spatial` obtains the HCA files described in `provenance/hca_file_catalog.json`.
It verifies SHA-256, including files already present, and exits with failure if
any download fails. These raw downloads are several GB. Then `spatial` performs
normalization, all-section models and paired-region tests. The fixed independent
TF targets come from `resources/independent_TF_target_sets.csv`, exactly matching
the completed revision. All 144 planned associations remain in the BH family.

`prepare-bulk` requires `--bulk-inputs` containing the original probe-level files:

```
GSE17048_matrix_exp.txt
GSE17048_series_matrix.txt
GSE17048_GPL6947.txt
GSE32915_matrix_exp.txt
GSE32915_series_matrix.txt
GSE32915_GPL4133.txt
```

Probe mappings must contain `ID` and `gene_name` columns. The data archive supplies
the GSE203241 series matrix and GPL571 annotation under `runtime/`, along with
the original candidate-gene and mRNA universe inputs. Existing reconstructed
gene/feature matrices allow this upstream step to be skipped. GEO source records:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17048,
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32915,
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203241.

## 4. Figures and Illustrator

`python run.py run figures` produces R statistical plots in the workspace.
Optional finishing: use Illustrator's File > Scripts > Other Script to select
`code/11_illustrator_finish.jsx`, then select the prepared workspace when prompted.
It imports PDF vectors and exports editable AI, EPS (main figures), PDF, PNG and
TIFF. Adobe Illustrator 2026, Arial/Arial Bold and local file access are required.
The statistical results do not depend on Illustrator.

## 5. Historical and manuscript tools

`archive/legacy` retains historical analysis behavior for provenance. It includes
the previous cell-level tests, slice screening and high/low expression signature
explorations; these are superseded by `code/` for the revised scientific claims.
Some legacy stages need old RData/models and full original inputs not shipped here.
`08a`/`08b` are split from one mixed-language .sh source and are not a revised KO.

`archive/revision_audit/06b_fixed_module_denominators.R` was a one-time migration,
not a pipeline stage; never run it on newly generated fixed-denominator results.
The separate local manuscript-tools companion requires original Word files, Microsoft Word COM,
and optional Poppler paths (`POPPLER_PATH`). It is not needed for any analysis and is
excluded from the public upload because it contains review-specific writing payloads.

## 6. Verification scope

Offline checks cover complete source mapping, content hashes, Python/R parsing,
configuration parsing, path removal, archive path safety, and real launcher stages.
Recorded local verification additionally extracts the exact data archive, recomputes
key statistics and rebuilds figures from the portable repository. It does not assert
an independent rerun of all upstream raw-data stages or a clean-machine dependency
installation. Package versions are observed direct dependencies, not a complete renv lock.
