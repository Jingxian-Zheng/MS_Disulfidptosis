"""Portable command-line entry point for the revised MS study."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
STAGES = {
    "download-spatial": ("python", "02_download_spatial.py", "Download checksum-verified public HCA objects"),
    "single-nucleus": ("r", "03_single_nucleus.R", "Normalize the complete curated atlas; donor pseudobulk"),
    "prepare-bulk": ("python", "04_prepare_bulk.py", "Reconstruct bulk features from original probe inputs"),
    "nested-rank": ("r", "05_nested_machine_learning.R", "Separate rank-harmonized nested/external workflow"),
    "nested-original": ("r", "05b_original_expression_sensitivity.R", "Primary original-expression nested workflow"),
    "spatial": ("r", "06_all_section_spatial.R", "All-section spatial models from public raw h5ad files"),
    "communication": ("r", "07_communication_sensitivity.R", "Matched B-cell inclusion analysis; full atlas needed"),
    "virtual-ko": ("r", "08_virtual_KO.R", "Three official DSE virtual-network perturbations"),
    "case-only": ("r", "09_case_only_immune_drug.R", "99-case immune and all-198-drug reassessment"),
    "figures": ("r", "10_figures.R", "Six main and four supplementary statistical figures"),
    "minimal-data": ("r", "13_minimal_data.R", "Minimal lineage/condition statistics; full atlas needed"),
    "spatial-fdr": ("r", "14_finalize_FDR.R", "144 planned spatial tests, retaining failed fits"),
    "validate": ("r", "15_validate_processed.R", "Recompute key statistics from the supplied processed inputs"),
}

def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workspace', type=Path, default=ROOT/'work'/'reproduction')
    p.add_argument('--rscript', default='Rscript')
    p.add_argument('--r-library', type=Path)
    p.add_argument('--atlas', type=Path)
    p.add_argument('--bulk-inputs', type=Path)
    sub = p.add_subparsers(dest='command', required=True)
    sub.add_parser('list', help='List analysis stages without running them')
    sub.add_parser('check', help='Check source integrity and release structure without scientific dependencies')
    prep = sub.add_parser('prepare', help='Verify and extract the companion S1_Data.zip into an ignored workspace')
    prep.add_argument('--archive', type=Path, required=True)
    run = sub.add_parser('run', help='Run exactly one selected stage')
    run.add_argument('stage', choices=STAGES)
    run.add_argument('--dry-run', action='store_true')
    return p

def main(argv=None):
    a = parser().parse_args(argv)
    if a.command == 'list':
        for name, (_, script, description) in STAGES.items():
            print(f'{name:18s} code/{script}\n  {description}')
        return 0
    if a.command == 'check':
        return subprocess.call([sys.executable, str(ROOT/'tools/check_release.py')])
    if a.command == 'prepare':
        return subprocess.call([sys.executable, str(ROOT/'tools/prepare_data.py'),
                                '--archive', str(a.archive), '--workspace', str(a.workspace)])
    workspace = a.workspace.resolve()
    if not workspace.is_dir():
        raise SystemExit('Workspace does not exist. Run prepare --archive PATH_TO_S1_Data.zip first.')
    if a.stage in {'single-nucleus', 'communication', 'minimal-data'} and not (a.atlas or os.getenv('ATLAS_RDS')):
        raise SystemExit('This upstream stage requires --atlas pointing to the complete curated Seurat atlas.')
    env = os.environ.copy()
    env.update(REVISION_ROOT=str(workspace), LEGACY_STUDY_ROOT=str(workspace/'legacy_inputs'), MS_CODE_ROOT=str(ROOT))
    if a.atlas:
        env['ATLAS_RDS'] = str(a.atlas.resolve())
    if a.bulk_inputs:
        env['BULK_INPUT_ROOT'] = str(a.bulk_inputs.resolve())
    if a.r_library:
        env['R_LIBS_USER'] = str(a.r_library.resolve())
    lang, name, _ = STAGES[a.stage]
    exe = sys.executable if lang == 'python' else shutil.which(a.rscript)
    if not exe:
        raise SystemExit('Rscript was not found. Supply --rscript or add R to PATH.')
    command = [exe, str(ROOT/'code'/name)]
    if a.dry_run:
        print(json.dumps({'stage':a.stage,'command':command,'workspace':str(workspace)}, indent=2))
        return 0
    return subprocess.call(command, cwd=workspace, env=env)

if __name__ == '__main__':
    raise SystemExit(main())
