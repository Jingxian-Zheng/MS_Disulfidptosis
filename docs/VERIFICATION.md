# Local verification — 2026-10-27

| Check | Result |
|---|---|
| Original source-file inventory | 82 records; original SHA-256 hashes unchanged |
| Public scientific/source records | 76; duplicate versions map to maintained scripts |
| Separate local manuscript-production tools | 6; excluded from public source ZIP |
| R source parsing | All 47 files passed, including historical and vendored source |
| Python source parsing | All 21 files passed |
| Data package preparation | Archive SHA-256 and all 144 manifest entries verified |
| Processed-data statistical validation | All 17 checks passed |
| Figure rebuilding | All 6 main and 4 supplementary plots exported as PDF/SVG/PNG |
| Illustrator finishing source | JavaScript syntax checked without invoking Adobe |
| Offline interface tests | Archive traversal handling and real launcher stages passed |

Detailed statistical checks: `results/reference/processed_validation.csv`.
Source identities and coverage: `provenance/source_inventory.csv`.
The release content manifest can be checked by `python run.py check`.

The source-only GitHub workflow runs offline integrity checks, interface tests and
R parsing. It does not download patient datasets or run expensive upstream analyses.
Local validation used the recorded R environment and existing processed data. It
does not assert a clean-machine installation, all raw-data reruns, or remote upload.
