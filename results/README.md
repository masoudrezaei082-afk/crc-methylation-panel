# Results of the reported run

These are the outputs of the run reported in the manuscript. The discovery stages ran on 2026-10-24; the stages from `network` onward were rerun on 2026-11-06 (`logs/run_2026-11-06_from_network.log`).

| Path | Content |
|---|---|
| `metrics.json` | Every number reported in the manuscript, by stage. |
| `SUMMARY.md` | Readable summary of the main tables. |
| `tables/` | One CSV per table; the supplementary workbook is built from these. |
| `figures/` | Figures at 300 dpi. |
| `logs/` | Per-stage logs (`stage_*.json`), run logs and `pip_freeze.txt`. |
| `config_used.yaml`, `environment.json` | Configuration and package versions of the run. |
| `provenance_downloads.json` | URL and checksum of every downloaded input. |
| `paired_sensitivity/` | Paired tumour–normal sensitivity analysis (`sensitivity/paired_wilcoxon.py`); `NOTE.md` explains each file. |
| `MD5SUMS` | Checksums of every file in this folder. |

Panel identifiers follow the code (`blood_compatible_11`, `screened_11`, `network_prioritized_11`); the main README maps them to the names used in the manuscript.
