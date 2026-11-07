# Paired tumour–normal sensitivity analysis (Design B)

All 38 TCGA colon normals come from patients who also have a tumour in the data. The main analysis compares 295 tumours with 38 normals using the Mann–Whitney U test, which treats the two groups as independent. `sensitivity/paired_wilcoxon.py` repeats the comparison as a Wilcoxon signed-rank test on the 38 within-patient differences, on the same matrix and with the same thresholds (FDR < 0.05, |Δβ| > 0.2).

## Results

- **Panel CpGs.** All 21 distinct CpGs of the blood-compatible and network-prioritized panels remain hypermethylated. Their paired FDR is at most 2.24 × 10⁻⁵, and the paired Δβ differs from the main Δβ by at most 0.07.
- **Genome-wide.** There are 38,947 DMCs in the main analysis and 41,103 in the paired analysis. No main DMC changes sign. Of the 2,734 that drop out, 2,730 fall just below |Δβ| 0.2. Across all probes, the Pearson r of the two Δβ is 0.989.
- **Selection funnel.** The funnel was rerun under six variants (`paired_funnel.csv`):

| Variant | DMC call | Ten-cancer screen |
|---|---|---|
| A | main | main (rank-sum, all samples) |
| B | signed-rank FDR, main Δβ | main |
| C | signed-rank FDR, paired Δβ | main |
| D | as C | signed-rank on each cancer's pairs, BH over the screened probes |
| E | as D | genome-wide BH upper bound |
| F | as C | rank-sum on the same paired samples, pairing ignored |

  - B gives the same panels as A.
  - C keeps 7 of the 11 blood-compatible CpGs; the change comes from using only the 38 paired tumours.
  - D–F keep far fewer colon-specific CpGs, almost entirely because of ESCA. The 16 ESCA tumours that have a matched normal are much more methylated at the candidate CpGs than the 169 without one (median β 0.53 vs 0.24; `screen_composition.csv`). F, which ignores the pairing, gives the same panel as D, so the drop reflects which tumours have a normal, not the pairing.

## Files

| File | Content |
|---|---|
| `S_paired_all_cpgs_designB.csv.gz` | All 397,605 probes: main and paired Δβ, p, FDR and calls, plus the rank-sum on the same 38 pairs. |
| `paired_panels.csv` | The 22 panel CpGs: main and paired statistics, per-cancer paired screen and status in each variant. |
| `paired_funnel.csv` | Counts, exclusions per cancer and panels for variants A–F. |
| `paired_summary.json` | Genome-wide concordance and the pairs available per screened cancer. |
| `screen_composition.csv` | Screened cancers: tumours with and without a matched normal. |
| `S_paired_screen_candidates.csv.gz` | Per-cancer paired screen statistics of every candidate CpG. |
| `screen_variant_*.csv.gz` | The ten-cancer screen of each variant. |
| `run.log`, `pip_freeze.txt` | Log and environment of the run (Python 3.13.16). |

## Notes

- The signed-rank p is exact (no ties or zero differences, n ≤ 50); otherwise it uses the normal approximation with tie and continuity correction.
- The pan-cancer stage keeps only the probes with COAD Δβ > 0.15, plus *SEPT9* (30,763 probes). The paired screen tests these rows. BH over them is valid because they were selected on COAD data only.
- CESC (3 pairs) and SARC (4 pairs) cannot reach p < 0.05 with a signed-rank test, so they keep their main-screen call.
- Candidates without blood data were treated once as not compatible and once as compatible. The overlap with the main panels is the same either way.
- RNA-seq normals are also paired with tumours. The differential expression model has no patient term, and this was not tested.
