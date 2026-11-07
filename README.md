# crc-methylation-panel

Code and results for the manuscript

> **Genome-Wide Methylation Profiling Uncovers Specific Biomarkers for Colon Cancer Diagnosis**
> Ali Mahmoudabadi, Masoud Rezaei, Asma Askari

The pipeline downloads public data, selects CpG markers for colon adenocarcinoma, and tests them in TCGA, in two GEO validation cohorts, against blood and in plasma cell-free DNA. A single command reproduces every table and figure.

## Data

All inputs are public, and the pipeline downloads them itself:

- TCGA methylation (HumanMethylation450) and RNA-seq (STAR counts) from the GDC hub of [UCSC Xena](https://xenabrowser.net): COAD plus the ten cancers of the specificity screen;
- GEO GSE48684 and GSE131013 (validation), GSE40279 (whole blood), GSE35069 (sorted leukocytes) and GSE122126 (plasma cfDNA);
- platform GPL13534 and the [Zhou et al. probe mask](https://zwdzwd.github.io/InfiniumAnnotation);
- STRING v12.5 and JASPAR 2024 for the network prioritization.

The URL, date and MD5 checksum of the downloaded files are recorded in `results/provenance_downloads.json`.

## Layout

| Path | Content |
|---|---|
| `run_all.py` | Cached end-to-end pipeline. Stages, in order: `annot coad pancan dmc rna integrate specificity network geo classify blood report`. |
| `config.yaml` | Every threshold, seed and parameter used in the manuscript. |
| `src/crcmeth/` | One module per stage. |
| `tests/` | Unit tests on synthetic data (`test_core.py`) and an end-to-end smoke test (`smoke_pipeline.py`), both offline. |
| `notebooks/crc_meth_repro.ipynb` | Google Colab driver for the free CPU runtime. |
| `sensitivity/paired_wilcoxon.py` | Paired tumour–normal sensitivity analysis (see below). |
| `results/` | Outputs of the run reported in the manuscript (see `results/README.md`). |

## Running

The reported run used Python 3.13. The pipeline needs about 12 GB of memory and no GPU. A full run takes several hours, mostly to download and stream the TCGA matrices.

```bash
git clone https://github.com/masoudrezaei082-afk/crc-methylation-panel.git
cd crc-methylation-panel
pip install numpy pandas scipy scikit-learn statsmodels matplotlib pyarrow pyyaml requests openpyxl
pip install -r requirements.txt
python -m pytest -q tests
python run_all.py --workdir WORK --resultsdir RESULTS                   # all stages
python run_all.py --workdir WORK --resultsdir RESULTS --stage blood     # one stage
python run_all.py --workdir WORK --resultsdir RESULTS --from network    # from a stage onward
```

`requirements.txt` lists only the packages that Colab does not ship. The other versions used are in `results/environment.json`, and the full `pip freeze` is in `results/logs/pip_freeze.txt`.

Each stage caches its output in `WORK/cache`, so an interrupted run resumes where it stopped. `--force` recomputes. On Colab, open `notebooks/crc_meth_repro.ipynb` and run the cells in order.

## Names used in the code

| Code | Manuscript |
|---|---|
| `blood_compatible_11` | blood-compatible panel |
| `screened_11` | screen-derived panel |
| `network_prioritized_11` | network-prioritized panel |
| `SEPT9_promoter` | *SEPT9* promoter probes |
| Design B | primary analysis: COAD tumours vs colon normals |
| Design A | sensitivity analysis: colon normals plus 60 normals from six other tissues |

## Paired tumour–normal sensitivity analysis

In TCGA, every colon normal comes from a patient who also has a tumour in the data. The main analysis uses the Mann–Whitney U test, which treats the two groups as independent. `sensitivity/paired_wilcoxon.py` repeats the Design B comparison as a Wilcoxon signed-rank test on the 38 tumour–normal pairs. It also reruns the selection funnel, including a paired ten-cancer specificity screen. It reads the cache of a full pipeline run:

```bash
python sensitivity/paired_wilcoxon.py WORK results/paired_sensitivity
```

Its outputs are in `results/paired_sensitivity/`.

## License

MIT; see `LICENSE`. If you use this code, please cite the article (see `CITATION.cff`).
