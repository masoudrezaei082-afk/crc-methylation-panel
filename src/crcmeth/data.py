"""Configuration, downloads with provenance, TCGA barcodes and matrix readers."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

TUMOR_CODE = 1          # TCGA sample type 01 = primary solid tumour
NORMAL_CODE = 11        # TCGA sample type 11 = solid tissue normal


# --------------------------------------------------------------------------- config
def load_config(path: str | Path, workdir: str | None = None,
                resultsdir: str | None = None) -> dict:
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    if workdir:
        cfg["workdir"] = workdir
    if resultsdir:
        cfg["resultsdir"] = resultsdir
    for key in ("workdir", "resultsdir"):
        Path(cfg[key]).mkdir(parents=True, exist_ok=True)
    return cfg


def rng(cfg: dict, offset: int = 0) -> np.random.Generator:
    return np.random.default_rng(cfg["seed"] + offset)


# --------------------------------------------------------------------------- downloads
def _provenance_path(cfg: dict) -> Path:
    return Path(cfg["resultsdir"]) / "provenance_downloads.json"


def record_provenance(cfg: dict, url: str, path: Path) -> dict:
    md5 = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            md5.update(block)
    entry = {"url": url, "file": path.name, "bytes": path.stat().st_size,
             "md5": md5.hexdigest(),
             "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    prov_file = _provenance_path(cfg)
    prov = json.loads(prov_file.read_text()) if prov_file.exists() else {}
    prov[path.name] = entry
    prov_file.write_text(json.dumps(prov, indent=1, sort_keys=True))
    return entry


def download(cfg: dict, url: str, name: str | None = None, retries: int = 5) -> Path:
    """Download once into workdir/raw, resuming partial files; record md5."""
    import requests

    raw = Path(cfg["workdir"]) / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    dest = raw / (name or url.rsplit("/", 1)[-1])
    if dest.exists():
        return dest
    part = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(retries):
        try:
            have = part.stat().st_size if part.exists() else 0
            headers = {"Range": f"bytes={have}-"} if have else {}
            with requests.get(url, stream=True, timeout=120, headers=headers) as r:
                if r.status_code == 416:          # already complete
                    break
                r.raise_for_status()
                mode = "ab" if r.status_code == 206 else "wb"
                with open(part, mode) as fh:
                    for block in r.iter_content(1 << 22):
                        fh.write(block)
            break
        except Exception as exc:  # network hiccup: back off and resume
            if attempt == retries - 1:
                raise
            print(f"  retry {attempt + 1} for {url}: {exc}")
            time.sleep(10 * (attempt + 1))
    part.rename(dest)
    record_provenance(cfg, url, dest)
    return dest


def gdc_url(cfg: dict, project: str, kind: str) -> str:
    return f"{cfg['urls']['gdc_hub']}/TCGA-{project}.{kind}.tsv.gz"


# --------------------------------------------------------------------------- barcodes
def parse_barcode(sample: str) -> tuple[str, int, str]:
    """TCGA-AA-3517-01A -> (patient 'TCGA-AA-3517', type 1, vial 'A')."""
    parts = sample.split("-")
    code = parts[3]
    return "-".join(parts[:3]), int(code[:2]), code[2:3] or "A"


def patient_of(sample: str) -> str:
    """TCGA-AA-3517-01A -> 'TCGA-AA-3517': the grouping key of every TCGA split."""
    return parse_barcode(sample)[0]


def select_samples(columns, code: int) -> list[str]:
    """Samples of one TCGA sample type, one per patient (first vial)."""
    best: dict[str, tuple[str, str]] = {}
    for col in columns:
        try:
            patient, stype, vial = parse_barcode(col)
        except (IndexError, ValueError):
            continue
        if stype != code:
            continue
        if patient not in best or (vial, col) < best[patient]:
            best[patient] = (vial, col)
    return sorted(v[1] for v in best.values())


def header_columns(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n").split("\t")


# --------------------------------------------------------------------------- readers
def iter_matrix(path: Path, usecols: list[str], chunk_rows: int = 40000):
    """Yield float32 chunks (probes x samples) of a Xena tsv.gz matrix."""
    first = header_columns(path)[0]
    reader = pd.read_csv(path, sep="\t", usecols=[first] + list(usecols),
                         index_col=0, chunksize=chunk_rows,
                         dtype={c: np.float32 for c in usecols}, engine="c")
    for chunk in reader:
        chunk.index.name = "probe"
        yield chunk[list(usecols)]


def read_matrix(path: Path, usecols: list[str], chunk_rows: int = 40000) -> pd.DataFrame:
    return pd.concat(list(iter_matrix(path, usecols, chunk_rows)))


def parse_series_matrix(path: Path, with_table: bool = True):
    """GEO series matrix -> (sample metadata DataFrame, data table or None)."""
    meta: dict[str, list] = {}
    char_rows: list[list[str]] = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            key, *vals = line.rstrip("\n").split("\t")
            vals = [v.strip('"') for v in vals]
            if key == "!Sample_characteristics_ch1":
                char_rows.append(vals)
            else:
                meta.setdefault(key[8:], vals)
    md = pd.DataFrame(meta)
    for row in char_rows:                        # "key: value" cells, keys vary per row
        for i, cell in enumerate(row):
            if ":" in cell:
                k, v = cell.split(":", 1)
                md.loc[i, k.strip()] = v.strip()
    md = md.set_index("geo_accession")
    table = None
    if with_table:
        table = pd.read_csv(path, sep="\t", comment="!", index_col=0, engine="c")
        table = table.astype(np.float32)
        table.index.name = "probe"
        if table.shape[0] == 0:
            table = None
    return md, table


def series_table(path: Path, probes, chunk_rows: int = 50000) -> pd.DataFrame:
    """Rows `probes` of a series-matrix data table, streamed (whole-blood sets do not fit as float64)."""
    probes = set(probes)
    out = []
    reader = pd.read_csv(path, sep="\t", comment="!", index_col=0, chunksize=chunk_rows, engine="c")
    for chunk in reader:
        out.append(chunk.loc[chunk.index.isin(probes)].astype(np.float32))
    table = pd.concat(out)
    table.index.name = "probe"
    return table


def save_frame(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    out = df.copy()
    out.columns = [str(c) for c in out.columns]
    out.to_parquet(path)


def load_frame(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path)


def write_json(obj, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, sort_keys=True, default=_json_default))


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.ndarray, set, tuple)):
        return list(o)
    raise TypeError(type(o))


def env_versions() -> dict:
    import importlib
    import platform
    out = {"python": platform.python_version()}
    for m in ("numpy", "pandas", "scipy", "statsmodels", "sklearn",
              "pydeseq2", "gseapy", "pyarrow", "matplotlib"):
        try:
            out[m] = importlib.import_module(m).__version__
        except Exception:
            out[m] = None
    out["cpu_count"] = os.cpu_count()
    return out
