"""450K probe annotation (Illumina GPL13534, hg19) and Zhou et al. probe masks."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

ISLAND_ORDER = ["Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"]
REGION_ORDER = ["TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR"]


def load_gpl13534(path: Path) -> pd.DataFrame:
    """One row per probe: chrom, pos, island relation, gene/group lists."""
    with gzip.open(path, "rt", errors="replace") as fh:
        skip = next(i for i, line in enumerate(fh) if line.startswith("IlmnID"))
    cols = ["IlmnID", "CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group",
            "Relation_to_UCSC_CpG_Island"]
    a = pd.read_csv(path, skiprows=skip, usecols=cols, low_memory=False)
    a = a[a["IlmnID"].astype(str).str.startswith(("cg", "ch"))]
    a = a.rename(columns={"IlmnID": "probe", "CHR": "chrom", "MAPINFO": "pos",
                          "Relation_to_UCSC_CpG_Island": "island"})
    a["island"] = a["island"].fillna("OpenSea")
    a["chrom"] = a["chrom"].astype(str)
    names = a["UCSC_RefGene_Name"].fillna("").str.split(";")
    groups = a["UCSC_RefGene_Group"].fillna("").str.split(";")
    a["gene_group"] = [
        sorted({(g, r) for g, r in zip(n, gr) if g and r}) for n, gr in zip(names, groups)
    ]
    a = a.drop(columns=["UCSC_RefGene_Name", "UCSC_RefGene_Group"])
    return a.set_index("probe")


def genes_in(annot: pd.DataFrame, groups) -> pd.Series:
    """probe -> sorted gene symbols whose annotated region is in `groups`."""
    groups = set(groups)
    return annot["gene_group"].map(lambda gg: sorted({g for g, r in gg if r in groups}))


def regions_of(annot: pd.DataFrame) -> pd.Series:
    return annot["gene_group"].map(lambda gg: sorted({r for _, r in gg}))


def load_zhou_mask(path: Path) -> set[str]:
    """Probes with M_general == TRUE (Zhou et al. 2017 general-purpose mask)."""
    m = pd.read_csv(path, sep="\t", usecols=["Probe_ID", "M_general"], low_memory=False)
    return set(m.loc[m["M_general"].astype(str).str.upper() == "TRUE", "Probe_ID"])


def gene_probemap(path: Path) -> pd.DataFrame:
    """GENCODE v36 (hg38) gene map: Ensembl id (versioned) -> symbol, TSS."""
    g = pd.read_csv(path, sep="\t")
    g.columns = ["id", "gene", "chrom", "start", "end", "strand"]
    g["tss"] = g["start"].where(g["strand"] == "+", g["end"])
    return g.set_index("id")
