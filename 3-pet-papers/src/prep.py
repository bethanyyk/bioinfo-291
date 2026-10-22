"""Reproduce the authors' preprocessing from Rscript.R and expose the review's inputs."""
import hashlib
import numpy as np
import pandas as pd

import json
import os

# Artifact version ids of the four supplied inputs.  Paths are resolved at run time
# (goold_2016/paths.json is written by the review script via host.artifact_path).
ARTIFACTS = {
    "rawData.csv": "094c1619-8c9c-43a5-be47-6e3c4f25e7ed",
    "finalMIdata_averaged_rounded.csv": "2224cbbb-dab3-409a-bc95-bd8658466fb9",
    "Rscript.R": "465b1417-e091-417d-9006-e8fef5451ece",
    "paper_text": "fc638965-1d60-4359-9592-07b6ca90bcfe",
}
_P = json.load(open(os.path.join(os.path.dirname(__file__), "paths.json")))
RAW_CSV = _P["rawData.csv"]
MI_CSV = _P["finalMIdata_averaged_rounded.csv"]
RSCRIPT = _P["Rscript.R"]
PAPER_TXT = _P["paper_text"]

# Rscript.R line 250 / 285: descriptors dropped before network estimation
DROP_HANDLER_LRT = ["FHO", "OCD", "OCN", "ANX", "QU", "TC", "CC", "TUG"]  # clmm vs clm LRT
DROP_LOW_VARIATION = ["CB"]

# Table 1 of the paper: a = desirable, b = undesirable
DESIRABLE = ["ACT", "ADP", "CUR", "FIT", "FL", "FSH", "PLA", "PS", "PSV",
             "REC", "SLP", "SOC", "TOY", "WIL"]
UNDESIRABLE = ["DA", "FDA", "FoH", "GUS", "GWL", "STR"]
NODES20 = sorted(DESIRABLE + UNDESIRABLE, key=str.upper)

LABELS = {
    "ACT": "active and nimble", "ADP": "adapts to new situations quickly",
    "CUR": "curious", "DA": "dog aggressive", "FDA": "food aggressive",
    "FIT": "physically fit", "FL": "fearless", "FoH": "fear of heights",
    "FSH": "stays focused during searches", "GUS": "gives up searches quickly",
    "GWL": "growls at strangers", "PLA": "playful", "PS": "problem solving",
    "PSV": "persevering", "REC": "recalls", "SLP": "walks on slippery surfaces",
    "SOC": "socially attached to you", "STR": "nervous/tense when startled",
    "TOY": "willing to give you a toy", "WIL": "willing to please",
}


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 16), b""):
            h.update(blk)
    return h.hexdigest()[:16]


def load_all():
    raw = pd.read_csv(RAW_CSV)
    mi = pd.read_csv(MI_CSV)
    rscript = open(RSCRIPT, encoding="utf-8", errors="replace").read()
    paper = open(PAPER_TXT, encoding="utf-8", errors="replace").read()
    return raw, mi, rscript, paper


def rebuild_preprocessing(raw):
    """Follow Rscript.R lines 76-186 on rawData.csv; returns a dict of audit numbers."""
    desc_all = [c for c in raw.columns if c not in ("workType", "handlerType", "breed", "newID")]
    z = raw[desc_all].copy()
    zero_frac = (z == 0).sum() / len(z)                     # NAs are not zeros
    drop_zero10 = sorted(zero_frac[zero_frac > 0.10].index)  # line 86-88
    kept38 = [c for c in desc_all if c not in drop_zero10]
    per_dog_zero = (raw[kept38] == 0).sum(axis=1) / len(kept38)
    worst_dog = int(per_dog_zero.values.argmax())
    rm170 = raw.drop(index=raw.index[169]).reset_index(drop=True)   # line 96
    s = rm170[kept38].replace(0, np.nan)
    na_frac = s.isna().sum() / len(s)
    drop_na5 = sorted(na_frac[na_frac > 0.05].index)                # line 149-151
    kept34 = [c for c in kept38 if c not in drop_na5]
    return {
        "zero_frac": zero_frac, "drop_zero10": drop_zero10, "kept38": kept38,
        "per_dog_zero": per_dog_zero, "worst_dog_idx0": worst_dog,
        "worst_dog_frac": float(per_dog_zero.max()),
        "na_frac": na_frac, "drop_na5": drop_na5, "kept34": kept34,
        "rm170": rm170,
    }


def network_inputs(mi):
    """The two 20-node ordinal matrices the paper's networks were estimated from."""
    patrol = mi.loc[mi.workType == 1, NODES20].to_numpy(int)
    detect = mi.loc[mi.workType == 0, NODES20].to_numpy(int)
    return patrol, detect
