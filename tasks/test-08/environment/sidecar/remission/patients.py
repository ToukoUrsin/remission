"""Hidden side of Remission: one virtual patient per seed.

Never shipped to the agent's container. A patient is a set of parameter values
for the public model in model.py: drug potencies, one healthy cell population,
a dominant tumour clone and zero to two minor clones that are already resistant
to something.

The distributions are calibrated against public cell-line data (GDSC2) and
checked against published trial response rates; see science/ and docs/MODEL.md.
"""
from __future__ import annotations

import hashlib

import numpy as np

from . import model as M

# Typical target potencies (micromolar), their spread (sd of ln K) across patients, and the off-target
# toxicity to normal cells: public knowledge, kept in model.py.
K_MED, K_SD = M.K_TYPICAL, M.K_LN_SD
TOX_MED, TOX_SD = M.TOX_TYPICAL, M.TOX_LN_SD

ARCHETYPES = ("BRAF_CRC", "BRAF_MEL", "KRAS", "EGFR")
ARCH_P = (0.28, 0.24, 0.22, 0.26)
TISSUE = {"BRAF_CRC": "colorectal", "BRAF_MEL": "skin (melanoma)", "KRAS": "colorectal", "EGFR": "lung"}
MINORS = {
    "BRAF_CRC": ("ras", "alt", "bclx", "mcl1", "pten"),
    "BRAF_MEL": ("ras", "mek1", "alt", "bclx", "mcl1", "pten"),
    "KRAS": ("alt", "bclx", "mcl1", "pten"),
    "EGFR": ("t790m", "ras", "alt", "bclx", "mcl1", "pten"),
}

# What sequencing calls each alteration.
ALTERATION_LABEL = {
    "braf": "BRAF V600E",
    "egfr": "EGFR L858R",
    "ras": "KRAS G12C",
    "pik3ca": "PIK3CA H1047R",
    "pten": "PTEN loss",
    "mek1": "MAP2K1 (MEK1) K57N",
    "alt": "MET amplification",
    "t790m": "EGFR T790M",
    "mcl1": "MCL1 amplification",
    "bclx": "BCL2L1 (BCL-xL) amplification",
}


def seed_from(salt: str, index: int) -> int:
    """Stable 31-bit seed from a salt and an index (so a suite is one string)."""
    h = hashlib.sha256(f"{salt}:{index}".encode()).digest()
    return int.from_bytes(h[:4], "big") & 0x7FFFFFFF


def _lognorm(rng, med, sd):
    return float(med * np.exp(rng.normal(0.0, sd)))


death_drive0 = M.death_drive0


def make_patient(seed: int) -> dict:
    rng = np.random.default_rng([int(seed), 7919])
    arch = str(rng.choice(ARCHETYPES, p=ARCH_P))
    K = np.array([_lognorm(rng, K_MED[d], K_SD) for d in M.DRUGS])

    # Healthy cells are not primed for apoptosis: their buffer is 1.6 to 2.6 times their untreated death drive.
    healthy = M.make_clone(
        name="healthy", healthy=True, egfr=1.0, alt=0.3,
        fb=rng.uniform(1.0, 5.0), kp=rng.uniform(0.24, 0.50), kdmax=rng.uniform(0.6, 0.9),
        bfrac=rng.uniform(0.15, 0.60), mfrac=rng.uniform(0.15, 0.60),
        sE=rng.uniform(0.45, 0.65), sA=rng.uniform(0.45, 0.65))
    healthy["bcl"] = death_drive0(healthy) * rng.uniform(1.6, 2.6)
    healthy["tox"] = np.array([_lognorm(rng, TOX_MED, TOX_SD) for _ in M.DRUGS])

    alts = []
    if arch == "BRAF_CRC":
        major = dict(braf=1.0, ras=0.0, egfr=_lognorm(rng, 1.0, 0.45), alt=_lognorm(rng, 0.12, 0.7),
                     fb=rng.uniform(2.5, 6.0), sE=rng.uniform(0.80, 0.95), sA=rng.uniform(0.25, 0.50))
        alts.append("braf")
    elif arch == "BRAF_MEL":
        major = dict(braf=1.0, ras=0.0, egfr=_lognorm(rng, 0.08, 0.6), alt=_lognorm(rng, 0.10, 0.8),
                     fb=rng.uniform(1.0, 4.0), sE=rng.uniform(0.80, 0.95), sA=rng.uniform(0.25, 0.50))
        alts.append("braf")
    elif arch == "KRAS":
        major = dict(braf=0.0, ras=rng.uniform(0.7, 1.1), egfr=_lognorm(rng, 0.8, 0.45), alt=_lognorm(rng, 0.15, 0.7),
                     fb=rng.uniform(2.0, 6.0), sE=rng.uniform(0.75, 0.92), sA=rng.uniform(0.35, 0.60))
        alts.append("ras")
    else:  # EGFR-mutant, RAS/RAF wild type
        major = dict(braf=0.0, ras=0.0, egfr=rng.uniform(2.5, 5.0), alt=_lognorm(rng, 0.10, 0.6),
                     fb=rng.uniform(0.5, 2.0), sE=rng.uniform(0.70, 0.90), sA=rng.uniform(0.45, 0.70))
        alts.append("egfr")
    major.update(kp=rng.uniform(0.50, 0.90), kdmax=rng.uniform(0.9, 1.5),
                 bfrac=rng.uniform(0.10, 0.40), mfrac=rng.uniform(0.10, 0.40),
                 pi3k=0.0, pten=1.0, kshift=np.ones(M.ND))
    if arch == "EGFR":
        major["kshift"][M.iEGFR] = 0.1      # the EGFR inhibitor binds the mutant receptor ten times better
    u = rng.random()
    if u < 0.20:
        major["pi3k"] = rng.uniform(0.4, 0.8)
        major["sA"] = min(major["sA"] + 0.2, 0.9)
        alts.append("pik3ca")
    elif u < 0.30:
        major["pten"] = rng.uniform(0.20, 0.35)
        major["sA"] = min(major["sA"] + 0.2, 0.9)
        alts.append("pten")
    # Tumour cells are primed: their buffer barely covers their untreated death drive.
    # BRAF-mutant colorectal tumours are the least primed (their response rate to RAF/EGFR blockade is the
    # lowest of the four tumour types: 20 to 27 % in BEACON CRC).
    major["bcl"] = death_drive0(M.make_clone(**major)) * rng.uniform(1.30, 2.00) * (1.3 if arch == "BRAF_CRC" else 1.0)
    tumour = [M.make_clone(name="major", **major)]
    tumour[0]["alterations"] = list(alts)

    # minor clones that are already resistant to something
    options = [k for k in MINORS[arch] if k not in alts]
    n_minor = int(rng.choice([0, 1, 2], p=[0.08, 0.52, 0.40]))
    kinds = [str(k) for k in rng.choice(options, size=n_minor, replace=False)] if n_minor else []
    fracs = []
    for kind in kinds:
        c = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in major.items()}
        if kind == "ras":
            c["ras"] = rng.uniform(0.7, 1.0)
        elif kind == "mek1":
            c["kshift"][M.iMEK] = rng.uniform(15.0, 40.0)
        elif kind == "t790m":
            c["kshift"][M.iEGFR] = rng.uniform(3.0, 10.0)
        elif kind == "pten":
            c["pten"] = rng.uniform(0.20, 0.35)
            c["sA"] = min(c["sA"] + 0.2, 0.9)
        elif kind == "alt":
            c["alt"] = rng.uniform(0.6, 1.5)
        elif kind == "mcl1":
            c["bcl"] = c["bcl"] * rng.uniform(1.6, 2.4)
            c["bfrac"], c["mfrac"] = 0.05, rng.uniform(0.80, 0.90)
        elif kind == "bclx":
            c["bcl"] = c["bcl"] * rng.uniform(1.6, 2.4)
            c["bfrac"], c["mfrac"] = rng.uniform(0.80, 0.90), 0.05
        clone = M.make_clone(name=f"minor-{kind}", **c)
        clone["alterations"] = list(alts) + [kind]
        tumour.append(clone)
        fracs.append(float(np.exp(rng.uniform(np.log(0.004), np.log(0.06)))))
    frac = np.array([1.0 - sum(fracs)] + fracs)

    return {"seed": int(seed), "archetype": arch, "tissue": TISSUE[arch], "K": K,
            "healthy": healthy, "tumour": tumour, "frac": frac}


def chart(patient: dict) -> dict:
    """The bulk sequencing report in the patient's chart (clones above 10 % only)."""
    major = patient["tumour"][0]
    found = [ALTERATION_LABEL[a] for a in major["alterations"]]
    wt = [g for g, key in (("BRAF", "braf"), ("KRAS", "ras"), ("EGFR", "egfr"), ("PIK3CA", "pik3ca"), ("PTEN", "pten"))
          if key not in major["alterations"]]
    return {"tissue": patient["tissue"],
            "bulk_sequencing": {"alterations": found, "wild_type": wt,
                                "note": "bulk panel; clones under about 10 % of the tumour are not reported"}}
