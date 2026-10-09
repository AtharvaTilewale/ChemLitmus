#!/usr/bin/env python
"""Benchmark: how much do structural-alert verdicts depend on molecule preparation?

Reproducible recipe for the headline numbers quoted in the README and documentation.

Inputs
------
* Alert catalogue: the ChEMBL structural-alert table as redistributed by ``rd_filters``
  (CC BY-SA 3.0; https://github.com/PatWalters/rd_filters). Downloaded on first run to
  ``benchmarks/data/alert_collection.csv``; the SHA-256 of the file used is printed and written
  into the result so a later run can be compared against it.
* Molecules: the reference library bundled with ChemLitmus (9,272 ChEMBL_37 structures), or any
  ``.smi``/``.csv`` passed with ``--library``.

Every number is reported **per rule set and for the union**, with its denominator, because the
catalogue contains eight published sets and a measurement over the union is not a measurement of
PAINS.

Usage
-----
    python benchmarks/alert_preparation_sensitivity.py --output benchmarks/results.json
    python benchmarks/alert_preparation_sensitivity.py --library my_library.smi --max-molecules 2000
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np

from chemlitmus.core.smartsaudit import PREPARATIONS, _build_library, load_patterns, load_reference_library

ALERTS_URL = "https://raw.githubusercontent.com/PatWalters/rd_filters/master/rd_filters/data/alert_collection.csv"
DATA_DIR = Path(__file__).parent / "data"


def fetch_alerts() -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    p = DATA_DIR / "alert_collection.csv"
    if not p.exists():
        print(f"downloading {ALERTS_URL}", file=sys.stderr)
        urllib.request.urlretrieve(ALERTS_URL, p)
    return p


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alerts", type=Path, default=None, help="Alert CSV (default: fetch the rd_filters copy).")
    ap.add_argument("--library", type=Path, default=None, help="Molecule library (default: the bundled reference set).")
    ap.add_argument("--max-molecules", type=int, default=None)
    ap.add_argument("--output", type=Path, default=Path(__file__).parent / "results.json")
    args = ap.parse_args()

    alerts_path = args.alerts or fetch_alerts()
    patterns = load_patterns(alerts_path)
    mols, lib_source = load_reference_library(args.library, max_molecules=args.max_molecules)
    print(f"{len(patterns)} patterns · {len(mols)} molecules ({lib_source})", file=sys.stderr)

    from rdkit import Chem

    queries, sets = [], []
    n_unparseable = 0
    for smarts, _name, rule_set in patterns:
        q = Chem.MolFromSmarts(smarts)
        if q is None:
            n_unparseable += 1
        queries.append(q)
        sets.append(rule_set or "(unnamed)")

    matrices = {}
    evaluated = {}
    for prep in PREPARATIONS:
        lib = _build_library(mols, prep)
        M = np.zeros((len(queries), len(mols)), dtype=bool)
        for i, q in enumerate(queries):
            if q is not None:
                M[i], err = lib.match(q)
                if err:
                    M[i] = False
        matrices[prep] = M
        evaluated[prep] = lib.evaluated
        print(f"  {prep}: {lib.n_evaluated} molecules evaluated, {lib.n_failed} preparation failures", file=sys.stderr)

    by_set = defaultdict(list)
    for i, s in enumerate(sets):
        by_set[s].append(i)
    default = PREPARATIONS[0]

    result = {
        "alerts_file": str(alerts_path.name), "alerts_sha256": sha256(alerts_path), "n_patterns": len(patterns),
        "n_unparseable_patterns": n_unparseable, "library_source": lib_source, "n_molecules": len(mols),
        "preparations": PREPARATIONS, "default_preparation": default, "rule_sets": {}, "union": {},
        "note": ("Each figure is for the rule set named; the 'union' entry covers every set in the file at once and "
                 "must not be attributed to any single published set."),
    }
    for name, idx in list(by_set.items()) + [("union", list(range(len(queries))))]:
        sub = {}
        flagged = {}
        for prep in PREPARATIONS:
            M = matrices[prep][idx]
            f = M.any(axis=0)
            flagged[prep] = f
            sub[f"flagged_{prep}"] = int(f.sum())
            sub[f"flagged_fraction_{prep}"] = round(float(f.mean()), 4)
            sub[f"patterns_firing_{prep}"] = int((M.sum(axis=1) > 0).sum())
        for prep in PREPARATIONS[1:]:
            both = evaluated[prep] & evaluated[default]
            sub[f"verdict_flips_vs_{default}_{prep}"] = int(((flagged[prep] != flagged[default]) & both).sum())
            sub[f"verdict_flip_fraction_{prep}"] = round(float(((flagged[prep] != flagged[default]) & both).mean()), 4)
        sub["n_patterns"] = len(idx)
        sub["n_molecules"] = len(mols)
        (result["rule_sets"] if name != "union" else result)[name if name != "union" else "union"] = sub

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2))
    u = result["union"]
    print(json.dumps({k: v for k, v in u.items() if "fraction" in k or k.startswith("flagged_")}, indent=2))
    print(f"written to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
