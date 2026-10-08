"""Audit a set of SMARTS patterns (structural alerts, substructure filters) for correctness,
redundancy, breadth and reproducibility.

Structural-alert catalogues (PAINS, Brenk, Glaxo, in-house filters) are redistributed as plain
SMARTS text and applied with undocumented molecule-preparation settings. This module evaluates a
pattern set against a reference population of real molecules and reports:

* **compile**      - patterns RDKit cannot parse, and patterns that need explicit hydrogens.
* **breadth**      - patterns that match an unexpectedly large share of the reference set.
* **dead**         - patterns that match nothing, triaged into *rare combination* (every query
                     atom is individually realisable) versus *never-matching atom* (at least one
                     query atom matches no atom in any reference molecule).
* **redundancy**   - exact duplicate SMARTS, patterns that are indistinguishable on the reference
                     set, and patterns strictly subsumed by another pattern in the same set.
* **sensitivity**  - how hit counts and per-compound verdicts change when the reference molecules
                     are prepared differently (implicit H, explicit H, kekulized).

Everything is offline. Matching uses RDKit's ``SubstructLibrary`` and runs multithreaded in C++.
"""

from __future__ import annotations

import gzip
import importlib.resources
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from pydantic import BaseModel, Field

try:
    import numpy as np
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdSubstructLibrary as _ssl

    RDLogger.DisableLog("rdApp.*")
    _RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover
    _RDKIT_AVAILABLE = False

# --------------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------------

PREPARATIONS: List[str] = ["implicit-h", "explicit-h", "kekule"]
"""Molecule preparations recognised by the audit and by ``--prep`` on other commands."""

AUDIT_CHECKS: List[str] = ["compile", "breadth", "dead", "redundancy", "sensitivity"]
"""Audit checks, in the order they run."""

DEFAULT_BREADTH_THRESHOLD = 0.10
"""A pattern matching more than this fraction of the reference set is flagged over-broad."""

DEFAULT_DEAD_SAMPLE = 2500
"""Reference molecules sampled when triaging dead patterns atom by atom."""


# --------------------------------------------------------------------------------------------
# Result models
# --------------------------------------------------------------------------------------------


class PatternAudit(BaseModel):
    """Audit record for a single SMARTS pattern."""

    index: int
    smarts: str
    name: Optional[str] = None
    rule_set: Optional[str] = None

    # compile
    parses: bool = False
    parse_error: Optional[str] = None
    n_query_atoms: int = 0
    requires_explicit_h: bool = Field(
        False, description="Pattern contains a hydrogen query atom, so it can only match molecules prepared with explicit hydrogens."
    )
    has_recursive_smarts: bool = False

    # breadth (default preparation)
    n_hits: Optional[int] = None
    hit_fraction: Optional[float] = None
    over_broad: bool = False

    # dead
    never_fires: bool = False
    dead_verdict: Optional[str] = Field(
        None,
        description="'rare combination' | 'never-matching atom' | 'fires only with <preparation>' | None when the pattern fires.",
    )
    never_matching_atoms: List[int] = Field(default_factory=list, description="Query-atom indices that match no reference atom.")

    # redundancy
    duplicate_of: Optional[int] = Field(None, description="Index of an earlier pattern with the identical SMARTS string.")
    equivalent_to: List[int] = Field(default_factory=list, description="Other patterns with an identical hit set on the reference library.")
    subsumed_by: Optional[int] = Field(None, description="Index of a pattern whose hit set strictly contains this one's.")

    # sensitivity
    hits_by_preparation: Dict[str, int] = Field(default_factory=dict)
    preparation_sensitive: bool = Field(False, description="Hit count differs between preparations.")

    @property
    def flags(self) -> List[str]:
        out: List[str] = []
        if not self.parses:
            out.append("unparseable")
        if self.requires_explicit_h:
            out.append("needs-explicit-h")
        if self.over_broad:
            out.append("over-broad")
        if self.never_fires:
            out.append("dead:" + (self.dead_verdict or "unknown").replace(" ", "-"))
        if self.duplicate_of is not None:
            out.append("duplicate")
        if self.equivalent_to:
            out.append("equivalent")
        if self.subsumed_by is not None:
            out.append("subsumed")
        if self.preparation_sensitive:
            out.append("prep-sensitive")
        return out

    @property
    def clean(self) -> bool:
        return not self.flags


class SensitivitySummary(BaseModel):
    """Catalogue-level reproducibility summary across molecule preparations."""

    preparations: List[str]
    n_molecules: int
    compounds_flagged: Dict[str, int] = Field(default_factory=dict, description="Molecules hit by >=1 pattern, per preparation.")
    total_hits: Dict[str, int] = Field(default_factory=dict)
    patterns_firing: Dict[str, int] = Field(default_factory=dict)
    verdict_flips: Dict[str, int] = Field(
        default_factory=dict,
        description="Molecules whose pass/fail verdict differs between the default preparation and each other one.",
    )
    n_sensitive_patterns: int = 0


class SmartsAuditResult(BaseModel):
    """Complete result of auditing a SMARTS pattern set."""

    n_patterns: int
    n_molecules: int
    library_source: str
    checks_run: List[str]
    breadth_threshold: float
    patterns: List[PatternAudit]
    sensitivity: Optional[SensitivitySummary] = None
    elapsed_seconds: float = 0.0
    error: Optional[str] = None

    # ----- convenience aggregates -----
    @property
    def n_unparseable(self) -> int:
        return sum(1 for p in self.patterns if not p.parses)

    @property
    def n_over_broad(self) -> int:
        return sum(1 for p in self.patterns if p.over_broad)

    @property
    def n_dead(self) -> int:
        return sum(1 for p in self.patterns if p.never_fires)

    @property
    def n_dead_never_matching_atom(self) -> int:
        return sum(1 for p in self.patterns if p.dead_verdict == "never-matching atom")

    @property
    def n_duplicates(self) -> int:
        return sum(1 for p in self.patterns if p.duplicate_of is not None)

    @property
    def n_equivalent(self) -> int:
        return sum(1 for p in self.patterns if p.equivalent_to)

    @property
    def n_subsumed(self) -> int:
        return sum(1 for p in self.patterns if p.subsumed_by is not None)

    @property
    def n_needs_explicit_h(self) -> int:
        return sum(1 for p in self.patterns if p.requires_explicit_h)

    @property
    def n_clean(self) -> int:
        return sum(1 for p in self.patterns if p.clean)

    def to_rows(self) -> List[Dict[str, object]]:
        """Flat per-pattern rows suitable for CSV export."""
        rows = []
        for p in self.patterns:
            row: Dict[str, object] = {
                "index": p.index,
                "name": p.name or "",
                "rule_set": p.rule_set or "",
                "smarts": p.smarts,
                "parses": p.parses,
                "parse_error": p.parse_error or "",
                "n_query_atoms": p.n_query_atoms,
                "requires_explicit_h": p.requires_explicit_h,
                "has_recursive_smarts": p.has_recursive_smarts,
                "n_hits": p.n_hits if p.n_hits is not None else "",
                "hit_fraction": round(p.hit_fraction, 5) if p.hit_fraction is not None else "",
                "over_broad": p.over_broad,
                "never_fires": p.never_fires,
                "dead_verdict": p.dead_verdict or "",
                "never_matching_atoms": ";".join(map(str, p.never_matching_atoms)),
                "duplicate_of": p.duplicate_of if p.duplicate_of is not None else "",
                "equivalent_to": ";".join(map(str, p.equivalent_to)),
                "subsumed_by": p.subsumed_by if p.subsumed_by is not None else "",
                "preparation_sensitive": p.preparation_sensitive,
                "flags": ";".join(p.flags),
            }
            for prep in PREPARATIONS:
                row[f"hits_{prep}"] = p.hits_by_preparation.get(prep, "")
            rows.append(row)
        return rows


class AtomExplanation(BaseModel):
    atom_index: int
    query: str = Field(description="Compact description of the atom primitive.")
    n_matching_molecules: int = Field(description="Reference molecules containing at least one atom satisfying this primitive alone.")
    is_hydrogen: bool = False


class SmartsExplanation(BaseModel):
    """Decomposition of a single SMARTS pattern against the reference library."""

    smarts: str
    parses: bool
    parse_error: Optional[str] = None
    normalized_smarts: Optional[str] = None
    n_query_atoms: int = 0
    n_query_bonds: int = 0
    requires_explicit_h: bool = False
    has_recursive_smarts: bool = False
    atoms: List[AtomExplanation] = Field(default_factory=list)
    hits_by_preparation: Dict[str, int] = Field(default_factory=dict)
    n_molecules: int = 0
    example_matches: List[str] = Field(default_factory=list)
    never_matching_atoms: List[int] = Field(default_factory=list)

    @property
    def verdict(self) -> str:
        if not self.parses:
            return "unparseable"
        default_hits = self.hits_by_preparation.get("implicit-h", 0)
        if default_hits > 0:
            return "fires"
        if self.never_matching_atoms:
            return "dead: never-matching atom"
        for prep, n in self.hits_by_preparation.items():
            if prep != "implicit-h" and n > 0:
                return f"fires only with {prep}"
        if self.requires_explicit_h:
            return "dead: rare combination (even with explicit hydrogens)"
        return "dead: rare combination"


# --------------------------------------------------------------------------------------------
# Pattern-file parsing
# --------------------------------------------------------------------------------------------

_SMARTS_COLUMNS = ("smarts", "pattern", "query", "smarts_pattern")
_NAME_COLUMNS = ("name", "description", "rule_id", "id", "alert", "label")
_SET_COLUMNS = ("rule_set_name", "rule_set", "set", "catalog", "catalogue", "source")


def load_patterns(path: Path | str) -> List[Tuple[str, Optional[str], Optional[str]]]:
    """Read SMARTS patterns from a file.

    Supported layouts:

    * ``.csv`` / ``.tsv`` / ``.xlsx`` with a column named ``smarts`` (or ``pattern``/``query``);
      optional ``name``/``description`` and ``rule_set`` columns are carried through.
    * ``.smarts`` / ``.txt`` / ``.smi``: one pattern per line, optionally followed by whitespace
      and a name. Lines starting with ``#`` are ignored.

    Returns:
        List of ``(smarts, name, rule_set)`` tuples.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Pattern file not found: {path}")
    ext = path.suffix.lower()
    out: List[Tuple[str, Optional[str], Optional[str]]] = []

    if ext in (".csv", ".tsv", ".xlsx", ".xls"):
        import pandas as pd

        if ext == ".tsv":
            df = pd.read_csv(path, sep="\t")
        elif ext == ".csv":
            df = pd.read_csv(path)
        else:
            df = pd.read_excel(path)
        cols = {str(c).strip().lower(): c for c in df.columns}
        smarts_col = next((cols[c] for c in _SMARTS_COLUMNS if c in cols), None)
        if smarts_col is None:
            raise ValueError(
                f"No SMARTS column found in {path.name}. Expected one of {list(_SMARTS_COLUMNS)}; got {list(df.columns)}"
            )
        name_col = next((cols[c] for c in _NAME_COLUMNS if c in cols), None)
        set_col = next((cols[c] for c in _SET_COLUMNS if c in cols), None)
        for _, r in df.iterrows():
            s = r[smarts_col]
            if not isinstance(s, str) or not s.strip():
                continue
            name = str(r[name_col]) if name_col is not None and not pd.isna(r[name_col]) else None
            rs = str(r[set_col]) if set_col is not None and not pd.isna(r[set_col]) else None
            out.append((s.strip(), name, rs))
        return out

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            out.append((parts[0], parts[1].strip() if len(parts) > 1 else None, None))
    return out


# --------------------------------------------------------------------------------------------
# Reference library
# --------------------------------------------------------------------------------------------


def load_reference_library(path: Optional[Path | str] = None, max_molecules: Optional[int] = None) -> Tuple[List["Chem.Mol"], str]:
    """Load reference molecules.

    Args:
        path: A SMILES-bearing file (.smi/.csv/.tsv/.xlsx/.sdf). ``None`` loads the bundled
              ChEMBL-derived library (see ``chemlitmus/data/README.md`` for licence terms).
        max_molecules: Keep only the first N molecules.

    Returns:
        ``(molecules, source_description)``
    """
    if not _RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is required. Install with: pip install rdkit")

    if path is None:
        res = importlib.resources.files("chemlitmus.data").joinpath("reference_library.smi.gz")
        with res.open("rb") as raw, gzip.open(raw, "rt") as fh:
            smiles = [ln.split()[0] for ln in fh if ln.strip() and not ln.startswith("#")]
        source = "bundled (ChEMBL_37 subset)"
    else:
        from chemlitmus.utils.parsers import parse_compounds_file

        smiles = parse_compounds_file(Path(path))
        source = str(path)

    if max_molecules is not None:
        smiles = smiles[:max_molecules]
    mols = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        if m is not None:
            mols.append(m)
    if not mols:
        raise ValueError(f"No valid molecules loaded from {source}")
    return mols, source


def prepare_molecule(mol: "Chem.Mol", preparation: str) -> "Chem.Mol":
    """Return a copy of ``mol`` in the requested preparation state."""
    if preparation == "implicit-h":
        return Chem.Mol(mol)
    if preparation == "explicit-h":
        return Chem.AddHs(mol)
    if preparation == "kekule":
        k = Chem.Mol(mol)
        try:
            Chem.Kekulize(k, clearAromaticFlags=True)
        except Exception:
            return Chem.Mol(mol)
        return k
    raise ValueError(f"Unknown preparation {preparation!r}. Valid: {PREPARATIONS}")


def _build_library(mols: Sequence["Chem.Mol"], preparation: str) -> "_ssl.SubstructLibrary":
    # Pattern fingerprints pre-screen candidates before the full match. The fingerprint is
    # computed on the prepared molecule, so screening is consistent for every preparation.
    lib = _ssl.SubstructLibrary(_ssl.MolHolder(), _ssl.PatternHolder())
    for m in mols:
        lib.AddMol(prepare_molecule(m, preparation))
    return lib


def _match_vector(lib: "_ssl.SubstructLibrary", query: "Chem.Mol", n: int) -> "np.ndarray":
    v = np.zeros(n, dtype=bool)
    try:
        idx = lib.GetMatches(query, numThreads=-1, maxResults=-1)
    except Exception:
        return v
    if len(idx):
        v[np.fromiter(idx, dtype=np.int64)] = True
    return v


# --------------------------------------------------------------------------------------------
# Pattern introspection
# --------------------------------------------------------------------------------------------


def _requires_explicit_h(query: "Chem.Mol") -> bool:
    """True when some query atom can only be satisfied by a hydrogen atom.

    Negated primitives (``[!#1]``) and alternatives (``[#1,#6]``) do not count: they are
    satisfiable without hydrogens being present as graph atoms.
    """
    for a in query.GetAtoms():
        if a.GetAtomicNum() == 1:
            return True
        if not a.HasQuery():
            continue
        desc = a.DescribeQuery()
        head = desc.strip().splitlines()[0] if desc.strip() else ""
        if head.startswith("AtomOr"):
            continue
        if "AtomAtomicNum 1 = val" in desc:  # the negated form reads "AtomAtomicNum 1 != val"
            return True
    return False


def _has_recursive(query: "Chem.Mol") -> bool:
    return any(a.HasQuery() and "RecursiveStructure" in a.DescribeQuery() for a in query.GetAtoms())


def _atom_subquery(query: "Chem.Mol", keep: int) -> "Chem.Mol":
    rw = Chem.RWMol(query)
    for j in sorted((a.GetIdx() for a in query.GetAtoms() if a.GetIdx() != keep), reverse=True):
        rw.RemoveAtom(j)
    return rw.GetMol()


def _compact_atom_query(query: "Chem.Mol", idx: int) -> str:
    s = Chem.MolToSmarts(_atom_subquery(query, idx))
    return s if s else query.GetAtomWithIdx(idx).DescribeQuery().strip().replace("\n", " ")


# --------------------------------------------------------------------------------------------
# Main audit
# --------------------------------------------------------------------------------------------


def audit_smarts(
    patterns: Sequence[str] | Sequence[Tuple[str, Optional[str], Optional[str]]],
    library: Optional[Sequence["Chem.Mol"]] = None,
    library_source: str = "user-supplied",
    checks: Optional[Iterable[str]] = None,
    breadth_threshold: float = DEFAULT_BREADTH_THRESHOLD,
    preparations: Optional[Sequence[str]] = None,
    dead_sample: int = DEFAULT_DEAD_SAMPLE,
) -> SmartsAuditResult:
    """Audit a set of SMARTS patterns against a reference molecule library.

    Args:
        patterns: SMARTS strings, or ``(smarts, name, rule_set)`` tuples as returned by
                  :func:`load_patterns`.
        library: Reference molecules. ``None`` loads the bundled library.
        library_source: Label recorded in the result.
        checks: Subset of :data:`AUDIT_CHECKS` to run; ``None`` runs all. ``compile`` always runs.
        breadth_threshold: Hit fraction above which a pattern is flagged over-broad.
        preparations: Preparations for the sensitivity check; default all of :data:`PREPARATIONS`.
                      The first entry is the default preparation used for breadth/dead/redundancy.
        dead_sample: Molecules sampled for atom-level triage of dead patterns.

    Returns:
        :class:`SmartsAuditResult`.
    """
    t0 = time.time()
    if not _RDKIT_AVAILABLE:
        return SmartsAuditResult(
            n_patterns=0, n_molecules=0, library_source=library_source, checks_run=[],
            breadth_threshold=breadth_threshold, patterns=[], error="RDKit is not installed.",
        )

    requested = set(AUDIT_CHECKS) if checks is None else {c.strip().lower() for c in checks}
    bad = requested - set(AUDIT_CHECKS)
    if bad:
        raise ValueError(f"Unknown check(s): {sorted(bad)}. Valid: {AUDIT_CHECKS}")
    requested.add("compile")
    preps = list(preparations) if preparations else list(PREPARATIONS)
    for p in preps:
        if p not in PREPARATIONS:
            raise ValueError(f"Unknown preparation {p!r}. Valid: {PREPARATIONS}")
    default_prep = preps[0]

    if library is None:
        library, library_source = load_reference_library()
    mols = list(library)
    n_mol = len(mols)

    # ---- normalise input -------------------------------------------------------------------
    triples: List[Tuple[str, Optional[str], Optional[str]]] = []
    for item in patterns:
        if isinstance(item, str):
            triples.append((item, None, None))
        else:
            s, name, rs = (list(item) + [None, None])[:3]
            triples.append((str(s), name, rs))

    # ---- compile ---------------------------------------------------------------------------
    records: List[PatternAudit] = []
    queries: List[Optional["Chem.Mol"]] = []
    for i, (s, name, rs) in enumerate(triples):
        rec = PatternAudit(index=i, smarts=s, name=name, rule_set=rs)
        q = Chem.MolFromSmarts(s)
        if q is None or q.GetNumAtoms() == 0:
            rec.parses = False
            rec.parse_error = "RDKit could not parse this SMARTS" if q is None else "SMARTS has no atoms"
        else:
            rec.parses = True
            rec.n_query_atoms = q.GetNumAtoms()
            rec.requires_explicit_h = _requires_explicit_h(q)
            rec.has_recursive_smarts = _has_recursive(q)
        records.append(rec)
        queries.append(q if rec.parses else None)

    need_matrix = requested & {"breadth", "dead", "redundancy", "sensitivity"}
    matrices: Dict[str, "np.ndarray"] = {}
    if need_matrix and n_mol:
        preps_to_run = preps if "sensitivity" in requested else [default_prep]
        for prep in preps_to_run:
            lib = _build_library(mols, prep)
            M = np.zeros((len(records), n_mol), dtype=bool)
            for i, q in enumerate(queries):
                if q is not None:
                    M[i] = _match_vector(lib, q, n_mol)
            matrices[prep] = M
        M0 = matrices[default_prep]
        hits0 = M0.sum(axis=1)

        for rec, h in zip(records, hits0):
            if rec.parses:
                rec.n_hits = int(h)
                rec.hit_fraction = float(h) / n_mol
                rec.hits_by_preparation[default_prep] = int(h)

        # ---- breadth -----------------------------------------------------------------------
        if "breadth" in requested:
            for rec in records:
                if rec.parses and rec.hit_fraction is not None and rec.hit_fraction > breadth_threshold:
                    rec.over_broad = True

        # ---- dead --------------------------------------------------------------------------
        if "dead" in requested:
            sample_idx = np.linspace(0, n_mol - 1, num=min(dead_sample, n_mol), dtype=int)
            sample_mols = [mols[j] for j in sample_idx]
            sample_libs: Dict[str, "_ssl.SubstructLibrary"] = {}
            for rec, q in zip(records, queries):
                if not rec.parses or rec.n_hits:
                    continue
                rec.never_fires = True
                # dead under the default preparation but alive under another one: that is a
                # preparation requirement, not a dead rule
                alive_in = [p for p, M in matrices.items() if p != default_prep and M[rec.index].any()]
                if alive_in:
                    rec.dead_verdict = f"fires only with {alive_in[0]}"
                    continue
                # Hydrogen-bearing patterns are triaged against explicit-H molecules so the verdict
                # is about the rest of the pattern; their hydrogen dependence is already reported
                # through ``requires_explicit_h``.
                triage_prep = "explicit-h" if rec.requires_explicit_h else default_prep
                if triage_prep not in sample_libs:
                    sample_libs[triage_prep] = _build_library(sample_mols, triage_prep)
                sample_lib = sample_libs[triage_prep]
                bad_atoms = []
                for ai in range(q.GetNumAtoms()):
                    sub = _atom_subquery(q, ai)
                    try:
                        if len(sample_lib.GetMatches(sub, numThreads=-1, maxResults=1)) == 0:
                            bad_atoms.append(ai)
                    except Exception:
                        bad_atoms.append(ai)
                rec.never_matching_atoms = bad_atoms
                rec.dead_verdict = "never-matching atom" if bad_atoms else "rare combination"

        # ---- redundancy --------------------------------------------------------------------
        if "redundancy" in requested:
            seen: Dict[str, int] = {}
            for rec in records:
                if not rec.parses:
                    continue
                if rec.smarts in seen:
                    rec.duplicate_of = seen[rec.smarts]
                else:
                    seen[rec.smarts] = rec.index

            live = [i for i, r in enumerate(records) if r.parses and r.n_hits]
            groups: Dict[bytes, List[int]] = {}
            for i in live:
                groups.setdefault(np.packbits(M0[i]).tobytes(), []).append(i)
            for g in groups.values():
                if len(g) > 1:
                    for i in g:
                        records[i].equivalent_to = [j for j in g if j != i]

            if live:
                packed = np.packbits(M0[live], axis=1)
                counts = hits0[live]
                order = np.argsort(counts, kind="stable")
                for pos, a in enumerate(order):
                    ia = live[a]
                    candidates = order[pos + 1 :]
                    candidates = candidates[counts[candidates] > counts[a]]
                    if candidates.size == 0:
                        continue
                    # A ⊆ B  <=>  (A & ~B) == 0
                    contained = ~np.any(packed[a] & ~packed[candidates], axis=1)
                    if contained.any():
                        ib = live[int(candidates[np.argmax(contained)])]
                        records[ia].subsumed_by = ib

        # ---- sensitivity -------------------------------------------------------------------
        if "sensitivity" in requested and len(preps) > 1:
            summary = SensitivitySummary(preparations=preps, n_molecules=n_mol)
            flagged0 = M0.any(axis=0)
            for prep in preps:
                M = matrices[prep]
                hits = M.sum(axis=1)
                for rec, h in zip(records, hits):
                    if rec.parses:
                        rec.hits_by_preparation[prep] = int(h)
                summary.compounds_flagged[prep] = int(M.any(axis=0).sum())
                summary.total_hits[prep] = int(M.sum())
                summary.patterns_firing[prep] = int((hits > 0).sum())
                if prep != default_prep:
                    summary.verdict_flips[prep] = int((M.any(axis=0) != flagged0).sum())
            for rec in records:
                if rec.parses and len(set(rec.hits_by_preparation.values())) > 1:
                    rec.preparation_sensitive = True
            summary.n_sensitive_patterns = sum(1 for r in records if r.preparation_sensitive)
        else:
            summary = None
    else:
        summary = None

    return SmartsAuditResult(
        n_patterns=len(records),
        n_molecules=n_mol,
        library_source=library_source,
        checks_run=[c for c in AUDIT_CHECKS if c in requested],
        breadth_threshold=breadth_threshold,
        patterns=records,
        sensitivity=summary,
        elapsed_seconds=round(time.time() - t0, 3),
    )


# --------------------------------------------------------------------------------------------
# Single-pattern explanation
# --------------------------------------------------------------------------------------------


def explain_smarts(
    smarts: str,
    library: Optional[Sequence["Chem.Mol"]] = None,
    preparations: Optional[Sequence[str]] = None,
    n_examples: int = 5,
    sample: int = DEFAULT_DEAD_SAMPLE,
) -> SmartsExplanation:
    """Decompose a SMARTS pattern and report how each atom primitive, and the whole pattern,
    behaves against the reference library under each preparation."""
    if not _RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is required. Install with: pip install rdkit")
    q = Chem.MolFromSmarts(smarts)
    if q is None or q.GetNumAtoms() == 0:
        return SmartsExplanation(
            smarts=smarts, parses=False,
            parse_error="RDKit could not parse this SMARTS" if q is None else "SMARTS has no atoms",
        )
    if library is None:
        library, _ = load_reference_library()
    mols = list(library)
    preps = list(preparations) if preparations else list(PREPARATIONS)

    exp = SmartsExplanation(
        smarts=smarts, parses=True,
        normalized_smarts=Chem.MolToSmarts(q),
        n_query_atoms=q.GetNumAtoms(), n_query_bonds=q.GetNumBonds(),
        requires_explicit_h=_requires_explicit_h(q), has_recursive_smarts=_has_recursive(q),
        n_molecules=len(mols),
    )

    sample_idx = np.linspace(0, len(mols) - 1, num=min(sample, len(mols)), dtype=int)
    sample_mols = [mols[j] for j in sample_idx]
    atom_lib = _build_library(sample_mols, "explicit-h" if exp.requires_explicit_h else "implicit-h")
    for a in q.GetAtoms():
        sub = _atom_subquery(q, a.GetIdx())
        try:
            n = len(atom_lib.GetMatches(sub, numThreads=-1, maxResults=-1))
        except Exception:
            n = 0
        exp.atoms.append(AtomExplanation(
            atom_index=a.GetIdx(), query=_compact_atom_query(q, a.GetIdx()),
            n_matching_molecules=int(n), is_hydrogen=(a.GetAtomicNum() == 1),
        ))
        if n == 0:
            exp.never_matching_atoms.append(a.GetIdx())

    for prep in preps:
        lib = _build_library(mols, prep)
        try:
            idx = list(lib.GetMatches(q, numThreads=-1, maxResults=-1))
        except Exception:
            idx = []
        exp.hits_by_preparation[prep] = len(idx)
        if prep == preps[0] or (not exp.example_matches and idx):
            exp.example_matches = [Chem.MolToSmiles(mols[j]) for j in idx[:n_examples]]
    return exp


__all__ = [
    "PREPARATIONS",
    "AUDIT_CHECKS",
    "DEFAULT_BREADTH_THRESHOLD",
    "PatternAudit",
    "SensitivitySummary",
    "SmartsAuditResult",
    "AtomExplanation",
    "SmartsExplanation",
    "load_patterns",
    "load_reference_library",
    "prepare_molecule",
    "audit_smarts",
    "explain_smarts",
]
