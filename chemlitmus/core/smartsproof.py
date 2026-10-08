"""Static containment proofs for SMARTS patterns — no reference molecules needed.

``smartsaudit`` finds patterns that are *empirically* subsumed on a reference library. This
module *proves* containment: ``P ⊆ Q`` means every molecule matched by ``P`` is matched by
``Q``, for every molecule, by construction.

How
---
1. **Atom expressions** are read from RDKit's query tree (``QueryAtom.DescribeQuery``) and
   evaluated over a finite **universe of abstract atom states** — element × aromaticity ×
   formal charge × total H count × explicit degree × implicit H count × ring count × ring-bond
   count — built so that it contains every state any primitive in the patterns under study can
   distinguish, plus one unreferenced "other" value per dimension. Over that universe,
   "``A`` implies ``B``" is decided exactly by ``bits(A) & ~bits(B) == 0``.
2. **Bond expressions** are evaluated the same way over bond order × ring membership.
3. **Containment.** ``P ⊆ Q`` is proven by an injective monomorphism ``f`` from the atoms of ``Q``
   into the atoms of ``P`` such that every ``P`` atom expression implies the ``Q`` expression it is
   mapped from, and every ``Q`` bond maps onto a ``P`` bond whose expression implies it. This is the
   Chandra–Merlin homomorphism witness under RDKit's injective, non-induced match semantics. The
   witness is returned so the claim can be checked.

What it can and cannot say
--------------------------
* **Sound, not complete.** A witness proves containment. *No witness* proves nothing: the pattern
  pair is **undecided**, not refuted. Counts of proven containment are lower bounds.
* **Relative to the universe.** Isotopes (``[2H]``), chirality, total valence (``v``), ring size
  (``r``), recursive SMARTS (``$(...)``) and component-level grouping are outside it; a pattern
  that uses them is reported **not analysable**, with the reason.
* Default RDKit matching ignores chirality and bond direction, and so does the prover.
"""

from __future__ import annotations

import itertools
import re
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from pydantic import BaseModel, Field

try:
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    _RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover
    _RDKIT_AVAILABLE = False


class Unsupported(Exception):
    """A primitive outside the atom/bond universe."""


# ------------------------------------------------------------------------------------ query trees

@dataclass
class Node:
    op: str                      # 'and' | 'or' | 'leaf'
    name: str = ""
    value: int = 0
    negate: bool = False
    children: List["Node"] = field(default_factory=list)


def parse_describe(text: str) -> Node:
    """Parse the indented tree printed by ``DescribeQuery()``."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    stack: List[Tuple[int, Node]] = []
    root: Optional[Node] = None
    for ln in lines:
        indent = len(ln) - len(ln.lstrip())
        tok = ln.strip()
        head = tok.split()[0]
        if head in ("AtomAnd", "BondAnd"):
            node = Node("and")
        elif head in ("AtomOr", "BondOr"):
            node = Node("or")
        else:
            neg = ("!=" in tok) or (" not in " in tok)
            m = re.search(r"\s(-?\d+)\s*(?:=|!=)\s*val", tok)
            val = int(m.group(1)) if m else 0
            node = Node("leaf", name=head, value=val, negate=neg)
        while stack and stack[-1][0] >= indent:
            stack.pop()
        if stack:
            stack[-1][1].children.append(node)
        else:
            root = node
        stack.append((indent, node))
    assert root is not None
    return root


def _leaves(node: Node) -> Iterable[Node]:
    if node.op == "leaf":
        yield node
    for c in node.children:
        yield from _leaves(c)


UNSUPPORTED_ATOM = {
    "RecursiveStructure": "recursive SMARTS $(...)",
    "AtomIsotope": "isotope",
    "AtomTotalValence": "total valence v<n>",
    "AtomMinRingSize": "ring size r<n>",
    "AtomHasChiralTag": "chirality",
    "AtomNumRadicalElectrons": "radical count",
    "AtomHybridization": "hybridisation ^<n>",
    "AtomUnsaturated": "unsaturation",
    "AtomNumHeteroatomNeighbors": "heteroatom-neighbour count",
    "AtomNumAliphaticHeteroatomNeighbors": "heteroatom-neighbour count",
    "AtomHasHeteroatomNeighbors": "heteroatom-neighbour test",
    "AtomNonHydrogenDegree": "non-hydrogen degree",
}


# ------------------------------------------------------------------------------------ universe

AROMATIC_BOND = 12
_BASE_Z = [1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 53]


class Universe:
    """Every abstract atom state the primitives under study can tell apart (plus 'other' values)."""

    def __init__(self, leaves: Iterable[Node]):
        ref: Dict[str, set] = {k: set() for k in ("Z", "q", "H", "D", "h", "X", "R", "x")}
        for lf in leaves:
            n, v = lf.name, lf.value
            if n == "AtomType":
                ref["Z"].add(v % 1000)
            elif n == "AtomAtomicNum":
                ref["Z"].add(v)
            elif n == "AtomFormalCharge":
                ref["q"].add(v)
            elif n == "AtomHCount":
                ref["H"].add(v)
            elif n == "AtomImplicitHCount":
                ref["h"].add(v)
            elif n == "AtomExplicitDegree":
                ref["D"].add(v)
            elif n == "AtomTotalDegree":
                ref["X"].add(v)
            elif n == "AtomInNRings":
                if v >= 0:
                    ref["R"].add(v)
            elif n == "AtomRingBondCount":
                ref["x"].add(v)

        def dom(vals, base, sentinel):
            return sorted(set(vals) | set(base) | {sentinel})

        xmax = max(ref["X"] | {4})
        Z = dom(ref["Z"], _BASE_Z, 0)                               # 0 = any other element
        q = dom(ref["q"], [-1, 0, 1], 99)
        H = dom(ref["H"], range(0, 4), 9)
        D = dom(ref["D"], range(0, max(5, xmax + 1)), 9)
        h = dom(ref["h"], range(0, 4), 9)
        R = dom(ref["R"], [0, 1, 2], 9)
        x = dom(ref["x"], [0, 2, 3], 9)
        grid = np.array(list(itertools.product(Z, (0, 1), q, H, D, h, R, x)), dtype=np.int16)
        gZ, ga, gq, gH, gD, gh, gR, gx = (grid[:, i] for i in range(8))
        ok = (gx <= gD)                                              # ring bonds are bonds
        ok &= np.where(gR == 0, gx == 0, gx >= 2)                    # in a ring <=> >= 2 ring bonds
        ok &= (ga == 0) | (gR >= 1)                                  # aromatic atoms are in rings
        ok &= gH <= gh + gD                                          # total H = implicit + explicit-H neighbours
        self.Z, self.arom, self.q, self.H, self.D, self.h, self.R, self.x = (c[ok] for c in (gZ, ga, gq, gH, gD, gh, gR, gx))
        self.X = self.D + self.h
        self.size = int(ok.sum())
        self.domains = {"Z": Z, "q": q, "H": H, "D": D, "h": h, "R": R, "x": x}
        # bond universe: (order, in_ring)
        self.b_order = np.array([1, 2, 3, AROMATIC_BOND, 1, 2, 3], dtype=np.int16)
        self.b_ring = np.array([1, 1, 1, 1, 0, 0, 0], dtype=np.int16)

    # -- evaluation --------------------------------------------------------------------------
    def atom_bits(self, node: Node) -> np.ndarray:
        if node.op == "and":
            out = np.ones(self.size, bool)
            for c in node.children:
                out &= self.atom_bits(c)
            return out
        if node.op == "or":
            out = np.zeros(self.size, bool)
            for c in node.children:
                out |= self.atom_bits(c)
            return out
        n, v = node.name, node.value
        if n == "AtomNull":
            bits = np.ones(self.size, bool)
        elif n == "AtomType":
            bits = (self.Z == v % 1000) & (self.arom == (1 if v >= 1000 else 0))
        elif n == "AtomAtomicNum":
            bits = self.Z == v
        elif n == "AtomIsAromatic":
            bits = self.arom == 1
        elif n == "AtomIsAliphatic":
            bits = self.arom == 0
        elif n == "AtomFormalCharge":
            bits = self.q == v
        elif n == "AtomHCount":
            bits = self.H == v
        elif n == "AtomImplicitHCount":
            bits = self.h == v
        elif n == "AtomHasImplicitH":
            bits = self.h > 0
        elif n == "AtomExplicitDegree":
            bits = self.D == v
        elif n == "AtomTotalDegree":
            bits = self.X == v
        elif n == "AtomInNRings":
            bits = (self.R >= 1) if v < 0 else (self.R == v)
        elif n == "AtomInRing":
            bits = self.R >= 1
        elif n == "AtomRingBondCount":
            bits = self.x == v
        else:
            raise Unsupported(UNSUPPORTED_ATOM.get(n, n))
        return ~bits if node.negate else bits

    def bond_bits(self, node: Node) -> np.ndarray:
        if node.op == "and":
            out = np.ones(7, bool)
            for c in node.children:
                out &= self.bond_bits(c)
            return out
        if node.op == "or":
            out = np.zeros(7, bool)
            for c in node.children:
                out |= self.bond_bits(c)
            return out
        n, v = node.name, node.value
        if n == "BondNull":
            bits = np.ones(7, bool)
        elif n == "BondOrder":
            bits = self.b_order == v
        elif n == "SingleOrAromaticBond":
            bits = (self.b_order == 1) | (self.b_order == AROMATIC_BOND)
        elif n == "BondInRing":
            bits = self.b_ring == 1
        else:
            raise Unsupported(n)
        return ~bits if node.negate else bits


# ------------------------------------------------------------------------------------ compiled patterns

@dataclass
class Compiled:
    smarts: str
    name: Optional[str]
    atoms: List[int]                      # index into universe-wide unique atom expression table
    bonds: Dict[Tuple[int, int], int]     # (i, j) sorted -> unique bond expression index
    adjacency: Dict[int, List[int]]
    unanalysable: Optional[str] = None
    unsatisfiable_atom: Optional[int] = None


class Prover:
    """Holds the universe and the implication tables for a set of patterns."""

    def __init__(self, smarts_list: Sequence[str], names: Optional[Sequence[Optional[str]]] = None, max_container_atoms: int = 12, step_budget: int = 20000):
        if not _RDKIT_AVAILABLE:  # pragma: no cover
            raise ImportError("RDKit is required for smartsproof")
        self.max_container_atoms = max_container_atoms
        self.step_budget = step_budget
        names = list(names) if names is not None else [None] * len(smarts_list)
        mols, trees = [], []
        all_leaves: List[Node] = []
        self._atom_key: Dict[str, int] = {}
        self._bond_key: Dict[str, int] = {}
        atom_trees: List[Node] = []
        bond_trees: List[Node] = []
        for sm in smarts_list:
            q = Chem.MolFromSmarts(sm)
            mols.append(q)
            if q is None:
                trees.append(None)
                continue
            at, bt = [], []
            for a in q.GetAtoms():
                txt = a.DescribeQuery()
                if txt not in self._atom_key:
                    self._atom_key[txt] = len(atom_trees)
                    t = parse_describe(txt); atom_trees.append(t); all_leaves.extend(_leaves(t))
                at.append(self._atom_key[txt])
            for b in q.GetBonds():
                txt = b.DescribeQuery()
                if txt not in self._bond_key:
                    self._bond_key[txt] = len(bond_trees)
                    bond_trees.append(parse_describe(txt))
                bt.append(self._bond_key[txt])
            trees.append((at, bt))
        self.universe = Universe(all_leaves)
        # evaluate every unique expression once
        self.atom_bits: List[Optional[np.ndarray]] = []
        self.atom_unsupported: List[Optional[str]] = []
        for t in atom_trees:
            try:
                self.atom_bits.append(self.universe.atom_bits(t)); self.atom_unsupported.append(None)
            except Unsupported as exc:
                self.atom_bits.append(None); self.atom_unsupported.append(str(exc))
        self.bond_bits: List[Optional[np.ndarray]] = []
        self.bond_unsupported: List[Optional[str]] = []
        for t in bond_trees:
            try:
                self.bond_bits.append(self.universe.bond_bits(t)); self.bond_unsupported.append(None)
            except Unsupported as exc:
                self.bond_bits.append(None); self.bond_unsupported.append(str(exc))
        # implication tables: implies[i, j] == True  <=>  expression i implies expression j
        na, nb = len(atom_trees), len(bond_trees)
        self.atom_implies = np.zeros((na, na), bool)
        packed = [np.packbits(b) if b is not None else None for b in self.atom_bits]
        for i in range(na):
            if packed[i] is None:
                continue
            for j in range(na):
                if packed[j] is not None:
                    self.atom_implies[i, j] = not np.any(packed[i] & ~packed[j])
        self.bond_implies = np.zeros((nb, nb), bool)
        for i in range(nb):
            if self.bond_bits[i] is None:
                continue
            for j in range(nb):
                if self.bond_bits[j] is not None:
                    self.bond_implies[i, j] = not np.any(self.bond_bits[i] & ~self.bond_bits[j])
        # compile patterns
        self.patterns: List[Compiled] = []
        for sm, nm, q, tr in zip(smarts_list, names, mols, trees):
            if q is None:
                self.patterns.append(Compiled(sm, nm, [], {}, {}, unanalysable="does not parse")); continue
            at, bt = tr
            reasons = sorted({self.atom_unsupported[k] for k in at if self.atom_unsupported[k]} | {self.bond_unsupported[k] for k in bt if self.bond_unsupported[k]})
            if "." in sm:
                reasons.append("disconnected pattern")
            bonds: Dict[Tuple[int, int], int] = {}
            adj: Dict[int, List[int]] = {i: [] for i in range(q.GetNumAtoms())}
            for b, k in zip(q.GetBonds(), bt):
                i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                bonds[(min(i, j), max(i, j))] = k
                adj[i].append(j); adj[j].append(i)
            c = Compiled(sm, nm, at, bonds, adj, unanalysable=("; ".join(reasons) or None))
            if c.unanalysable is None:
                for idx, k in enumerate(at):
                    if not self.atom_bits[k].any():
                        c.unsatisfiable_atom = idx; break
            self.patterns.append(c)

    # -- containment ------------------------------------------------------------------------
    def witness(self, p: int, q: int) -> Tuple[Optional[Dict[int, int]], str]:
        """Try to prove patterns[p] ⊆ patterns[q]. Returns (mapping q_atom -> p_atom, status).

        status: 'proven' | 'no witness' | 'budget exceeded' | 'not analysable: <reason>'
        """
        P, Q = self.patterns[p], self.patterns[q]
        for c in (P, Q):
            if c.unanalysable:
                return None, f"not analysable: {c.unanalysable}"
        if len(Q.atoms) > len(P.atoms) or len(Q.bonds) > len(P.bonds):
            return None, "no witness"
        if len(Q.atoms) > self.max_container_atoms:
            return None, "budget exceeded"
        # candidates: P atom a may host Q atom b iff expr(P_a) implies expr(Q_b)
        cands = [[a for a in range(len(P.atoms)) if self.atom_implies[P.atoms[a], Q.atoms[b]]] for b in range(len(Q.atoms))]
        if any(not c for c in cands):
            return None, "no witness"
        order = sorted(range(len(Q.atoms)), key=lambda b: (len(cands[b]), -len(Q.adjacency[b])))
        # BFS-ish: after the first, prefer atoms adjacent to already placed ones
        placed_order: List[int] = []
        remaining = set(order)
        while remaining:
            nxt = None
            for b in order:
                if b in remaining and any(n in placed_order for n in Q.adjacency[b]):
                    nxt = b; break
            if nxt is None:
                nxt = next(b for b in order if b in remaining)
            placed_order.append(nxt); remaining.discard(nxt)
        mapping: Dict[int, int] = {}
        used: set = set()
        steps = [0]

        def ok_bonds(b: int, a: int) -> bool:
            for nb in Q.adjacency[b]:
                if nb in mapping:
                    pa = mapping[nb]
                    pk = P.bonds.get((min(a, pa), max(a, pa)))
                    if pk is None:
                        return False
                    qk = Q.bonds[(min(b, nb), max(b, nb))]
                    if not self.bond_implies[pk, qk]:
                        return False
            return True

        def rec(i: int) -> Optional[bool]:
            if i == len(placed_order):
                return True
            b = placed_order[i]
            for a in cands[b]:
                steps[0] += 1
                if steps[0] > self.step_budget:
                    return None
                if a in used or not ok_bonds(b, a):
                    continue
                mapping[b] = a; used.add(a)
                r = rec(i + 1)
                if r:
                    return True
                if r is None:
                    return None
                del mapping[b]; used.discard(a)
            return False

        r = rec(0)
        if r:
            return dict(mapping), "proven"
        return None, ("budget exceeded" if r is None else "no witness")


# ------------------------------------------------------------------------------------ public API

class ProofResult(BaseModel):
    relation: str = Field(description="'subsumes' (B contains A: every match of A matches B) or 'equivalent'")
    a: str
    b: str
    status: str = Field(description="proven | no witness | budget exceeded | not analysable: <reason>")
    proven: bool
    witness: Optional[Dict[int, int]] = Field(None, description="For 'subsumes': atom index of B -> atom index of A.")
    witness_reverse: Optional[Dict[int, int]] = None
    status_reverse: Optional[str] = None
    universe_size: int = 0


def subsumes(a: str, b: str) -> ProofResult:
    """Prove that every molecule matched by ``a`` is matched by ``b`` (``a ⊆ b``; ``b`` is the broader pattern)."""
    pr = Prover([a, b])
    w, st = pr.witness(0, 1)
    return ProofResult(relation="subsumes", a=a, b=b, status=st, proven=st == "proven", witness=w, universe_size=pr.universe.size)


def equivalent(a: str, b: str) -> ProofResult:
    """Prove that ``a`` and ``b`` match exactly the same molecules (containment both ways)."""
    pr = Prover([a, b])
    w1, s1 = pr.witness(0, 1)
    w2, s2 = pr.witness(1, 0)
    proven = s1 == "proven" and s2 == "proven"
    if proven:
        st = "proven"
    elif s1.startswith("not analysable"):
        st = s1
    else:
        st = s1 if s1 != "proven" else s2
    return ProofResult(relation="equivalent", a=a, b=b, status=st, proven=proven, witness=w1, witness_reverse=w2, status_reverse=s2, universe_size=pr.universe.size)


class SatisfiabilityResult(BaseModel):
    smarts: str
    status: str = Field(description="satisfiable at atom and bond level | unsatisfiable | not analysable: <reason>")
    satisfiable: Optional[bool]
    atom_index: Optional[int] = Field(None, description="The atom whose expression matches no atom state, when unsatisfiable.")
    atom_expression: Optional[str] = None
    universe_size: int = 0


def satisfiable(smarts: str) -> SatisfiabilityResult:
    """Can any atom state satisfy every atom expression, and any bond state every bond expression?

    This is a *local* check: it catches ``[C;N]``, ``[c;!a]``, ``[R0;x2]`` and the like. A pattern
    that passes may still be impossible to realise for global reasons (valence, ring geometry).
    """
    pr = Prover([smarts])
    c = pr.patterns[0]
    if c.unanalysable:
        return SatisfiabilityResult(smarts=smarts, status=f"not analysable: {c.unanalysable}", satisfiable=None, universe_size=pr.universe.size)
    if c.unsatisfiable_atom is not None:
        q = Chem.MolFromSmarts(smarts)
        return SatisfiabilityResult(smarts=smarts, status="unsatisfiable", satisfiable=False, atom_index=c.unsatisfiable_atom,
                                    atom_expression=Chem.MolFragmentToSmarts(q, [c.unsatisfiable_atom]) if q else None, universe_size=pr.universe.size)
    for k in set(c.bonds.values()):
        if not pr.bond_bits[k].any():
            return SatisfiabilityResult(smarts=smarts, status="unsatisfiable", satisfiable=False, atom_expression="bond expression matches no bond state", universe_size=pr.universe.size)
    return SatisfiabilityResult(smarts=smarts, status="satisfiable at atom and bond level", satisfiable=True, universe_size=pr.universe.size)


class PatternProof(BaseModel):
    index: int
    name: Optional[str] = None
    smarts: str
    status: str = Field(description="proven redundant | proven equivalent | no witness | not analysable | unsatisfiable")
    reason: Optional[str] = None
    proven_subsumed_by: List[str] = Field(default_factory=list, description="Names/SMARTS of catalogue patterns proven to contain this one (strictly or equivalently).")
    proven_equivalent_to: List[str] = Field(default_factory=list)
    witness: Optional[Dict[int, int]] = Field(None, description="For the first container: its atom index -> this pattern's atom index.")
    n_undecided_pairs: int = Field(0, description="Candidate containers (<= this pattern's size) for which no witness was found: undecided, not refuted.")
    n_budget_exceeded: int = 0


class CatalogueProof(BaseModel):
    n_patterns: int
    universe_size: int
    universe_domains: Dict[str, List[int]]
    n_analysable: int
    n_not_analysable: int
    not_analysable_reasons: Dict[str, int] = Field(default_factory=dict)
    n_unsatisfiable: int = 0
    n_proven_redundant: int = Field(0, description="Patterns with at least one proven container in the catalogue (lower bound).")
    n_proven_equivalent: int = 0
    n_no_witness: int = Field(0, description="Analysable patterns with no proven container: undecided, not refuted.")
    n_pairs_tested: int = 0
    n_pairs_budget_exceeded: int = 0
    patterns: List[PatternProof] = Field(default_factory=list)


def prove_catalogue(patterns: Sequence[Tuple[str, Optional[str]]], max_container_atoms: int = 12, step_budget: int = 20000,
                    progress_callback=None) -> CatalogueProof:
    """For every pattern, look for catalogue patterns proven to contain it."""
    smarts = [p[0] for p in patterns]; names = [p[1] for p in patterns]
    pr = Prover(smarts, names, max_container_atoms=max_container_atoms, step_budget=step_budget)
    label = [n or s for n, s in zip(names, smarts)]
    reasons = Counter()
    out: List[PatternProof] = []
    n_pairs = n_budget = 0
    analysable = [i for i, c in enumerate(pr.patterns) if not c.unanalysable]
    for i, c in enumerate(pr.patterns):
        pp = PatternProof(index=i, name=c.name, smarts=c.smarts, status="no witness")
        if c.unanalysable:
            pp.status = "not analysable"; pp.reason = c.unanalysable
            for r in c.unanalysable.split("; "):
                reasons[r] += 1
            out.append(pp); continue
        if c.unsatisfiable_atom is not None:
            pp.status = "unsatisfiable"; pp.reason = f"atom {c.unsatisfiable_atom} matches no atom state"
            out.append(pp); continue
        for j in analysable:
            if j == i or pr.patterns[j].unsatisfiable_atom is not None:
                continue
            Pj = pr.patterns[j]
            if len(Pj.atoms) > len(c.atoms) or len(Pj.bonds) > len(c.bonds):
                continue
            n_pairs += 1
            w, st = pr.witness(i, j)
            if st == "proven":
                pp.proven_subsumed_by.append(label[j])
                if pp.witness is None:
                    pp.witness = w
                if len(Pj.atoms) == len(c.atoms) and len(Pj.bonds) == len(c.bonds):
                    w2, st2 = pr.witness(j, i)
                    if st2 == "proven":
                        pp.proven_equivalent_to.append(label[j])
            elif st == "budget exceeded":
                pp.n_budget_exceeded += 1; n_budget += 1; pp.n_undecided_pairs += 1
            else:
                pp.n_undecided_pairs += 1
        if pp.proven_equivalent_to and len(pp.proven_equivalent_to) == len(pp.proven_subsumed_by):
            pp.status = "proven equivalent"
        elif pp.proven_subsumed_by:
            pp.status = "proven redundant"
        out.append(pp)
        if progress_callback:
            progress_callback(i + 1, len(pr.patterns))
    res = CatalogueProof(
        n_patterns=len(pr.patterns), universe_size=pr.universe.size, universe_domains=pr.universe.domains,
        n_analysable=len(analysable), n_not_analysable=len(pr.patterns) - len(analysable), not_analysable_reasons=dict(reasons),
        n_unsatisfiable=sum(1 for p in out if p.status == "unsatisfiable"),
        n_proven_redundant=sum(1 for p in out if p.status in ("proven redundant", "proven equivalent")),
        n_proven_equivalent=sum(1 for p in out if p.proven_equivalent_to),
        n_no_witness=sum(1 for p in out if p.status == "no witness"),
        n_pairs_tested=n_pairs, n_pairs_budget_exceeded=n_budget, patterns=out,
    )
    return res


__all__ = ["subsumes", "equivalent", "satisfiable", "prove_catalogue", "Prover", "Universe", "ProofResult",
           "SatisfiabilityResult", "PatternProof", "CatalogueProof", "UNSUPPORTED_ATOM"]
