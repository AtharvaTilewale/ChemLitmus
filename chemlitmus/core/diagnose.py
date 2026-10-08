"""Deterministic, explainable diagnosis of invalid SMILES.

RDKit tells you a SMILES failed. This module tells you *where* and *why*, in a fixed order of
checks, each one pointing at a character position or an atom:

1. **characters**   - whitespace, non-ASCII, or characters outside the SMILES alphabet
2. **brackets**     - malformed bracket atoms ``[...]`` (unknown element, bad charge/H-count)
3. **parentheses**  - unbalanced branches, with the position of the unmatched symbol
4. **rings**        - ring-closure digits that are never closed, or closed onto the same atom
5. **syntax**       - anything RDKit's parser still rejects, with its reported position
6. **valence**      - atoms whose explicit valence exceeds what the element permits
7. **aromaticity**  - rings that cannot be kekulized (classic: a pyrrole ``n`` missing ``[nH]``)

Where a safe, mechanical repair exists (strip whitespace, close a dangling branch, drop an unclosed
ring digit, add ``[nH]`` to an unkekulizable ring) it is attempted and the result re-validated.
Repairs are suggestions, not chemistry: the output records exactly what was changed.
"""

from __future__ import annotations

import io
import logging
import re
from typing import List, Optional, Tuple

from pydantic import BaseModel, Field

try:
    from rdkit import Chem, RDLogger, rdBase

    RDLogger.DisableLog("rdApp.*")
    _RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover
    _RDKIT_AVAILABLE = False


DIAGNOSTIC_CATEGORIES: List[str] = [
    "characters", "brackets", "parentheses", "rings", "syntax", "valence", "aromaticity",
]

_ORGANIC = {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I", "b", "c", "n", "o", "p", "s", "*"}
_ELEMENTS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te",
    "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
    "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    "b", "c", "n", "o", "p", "s", "se", "as", "te", "si", "*",
}
_BRACKET_RE = re.compile(
    r"^\[(?P<isotope>\d+)?(?P<symbol>[A-Z][a-z]?|[a-z]{1,2}|\*)"
    r"(?P<chiral>@@|@|@TH\d|@AL\d|@SP\d|@TB\d{1,2}|@OH\d{1,2})?"
    r"(?P<hcount>H\d*)?"
    r"(?P<charge>\+\+|--|[+-]\d*)?"
    r"(?P<map>:\d+)?\]$"
)
_ALLOWED_CHARS = set("BCNOPSFIKHbcnopsAlrQu*[]()=#$:/\\.%@+-0123456789")  # superset; checked per token
_RD_POS_RE = re.compile(r"around position (\d+)")


class SmilesProblem(BaseModel):
    """One located problem in a SMILES string."""

    category: str = Field(description="One of DIAGNOSTIC_CATEGORIES.")
    message: str
    position: Optional[int] = Field(None, description="0-based character offset in the input, when known.")
    length: int = Field(1, description="Length of the offending span.")
    atom_index: Optional[int] = Field(None, description="RDKit atom index, for valence/aromaticity problems.")
    atom_indices: List[int] = Field(default_factory=list)
    suggestion: Optional[str] = None


class SmilesDiagnosis(BaseModel):
    """Full diagnosis of one SMILES string."""

    input_smiles: str
    is_valid: bool
    canonical_smiles: Optional[str] = None
    problems: List[SmilesProblem] = Field(default_factory=list)
    repaired_smiles: Optional[str] = Field(None, description="Result of mechanical repair, if one was attempted.")
    repaired_is_valid: Optional[bool] = None
    repairs_applied: List[str] = Field(default_factory=list)

    @property
    def primary_category(self) -> Optional[str]:
        return self.problems[0].category if self.problems else None

    def caret_line(self) -> str:
        """An ASCII marker line aligned under the input, pointing at problem positions."""
        marks = [" "] * len(self.input_smiles)
        for p in self.problems:
            if p.position is not None:
                for i in range(p.position, min(len(marks), p.position + max(1, p.length))):
                    marks[i] = "^"
        return "".join(marks).rstrip()


# --------------------------------------------------------------------------------------------
# Tokeniser
# --------------------------------------------------------------------------------------------


def _tokenize(s: str) -> List[Tuple[str, str, int]]:
    """Split a SMILES into ``(kind, text, position)`` tokens.

    kind ∈ {'atom', 'bracket', 'bond', 'open', 'close', 'ring', 'dot', 'bad'}.
    Does not validate chemistry; only lexes.
    """
    toks: List[Tuple[str, str, int]] = []
    i, n = 0, len(s)
    while i < n:
        ch = s[i]
        if ch == "[":
            j = s.find("]", i)
            if j == -1:
                toks.append(("bad", s[i:], i))
                break
            toks.append(("bracket", s[i : j + 1], i))
            i = j + 1
        elif ch in "()":
            toks.append(("open" if ch == "(" else "close", ch, i)); i += 1
        elif ch in "-=#$:/\\~":
            toks.append(("bond", ch, i)); i += 1
        elif ch == ".":
            toks.append(("dot", ch, i)); i += 1
        elif ch == "%":
            m = re.match(r"%\d\d", s[i:])
            if m:
                toks.append(("ring", m.group(0), i)); i += 3
            else:
                toks.append(("bad", ch, i)); i += 1
        elif ch.isdigit():
            toks.append(("ring", ch, i)); i += 1
        elif s[i : i + 2] in ("Cl", "Br"):
            toks.append(("atom", s[i : i + 2], i)); i += 2
        elif ch in "BCNOPSFIbcnops*":
            toks.append(("atom", ch, i)); i += 1
        else:
            toks.append(("bad", ch, i)); i += 1
    return toks


# --------------------------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------------------------


def _check_characters(s: str) -> List[SmilesProblem]:
    out = []
    for i, ch in enumerate(s):
        if ch.isspace():
            out.append(SmilesProblem(category="characters", position=i, message="Whitespace inside SMILES",
                                     suggestion="Remove whitespace, or split into separate records"))
        elif ord(ch) > 127:
            out.append(SmilesProblem(category="characters", position=i,
                                     message=f"Non-ASCII character {ch!r} (U+{ord(ch):04X})",
                                     suggestion="Replace with the ASCII equivalent (e.g. '-' for an en dash)"))
    return out


def _check_tokens(toks: List[Tuple[str, str, int]]) -> List[SmilesProblem]:
    out = []
    for kind, text, pos in toks:
        if kind == "bad":
            if text.startswith("["):
                out.append(SmilesProblem(category="brackets", position=pos, length=len(text),
                                         message="Bracket atom is never closed", suggestion="Add the missing ']'"))
            else:
                out.append(SmilesProblem(category="characters", position=pos,
                                         message=f"Character {text!r} is not valid in SMILES"))
        elif kind == "bracket":
            m = _BRACKET_RE.match(text)
            if not m:
                out.append(SmilesProblem(category="brackets", position=pos, length=len(text),
                                         message=f"Malformed bracket atom {text}",
                                         suggestion="Expected [isotope?symbol chirality?H-count?charge?:map?]"))
            elif m.group("symbol") not in _ELEMENTS:
                out.append(SmilesProblem(category="brackets", position=pos, length=len(text),
                                         message=f"Unknown element symbol {m.group('symbol')!r} in {text}"))
    return out


def _check_parentheses(toks: List[Tuple[str, str, int]]) -> List[SmilesProblem]:
    out, stack = [], []
    for kind, text, pos in toks:
        if kind == "open":
            stack.append(pos)
        elif kind == "close":
            if not stack:
                out.append(SmilesProblem(category="parentheses", position=pos,
                                         message="')' closes a branch that was never opened",
                                         suggestion="Remove this ')' or add a matching '('"))
            else:
                stack.pop()
    for pos in stack:
        out.append(SmilesProblem(category="parentheses", position=pos,
                                 message="'(' opens a branch that is never closed",
                                 suggestion="Add a matching ')'"))
    return out


def _check_rings(toks: List[Tuple[str, str, int]]) -> List[SmilesProblem]:
    out = []
    open_rings: dict = {}      # digit -> (position, atom ordinal)
    atom_ordinal = -1
    for kind, text, pos in toks:
        if kind in ("atom", "bracket"):
            atom_ordinal += 1
        elif kind == "ring":
            if text in open_rings:
                opos, oatom = open_rings.pop(text)
                if oatom == atom_ordinal:
                    out.append(SmilesProblem(category="rings", position=pos, length=len(text),
                                             message=f"Ring closure {text} opens and closes on the same atom"))
            else:
                if atom_ordinal < 0:
                    out.append(SmilesProblem(category="rings", position=pos, length=len(text),
                                             message=f"Ring closure {text} appears before any atom"))
                open_rings[text] = (pos, atom_ordinal)
    for digit, (pos, _) in open_rings.items():
        out.append(SmilesProblem(category="rings", position=pos, length=len(digit),
                                 message=f"Ring closure {digit} is opened but never closed",
                                 suggestion="Add the matching closure digit, or remove this one"))
    return out


_RDKIT_LOGGER: Optional[logging.Logger] = None


def _rdkit_logger() -> logging.Logger:
    """Route RDKit's C++ log stream into a dedicated Python logger, once.

    Switching the sink back and forth with ``LogToCppStreams`` re-enables stderr output that
    ``DisableLog`` can no longer silence, so the sink is set once and visibility is controlled
    entirely through this logger's level and handlers.
    """
    global _RDKIT_LOGGER
    if _RDKIT_LOGGER is None:
        rdBase.LogToPythonLogger()
        lg = logging.getLogger("rdkit")
        lg.propagate = False
        lg.handlers.clear()
        lg.setLevel(logging.CRITICAL)
        rdBase.DisableLog("rdApp.*")
        _RDKIT_LOGGER = lg
    return _RDKIT_LOGGER


def _rdkit_parse_error(s: str) -> Tuple[Optional[int], str]:
    """Parse with RDKit (no sanitization) and capture its syntax message and position."""
    lg = _rdkit_logger()
    buf = io.StringIO()
    handler = logging.StreamHandler(buf)
    lg.addHandler(handler)
    lg.setLevel(logging.DEBUG)
    rdBase.EnableLog("rdApp.error")
    try:
        Chem.MolFromSmiles(s, sanitize=False)
    finally:
        rdBase.DisableLog("rdApp.error")
        lg.removeHandler(handler)
        lg.setLevel(logging.CRITICAL)
    text = buf.getvalue()
    m = _RD_POS_RE.search(text)
    pos = int(m.group(1)) - 1 if m else None
    first = next((ln.split("SMILES Parse Error:", 1)[1].strip() for ln in text.splitlines() if "SMILES Parse Error:" in ln), "")
    first = first.split(" while parsing:")[0].strip() or "RDKit could not parse this SMILES"
    return pos, first


def _atom_positions(toks: List[Tuple[str, str, int]]) -> List[Tuple[int, int]]:
    return [(pos, len(text)) for kind, text, pos in toks if kind in ("atom", "bracket")]


def _check_chemistry(s: str, toks: List[Tuple[str, str, int]]) -> Tuple[Optional["Chem.Mol"], List[SmilesProblem]]:
    mol = Chem.MolFromSmiles(s, sanitize=False)
    if mol is None:
        pos, msg = _rdkit_parse_error(s)
        return None, [SmilesProblem(category="syntax", position=pos, message=msg)]
    apos = _atom_positions(toks)
    out = []
    for prob in Chem.DetectChemistryProblems(mol):
        ptype = prob.GetType()
        if ptype == "AtomValenceException":
            idx = prob.GetAtomIdx()
            atom = mol.GetAtomWithIdx(idx)
            pos, ln = apos[idx] if idx < len(apos) else (None, 1)
            sugg = None
            if atom.GetSymbol() == "N" and atom.GetFormalCharge() == 0:
                sugg = "A four-connected nitrogen must be written as [N+]"
            elif atom.GetSymbol() == "C":
                sugg = "Carbon cannot exceed four bonds; check for a missing branch or ring closure"
            out.append(SmilesProblem(category="valence", position=pos, length=ln, atom_index=idx,
                                     message=prob.Message(), suggestion=sugg))
        elif ptype == "KekulizeException":
            idxs = list(prob.GetAtomIndices())
            pos = apos[idxs[0]][0] if idxs and idxs[0] < len(apos) else None
            symbols = {mol.GetAtomWithIdx(i).GetSymbol() for i in idxs}
            if "N" in symbols:
                sugg = "A pyrrole-type nitrogen must carry its hydrogen: write it as [nH]"
            elif symbols <= {"C"} and len(idxs) % 2 == 1:
                sugg = "An odd all-carbon aromatic ring needs a charge ([cH-] / [c+]) or should be written with explicit double bonds"
            else:
                sugg = "This ring has no valid Kekulé form; check heteroatom hydrogens and charges, or write explicit double bonds"
            out.append(SmilesProblem(category="aromaticity", position=pos, atom_indices=idxs,
                                     message=prob.Message(), suggestion=sugg))
        elif ptype == "AtomKekulizeException":
            idx = prob.GetAtomIdx()
            pos, ln = apos[idx] if idx < len(apos) else (None, 1)
            out.append(SmilesProblem(category="aromaticity", position=pos, length=ln, atom_index=idx,
                                     message=prob.Message(),
                                     suggestion="Aromatic (lower-case) atoms must be in a ring; write this atom in upper case"))
        else:
            out.append(SmilesProblem(category="syntax", message=f"{ptype}: {prob.Message()}"))
    return mol, out


# --------------------------------------------------------------------------------------------
# Repairs
# --------------------------------------------------------------------------------------------


def _try_repairs(s: str, problems: List[SmilesProblem]) -> Tuple[Optional[str], List[str]]:
    """Apply mechanical repairs implied by ``problems``. Returns (repaired, descriptions)."""
    applied: List[str] = []
    cur = s
    cats = {p.category for p in problems}

    if "characters" in cats:
        new = "".join(ch for ch in cur if not ch.isspace())
        new = new.replace("\u2013", "-").replace("\u2014", "-").replace("\u2212", "-")
        if new != cur:
            applied.append("removed whitespace / normalised dashes"); cur = new

    toks = _tokenize(cur)
    paren_probs = _check_parentheses(toks)
    if paren_probs:
        # drop unmatched ')' and append ')' for unclosed '('
        drop = {p.position for p in paren_probs if p.message.startswith("')'")}
        unclosed = [p.position for p in paren_probs if p.message.startswith("'('")]
        # a '(' with nothing substantive after it is a dangling branch: remove it rather than close it
        n_open = 0
        for pos in unclosed:
            if cur[pos + 1 :].strip("()") == "":
                drop.add(pos)
            else:
                n_open += 1
        new = "".join(ch for i, ch in enumerate(cur) if i not in drop) + ")" * n_open
        if new != cur:
            applied.append(f"balanced parentheses (removed {len(drop)}, appended {n_open})"); cur = new

    toks = _tokenize(cur)
    ring_probs = [p for p in _check_rings(toks) if "never closed" in p.message]
    if ring_probs:
        drop = set()
        for p in ring_probs:
            drop.update(range(p.position, p.position + p.length))
        new = "".join(ch for i, ch in enumerate(cur) if i not in drop)
        if new != cur:
            applied.append(f"removed {len(ring_probs)} unclosed ring-closure digit(s)"); cur = new

    if Chem.MolFromSmiles(cur) is not None:
        return cur, applied

    # aromaticity: non-ring atoms marked aromatic -> upper-case them
    mol = Chem.MolFromSmiles(cur, sanitize=False)
    if mol is not None:
        bad = [q.GetAtomIdx() for q in Chem.DetectChemistryProblems(mol) if q.GetType() == "AtomKekulizeException"]
        if bad:
            apos = _atom_positions(_tokenize(cur))
            chars = list(cur)
            for idx in bad:
                if idx < len(apos):
                    pos, ln = apos[idx]
                    if ln == 1 and chars[pos].islower():
                        chars[pos] = chars[pos].upper()
            new = "".join(chars)
            if new != cur:
                applied.append(f"wrote {len(bad)} non-ring aromatic atom(s) in upper case"); cur = new
                if Chem.MolFromSmiles(cur) is not None:
                    return cur, applied

    # aromaticity: try [nH] on one aromatic n at a time (outside brackets)
    mol = Chem.MolFromSmiles(cur, sanitize=False)
    if mol is not None and any(p.GetType() == "KekulizeException" for p in Chem.DetectChemistryProblems(mol)):
        toks = _tokenize(cur)
        n_positions = [pos for kind, text, pos in toks if kind == "atom" and text == "n"]
        for pos in n_positions:
            cand = cur[:pos] + "[nH]" + cur[pos + 1 :]
            if Chem.MolFromSmiles(cand) is not None:
                applied.append(f"wrote aromatic nitrogen at position {pos} as [nH]")
                return cand, applied

    return (cur if applied else None), applied


# --------------------------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------------------------


def diagnose_smiles(smiles: str, try_repair: bool = True) -> SmilesDiagnosis:
    """Diagnose why a SMILES string is invalid, locating each problem.

    Args:
        smiles: Input string.
        try_repair: Attempt mechanical repairs and re-validate.

    Returns:
        :class:`SmilesDiagnosis`. For a valid input ``is_valid`` is True and ``problems`` is empty.
    """
    if not _RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is required. Install with: pip install rdkit")
    s = smiles if isinstance(smiles, str) else str(smiles)
    if not s.strip():
        return SmilesDiagnosis(input_smiles=s, is_valid=False,
                               problems=[SmilesProblem(category="characters", position=0, message="Empty SMILES")])

    mol = Chem.MolFromSmiles(s)
    inner = s.strip()
    if mol is not None and not any(ch.isspace() for ch in inner):
        return SmilesDiagnosis(input_smiles=s, is_valid=True, canonical_smiles=Chem.MolToSmiles(mol))

    problems: List[SmilesProblem] = []
    if mol is not None:
        # RDKit stops at the first whitespace and treats the remainder as a title, so the
        # record parsed, but not as the molecule the string describes.
        truncated = Chem.MolToSmiles(mol)
        for i, ch in enumerate(s):
            if ch.isspace() and i >= len(s) - len(s.lstrip()) and i < len(s.rstrip()):
                problems.append(SmilesProblem(
                    category="characters", position=i,
                    message=f"Whitespace splits the record: RDKit silently parses only '{s[:i].strip()}' (→ {truncated}) and treats the rest as a title",
                    suggestion="Remove the whitespace, or split into separate records",
                ))
                break
        diag = SmilesDiagnosis(input_smiles=s, is_valid=False, problems=problems)
        if try_repair:
            repaired, applied = _try_repairs(s, problems)
            if repaired is not None:
                rm = Chem.MolFromSmiles(repaired)
                diag.repaired_smiles = Chem.MolToSmiles(rm) if rm is not None else repaired
                diag.repaired_is_valid = rm is not None
                diag.repairs_applied = applied
        return diag

    problems += _check_characters(s)
    seen_positions = {p.position for p in problems}
    toks = _tokenize(s)
    problems += [p for p in _check_tokens(toks) if p.position not in seen_positions]
    problems += _check_parentheses(toks)
    problems += _check_rings(toks)
    lexical = bool(problems)
    _, chem_probs = _check_chemistry(s, toks)
    # An RDKit syntax error that we already explained lexically adds nothing.
    if lexical:
        chem_probs = [p for p in chem_probs if p.category != "syntax"]
    problems += chem_probs
    if not problems:
        pos, msg = _rdkit_parse_error(s)
        problems.append(SmilesProblem(category="syntax", position=pos, message=msg))

    order = {c: i for i, c in enumerate(DIAGNOSTIC_CATEGORIES)}
    problems.sort(key=lambda p: (order[p.category], p.position if p.position is not None else 10**9))

    diag = SmilesDiagnosis(input_smiles=s, is_valid=False, problems=problems)
    if try_repair:
        repaired, applied = _try_repairs(s, problems)
        if repaired is not None:
            rm = Chem.MolFromSmiles(repaired)
            diag.repaired_smiles = Chem.MolToSmiles(rm) if rm is not None else repaired
            diag.repaired_is_valid = rm is not None
            diag.repairs_applied = applied
    return diag


__all__ = ["DIAGNOSTIC_CATEGORIES", "SmilesProblem", "SmilesDiagnosis", "diagnose_smiles"]
