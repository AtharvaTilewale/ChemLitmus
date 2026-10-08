# SMARTS auditing

A structural-alert catalogue is a list of SMARTS patterns and a rule: reject any compound that matches one. The patterns were written by experts, published, converted between toolkits, patched by hand, and redistributed — and at no point validated as a *set*. `smartsaudit` treats the catalogue as data and measures it.

## What is measured

Everything is computed against a **reference population** of real molecules. The bundled default is 9,272 ChEMBL molecules (3,417 approved drugs plus a small-molecule sample); `--library` substitutes any file. The reference set is what makes "dead", "over-broad", "equivalent" and "subsumed" mean something concrete: they are facts about the catalogue *on this chemistry*.

### compile

Can RDKit parse the SMARTS at all? Does it contain a hydrogen query atom (`[H]`, `[#1]`) that can only match explicit-H molecules? Does it use recursive `$(...)` SMARTS? Pure syntax; no molecules involved.

Hydrogen detection is **positive-only**: `[!#1]` ("not hydrogen") and `[#1,#6]` ("hydrogen or carbon") do not require explicit hydrogens and are not flagged.

### breadth

Hit fraction on the reference set. Above `--breadth-threshold` (default 10%) a pattern is *over-broad*. In the ChEMBL alerts, MLSMR "Long aliphatic chain" matches 22% of the set and 24% of approved drugs — used as a rejection rule it would discard a quarter of marketed drug space.

### dead

Zero hits on the reference set, triaged into three verdicts:

- **`never-matching atom`** — some query atom, tested *alone* against a sample of reference molecules, matches no atom at all. `[N+]#[C-]` (isonitrile written with formal charges) is the canonical case: RDKit's sanitiser never produces that charge pattern, so the alert is unreachable. Almost certainly a defect.
- **`rare combination`** — every query atom is individually realisable, the combination just does not occur. Probably a working alert for rare chemistry.
- **`fires only with <prep>`** — dead under the default preparation but alive under another. A preparation requirement, not a dead rule.

Query atoms are tested in molecule context (a single-atom sub-query run through the substructure library), because recursive SMARTS cannot be evaluated against a detached atom. Hydrogen-bearing patterns are triaged against explicit-H molecules so the verdict is about the rest of the pattern.

### redundancy

Three kinds, each over the reference set:

- **duplicate** — the identical SMARTS string appears earlier in the file.
- **equivalent** — a different string with the identical hit set.
- **subsumed** — another pattern's hit set strictly contains this one's. As a *rejection* rule the subsumed pattern can never change a decision: everything it would reject is already rejected.

### sensitivity

The reference set is prepared three ways — implicit H, explicit H, kekulised — and the whole catalogue is re-matched against each. Per pattern: does the hit count change? Per catalogue: how many compounds are flagged under each preparation, and how many change pass/fail verdict relative to the default? See [Molecule preparation](molecule-preparation.md).

## How it runs

Matching uses RDKit's `SubstructLibrary` with pattern-fingerprint pre-screening, multithreaded in C++. The result is a boolean matrix (patterns × molecules) per preparation. Everything else — breadth, dead, equivalence, subsumption, sensitivity — is bit arithmetic on that matrix: equivalence is a hash of the match column, subsumption is a containment test over columns sorted by popcount. A thousand patterns against 9,272 molecules under three preparations takes about 80 s wall on a laptop, dominated by the recursive PAINS patterns.

## What the audit cannot tell you

These limits are not footnotes; they decide what the numbers mean.

**"Never fires" is an upper bound on dead rules, not a count.** Alert sets deliberately target liabilities that medicinal chemists have learned to avoid, so a well-written alert for a reactive group will fire rarely or never on drug-like chemistry. The `dead:*` triage separates the likely-defective from the merely rare; only `never-matching atom` is a strong signal of a bug.

**Equivalence and subsumption are empirical, not logical.** They are statements about hit sets on a finite library. Two patterns indistinguishable on 9,272 drug-like molecules may differ on natural products or agrochemicals. Nothing in `smartsaudit` *proves* that one SMARTS logically contains another — that would require static analysis of the query language, which is a different (and open) problem. Treat `subsumed` as "redundant on this reference set", and re-run with `--library` on your own chemistry before acting on it.

**Redundancy is not always a defect.** An alert duplicated across published sets records that two independent groups flagged the same substructure — provenance a user may want. The audit reports; it never silently deduplicates.

**Preparation sensitivity is a property of the catalogue, not an error in any one answer.** A pattern written for explicit-H molecules is not wrong; it is wrong to run it without saying so. The defect the audit exposes is the missing declaration.

**One reference set.** All empirical verdicts shift with the reference population. The bundled set is drug-like by construction; its percentages are not universal constants.

## Why this matters

Structural-alert filtering gates nearly every virtual-screening campaign and library purchase, and PAINS flags are a stated review criterion at major medicinal-chemistry journals. The filters themselves are distributed as unvalidated text and applied with undocumented settings. `smartsaudit` exists so that "we filtered for PAINS" can mean something reproducible.
