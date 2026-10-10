# Citing ChemLitmus

If ChemLitmus contributed to published work, please cite the software. A `CITATION.cff` file is included in the repository so GitHub and reference managers can generate a citation automatically.

```bibtex
@software{chemlitmus,
  author  = {Tilewale, Atharva},
  title   = {ChemLitmus: quality control for small-molecule data},
  year    = {2026},
  version = {1.0.0},
  url     = {https://github.com/AtharvaTilewale/ChemLitmus},
  license = {MIT}
}
```

A DOI will be minted on first release through Zenodo and added here.

## Please also cite

**RDKit**, which performs every chemical computation in ChemLitmus:

> RDKit: Open-source cheminformatics. https://www.rdkit.org

**ChEMBL**, if you used the bundled reference library (the default for `smartsaudit`):

> Zdrazil B, et al. The ChEMBL Database in 2023. *Nucleic Acids Res.* 2024;52(D1):D1180–D1192. doi:10.1093/nar/gkad1004

**PubChem**, if you used `lookup`, `batch`, `download` or `iupacname --online`:

> Kim S, et al. PubChem 2025 update. *Nucleic Acids Res.* 2025;53(D1):D1516–D1525. doi:10.1093/nar/gkae1059

## Reporting results from `smartsaudit`

If you publish audit numbers, state the reference library (bundled ChEMBL subset, or your own), the ChemLitmus version, and — for any claim of redundancy or deadness — that the finding is empirical over that library. See [SMARTS auditing — what the audit cannot tell you](../concepts/smarts-auditing.md#what-the-audit-cannot-tell-you).
