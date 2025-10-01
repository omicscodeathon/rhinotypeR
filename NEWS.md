# rhinotypeR News

## 1.3.0 (devel) — 2025-10-01
### Added
- New function: `alignToRefs()`, which aligns user sequences to the packaged
  rhinovirus prototype references using `msa::msa()` (ClustalW, ClustalOmega, or MUSCLE), with an
  option to trim alignments to the non-gap span of a chosen reference.

### Improved
- **SNPeek()** and **plotAA()**: now support zooming to specific genomic regions
  and highlighting individual sequences of interest.
- **assignTypes()**: now reports the genetic distance to the nearest prototype
  even when the query is `"unassigned"`.
- **Distance calculation**: switched to using **MSA2dist** for pairwise distance
  computation. This provides broader model support — access to all substitution
  models in **ape::dist.dna**, in addition to the `"IUPAC"` model.
- **Data access**: prototype and example datasets are now exported as packaged
  objects (`rhinovirusVP4`, `rhinovirusPrototypesVP4`) instead of being accessed
  via `system.file()`, simplifying workflows.
- Improved documentation:
  - Improved explanations and runnable examples added where feasible.
  - Visualization functions return results invisibly to avoid console clutter,
    while still allowing advanced users to capture outputs.
- Code style and formatting improved (e.g.,consistent internal helpers).

---

## 0.1.0 — Initial release
- First public release of **rhinotypeR**.
- Implements VP4/2-based rhinovirus genotyping workflow.
- Provides bundled prototype reference sequences and example datasets.
