 GC-Content Spikiness as an Anomaly Signal for Repeat-Expansion Disease Genes

## Hypothesis
Repeat-expansion disorders (e.g., Huntington's disease, Fragile X syndrome, spinocerebellar ataxias)
are caused by abnormal expansion of short repeated sequences (e.g., CAG, CTG, GAA) within specific genes.
This project tested whether such genes show a detectable "spike" in local GC-content — measurable via
sliding-window GC% analysis — that distinguishes them from housekeeping/control genes, using unsupervised
anomaly detection.

## Method
1. Fetched RefSeq mRNA sequences for each gene via the NCBI Entrez API.
2. Computed GC% across sliding windows (500 bp window, 100 bp step) per gene.
3. Extracted per-gene features: mean, max, min, std of GC%, "spikiness" (max-min range), and max z-score
   (how extreme the single highest-GC window is relative to that gene's own mean/SD).
4. Applied `IsolationForest` (scikit-learn) to flag anomalous genes, using both a full 6-feature set and a
   targeted 2-feature set (spikiness, max z-score).
5. Ran a stability check across 5 random seeds to confirm results were not seed-dependent.
6. Repeated the full pipeline on an expanded gene panel (n=25) after initial testing on a small panel (n=11).

## Results

**Initial run (n=11: 4 disease-associated, 7 housekeeping genes):**
- Targeted feature set correctly flagged 2/4 disease genes (ATXN1, FMR1) as anomalous.
- Also incorrectly flagged 2 housekeeping genes (RPL13A, TUBB) — both at the *low* end of the
  z-score distribution, revealing that IsolationForest isolates statistical extremes in both
  directions, not specifically "high" values.
- Stability check across 5 random seeds showed 100% consistent labeling — the mixed result was
  not due to randomness.

**Expanded run (n=25: 13 disease-associated, 12 housekeeping genes):**
- Only 2/13 disease genes (ATXN1, FMR1) were correctly flagged — the same two as before.
- 3 housekeeping genes were incorrectly flagged (HPRT1, PPIA, UBC).
- The apparent separation seen at n=11 did **not** hold at n=25.

## Conclusion
GC-content spikiness showed apparent separation between disease-associated repeat-expansion genes
and housekeeping genes in a small initial sample, but this did not generalize to a larger, more
diverse gene panel. This suggests the original result was likely a small-sample artifact rather
than a robust biological signal. GC%-based sliding-window analysis, while easy to compute, is
likely too indirect a proxy for repeat-expansion detection — a more direct approach (e.g., counting
specific repeat motifs like CAG/CTG/GAA rather than overall GC%) would more accurately reflect the
actual disease mechanism and is a natural next step.

## Limitations
- Small sample size, even after expansion (n=25).
- GC-spikiness is a proxy signal, not a direct measurement of repeat expansions.
- `contamination` parameter was set using prior knowledge of group sizes, making this a
  semi-supervised rather than fully blind unsupervised test.
- IsolationForest has no directional prior — it flags statistical extremes in both directions,
  which produced false positives among low-GC-variability housekeeping genes.
- No correction for gene length, transcript variant, or genomic position of the repeat locus.

## What I'd do differently
- Use direct repeat-motif counting (CAG/CTG/GAA frequency in defined windows) instead of GC%.
- Cross-reference the position of detected "spikes" against known repeat-locus coordinates from
  OMIM/ClinVar to confirm biological relevance, not just statistical anomaly.
- Use a supervised approach (given labels are already known) rather than unsupervised anomaly
  detection, and report a proper accuracy/precision/recall rather than eyeballing flagged genes.
